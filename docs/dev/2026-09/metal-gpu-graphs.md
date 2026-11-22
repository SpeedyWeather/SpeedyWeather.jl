# Enable and test Metal `gpu_graphs` (MPSGraph-fused batched Fourier transform)

> Status: **planned**. A candidate fix and test suite were drafted and run locally on
> 2026-09-17, then reverted the same session pending Metal-hardware verification; nothing
> from that attempt is on `main` or pushed. This document picks the work back up for a
> future PR.

Date of initial draft: 2026-09-18

Base revision: `3119ee10` (`gg/metal-support`)

## Originating prompt

> please write a summary like the other plans in docs/dev/ so we can pick up enabling and
> testing metal graphs in a future PR

## Revision log

- **2026-09-18, initial draft.** Written retroactively to capture same-day local
  experimentation (commits `8a5acf3a`..`abaef016` on `gg/metal-support`) that was tried and
  then reverted, so the investigation isn't lost. No code changes accompany this draft.

## Problem description

Metal is the only GPU backend where `gpu_graphs` is forced off unconditionally:

```julia
# SpeedyTransforms/ext/SpeedyTransformsMetalExt.jl
default_gpu_graphs(::Metal.MetalBackend) = false
```

introduced in `459b4439` ("Disable gpu_graphs for Metal GPU CI debugging"). Unlike CUDA/HIP
graphs, Metal's accelerated path doesn't use graph capture/replay of kernel launches; it
batches each north/south ring pair's Fourier transform into a single `MPSGraph` (see the
module-level comment in `SpeedyTransformsMetalExt.jl`) to cut CPU/GPU round-trips. That
`MPSGraph` path exists and is exercised indirectly, but is never turned on by default and has
no dedicated equivalence test — so it's untested against the plain per-ring path and its
correctness on real hardware is unverified. `docs/src/architectures_gpu.md` currently
documents this state plainly: "On **Metal**, `gpu_graphs` is **disabled by default** and not
currently supported."

## Background

- `default_gpu_graphs` is backend-dispatched in `SpeedyTransforms`; CUDA defaults to `true`,
  AMDGPU (HIP graphs) defaults to `false` pending broader hardware verification (see
  `docs/src/architectures_gpu.md`'s "GPU Graphs" section), and Metal is hardcoded `false`.
- The Metal fused-graph implementation (`forward_loop_fused!` / `inverse_loop_fused!` in
  `SpeedyTransformsMetalExt.jl`) builds a cached `MPSGraph` per `(SpectralTransform, nlayers)`
  pair, encodes it onto an `MPSCommandBuffer`, and calls `Metal.commit!(cmdbuf)` to submit the
  work asynchronously. The buffers it just wrote (`ring_complex_both` / `ring_real_both`) are
  then read straight back out via `copyto!` with no explicit wait in between.
- Commit `a475134e` ("use Simulation(model) to reproduce Milan's segfault") swapped
  `initialize!(model)` for `Simulation(model)` in `SpeedyWeather/test/GPU/primitive_wet.jl`,
  suggesting a colleague (Milan) had hit a segfault on Metal that the author was trying to
  reproduce around the same time as this investigation — likely the same class of issue as
  below, though the connection isn't confirmed in the commit history.
- PR [#1217](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/1217) is referenced (in
  the now-reverted synchronization fix, see below) as the origin of observed forward-only
  `gpu_graphs` failures that motivated suspecting a missing synchronization.

### What was tried and reverted on 2026-09-17

Four commits landed and were reverted within the same session, in this order:

1. `8a5acf3a` — **"Wait for Metal MPSGraph command buffer to complete before reading FFT
   results."** Added `KernelAbstractions.synchronize(S.architecture)` right after each
   `Metal.commit!(cmdbuf)` call in `forward_loop_fused!` and `inverse_loop_fused!`. Rationale:
   `commit!` submits asynchronously, and nothing guaranteed the subsequent `copyto!` calls
   were ordered after the GPU finished writing — a plausible-but-wrong-values race, suspected
   cause of the PR #1217 failures. This same commit also, as a side effect, temporarily
   commented out most of `SpeedyWeather/test/GPU/runtests.jl` (`kernels_GPU.jl`,
   `broadcasting.jl`, `spectral_transform.jl`, `interpolate.jl`, `set.jl`,
   `vertical_integration.jl`, `barotropic.jl`, `shallowwater.jl`, `primitive_wet.jl`,
   `gpu_graphs_shared.jl`, and `MetalGPU/metal.jl`) "to save CI time while iterating" — marked
   in comments as needing to be restored before merging.
2. `92d1ba0b` — **"include metal graph tests."** Added
   `include("metal_graphs.jl")` to the `:Metal` branch of `runtests.jl`.
3. `ab081c7b` — **"add forgotten file for metal graph tests."** Added
   `SpeedyWeather/test/GPU/metal_graphs.jl` itself (the previous commit referenced it before
   it existed) with two testsets, deliberately *not* gated on `default_gpu_graphs` since
   that's exactly what's under test:
   - **Fourier equivalence**: for each grid in a local `grid_list` (`FullGaussianGrid`,
     `OctahedralGaussianGrid`, `OctahedralClenshawGrid`), build one `SpectralTransform` with
     `gpu_graphs = false` and one with `gpu_graphs = true`, transform the same random field/
     coefficients both ways, and check `≈` agreement (`rtol = sqrt(eps(Float32))`). Also checks
     that the graph/buffer caches (`FOURIER_GRAPH_CACHES` / `FOURIER_BUFFER_CACHES`) are
     actually populated, and that repeated `transform!` calls reuse the cached graph/buffers
     without growing the cache.
   - **Full model run**: construct a `PrimitiveWetModel` with `gpu_graphs = true` explicitly,
     run 5 steps, assert `feedback.nans_detected == false`.
4. `2938515d`, `3a11d1ff`, `dcd46f8a` — straight `git revert` of all three commits above, back
   to back, same session. **The commit messages give no reason.** Because nothing was pushed
   (`origin/gg/metal-support` still points at the much older `63820c86`), this reads as
   in-session experimentation that was backed out rather than a fix that was tested and found
   wrong — but that is inferred from the absence of a stated reason, not confirmed. Two
   plausible, non-exclusive explanations, neither verified:
   - The `synchronize` fix was never actually run against Metal hardware in this environment
     (no Apple Silicon available to the assistant), so it was reverted as unverified rather
     than disproven.
   - Shipping the fix would have required also restoring the CI-time-saving comment-outs in
     `runtests.jl` from commit `8a5acf3a`, which were explicitly marked not-for-merge; reverting
     everything together was the simplest way to leave the branch clean.

Only `abaef016` ("add documentaiton for metal support") was kept, which documents the
*current* (unfixed, `gpu_graphs` off) state in `README.md` and
`docs/src/architectures_gpu.md` — it does not depend on any of the reverted code.

## Summary of changes

None yet — this document is scoped to planning. The reverted diffs above (recoverable via
`git show 8a5acf3a`, `git show 92d1ba0b`, `git show ab081c7b`) are a reasonable starting point
for a future PR, not a final answer.

## Testing and verification

Nothing has been run against real Metal hardware as part of this investigation; all work so
far happened in an environment without Apple Silicon access. Before re-attempting:

1. Reproduce the PR #1217 forward-only `gpu_graphs` failure on Metal hardware first, to
   confirm there is in fact a live bug to fix (vs. one already fixed by unrelated changes
   since).
2. Re-apply the `KernelAbstractions.synchronize(S.architecture)` calls (or a cheaper
   equivalent — see Known limitations) and confirm the equivalence test in `metal_graphs.jl`
   passes reliably (not just once — race conditions can pass intermittently).
3. Restore `runtests.jl` to its un-commented state (all the includes commit `8a5acf3a`
   disabled) before merging; the CI-time-saving shortcut must not land on `main`.
4. Run the full `SpeedyWeather/test/GPU/runtests.jl` suite (`MetalGPU/metal.jl` plus the new
   `metal_graphs.jl`) on real Metal CI, not just the new testset in isolation, since
   synchronization bugs elsewhere in the Metal extension could interact.
5. Only flip `default_gpu_graphs(::Metal.MetalBackend)` to `true` once the above is green on
   Metal CI across multiple runs (per the AMDGPU precedent in `architectures_gpu.md`, which
   is still `false` "pending broader hardware verification" despite being stable on some
   hardware) — Metal's is a smaller, more homogeneous hardware population, but it still isn't
   a default to flip on a single passing run.

## Documentation changes

Once `gpu_graphs` is verified and enabled by default on Metal, update the "On **Metal**"
paragraph in `docs/src/architectures_gpu.md`'s "GPU Graphs" section (currently states it is
"disabled by default and not currently supported") to match, mirroring how the CUDA/AMDGPU
paragraphs are written.

## Known limitations

- `KernelAbstractions.synchronize(S.architecture)` is a full device sync — correct but
  possibly more heavyweight than necessary (e.g. waiting on the specific command buffer /
  queue rather than the whole device). Worth profiling once correctness is confirmed, since
  the whole point of the fused-graph path is reducing overhead.
- The suspected race (unsynchronized async `commit!` followed by an immediate read) was never
  confirmed as the actual root cause of the PR #1217 failures on real hardware — it's a
  plausible hypothesis backed by the async-submission semantics of `MPSCommandBuffer`, not a
  reproduced-and-fixed bug.

## Future work

- Carry out the verification steps above on Metal hardware (Apple Silicon CI or a
  contributor's machine) and land the fix + tests properly, including bumping
  `SpeedyTransforms`'s version per `CLAUDE.md`'s submodule-versioning rules.
- Consider whether the equivalence test belongs in `metal_graphs.jl` alongside
  `cuda_graphs.jl` / `hip_graphs.jl` conventions, or whether those three should eventually
  share more structure (the module comment in `SpeedyTransformsMetalExt.jl` already notes
  Metal doesn't share `gpu_graphs_common.jl`'s CUDA/HIP graph-capture machinery, so full
  unification may not be worthwhile).
