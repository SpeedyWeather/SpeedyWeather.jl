# Stale GPU-graph replays in the batched Fourier transform, and the benchmark's repeat loop

> Status: **in progress**. Both fixes implemented and verified on the A40 (Fix 1 on `mg/fix-gpu-graphs-issue`, Fix 2 uncommitted on `mg/benchmark-reverse`); the full H100 benchmark run is still to do.

Date of initial draft: 2026-11-05

Base revision: `2a4024ae0c72fa80c66ceaea6809f1613382e897` (`mg/benchmark-reverse`). The files Fix 1 touches are
identical on `origin/main` (`6875c835`, which includes #1380).

## Originating prompt

> The gpu-nvidia benchmark still failed. Here's the log

(`benchmarks-2921098.log`)

> I don't understand, are both fixes needed, or would only 1) suffice?

> why is the recording of the graph ever in the timed run, shouldn't it be in the warm up?

> okay write hand off markdown for both fix 1 and 2, I want to do them in another session and process

## Revision log

> Work on executing this plan, also keep in mind that we also might need to apply a similar fix for the HIP graphs of AMD.

- 2026-11-05: initial draft.
- 2026-11-05: implementation, with the recommended answer to each decision for sign-off (scratch size in the key,
  Fix 1 as its own PR from `origin/main`, a `reset_simulation!` hook).
  - **AMD:** Fix 1 covers the HIP graphs too: both extensions call the shared `graph_key`. That path is untested
    here (no AMD GPU, and the graph tests are gated on `default_gpu_graphs`, which is `false` on AMDGPU).
  - **Correction, gpu-amd:** the benchmark does *not* force HIP graphs on. `hasproperty(suite, :gpu_graphs)` in
    `manual_benchmarking.jl` is always `false`, because `BenchmarkSuite` has no `gpu_graphs` field (and never had,
    since the switch was added in #1288). So gpu-amd benchmarks have run without graphs, and only AMD users who
    pass `gpu_graphs = true` are affected.
  - **Test, changed from the sketch:**
    - The reference transforms use their own scratch memory. The sketch's reference wrote into the scratch
      memory under test, and that can hide a stale read.
    - Every call gets new data.
    - A second round exercises the replayed graphs: the first call for a key returns the warm-up's result, not
      the graph's.
    - Before the fix, 6 of its 10 tests fail: the four calls with the second scratch memory, and both graph counts.
  - **Fix 2, changed from the sketch:** `GC.gc()` now runs after `reset_simulation!` instead of before it. Reactant's
    previous deepcopy is then freed before, not during, the timed run.
  - **Reproducers:** copied next to this plan on `mg/fix-gpu-graphs-issue`, and run with the CUDA test environment.
    - `gpu-graph-two-scratch.jl` uses a reference with its own scratch memory and different coefficients per call.
    - `gpu-graph-stale-repeats.jl` understands both key types and has the new `MODE=copy`.
  - **Correction, who is affected:** `initialize!(model)` doesn't create new scratch memory. `TransformScratchMemory`
    is allocated as `model.spectral_transform.scratch_memory`, so every `initialize!(model)` shares the transform's
    scratch memory. Only copies of the `Variables` (`deepcopy`, as in the benchmark's repeat loop) bring new scratch
    memory and expose the bug. The fix is unchanged.
  - **Reproducer bug:** the original `device_pointers!` in `gpu-graph-stale-repeats.jl` didn't recurse into `Field`s
    and `LowerTriangularArray`s, which are `AbstractArray`s with a bits eltype. So most field addresses were missing
    from a copy's address ranges, and the counts of replayed and stale graphs were undercounts. The original's
    conclusion holds (stale replays did happen). Fixed in the copy here.

## Problem description

The gpu-nvidia benchmark (job 2921098) failed after #1380 had fixed the `deepcopy` blow-up. It ran on an H100 on
2026-11-05 from 13:26 to 13:45, on branch `mg/benchmark-reverse` at `2a4024ae`, using `submit_bench.sh`
(`gpu --length 2 --repeats 10 --hi-res`). The Slurm state is FAILED.

- **benchmark100 is clean**, including PrimitiveDry and PrimitiveWet. So #1380 worked.
- **benchmark200 goes bad mid-configuration.** In ShallowWater T256, timed repeats 1–5 run normally (about 93 m/s).
  Repeats 6–10 are NaN from their first progress update.
- **benchmark201 configs 1–4 are fine; configs 5–9 are NaN**, with "NaN or Inf detected at time step 51". Config 9
  shows the pattern clearly: warm-up and repeats 1–2 are fine, then repeat 3 is NaN from 3 %. It then fails with
  `CUFFTError … CUFFT_EXEC_FAILED`.
- **The job then dies.** `free_memory!` calls `clear_fourier_graph_cache!` → `CUDA.reclaim()` →
  `device_synchronize`, which runs outside the suite's `try`. It throws
  `CUDA error: an illegal memory access was encountered (code 700, ERROR_ILLEGAL_ADDRESS)`. Most of the
  log's 339k lines are about 6,700 warnings of the form "Error while freeing … of GPU memory: CuError(code=0x2bc)",
  where 0x2bc is 700.

A run that is NaN from its very first update isn't blowing up numerically; it reads the wrong memory. The cause is
in SpeedyTransforms, not in the model.

## Background

### How the GPU graphs are keyed

Since #1147, the batched Fourier transform on CUDA and AMDGPU can record its ~2·nlat_half small cuFFTs and
gather/scatter kernels once as a GPU graph, and replay it afterwards
([gpu_graphs_common.jl](../../../SpeedyTransforms/src/gpu_graphs_common.jl), the CUDA and AMDGPU extensions in
`SpeedyTransforms/ext/`).

A captured graph bakes in the device address of every buffer it touches:

| Buffer | Owner | Stable for the graph's lifetime? |
|---|---|---|
| `field.data`, the grid field | caller (`Variables`) | no |
| `f_north`, `f_south` (forward) or `g_north`, `g_south` (inverse) | caller: `scratch_memory.north/.south` | no |
| packed work buffers, index arrays (`nlons`, offsets, …) | `GPUFourierGraphCache` | yes |
| FFT plans | `SpectralTransform` | yes |

The graphs are cached per `GPUFourierGraphCache`, a global `GRAPH_CACHES` with one entry per transform size. Within
a cache they are keyed **only by the field's address**, `graph_key(field.data) = reinterpret(UInt, pointer(data))`.
The comment above `graph_key` justifies this: "Safe because model buffers live for the whole run, so an address is
never freed and recycled under a stale graph."

That assumption is wrong. In a model the scratch memory is `vars.scratch.transform_memory`
([time_stepping/transform.jl](../../../SpeedyWeather/src/time_stepping/transform.jl)), which belongs to the
`Variables`, not the model. The `SpectralTransform`, and with it the graph cache, outlives the `Variables`:

- every `initialize!(model)` creates new `Variables`;
- so does `deepcopy(vars)`.

(Correction, see the revision log: `initialize!(model)` gives the new `Variables` the transform's own scratch memory,
`allocate(::AbstractVariable{TransformScratchMemory}, model) = model.spectral_transform.scratch_memory`, so only copies
of the `Variables`, e.g. `deepcopy(vars)`, come with new scratch memory.)

When a new copy's field lands at an address that the memory pool recycled from a freed copy, the cache finds the old
graph under that key and replays it. The replayed graph writes the Fourier output into the old scratch memory,
which has been freed, while the Legendre transform reads the new scratch memory, which wasn't written. The results
are wrong: often still finite, sometimes NaN. If the freed memory has been reused or released, the result is a
`DomainError` in a later kernel, `CUFFT_EXEC_FAILED` or an illegal address.

### Why the benchmark triggers it

Since `6a240474` (2026-11-03), `run_benchmark!`
([benchmark_suite.jl:209-258](../../../SpeedyWeather/benchmark/benchmark_suite.jl#L209-L258)) works like this:

1. It warms up on the variables of a first `initialize!(model)`.
2. It creates `initial_variables` with a second `initialize!(model)`. These are new buffers.
3. Every timed repeat runs `GC.gc()` and then `Simulation(deepcopy(initial_variables), model)`, so each repeat gets
   new buffers again, and the previous copy's GPU memory is freed. The graph cache is deliberately kept. The comment
   says that `free_memory!` "would also clear the cached Fourier graphs and so recapture them in the timed run".

Graphs are per buffer, so the warm-up's graphs never apply to the timed copies. Each timed repeat therefore either:

- captures new graphs inside the timed region (on the A40 the graph count grew 26 → 32 → 32 → 38 → … → 44 over 10
  ShallowWater T256 repeats);
- or replays stale graphs at recycled addresses.

### Evidence (login02 A40, 2026-11-05)

[gpu-graph-stale-repeats.jl](gpu-graph-stale-repeats.jl) runs the benchmark's exact sequence. It reports NaN in the
final state per repeat, and the number of graph keys inside the current copy that were captured while
`transform_memory` was at another address ("stale"):

| 10 timed repeats | Graph cache kept (`MODE=keep`, as in the benchmark) | Cache emptied before each repeat (`MODE=clear`) |
|---|---|---|
| ShallowWater T256 L1, 238 steps | 1 and 2 NaN in two runs (repeat 6 both times, as on the H100); 9 of 10 repeats replay 6 stale graphs | 0 NaN |
| PrimitiveWet T128 L8, 100 steps | 7 NaN in one run; `DomainError` in `reinitialize!(::ImplicitPrimitiveEquation)` at repeat 3 in another | 0 NaN |

- **Wrong results are routine:** almost every repeat replays stale graphs, and most of those repeats are not NaN.
  Their results are wrong but finite.
- **Damage persists:** PrimitiveWet is also NaN in repeats without a stale replay of their own. A stale graph wrote
  into freed memory that the pool had handed to an array that is still live.

[gpu-graph-two-scratch.jl](gpu-graph-two-scratch.jl) reproduces the bug deterministically, without relying on the
allocator. It transforms the same field with two `ScratchMemory`s and compares against `gpu_graphs = false`:

```
grid → spectral, graphs vs plain, scratch 1: 0.0
grid → spectral, graphs vs plain, scratch 2: 1.7944056
spectral → grid, graphs vs plain, scratch 1: 0.0
spectral → grid, graphs vs plain, scratch 2: 0.0
```

The inverse only passes because both calls use the same coefficients, so the stale scratch holds identical data. A
test must use different coefficients for the second call.

### Relation to the earlier blow-up

[2026-11-04-gpu-primitive-wet-blowup.md](2026-11-04-gpu-primitive-wet-blowup.md) found the first bug: `deepcopy`
detaching the fused `:uv_dpres_grad_spectral` views on CUDA, fixed in #1380. This is a second, independent bug,
exposed by the same repeat loop.

### Who else is affected

- **Copies of the variables on GPU:** any GPU run with `gpu_graphs = true` (the default on CUDA) that transforms with
  `Variables` whose scratch memory isn't the transform's own, i.e. copies of `Variables` such as `deepcopy(vars)`, or a
  `ScratchMemory` passed to `transform!` explicitly. `initialize!(model)` twice, or several simulations from one model,
  is safe on its own: all of them use the transform's scratch memory, so a replay at a recycled field address writes
  into the scratch memory that is in use.
- **AMDGPU:** the AMDGPU extension uses the same key. HIP graphs are off by default, so only runs with an explicit
  `gpu_graphs = true` are affected. The benchmark meant to force them on for gpu-amd (`manual_benchmarking.jl:257`),
  but that switch has no effect (see the revision log).
- **Reactant and CPU** have no graphs, so they are not affected.

## Summary of changes

Fix 1 is needed for correctness. Fix 2 is needed so the benchmark measures graph replay instead of captures or the
non-graph fallback. Fix 2 alone would also stop the benchmark failures, but users would keep the bug. Do Fix 1 first.

### Fix 1: key the graphs on every buffer the caller passes in (SpeedyTransforms)

A graph may only replay when every buffer address it baked in belongs to the current call. Buffers owned by the
cache or by the transform stay as they are.

**[gpu_graphs_common.jl](../../../SpeedyTransforms/src/gpu_graphs_common.jl)**

- Replace `graph_key(data)` (line 356) with a key over the field and both scratch arrays:

  ```julia
  # address of a device array (or of a view into one, e.g. the per-step views of the time stepping)
  @inline device_address(data) = reinterpret(UInt, pointer(data))

  # Graph-cache key: everything a captured graph bakes in that the caller provides. The scratch size is
  # included because the scatter/gather kernels index the scratch arrays with their own size, and a
  # ScratchMemory with another number of layers could be allocated at a recycled address.
  const GraphKey = Tuple{UInt, UInt, UInt, NTuple{3, Int}}
  @inline graph_key(field::AbstractField, north::AbstractArray, south::AbstractArray)::GraphKey =
      (device_address(field.data), device_address(north), device_address(south), size(north))
  ```

- Rewrite the comment above `graph_key`. Keep the explanation of why the key is an address and not the wrapper's
  identity (per-step views are new `CuArray` wrappers over the same memory). Replace the claim "Safe because model
  buffers live for the whole run" with why the scratch addresses are part of the key.
- Struct fields (lines 110 and 113):
  - `forward_execs::Dict{GraphKey, Union{Nothing, E}}`
  - `inverse_execs::Dict{Tuple{GraphKey, Bool}, Union{Nothing, E}}`
- `build_cache` (lines 169–170): construct both `Dict`s accordingly. `GraphKey` must be defined before the struct.
- Docstring of `MAX_GRAPHS` and of `GPUFourierGraphCache`: "one per distinct `field` buffer" becomes "one per
  distinct (field, scratch) buffer set".

**[SpeedyTransformsCUDAExt.jl](../../../SpeedyTransforms/ext/SpeedyTransformsCUDAExt.jl)** (lines 90 and 114) and
**[SpeedyTransformsAMDGPUExt.jl](../../../SpeedyTransforms/ext/SpeedyTransformsAMDGPUExt.jl)** (lines 102 and 127):

```julia
run_graph!(cache.forward_execs, graph_key(field, f_north, f_south), () -> forward_loop!(...), CUDA_GRAPH_BACKEND)
run_graph!(cache.inverse_execs, (graph_key(field, g_north, g_south), add), () -> inverse_loop!(...), CUDA_GRAPH_BACKEND)
```

Use `HIP_GRAPH_BACKEND` in the AMDGPU extension. Also update the CUDA extension's header comment (lines 36–43, "graphs
are keyed by the device pointer of `field.data`"), and check the AMDGPU header for a description of the key.

**Bookkeeping**

- SpeedyTransforms is already `1.0.0+DEV`, so no version bump is needed.
- `CHANGELOG.md`, first bullet under `## Unreleased`, for example: "Fix GPU graphs of the Fourier transform replaying
  with the scratch memory of freed `Variables` (e.g. after a second `initialize!(model)`), giving wrong results or
  NaNs on CUDA [#NNNN](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/NNNN)".
- Suggested as its own PR from `origin/main`, because it's a library bug that affects users. Merge `main` into
  `mg/benchmark-reverse` afterwards, or cherry-pick, to benchmark with it.

### Fix 2: run every timed repeat on the warm-up's buffers (benchmark)

Keep one working simulation, the warm-up's, and reset it to the initial state with `copy!` before each repeat.
`copy!(::Variables, ::Variables)` is correct on GPUArrays since #1380. It copies the clock and `Ref` values and skips
the write-before-read `scratch` group. The buffers keep their addresses, so the timed runs replay exactly the graphs
captured during the warm-up, with nothing captured in the timed region.

Reactant keeps the `deepcopy`: the existing comment says `copy!` would share the clock's device numbers there, and
Reactant has no GPU graphs.

**[benchmark_suite.jl](../../../SpeedyWeather/benchmark/benchmark_suite.jl)**, next to the other hooks
(`build_runner`, `synchronize_device`, …):

```julia
"""Reset `simulation` to the initial state `initial_variables` for a timed run. By default copied into the
variables of `simulation`, so every timed run steps the same buffers as the warm-up and replays the GPU graphs
captured there (they are per buffer, see SpeedyTransforms' `graph_key`)."""
function reset_simulation!(arch, simulation, initial_variables, model)
    copy!(simulation.variables, initial_variables)
    return simulation
end
```

**[reactant_runner.jl](../../../SpeedyWeather/benchmark/reactant_runner.jl)**:

```julia
# `copy!` would share the clock's device numbers, so run a deepcopy (no GPU graphs on Reactant)
reset_simulation!(::ReactantDevice, simulation, initial_variables, model) = Simulation(deepcopy(initial_variables), model)
```

**In `run_benchmark!`** (lines 242–258), keep the warm-up `simulation` instead of overwriting it:

```julia
initial_simulation = initialize!(model)     # initial state of every timed run, never run itself
initialize!(initial_simulation; period)
initial_variables = initial_simulation.variables
initial_simulation = nothing
times = Float64[]
while length(times) < suite.repeats && (isempty(times) || sum(times) + times[end] <= suite.max_repeat_time)
    GC.gc()
    simulation = reset_simulation!(architecture, simulation, initial_variables, model)
    initialize!(model.feedback, simulation.variables, model)    # a new progress meter
    synchronize_device(architecture, simulation)
    t0 = time()
    runner.time_stepping(simulation)
    synchronize_device(architecture, simulation)
    push!(times, time() - t0)
    finalize!(simulation)       # no `simulation = nothing`: the next repeat resets it
end
```

Also rewrite the comment above it (lines 233–241): the `deepcopy` reasoning only applies to Reactant now. Remove the
note about `free_memory!` and recapturing, or turn it into the explanation of why the warm-up's buffers are reused.

**Bookkeeping:** #1337 ("Benchmarks: completely revised …") is still open (not on `origin/main`), so Fix 2 is a commit
on `mg/benchmark-reverse` and needs no CHANGELOG line of its own. The benchmark directory is not package source, so
no version bump.

### Decisions for sign-off

1. Fix 1: include the scratch size in the key, as above? Recommended: it costs nothing and closes the case of a
   differently sized scratch memory at a recycled address.
2. Fix 1 as its own PR from `origin/main`, separate from #1337? Recommended.
3. Fix 2: a `reset_simulation!` hook with a Reactant method, matching how the benchmark already splits Reactant
   specifics into `reactant_runner.jl`?

## Testing and verification

Run everything from the repository root on login02, which has an A40. In the Bash tool, put Julia on the path with
`export PATH=/p/system/packages_rhel9/tools/julia/1.12.2/bin:$PATH` (`module load` doesn't work in a non-interactive
shell). Load `using CUDA` before SpeedyWeather, and don't load the system `cuda` modules.

### Fix 1

- **New test.** Add a testset to [gpu_graphs_shared.jl](../../../SpeedyWeather/test/GPU/gpu_graphs_shared.jl), gated
  like the others on `ext !== nothing && SpeedyTransforms.default_gpu_graphs(SpeedyWeather.GPU())`. This sketch
  hasn't been run yet:

  ```julia
  # Graphs bake in the scratch memory the caller passes (in a model: `vars.scratch.transform_memory`, part of
  # the Variables, not the model). The same field with another scratch memory must capture a new graph, not
  # replay the one recorded with the first.
  @testset "$prefix Graphs: another scratch memory doesn't replay a stale graph" begin
      if ext !== nothing && SpeedyTransforms.default_gpu_graphs(SpeedyWeather.GPU())
          architecture = SpeedyWeather.GPU()
          spectral_grid = SpectralGrid(; truncation = 32, nlayers = 8, architecture)
          (; NF, grid, spectrum, nlayers) = spectral_grid
          S = SpectralTransform(spectral_grid; gpu_graphs = true)
          S_reference = SpectralTransform(spectral_grid; gpu_graphs = false)
          scratches = [SpeedyTransforms.ScratchMemory(NF, architecture, grid, nlayers) for _ in 1:2]
          SpeedyTransforms.clear_fourier_graph_cache!()

          # grid → spectral, the same field
          field = rand(NF, grid, nlayers)
          reference = transform!(zeros(Complex{NF}, spectrum, nlayers), field, scratches[1], S_reference)
          for scratch in scratches
              coeffs = transform!(zeros(Complex{NF}, spectrum, nlayers), field, scratch, S)
              @test Array(coeffs.data) ≈ Array(reference.data)
          end

          # spectral → grid into the same field; different coefficients per call, as the first
          # scratch still holds the Legendre output of the first call
          out = zeros(NF, grid, nlayers)
          for scratch in scratches
              coeffs = rand(Complex{NF}, spectrum, nlayers)
              reference = transform!(zeros(NF, grid, nlayers), coeffs, scratch, S_reference)
              transform!(out, coeffs, scratch, S)
              @test Array(out.data) ≈ Array(reference.data)
          end

          n_graphs(execs) = sum(length(getfield(c, execs)) for c in values(SpeedyTransforms.GRAPH_CACHES); init = 0)
          @test n_graphs(:forward_execs) == 2     # one graph per scratch memory
          @test n_graphs(:inverse_execs) == 2
          SpeedyTransforms.clear_fourier_graph_cache!()
      end
  end
  ```

  It must **fail before the fix**: the forward comparison with the second scratch is off by about 1.8, and the
  inverse uses the first scratch's data. It must pass after the fix.
- **Existing graph tests must still pass,** in particular:
  - "per-step views of one buffer reuse a single graph", where the key must stay stable across new view wrappers;
  - "bounded graph cache over a GPU model run";
  - "add=true accumulates exactly once".
- **Commands:**

  ```bash
  julia --project=SpeedyWeather/test/GPU/CUDA -e 'using Pkg; Pkg.instantiate()'
  # graph tests only (not tried yet in this form):
  julia --project=SpeedyWeather/test/GPU/CUDA --check-bounds=yes -e 'using CUDA, SpeedyWeather, Test; include("SpeedyWeather/test/GPU/gpu_graphs_shared.jl"); include("SpeedyWeather/test/GPU/cuda_graphs.jl")'
  # full GPU suite
  julia --project=SpeedyWeather/test/GPU/CUDA --check-bounds=yes SpeedyWeather/test/GPU/runtests.jl
  ```

- **Reproducers:**
  - `julia --project=SpeedyWeather/test/GPU/CUDA docs/dev/2026-11/gpu-graph-two-scratch.jl` must print 0 (or ≈ 0)
    for all four lines.
  - `MODE=keep julia --project=SpeedyWeather/test/GPU/CUDA docs/dev/2026-11/gpu-graph-stale-repeats.jl` must show 0
    NaN repeats and 0 stale replays.
- **SpeedyTransforms CPU tests are unaffected** (graphs are GPU only), but run them as a check:
  `julia --project=SpeedyTransforms --check-bounds=yes -e 'using Pkg; Pkg.test("SpeedyTransforms")'`.

### Fix 2

- **A40 check of the mechanism.** Add a `MODE=copy` to `gpu-graph-stale-repeats.jl` that mirrors the new loop: keep
  the warm-up simulation and `copy!(simulation.variables, initial_variables)` per repeat. Then check:
  - the total number of graphs doesn't change over the 10 timed repeats;
  - 0 NaN, including without Fix 1.
- **Quick benchmark runs on login02.** These timings aren't representative, since the node is shared.
  - `julia --project=SpeedyWeather/benchmark SpeedyWeather/benchmark/manual_benchmarking.jl gpu --suites benchmark100,benchmark200 --repeats 10`
  - Check `reactant-gpu --suites benchmark100 --repeats 2` and `cpu --suites benchmark100 --repeats 2`, as both go
    through the new hook.
  - A partial run replaces that architecture's stored results with just the suites it ran. Restore
    `SpeedyWeather/benchmark/README.md` and `SpeedyWeather/benchmark/assets/benchmark_results.json` with
    `git checkout` afterwards unless the results are meant to be kept.
- **Full H100 run** with both fixes checked out: `sbatch submit_bench.sh` (QOS `gpushort`, at most 2 running jobs per
  user). The job runs whatever is checked out when it starts, so don't switch branches while it is pending. Check its
  log for:
  - 0 × "NaN or Inf detected";
  - 0 × "failed on";
  - 0 × `ILLEGAL_ADDRESS` or `0x000002bc`;
  - a small spread between fastest and slowest in the "timed run(s)" lines. Job 2921098 had ShallowWater T256 at
    median 3.82 s, fastest 0.814 s.

### Results (login03 A40, 2026-11-05)

Julia 1.12.2, CUDA test environment. Fix 1 was checked on `mg/fix-gpu-graphs-issue` and Fix 2 in a worktree of
`mg/benchmark-reverse`.

**Fix 1**

- **New testset** "another scratch memory doesn't replay a stale graph": 6 of 10 tests fail before the fix and all 10
  pass after it. The other four CUDA graph testsets pass.
- **`gpu-graph-two-scratch.jl`:** before the fix the scratch 2 line was 1.79; after it all four lines are `0.0`.
- **`gpu-graph-stale-repeats.jl`,** with the benchmark's old `deepcopy` loop (`MODE=keep`): 0 of 10 repeats NaN and 0
  stale replays for both ShallowWater T256 and PrimitiveWet T128 L8. New graphs are still captured in every repeat
  (+6 and +3), which Fix 2 addresses.
- **Full GPU suite:** all testsets except `reactant.jl` pass. That file failed with a Reactant MLIR `CompilationError`
  because the gitignored Manifest of `SpeedyWeather/test/GPU/CUDA` pins Reactant 0.2.291, while
  `SpeedyWeather/Project.toml` requires 0.2.305. In a freshly resolved environment (Reactant 0.2.305),
  `copy_variables.jl` and `reactant.jl` pass, so the whole suite passes.
- **SpeedyTransforms CPU tests:** 4118 passed and 4 errored. The errors are Enzyme `test_forward` assertions on Julia
  1.12, and the base commit `6875c835` gives the identical result.
- **Performance, base `6875c835` against the fix:** same package environment, interleaved, 2 rounds each.
  - Fourier transforms (graph replay) at T32, T64 and T128 L8: the same within 0.3%, e.g. 762.6 µs and 762.6 µs for a
    forward Fourier transform at T128 L8.
  - Full transforms: the same within 1%.
  - Time stepping (200 steps, median of 5): ShallowWater T64 0.555 s and 0.555 s, T128 1.294 s and 1.294 s;
    PrimitiveWet T32 0.252 s and 0.244–0.248 s, T64 0.391 s and 0.390–0.391 s.
  - No regression: the key computes two more pointers and hashes a longer tuple, which is negligible next to a graph
    launch.

**Fix 2**

- **`MODE=copy` without Fix 1:** 0 NaN and 0 graphs captured during the 10 repeats. Every repeat replays the warm-up's
  graphs: 6 for ShallowWater T256, 3 for PrimitiveWet's time loop. The same holds with Fix 1.
- **`gpu --suites benchmark100,benchmark200 --repeats 10`:**
  - 0 × "NaN or Inf detected", 0 × "failed on", 0 × illegal address or `CUFFT_EXEC_FAILED`.
  - Small spreads between repeats, e.g. median 2.69 s, fastest 2.69 s, slowest 2.70 s.
- **`reactant-gpu` and `cpu`** with `--suites benchmark100 --repeats 2`: all configurations ok, no NaN.
- The README and `benchmark_results.json` were restored afterwards.
- The full H100 run is still to do.

## Documentation changes

- Code comments as listed above (`graph_key`, both extension headers, the benchmark repeat loop).
- Optional: a sentence in `docs/src/architectures_gpu.md`, section "GPU Graphs". Graphs are cached per set of
  buffers, so each new `Variables` for the same model captures new ones. After many re-initialisations the cache is
  full (`MAX_GRAPHS`) and transforms fall back to the slower path without graphs. `clear_fourier_graph_cache!()`
  resets it.

## Known limitations

- **Dead graphs stay cached.** Graphs for buffers that were freed stay in the cache until `clear_fourier_graph_cache!`.
  Each new set of `Variables` for the same model adds about 6 (ShallowWater T256) to about 25 (PrimitiveWet T128 L8)
  graphs. After `MAX_GRAPHS = 64` per direction, transforms on new buffers run without graphs. That is correct but
  slower, e.g. for an ensemble of many simulations from one model.
- **Invalid gpu-nvidia results.** gpu-nvidia results from runs after `6a240474` (Nov 3) are invalid, even where no NaN
  was reported. None are committed: the stored gpu-nvidia entry in `benchmark_results.json` is from 2026-10-20.
- **The Oct 20 baseline is not affected by this bug.** It used a single timed run after a second `initialize!(model)`.
  The benchmark code of that time (last changed in `e410eed05`, Oct 17) has no `deepcopy`. A second `initialize!`
  shares the transform's scratch memory (see the revision log), so a replay at a recycled field address there was
  correct. The 12 NaN warnings in its log (`benchmarks-2618118.log`) are therefore not stale replays. Whether they are
  genuine blow-ups of high-resolution L16/L24 configurations is still unverified.

## Future work

- **Evict dead graphs:** for example, tie cache entries to the scratch memory's lifetime with a finalizer, evict the
  least recently used entry when full, or clear a transform's graphs when `Variables` are created for its model.
- **Make `free_memory!` failures non-fatal** in `run_benchmark_suite!`. It runs outside the `try`, although after a
  sticky CUDA error the context is unusable anyway.
- **Rerun benchmarks:** rerun gpu-nvidia and, if available, gpu-amd with both fixes. Then decide whether the Oct 20
  gpu-nvidia baseline needs a rerun.
