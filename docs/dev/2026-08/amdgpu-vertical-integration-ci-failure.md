# AMDGPU CI failure in `vertical_integration.jl` — investigation notes

> Status: **investigating**. Root cause not yet confirmed; likely pre-existing/environment, unrelated to PR #1193's own diff. Minimal reproducers written 2026-08-17, not yet run on AMD hardware (none available in this session).

Date of note: 2026-08-14

Written to resume the investigation quickly on Monday (2026-08-17).

## Update 2026-08-17: minimal reproducers + `clamp` hypothesis ruled out

Three tiered, standalone reproducer scripts were written in `docs/dev/2026-08/repro/`:

- `tier1_initialize_model.jl` — just `initialize!(model)` on a default `PrimitiveWetModel`,
  no `run!`, no CPU comparison, no vertical-integration checks. Smallest SpeedyWeather-level repro.
- `tier2_set_kernel.jl` — isolates the exact failing kernel launch: `set!` on a bare 3D `Field`
  with a `JablonowskiVorticity` functor, skipping model/ocean/land/greenhouse-gas init entirely.
  Distinguishes "bug is in `set!`'s kernel path" vs. "bug is elsewhere in `initialize!`'s chain
  (orography, `JablonowskiTemperature`, `ConstantRelativeHumidity`, ...)".
- `tier3_standalone_kernel.jl` — **no SpeedyWeather dependency at all**. Verbatim copy of the
  `JablonowskiVorticity` functor and `set_field_3d_kernel!` kernel, launched directly via
  `KernelAbstractions` + `AMDGPU`. If this alone crashes, it's a pure AMDGPU.jl/GPUCompiler bug
  reproducible without building SpeedyWeather at all — worth filing upstream against
  JuliaGPU/AMDGPU.jl directly. Verified the kernel/functor logic is correct by running tier3's
  body against the KernelAbstractions **CPU** backend (`Array` instead of `ROCArray`) in this
  session — it runs and produces a sane value. The AMDGPU-specific crash itself could **not** be
  verified in this session (no AMD GPU hardware available locally) — **run these three scripts on
  the buildkite AMDGPU runner next** and note which ones crash vs. pass.

**`clamp(X, 0, 1)` Float64-promotion hypothesis (next-steps item 3 from 2026-08-14) is ruled out.**
Checked directly on Julia 1.12: `clamp(::Float32, ::Int, ::Int)` stays `Float32`
(`promote_type(Float32, Int64) == Float32` on this Julia version, unlike the older
`Float32`/`Int64` → `Float64` promotion rule assumed in the original hypothesis). Also ran
`@code_warntype` on `JablonowskiVorticity(...)(λ, φ, η)` with `Float32` arguments end-to-end:
every intermediate value is concretely `Float32`, no `Union`/`Any`, no dynamic dispatch, `x^2`
uses `Base.literal_pow` throughout. The functor itself is fully type-stable on the CPU front end,
so there is no obvious Float64-promotion or type-instability bug left to find by further staring
at `initial_conditions.jl:457-511` — a line-by-line audit (next-steps item 3) is likely a dead end
unless AMDGPU's GPU-specific codegen diverges from CPU inference in a way `@code_warntype` can't
show (possible, but would need `@device_code_llvm`/`AMDGPU.@device_code_warntype` on real hardware
to confirm).

This shifts weight toward **hypothesis 2 (toolchain/version flakiness)** from the 2026-08-14 notes:
the `jl_f_throw_methoderror` / `hostcall_device_write_args!` signature in the original error looks
like GPUCompiler failing to prove a runtime-dispatch/exception path unreachable — something that's
sensitive to exact ROCm/LLVM/AMDGPU.jl versions and GPU architecture (this runner vs. LUMI's
MI250x), not necessarily a bug in SpeedyWeather's source at all. The tier3 standalone script is the
fastest way to confirm this: if it reproduces with only AMDGPU+KernelAbstractions loaded, the fix
(if any) belongs upstream, not in this repo.

## Update 2026-08-17 (later): found the likely unreachable-throw-path source, applied a fix

The MRE tests were wired into `.buildkite/pipeline.yml`/`SpeedyWeather/test/GPU/runtests.jl` (only
`mre_amdgpu_crash.jl` runs) and pushed via the `gd/test-buildkite-on-main` scratch branch (`DO NOT
MERGE`, pinned to Julia 1.12 for AMDGPU to match the original failure environment). A run against
the **full**, non-restricted test suite on that AMDGPU runner (before the restriction landed, or a
related run) surfaced this, which did not appear on the CUDA runner:

```
┌ Info: Global hostcalls detected!
│ - Source: MethodInstance for SpeedyWeather.gpu_set_field_3d_kernel!(..., JablonowskiVorticity{Float32}, ...)
│ - Hostcalls: [:malloc_hostcall, :malloc_hostcall]
│ Use `AMDGPU.synchronize(; stop_hostcalls=true)` to synchronize and stop them.
└ Otherwise, performance might degrade if they keep running in the background.
```

...also seen for `JablonowskiDivergence`, and separately (already-`inbounds=true`-annotated)
`_initialize_particles_kernel!` and the very large `column_parameterizations_kernel!`. These are
not test failures (tests still pass) but confirm these kernels retain a dynamic/exception code path
that AMDGPU has to keep a resident "hostcall" listener thread around for — the same general failure
class (GPUCompiler being unable to eliminate a supposedly-unreachable dynamic-dispatch/exception
path) as the `InvalidIRError` crash, just not (on this runner, this run) fatal.

**Likely source, specific to `set_field_3d_kernel!`/`JablonowskiVorticity`/`JablonowskiDivergence`:**
`acos(X)` and `sqrt(...)` are domain-checked in Base (`DomainError` for `|x|>1` / `x<0`). `X` is
`clamp`-ed to `[0,1]` and `cosηᵥ = cos(...)` is asserted positive by the existing comment, so both
are provably in-domain at every actual call site — but the compiler can't derive that from `clamp`
or `cos`'s general signature, so it keeps the `throw(DomainError(...))` branch (and everything that
branch needs: GC frame, exception object allocation) in the compiled kernel. That matches the error
signature from the original `InvalidIRError` almost exactly (`jl_f_throw_methoderror`,
`julia.new_gc_frame`, `hostcall_device_write_args!`, ...).

**Fix applied** (`SpeedyWeather/src/dynamics/initial_conditions.jl`,
`SpeedyWeather/src/variables/set.jl`):
- Wrapped the `acos`/`sqrt` expressions in `JablonowskiVorticity`/`JablonowskiDivergence` in
  `@fastmath`, which uses the domain-check-free fast-math variants and removes the throw path.
- Added the missing `@kernel inbounds = true` to `set_field_3d_kernel!` (every other performance
  kernel in the codebase has this per `CLAUDE.md`; this one didn't).
- Bumped `SpeedyWeather` to `0.22.0+DEV` and added a `CHANGELOG.md` entry per convention.
- All existing CPU unit tests re-run and pass, including `test/dynamics/initial_conditions.jl`'s
  `ZonalWind` testset (which exercises both functors) — `@fastmath` doesn't change results for
  in-domain inputs.

**Not fixed / still open:** `column_parameterizations_kernel!` (`parameterizations/tendencies.jl`)
and `_initialize_particles_kernel!` (`dynamics/particle_advection.jl`) *already* have
`inbounds=true` and still triggered hostcalls, proving `inbounds=true` alone isn't sufficient in
general — there's almost certainly other domain-checked math (`sqrt`/`log`/`acos`/`asin`/`^`) buried
somewhere in the many parameterization components `column_parameterizations_kernel!` dispatches
into (radiation, convection, boundary layer, ...). That's a much bigger audit across many files, not
attempted here. `set_field_3d_kernel!` was the one actually implicated in the `InvalidIRError`
stacktrace, so it was the priority.

**Not yet verified**: none of this has run on real AMDGPU hardware yet (no AMD GPU available in this
session). Next step is pushing this fix on `gd/test-buildkite-on-main` (or #1193) and checking
whether (a) the "Global hostcalls detected" messages for `set_field_3d_kernel!` disappear, and (b)
the original `vertical_integration.jl` `InvalidIRError` is gone once the full test suite is
restored (see the `TEMPORARY` marker in `SpeedyWeather/test/GPU/runtests.jl`).

## Originating context

Two separate threads from a 2026-08-14 session, both worth resuming:

1. **PR #1147 ("Towards HIP graphs")**, branch `gd/hip-graphs` — mostly wrapped up.
2. **PR #1193 ("test AMD CI in 1.12")**, branch `gd/fix-amdgpu-triangular-loop` — CI failure investigation, **not yet resolved**. This is the main thing to pick back up.

---

## 1. PR #1147 status (branch `gd/hip-graphs`)

Two outstanding review comments were discussed:

- **r3783127046** (@maximilian-gelbrecht): `GraphBackend`'s `synchronize` field should reuse the
  codebase's existing `synchronize(architecture)` instead of a bespoke per-backend closure.
  → **Done and committed**: commit `e8bfc4ee` "use SpeedyWeather's synchronize(architecture) in
  GraphBackend instead of a custom closure", already pushed to `origin/gd/hip-graphs`. Changed
  `gpu_graphs_common.jl`, `SpeedyTransformsCUDAExt.jl`, `SpeedyTransformsAMDGPUExt.jl`. No GPU
  hardware was available in the session to actually run the GPU test suite against this — only
  verified the three files parse. **TODO: run `SpeedyWeather/test/GPU/runtests.jl` (or at least
  the `gpu_graphs` tests) for real on CUDA/AMDGPU hardware before considering this closed.**

- **r3499445082** (@milankl): can `_fourier_batched!` be moved out of the extensions into
  `src/fourier.jl` for `GPUArray`? → Already answered in the PR thread: tried, reverted because it
  broke Enzyme's CPU-only autodiff tests (a GPU-dispatched method in `src/` loads unconditionally
  even with no GPU backend, which confuses Enzyme's static activity analysis). No code change
  needed here, just reply/resolve the thread.

  **Snag found while checking this**: the header comment in
  `SpeedyTransforms/ext/gpu_graphs_common.jl` (lines 8–13) cites
  `docs/dev/2026-08/gpu-graphs-common-interface.md` as the write-up for this Enzyme rationale —
  **that file does not exist in the repo**. Either write it (per the CLAUDE.md implementation-plan
  convention) or fix the comment to point somewhere real. This overlaps with the "leftover
  documentation references" cleanup a reviewer already flagged on the PR.

- **Not yet done**: `SpeedyTransforms/Project.toml` version hasn't been bumped at all on this
  branch (`0.2.0`, identical to `main`) despite `src/spectral_transform.jl` and the ext files
  having changed relative to `main`. Per CLAUDE.md convention this needs a `+DEV` (or `-DEV` if
  the `GraphBackend` field rename counts as breaking — it's an internal/extension-only API so
  probably just `+DEV`) bump before merge.

---

## 2. PR #1193 CI failure — main thing to resume

### Reproduction

Buildkite CI on `amdgpu2.luraess.com` (Julia 1.12, AMDGPU.jl v2.7.0, commit `b444215e` = tip of
`gd/fix-amdgpu-triangular-loop` / PR #1193) fails in `SpeedyWeather/test/GPU/vertical_integration.jl`
with an `InvalidIRError` while compiling a GPU kernel:

```
InvalidIRError: compiling MethodInstance for SpeedyWeather.gpu_set_field_3d_kernel!(...) resulted
in invalid LLVM IR
Reason: unsupported call to an unknown function (call to julia.new_gc_frame / julia.get_pgcstack /
julia.push_gc_frame / julia.gc_alloc_bytes / julia.pop_gc_frame / julia.get_gc_frame_slot)
Reason: unsupported dynamic function invocation (call to kernel_state(), convert, getproperty,
reinterpret, unsafe_load, hostcall_device_write_args!, hostcall_device_trigger_and_return!)
Reason: unsupported call to an unknown function (call to jl_f_throw_methoderror)
```

Full log was pasted into the session transcript (2026-08-14) if the raw text is needed again.

### What's actually failing (confirmed by reading source)

The kernel is `set_field_3d_kernel!` in `SpeedyWeather/src/variables/set.jl:233` (KernelAbstractions
apparently reports it internally as `gpu_set_field_3d_kernel!`, not a name that exists in our
source — just KA's generated binding for the GPU-specialized variant, not a mismatch worth
chasing). Call chain, confirmed via the stacktrace:

```
initialize!(model)                              # PrimitiveWetModel
 → initialize!(vars, initial_conditions, model)  # initial_conditions.jl:52 / :60 / :451
 → set!(vars, geometry; vorticity=JablonowskiVorticity(...), divergence=JablonowskiDivergence(...), static_func=true)
 → set!(var::AbstractField3D, f::Function, geometry, S; static_func=true)   # set.jl:198
 → _set_function_3d!(var, f, londs, latds, σ_levels_full)                  # set.jl:216 (static_func==true branch → runs directly on GPU, no CPU fallback)
 → launch!(..., set_field_3d_kernel!, var, londs, latds, σ_levels_full, f, kernel_func)
 → kernel body: var[ij,k] = kernel_func(var[ij,k], f(londs[ij], latds[ij], σ_levels_full[k]))
     where f::JablonowskiVorticity{Float32}, kernel_func = (a,b)->b  (add=false)
```

The crashing kernel instantiation's argument types match this exactly: `ROCDeviceMatrix{Float32,1}`
(var), 3×`ROCDeviceVector{Float32,1}` (londs/latds/σ_levels_full), `JablonowskiVorticity{Float32}`
(f), `var"#220#221"` (the `kernel_func` closure).

### Key finding: this is NOT caused by PR #1193's diff

```
git diff c343adf6 pr-1193 --stat
 .buildkite/pipeline.yml                                            | 2 +-
 CHANGELOG.md                                                       | 2 ++
 SpeedyWeather/src/time_stepping/implicit/implicit_primitive_equations.jl | 7 ++++++-
```

`SpeedyWeather/src/variables/set.jl` and `SpeedyWeather/src/dynamics/initial_conditions.jl` are
**byte-identical** between `main` (`c343adf6`) and the PR #1193 tip — PR #1193 only touches the
triangular-loop implicit-kernel workaround, `.buildkite/pipeline.yml` (added Julia 1.12 GPU CI),
and the changelog. So this failure would reproduce on plain `main` on this same runner/toolchain —
**it is pre-existing/environment, not a regression from #1193's actual change.**

### Why it's plausible this is the *same bug class* as an earlier, already-fixed issue

`JablonowskiVorticity`'s functor (`initial_conditions.jl:470-489`) has this comment right above the
suspicious line, added in PR #1171 ("Run GPU CI with 1.12", commit `f947d9c6`, authored by
Maximilian Gelbrecht):

```julia
# NOTE: `x^(3//2)` does NOT stay in Float32. Base computes a fractional power as
# exp2(e * log2(x)) and promotes to Float64 internally, so the kernel ends up needing
# Float64 log2/exp2, which AMDGPU cannot compile (GPUCompiler segfaults in check_ir!).
# cos is strictly positive here (η in [0, 1] gives cos in [0.39, 1]), so x*sqrt(x) is equivalent.
cosηᵥ = cos((η - η₀) * π / 2)
ζ = -4 * u₀ / radius * (cosηᵥ * sqrt(cosηᵥ)) * sind(φ) * cosd(φ) * (2 - 5sind(φ)^2)
```

That's the *exact* symptom class we're seeing now (dynamic dispatch / GPU-incompatible IR that
AMDGPU's compiler can't eliminate), already fought once in this very function. Two live
hypotheses, not yet distinguished:

1. **A second instance of the same pattern** still lurking in `JablonowskiVorticity` or
   `JablonowskiDivergence` (`initial_conditions.jl:457-511`) that the `x^(3//2)` fix didn't cover —
   worth auditing every `^`, `sqrt`, `acos`, `clamp`, `exp` call in both functors for anything that
   could promote to Float64 or hit a non-`@inline`-able generic method. `clamp(X, 0, 1)` with `Int`
   literals `0`/`1` against a `Float32` `X` is a mild suspect (goes through `promote()`) but is a
   very common, normally-safe pattern.
2. **Toolchain/version flakiness**: this runner has AMDGPU.jl v2.7.0 + Julia 1.12 + whatever
   ROCm/LLVM + GPU it has; LUMI (where this reportedly passed) likely has a different AMDGPU.jl
   version, ROCm version, and GPU architecture (MI250x). GPUCompiler's ability to dead-code-
   eliminate an unreachable dynamic-dispatch path is known to be sensitive to exact ROCm/LLVM
   versions — matches the general AMDGPU-compiler-fragility theme already documented in PR #1147
   (ROCm's stream-capture validator "doesn't reliably reject illegal-to-capture operations").

### Next steps (start here Monday)

1. **Isolate**: run `vertical_integration.jl` (or just `initialize!` on a tiny `PrimitiveWetModel`)
   against plain `main` on the *same* CI runner/toolchain. If it also fails, this is fully
   pre-existing/environment and should be reported/tracked as its own issue rather than blocking
   #1193 — probably worth flagging to Maximilian since he wrote the earlier fix and best knows this
   runner's quirks.
2. If it fails only on this branch and not `main` (unexpected given the diff above, but worth
   double-checking rather than trusting `git diff` blindly — e.g. check `Manifest.toml`/dependency
   version drift between the two CI runs too, not just source diff).
3. Audit `JablonowskiVorticity`/`JablonowskiDivergence` (`initial_conditions.jl:457-511`) line by
   line for any remaining fractional-power / Float64-promotion / type-instability pattern, applying
   the same `x*sqrt(x)`-style rewrite used for `cosηᵥ`.
4. If no source-level culprit is found, check whether pinning/bumping AMDGPU.jl to a different
   version on the CI runner changes the outcome (the Buildkite log shows AMDGPU v2.7.0 marked `⌃`
   i.e. upgradable).
5. Once root-caused, remember: CLAUDE.md requires a `CHANGELOG.md` entry and a submodule version
   bump for whatever gets touched to fix this.

### Useful local refs from this session

- `pr-1193` — local ref for `refs/pull/1193/head` (fetched from `upstream`), tip commit `b444215e`.
- Base commit for comparison: `c343adf6` (`v0.22.0`, tag on `upstream/main`).
