# Handoff: `full_diff_CI.jl` Enzyme codegen assertion on Julia 1.12 (x86)

**Goal:** make `SpeedyWeather/test/differentiability/full_diff_CI.jl` pass on Julia 1.12.
It **passes on 1.10** (≈13 min compile) and **aborts on 1.12** with an Enzyme codegen assertion.

## ⚠️ Platform note (why this is an x86 task)
- **x86 (CI, assertions-enabled Enzyme_jll):** aborts **fast** with the assertion below (exit 134).
  This is the fast, debuggable signal — iterate here. It took the GitHub Action CI 6 minutes to hit this error. 
- **ARM / Apple Silicon (local dev):** release Enzyme_jll hits the SAME bug but **HANGS** (compile
  loops forever, RSS ~140 MB, no assertion). You cannot see the assertion or iterate on ARM.
- ⇒ **Do all reproduction/localization/verification on an x86 Julia 1.12 node.**

## The failure (verbatim from CI)
```
Run: julia --project=SpeedyWeather/test/differentiability SpeedyWeather/test/differentiability/full_diff_CI.jl
[ Info: Time step changed ... (setup + run! Hour(6) succeed) ...
 VT: i64 idxs:{} start=0 size=4 storeSize=8 val=  %.fca.1.extract = extractvalue [3 x i64] %2299, 1
julia: .../Enzyme/DiffeGradientUtils.cpp:433: ... addToDiffe(...):
       Assertion `0 && "unhandled accumulate with partial sizes"' failed.
signal 6 (Aborted), exit 134
  addToDiffe            (DiffeGradientUtils.cpp:433)
  visitCommonStore      (AdjointGenerator.h:1219)
  visitStoreInst        (AdjointGenerator.h:899)
  CreatePrimalAndGradient (EnzymeLogic.cpp:4521)
```
It aborts during `CreatePrimalAndGradient` — i.e. **compiling the reverse of a `store`**. Enzyme
cannot accumulate a gradient into a value with "partial sizes": a `[3 x i64]` (24-byte) aggregate,
storing a **4-byte** piece into an **8-byte** slot. `size=4 storeSize=8` looks like a `Float32`
(4 B) stored into a `ComplexF32` (8 B) slot — i.e. a real→complex *partial* store. The spectral
prognostic is `ComplexF32`, so this smells like a complex-valued store in the reverse.

## What `full_diff_CI` does
`PrimitiveWetModel(trunc=5, nlayers=1)` (FULL model — physics included), `run!(Hour(6))`, then
```julia
dvars  = make_zero(vars);  dmodel = make_zero(model)
dvars.prognostic.{vorticity,divergence,humidity,temperature,pressure} .= 1 + im
autodiff(set_runtime_activity(Reverse), time_step!, Const,
    Duplicated(vars, dvars),
    Duplicated(model.time_stepping, dmodel.time_stepping),
    Duplicated(model, dmodel))          # <-- PARAMETER AD (Duplicated model)
@test sum(to_vec(dvars)[1]) != 0
```
Key: it is **parameter AD** (`Duplicated(model)`) of the **whole** `time_step!` (dynamics + physics
+ implicit + leapfrog). The file's own top comment already says it was "broken in Julia 1.11" — so
this is a long-standing 1.12 problem, not introduced by the recent Enzyme-rule work.

## Diagnosis so far (done on ARM via hang-vs-complete + one x86 CI datapoint)
- **Not a type-analysis issue:** `maxtypeoffset!(4096)` does NOT help (verified — identical
  hang/abort). So restoring the maxtypeoffset workaround is not the fix.
- **Not the analytic-adjoint `transform!` rule:** that rule is an AD boundary; `transform!(::Variables)`
  reverse compiles+runs on 1.12 in isolation (`_mwe5_pwet_fast.jl`), and full_diff_CI passes on 1.10
  with the identical rule.
- **The fused view-backed prognostic + leapfrog `get_step` is ONE confirmed trigger, and it is now
  FIXED for state AD:** `vars.prognostic.vorticity` is an `LTA{ComplexF32,3,SubArray}` (a view of the
  3-D `vars.fused.prognostic.parent`). `update_prognostic!` (leapfrog) does `get_step(var,1/2/lf)` →
  view-of-a-view, and `leapfrog_kernel!` stores into it. Enzyme's default `make_zero(vars)`
  **materialised the shadow to a plain `Array`**, mismatching the view-backed primal → the store
  reverse over view/Array types failed (`Duplicated(view-LTA, Array-LTA)` on barotropic;
  "partial sizes" on 1.12). **Fix committed:** view-preserving `make_zero(::Variables)` in
  `SpeedyWeather/ext/SpeedyWeatherEnzymeExt.jl` (deepcopy + zero via fuse parents; shadow stays
  `LTA{SubArray}`). Verified on 1.12 (ARM): `update_prognostic!` with **`Const(model)` (state AD)**
  compiles + autodiff-returns (was hang). Also fixes the barotropic `Duplicated(view,Array)` error.
- **A SECOND trigger remains in `full_diff_CI` and is NOT fixed by make_zero:** CI ran commit
  `634df776` (which HAS the make_zero fix) and still aborts identically. The remaining trigger is in
  the **parameter-AD** (`Duplicated(model)`) reverse of the full `time_step!` — my make_zero test
  only exercised **state AD** (`Const(model)`). This second store is the one to localize on x86.

## Suggested localization on x86 (fast aborts → quick iteration)
Reduce `full_diff_CI`'s single `autodiff` call and re-run on x86; each variant aborts in seconds if
it hits the bug, else compiles. Bisect what makes the abort appear/disappear:
1. `Const(model)` (state AD) instead of `Duplicated(model)` — does the abort disappear? If yes, the
   trigger is specifically the **parameter-AD (model shadow)** path (most likely, given make_zero
   fixed state-AD `update_prognostic!`). Ready-made script: `_ts_const_model.jl`.
2. `dynamics_only=true` model (no physics) with `Duplicated(model)` — abort or not? Isolates physics
   vs dyn-core. Ready-made: `_ts_dynonly.jl`.
3. Differentiate individual components with `Duplicated(model)` on 1.12 and see which aborts:
   `dynamics_tendencies!`, `implicit_correction!`, `parameterization_tendencies!`,
   `update_prognostic!`. (The pwet dyn-core tests in `primitivewet.jl` only did `Const(model)`, so
   parameter AD of these is untested.) `_driver_pwet.jl` is a good harness to adapt.
4. Once the component is known, dump the offending LLVM store: run with
   `ENV["ENZYME_..."]`/`Enzyme.API.printall!(true)` or inspect around the `[3 x i64]` /
   `size=4 storeSize=8` store — expect a `Float32`→`ComplexF32` partial store (real part of a complex
   spectral field, or a real model-parameter gradient flowing into a complex field).

## Candidate restructures (once localized)
- If it is a real→complex partial store: restructure that op so the reverse accumulates the full
  `ComplexF32` (e.g. avoid storing a bare real component into a complex slot under AD; combine
  real/imag; or route through a helper Enzyme can reverse).
- If a specific model component's parameter gradient is the culprit and is not needed, consider
  `EnzymeRules.inactive` / `inactive_type` for that component (as done for `SpectralTransform`) —
  only if that parameter is genuinely not a differentiation target.
- If it is intrinsic to differentiating the fused 3-D view stores under parameter AD, a custom
  reverse rule for the offending kernel (like the `transform!` analytic-adjoint rule) may be needed.

## Fallback if not fixable SpeedyWeather-side
Reduce a minimal Enzyme+SpeedyWeather (or Enzyme-only) reproducer of the "partial sizes" store and
file it on Enzyme.jl (it is an Enzyme codegen bug); meanwhile `@test_broken`/skip `full_diff_CI` on
Julia ≥ 1.11 with a comment linking the issue, so CI is green with the genuine fixes below.

## Fixes already committed on `mg/enzyme-1-12` (all validated, keep these)
- Analytic-adjoint Enzyme reverse rules for the whole `transform!` (chunked + batched) — fixes the
  chunked/fused transform reverse (error #3). `SpeedyTransformsEnzymeExt.jl`. (#1151)
- `inactive_type(::Type{<:SpectralTransform})` — fixes the large-aggregate type-analysis failure
  (error #2) exactly; removed the `maxtypeoffset!(4096)` workaround. Also fixed the 1.10
  `UV_from_vordiv!` KA GC bug as a bonus.
- `reset`/`add` kwargs on `_legendre!`/`_fourier!` (unified) enabling allocation-free adjoint pullbacks.
- `to_vec(::Variables)` (fused view-aware) in `SpeedyWeatherFiniteDifferencesExt.jl` — fixes the FD
  `Array→SubArray` harness error (barotropic went 24 to_vec errors → 0). Validated `_tovec_variables_check.jl` 8/8.
- View-preserving `make_zero(::Variables)` in `SpeedyWeatherEnzymeExt.jl` — fixes the leapfrog
  state-AD trigger + the barotropic `Duplicated(view,Array)` error.
- AD-rules unit tests un-gated on 1.12 (`SpeedyTransforms/test/runtests.jl`); `test_reverse` +
  chunked-vs-batched + accumulation tests added; full SpeedyTransforms suite green on 1.12.

## Ready-made scripts in this directory (for the x86 session)
- `full_diff_CI.jl` — the failing test itself.
- `_ts_const_model.jl` — full `time_step!`, `Const(model)` (state AD). (ARM: hung; x86: should show whether state AD aborts.)
- `_ts_dynonly.jl` — full `time_step!`, `dynamics_only=true`, `Duplicated(model)`.
- `_full_ci_offset_test.jl` — full_diff_CI body + `Enzyme.API.maxtypeoffset!(ARGS[1])` (confirmed offset doesn't help).
- `_viewshadow_test.jl` — `update_prognostic!` with a view-preserving shadow (the make_zero fix, state AD) — compiles on 1.12.
- `_tovec_variables_check.jl` — validates the `to_vec(::Variables)` fix (8/8).
- `_adjoint_check.jl`, `_chunked_rules_test.jl`, `_test_reverse_check.jl` — transform-rule validation.
- Full running notes / history: `MWE_REDUCTION_STATE.md` (see the L30 section at the top).

Reproduce on x86:
```
julia +1.12 --project=SpeedyWeather/test/differentiability SpeedyWeather/test/differentiability/full_diff_CI.jl
# expect: fast Aborted (exit 134) with the "unhandled accumulate with partial sizes" assertion
```
