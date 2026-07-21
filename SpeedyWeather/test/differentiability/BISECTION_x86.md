# Take-over guide (x86): localizing the Julia 1.12 "partial sizes" Enzyme abort

**Read first:** [`HANDOFF_full_diff_CI_x86.md`](HANDOFF_full_diff_CI_x86.md) (full diagnosis) and the
`L30` section at the top of [`MWE_REDUCTION_STATE.md`](MWE_REDUCTION_STATE.md) (history).

## The one open bug
`full_diff_CI.jl` — full `PrimitiveWetModel` `time_step!` (physics included) under **parameter AD**
(`Duplicated(model)`) — **passes on Julia 1.10**, **aborts on Julia 1.12 x86** with:
```
DiffeGradientUtils.cpp:433: addToDiffe: Assertion `0 && "unhandled accumulate with partial sizes"'
VT: i64 start=0 size=4 storeSize=8  ... extractvalue [3 x i64] ...     signal 6 (Aborted), exit 134
```
i.e. Enzyme compiling the **reverse of a store** cannot accumulate a 4-byte piece into an 8-byte
slot — smells like a `Float32`→`ComplexF32` partial store in the parameter-AD reverse.

**Platform:** aborts fast only on **x86** (assertions-enabled Enzyme_jll). On **ARM it HANGS** — do
all of this on an x86 Julia 1.12 node (or the CI below).

Already fixed and NOT the target: the **state-AD** (`Const(model)`) leapfrog/fused-view trigger
(fixed by view-preserving `make_zero(::Variables)`). This is a *second, distinct* trigger on the
`Duplicated(model)` path.

## How to run the bisection (x86 CPU node — NOT via CI)

Everything runs directly on an x86 Julia 1.12 node. `bisect_pwet_param.jl` runs exactly ONE
`autodiff` call, reusing the `test_utils.jl` harness (`ADSimulation`/`ADseed`) so each variant is the
`Duplicated(model)` twin of a `Const(model)` test in `primitivewet.jl`. It prints
`RESULT: SUCCESS ...` when the reverse compiles; otherwise it aborts (exit 134).

**Step 0 — set up the env ONCE** (do NOT let parallel SLURM tasks instantiate concurrently):
```
julia --project=SpeedyWeather/test/differentiability \
      SpeedyWeather/test/differentiability/setup_diff_env.jl
```

**Interactive node — sequential sweep** (each variant isolated; prints a PASS/ABORT table):
```
SpeedyWeather/test/differentiability/run_bisect_local.sh            # all 11 variants
SpeedyWeather/test/differentiability/run_bisect_local.sh A2_param C2_paramtend   # a subset
```

**Batch / extensive — SLURM array** (one task per variant; edit partition/account/Julia first):
```
sbatch SpeedyWeather/test/differentiability/bisect_pwet_param.slurm
# logs: SpeedyWeather/test/differentiability/_logs/bisect_<jobid>_<taskid>.log
# exit 134 in a task = that variant's component is the culprit
```

**Single variant, directly:**
```
julia --project=SpeedyWeather/test/differentiability \
      SpeedyWeather/test/differentiability/bisect_pwet_param.jl <variant>
```

(Deliberately no GitHub-Actions job for this — run it on the x86 node instead.)

## Decision tree (what each variant tells you)

The abort is a **compile-time** event, so it depends only on *which function* and *which args are
`Duplicated`* — not on resolution or seed values (trunc=5,nlayers=1 == full_diff_CI reproduces it).

```
A1_state (full time_step!, Const model) ......... expect PASS  (confirms make_zero fixed state AD)
A2_param (full time_step!, Duplicated model) .... expect ABORT (== full_diff_CI, the target)

If A2 aborts, split on physics:
  B1_dynonly (dynamics_only, Duplicated model)
    ├─ ABORT  → trigger is in the DYN-CORE parameter-AD reverse → look at C1/C3/C4/C5/C6
    └─ PASS   → trigger REQUIRES PHYSICS → look at C2/C7/C8 (and greenhouse/sea_ice)

Per-component (all Duplicated(model)); the RED one(s) localize the store:
  C1_dyntend     dynamics_tendencies!          (Const(model) twin PASSES in primitivewet.jl)
  C2_paramtend   parameterization_tendencies!  (never tested under param AD — prime physics suspect)
  C3_implicit    implicit_correction!
  C4_diffusion   horizontal_diffusion!
  C5_updateprog  update_prognostic!            (state-AD twin fixed by make_zero; param AD untested)
  C6_transform   transform!                    (analytic-adjoint rule; should PASS)
  C7_ocean       ocean_timestep!               (physics)
  C8_land        land_timestep!                (physics)
```
Note `dynamics_tendencies!` / `implicit_correction!` / `transform!` already **pass with `Const(model)`**
in `primitivewet.jl` on 1.12 — so if their `C*` (Duplicated) twin aborts, the trigger is
specifically the **model-shadow (parameter-AD) reverse** of that component, not its state AD.

## Once the culprit component is known — narrow to the store
1. Re-run just that component under `Enzyme.API.printall!(true)` (or dump the LLVM around the
   `[3 x i64]` / `size=4 storeSize=8` store) to find the offending `store`. Expect a `Float32`
   written into a `ComplexF32` slot (real part of a complex spectral field, or a real model-parameter
   gradient flowing into a complex field).
2. Bisect *within* the component: copy it locally (the method that cracked errors #1/#2 — see
   `MWE_REDUCTION_STATE.md`) and remove sub-calls until the abort disappears.

## Candidate fixes (once localized)
- **Real→complex partial store:** restructure that op so the reverse accumulates the full
  `ComplexF32` (combine real/imag; avoid storing a bare real into a complex slot under AD; or route
  through a helper Enzyme can reverse). This is the most likely real fix.
- **Unneeded parameter gradient:** if a specific component's parameter is not a differentiation
  target, mark it `EnzymeRules.inactive` / `inactive_type` (as done for `SpectralTransform`).
- **Intrinsic to the fused 3-D view stores under param AD:** a custom reverse rule for the offending
  kernel (like the `transform!` analytic-adjoint rule) may be needed.

## Fallback
If not fixable SpeedyWeather-side, reduce a minimal Enzyme(+SpeedyWeather) reproducer of the
"partial sizes" store, file it on Enzyme.jl, and `@test_broken`/skip `full_diff_CI` on Julia ≥ 1.11
with the issue link so CI stays green on the genuine fixes. Then wire the real
`full_diff_CI`/Enzyme CI to also run 1.12 (currently `CI_Enzyme.yml` is 1.10-only).

## Files for this task (all on branch `mg/enzyme-1-12-ci-x86`)
- `bisect_pwet_param.jl` — the parametrized one-autodiff-per-variant driver (this guide's variants).
- `setup_diff_env.jl` — one-time env setup (Pkg.develop + instantiate + precompile). Run first.
- `run_bisect_local.sh` — sequential sweep on an interactive node; prints a PASS/ABORT table.
- `bisect_pwet_param.slurm` — SLURM array (one task per variant) for extensive/batch runs.
- `full_diff_CI.jl` — the actual failing test (== variant A2_param).
- `_ts_const_model.jl` / `_ts_dynonly.jl` — standalone A1/B1 (kept; referenced by the handoff).
- `primitivewet.jl` — the existing `Const(model)` per-component tests these variants mirror.
- `test_utils.jl` — `ADSimulation` / `ADseed` / `initialize_with_spinup!` harness.
- `HANDOFF_full_diff_CI_x86.md`, `MWE_REDUCTION_STATE.md`, `ENZYME_DIFF_PROTOCOL.md` — diagnosis/history.
- `MWE_enzyme_runtime_ntuple.jl`, `MWE_enzyme_fat_aggregate.jl`, `MWE_enzyme_mutable_transform.jl`
  — dependency-free reproducers of the *already-fixed* errors #1/#2 (for the upstream issue #3275).
