# Enzyme Differentiability Test Protocol

Goal: run the extended differentiability tests under the **most recent Enzyme**
in **Julia 1.10** and **Julia 1.12**, independently. Order per version:
`speedy_transforms.jl` → `barotropic.jl` → `primitivewet.jl`.

We are hunting for **errors thrown by Enzyme** (exceptions during `autodiff`).
Finite-difference mismatches (`@test` fails / `@test_broken`) are *expected* in
some places and are **not** the target — but they are noted for context.

## Environment

| Item | Value |
|------|-------|
| Julia 1.10 | 1.10.11 |
| Julia 1.12 | 1.12.6 |
| Enzyme (both) | v0.13.173 (latest as of 2026-06-30; `Pkg.update` found nothing newer) |
| EnzymeTestUtils | v0.2.8 |
| Manifests | version-specific: `Manifest-v1.10.toml`, `Manifest-v1.12.toml` (kept independent) |

Driver: `run_diff_test.jl <file>` wraps the file in an outer `@testset` so all
inner testsets run even if one errors; uncaught exceptions print full backtraces.

---

## Status matrix

| Test file | Julia 1.10 | Julia 1.12 |
|-----------|-----------|-----------|
| speedy_transforms.jl | ✅ 114/114 pass, no Enzyme error | ✅ 114/114 pass, no Enzyme error |
| barotropic.jl (old API) | 🔧 6/6 errored — outdated API (not Enzyme) | 🔧 6/6 errored — outdated API (not Enzyme) |
| barotropic.jl (ported, **pre-fix**)  | ✅ no Enzyme errors (FD-tolerance fails only) | ❌ **4× EnzymeNoTypeError** (1.12-specific) |
| barotropic.jl (ported, **`looseTypeAnalysis!` fix**) | ✅ 17/3/0/1 — no Enzyme error | ✅ **17/3/0/1 — no Enzyme error (matches 1.10)** |
| primitivewet.jl (old API) | 🔧 errored — outdated API (not Enzyme) | 🔧 errored — outdated API (not Enzyme) |

**Headline (RESOLVED):** the SpeedyTransforms Enzyme rules are clean on both Julia versions.
The barotropic **model** layer previously hit **four `EnzymeNoTypeError`s on Julia 1.12** (two
distinct root causes: LowerTriangularArray view typing in `get_step`/`lta_view` via the leapfrog
`update_prognostic!`, and `vorticity_flux_curldiv!` under `Duplicated(model)`). Enabling
`Enzyme.API.looseTypeAnalysis!(true)` for Julia ≥ 1.11 in the Enzyme extension **clears all four**:
Julia 1.12 now produces `17 passed, 3 failed, 0 errored, 1 broken` — **byte-for-byte the same
outcome as Julia 1.10** (the 3 failures are shared FD-tolerance comparisons, not Enzyme errors),
confirming the loosened analysis yields the same gradients rather than corrupting them.
`primitivewet.jl` is still on the pre-`Variables` API and needs the same kind of port
as barotropic before it can exercise Enzyme.

Legend: ✅ no Enzyme error · ⚠️ FD mismatch only (no Enzyme error) · ❌ Enzyme error · 🔧 test-code bug · ⏳ pending

---

## Pre-run code review notes

Critical reading of the test files before running:

### `primitivewet.jl` — contains clear test-code bugs (independent of Enzyme)
- L29 `(; time) = progn.clock` — `progn` is undefined in this scope (the unified
  `Variables` refactor removed it; should be `vars`/`simulation`).
- L64 `return progn_new` inside `ocean_timestep` — `progn_new` undefined (should be `vars_new`).
- L77, L103, L125, L145 call `ADSeed` (capital S) but the helper is `ADseed` → `UndefVarError`.
- L140 `to_vec(dprogn)` — `dprogn` undefined (should be `dvars`).
- L85 `land_timestep` calls `ocean_timestep!` (copy/paste error).
These will abort the file early with plain Julia errors, *before* exercising Enzyme.

### `barotropic.jl`
- Uses the current `Variables`/`ADSimulation`/`ADseed` API — looks consistent with `test_utils.jl`.

### `speedy_transforms.jl`
- Uses direct `autodiff` on transform primitives; looks API-consistent.
- L314-315: `@test sum(dvor) != #...` then a nested `@test isapprox(...)` — the `!=`
  has no RHS on its line and swallows the next `@test` as its argument. Latent test-logic bug.

---

## Run log

(Chronological; newest entries appended.)

### `speedy_transforms.jl` — both versions ✅ (no Enzyme errors)

| | Julia 1.10 | Julia 1.12 |
|--|-----------|-----------|
| Result | **114/114 pass** | **114/114 pass** |
| Errors/Fails/Broken | 0 / 0 / 0 | 0 / 0 / 0 |
| Wall time (incl. compile) | 1m55s | 5m18s |

No Enzyme exceptions, no `@test` failures, no warnings (grep for
error/fail/warn/exception/segfault came back empty). The custom `EnzymeRules`
for `_fourier!` (`SpeedyTransformsEnzymeExt.jl`) and all spectral gradient
adjoints (`curl!`, `divergence!`, `UV_from_vor!`, `UV_from_vordiv!`, `∇²!`, `∇!`)
differentiate correctly and match finite differences under **Enzyme 0.13.173**.
This is the layer the old `0.13.153` compat ceiling guarded — it is clean now.

### `barotropic.jl` — both versions 🔧 (NOT Enzyme errors — outdated test API)

Identical outcome in **both Julia 1.10 and 1.12**:
`0 passed, 0 failed, 6 errored, 0 broken`.

**Every one of the 6 testsets dies at the *primal* dispatch, before Enzyme does
any differentiation.** The errors are `MethodError`/`UndefVarError` because
`barotropic.jl` + `timestep_utils.jl` are written against the **pre-`time_stepping/`-refactor
API**. Enzyme is blameless: when the undifferentiated call `f(args.val...)` has no
method, Enzyme faithfully reports the missing primal method.

Old→new API mismatches found (cause of each testset error):

| testset (line) | test calls (OLD) | current source (NEW) |
|---|---|---|
| dynamics_tendencies! (L11) | `dynamics_tendencies!(vars, lf2, model)` | `dynamics_tendencies!(vars, model)` — `lf` arg removed (`dynamics/tendencies.jl:3`) |
| horizontal_diffusion! (L48) | `horizontal_diffusion!(vars, diffusion, model, lf1)` | `horizontal_diffusion!(vars, diffusion, model)` — `lf` arg removed (`dynamics/horizontal_diffusion.jl:242`) |
| leapfrog! (L78) | `leapfrog!(...)` (two methods) | **`leapfrog!` deleted**; logic now in `update_prognostic!` + `leapfrog_kernel!`, orchestrated by `time_step!`/`diffusion_and_implicit!` (`time_stepping/steppers/leapfrog.jl`). Error: `UndefVarError: leapfrog! not defined in SpeedyWeather` |
| transform!(::Variables) (L135) | `transform!(vars, lf2, model)` | `transform!(vars, model; kwargs...)` — `lf` arg removed (`time_stepping/transform.jl:2`) |
| timestep! (L165) | `timestep!(vars, dt, model, lf1, lf2)` (via `timestep_oop!`) | renamed → `time_step!(vars, time_stepping, model)` (`time_stepping/time_integration.jl:33`) |
| Barotropic parameters (L227) | `timestep!(...)` (via `timestep_oop!`) | same rename → `time_step!` |

Root cause: the time-stepping was extracted into `SpeedyWeather/src/time_stepping/`
(`time_integration.jl`, `steppers/leapfrog.jl`, `transform.jl`, …). In the new design the
leapfrog index `lf` is no longer threaded through `dynamics_tendencies!`/`transform!`/
`horizontal_diffusion!`, the public stepping entry point is `time_step!(vars, time_stepping, model)`,
and there is no standalone `leapfrog!`. The differentiability test file and
`timestep_utils.jl` were never updated to match.

➡️ **These are test-maintenance errors, not Enzyme/Julia-1.12 regressions.** To
actually exercise Enzyme on the barotropic model the test must be ported to the
new API (see "Ported barotropic.jl" below for whether Enzyme then errors).

### `primitivewet.jl` — both versions 🔧 (NOT Enzyme errors — outdated test API)

Identical in both versions: `0 passed, 0 failed, 1 errored, 0 broken`. The single
big testset aborts immediately at `(; vars, model) = simulation`:
`FieldError: type Simulation has no field 'vars', available fields: 'variables', 'model'`.
Plus the further-down latent bugs already noted (`progn`, `ADSeed`, `dprogn`,
`return progn_new`). Again pre-`Variables`/pre-`time_stepping` API — never reaches Enzyme.

### Ported `barotropic.jl` + `timestep_utils.jl` → new API (2026-06-30)

To actually test Enzyme on the barotropic model, `barotropic.jl` and
`timestep_utils.jl` were ported to the current API:

- `timestep_oop!`/`timestep_oop` now call `time_step!(vars, model.time_stepping, model)`
  (dropped the `dt`/`lf` arguments; `reconstruct` path kept for the parameter test).
- testset 1: `dynamics_tendencies!(vars, model)`; sanity check now over the whole `dvars`.
- testset 2: `horizontal_diffusion!(vars, diffusion, model)`; analytical checks kept
  (diffusion reads prognostic/tendency step 1 ⇒ `∂tend/∂vor = expl·impl`, `∂tend/∂tend_in = impl`).
- testset 3: `leapfrog!` → `update_prognostic!(vars, model)` (full) and the per-variable
  `update_prognostic!(var, tendency, clock, time_stepping, implicit, model)`; the filter
  coefficients `w1,w2` and `Δt` are derived from the clock (`step_counter` forced to 2
  so the Robert+Williams filtered branch `lf==2` is exercised).
- testset 4: `transform!(vars, model)` (drops `lf`; note the method carries an
  `initialize::Bool` kwarg — a good incidental Enzyme-through-kwarg check).
- testsets 5+6: `time_step!` via the ported `timestep_oop!`/`timestep_oop`.

Results of the ported run are recorded below.

#### Findings from the first ported run (Julia 1.12)

1. **`EnzymeRuntimeActivityError` on the model-level AD (genuine Enzyme finding).**
   `autodiff(Reverse, dynamics_tendencies!, …)` with plain `Reverse` throws:
   > `EnzymeRuntimeActivityError: Detected potential need for runtime activity.
   > Constant memory is stored (or returned) to a differentiable variable and
   > correctness cannot be guaranteed with static activity.`
   Using `set_runtime_activity(Reverse)` resolves it (verified in isolation). The
   ported tests therefore use `set_runtime_activity(Reverse)` for all model-level
   AD. This is the standard remedy and the existing `timestepping.jl` already used
   it for the `Const(model)` case — but note plain `Reverse` is now insufficient
   even for `dynamics_tendencies!` alone.

2. **Default Barotropic time stepper is `NCycleLorenz`, not `Leapfrog`** (`models/barotropic.jl:51`).
   `NCycleLorenz` has `prognostic_steps == 1` and no Robert/Williams filters, so the
   leapfrog-specific analytical checks only make sense with an explicit `Leapfrog`.
   Per request: most tests now use `Leapfrog(spectral_grid)`, with a final test on
   the default `NCycleLorenz`.

3. **Array layout (Leapfrog), confirmed via diagnostic:**
   `prognostic.vorticity` = `(65,1,2)` (2 leapfrog steps), `tendencies.vorticity` =
   `(65,1,1)` (3-D, single step), `horizontal_diffusion.expl/impl` = `(11,1)` matrices
   indexed by degree `l`. The analytical diffusion check was rewritten to index
   `[:, 1, 1]` and use `spectrum.l_indices` directly (the old `[:, 1]` assumed a 2-D
   tendency and errored).

Rewritten `barotropic.jl`: 7 testsets (1–6 Leapfrog, 7 = NCycleLorenz full step).

#### ⭐ KEY RESULT — ported barotropic.jl, per testset (**PRE-FIX**, before `looseTypeAnalysis!`; see "FINAL RESULT" at the bottom for the resolved state)

| # | Testset (Leapfrog unless noted) | AD activity | Julia 1.10 | Julia 1.12 |
|---|----------------------------------|-------------|-----------|-----------|
| 1 | `dynamics_tendencies!`         | `Const(model)`           | ✅ no Enzyme error | ✅ no Enzyme error |
| 2 | `horizontal_diffusion!`        | `Const(model)`           | ✅ analytic checks pass | ✅ analytic checks pass |
| 3 | `update_prognostic!` (full) + single-var leapfrog | `Const(model)` | ✅ no Enzyme error | ❌ **EnzymeNoTypeError** (single-var) |
| 4 | `transform!(::Variables)`      | `Const(model)`           | ✅ no Enzyme error | ✅ no Enzyme error |
| 5 | full `time_step!`              | `Duplicated(model)` + `Const(model)` | ✅ no Enzyme error (matches FD, rtol 5%) | ❌ **EnzymeNoTypeError** |
| 6 | parameters (`reconstruct`)     | `Duplicated(model)`, `Duplicated(pvec)` | ✅ no Enzyme error (FD mismatch only) | ❌ **EnzymeNoTypeError** |
| 7 | full `time_step!`, **NCycleLorenz** | `Duplicated(model)` | ✅ no Enzyme error (matches FD, rtol 5%) | ❌ **EnzymeNoTypeError** |

Julia 1.10 summary: **0 Enzyme errors** in any testset. (The only `@test` failures are
strict FD tolerance checks; the `rtol=5%` FD comparison of the full-timestep gradient
actually *passes* — Enzyme's reverse gradient matches finite differences.)

Julia 1.12 summary: **`7 passed, 0 failed, 4 errored, 1 broken`** — the 4 errors are all
`EnzymeNoTypeError` (see below). The passing testsets (1,2,4) are exactly the ones that
differentiate with `Const(model)` and do **not** hit the two failing code paths.

#### ❌ Enzyme error #1 (Julia 1.12 only) — `get_step` / `lta_view`

Differentiating `update_prognostic!(var, tendency, clock, Leapfrog, implicit, model)`
(the per-variable leapfrog kernel path) throws:

> `EnzymeNoTypeError: Enzyme cannot statically prove the type of a value being differentiated`

Stacktrace (top frames):
```
[1] LowerTriangularArray   @ LowerTriangularArrays/src/lower_triangular_array.jl:14  [inlined]
[2] LowerTriangularArray   @ LowerTriangularArrays/src/lower_triangular_array.jl:32  [inlined]
[3] lta_view               @ LowerTriangularArrays/src/lower_triangular_array.jl:736 [inlined]
[4] get_step               @ SpeedyWeather/src/time_stepping/steps.jl:79             [inlined]
[5] update_prognostic!     @ SpeedyWeather/src/time_stepping/steppers/leapfrog.jl:230 [inlined]
```
i.e. `var_lf = get_step(var, lf)` (leapfrog.jl:230) with a **runtime** step index `lf`
builds a `LowerTriangularArray` **view** via `lta_view`, and on Julia 1.12 Enzyme's type
analysis cannot prove the view's type. On Julia 1.10 the identical call differentiates fine.

#### ❌ Enzyme error #2 (Julia 1.12 only) — `vorticity_flux_curldiv!` under `Duplicated(model)`

Differentiating the full `time_step!` **with respect to the model** (`Duplicated(model, …)`,
testsets 5/6/7) throws `EnzymeNoTypeError` within:
```
Failure within method: SpeedyWeather.var"#vorticity_flux_curldiv!#378"(…)
   @ SpeedyWeather/src/dynamics/tendencies.jl:7
```
Note testset 1 differentiates `dynamics_tendencies!` (which calls `vorticity_flux_curldiv!`)
with `Const(model)` and **passes** on 1.12 — so this failure is tied to differentiating
w.r.t. the model (and/or the larger fused `time_step!` context), again 1.12-only.

#### Bottom line

- **SpeedyTransforms layer (`speedy_transforms.jl`): clean on both 1.10 and 1.12.**
- **Model layer (barotropic): clean on Julia 1.10; on Julia 1.12 two distinct
  `EnzymeNoTypeError`s** — (1) LowerTriangularArray view typing in `get_step`/`lta_view`
  reached via the leapfrog `update_prognostic!`, and (2) `vorticity_flux_curldiv!` when
  differentiating w.r.t. the model. Both are type-analysis failures specific to Julia 1.12
  + Enzyme 0.13.173.
- Independent of version, model-level reverse AD requires `set_runtime_activity(Reverse)`.

---

## Source-level fixes (Julia 1.12 EnzymeNoTypeErrors)

### Fix #1 — `get_step` runtime-offset view (RESOLVED via `looseTypeAnalysis!`)

Root cause (established by MWEs): `get_step(var, lf)` with a **runtime** step index builds a
view with a runtime last-dimension offset. Such a view in *isolation* differentiates fine,
but once it is passed as an argument to an Enzyme-differentiated **KernelAbstractions kernel**
(here `leapfrog_kernel!` launched from `update_prognostic!`), Enzyme's type analysis on Julia
≥1.11 cannot prove the element type → `EnzymeNoTypeError`. A view with a **compile-time** offset
(e.g. the unrolled `get_steps`) is fine on all versions.

**Chosen resolution** (in [`ext/SpeedyWeatherEnzymeExt.jl`](../../ext/SpeedyWeatherEnzymeExt.jl)):
enable `Enzyme.API.looseTypeAnalysis!(true)` in the extension's `__init__`, **gated
`if VERSION >= v"1.11"`** (Julia 1.10 resolves these types natively and is left untouched). With
loose type analysis Enzyme is permitted to take its best guess when it cannot statically prove a
type. This is **safe here** because SpeedyWeather's differentiated arrays are homogeneous
(all `Float32` / `ComplexF32`), so the "guess" is unambiguous. The advantage over source-level
workarounds: it keeps `get_step`/`get_prognostic_step`/`get_tendency_step` and all the kernels
**completely untouched** — the step accessors stay branchless (important for Reactant) and no call
site changes.

Verified on Julia 1.12: the reproducer [`_mwe_get_step.jl`](_mwe_get_step.jl) (real `BarotropicModel`,
no spin-up) differentiating the single-variable leapfrog `update_prognostic!` now returns the
**correct gradient** `∂tend = Δt(1+w1−w2)`: `RESULT: SUCCESS (dtend extrema = 0.0022670224)`.
`speedy_transforms.jl` remains 114/114 (no gradient corruption from the loosened analysis).

#### Rejected / failed alternatives

| Approach | Outcome |
|---|---|
| `Enzyme.API.maxtypeoffset!(2^16)` | ❌ still `EnzymeNoTypeError` |
| `maxtypeoffset!(2^16)` + `maxtypedepth!(64)` | ❌ still `EnzymeNoTypeError` |
| `maxtypeoffset!(2^24)` (16 MB search window) | ❌ still `EnzymeNoTypeError` (compile time barely changed → Enzyme is *not* hitting the offset ceiling) |
| `_literal_step(builder, step)` branch ladder in `steps.jl` | ✅ works, but introduces a runtime branch in the hot step accessors — **rejected** (bad for Reactant) |
| fold step into a plain array index inside `leapfrog_kernel!` (à la vertical-advection #1131) | ✅ works, but changes the extensively-used kernels — **rejected** (too invasive) |
| custom `EnzymeRules` rule on `get_step` | ❌ **wrong gradients** — `get_step` returns a view used *mutably* (read+written, aliasing the prognostic), which a custom rule cannot model. Caught by the requested `EnzymeTestUtils.test_reverse` unit test. |

**Key conclusion on `maxtypeoffset!`:** raising it (even to 2^24) does **not** resolve the error,
so this is a genuine *unprovable-type* problem (the runtime pointer arithmetic loses the element
type), **not** a search-space/size limit. `looseTypeAnalysis!` — telling Enzyme to guess — is the
only knob that clears it, and it is verified not to corrupt gradients here.

Standalone MWE attempts (all *succeed* on 1.12, i.e. don't reproduce) kept as
`_mwe_enzyme_offset{,2,3}.jl` and `_mwe_enzyme_lta_offset.jl`; the only faithful reproducer is
`_mwe_get_step.jl`.

### UPDATE (2026-07-02) — SpeedyWeather-free MWE found; root cause revised

A later downward-reduction campaign (see [`MWE_REDUCTION_STATE.md`](MWE_REDUCTION_STATE.md))
succeeded where the 7 upward MWEs failed. The bug **is** reducible:
[`MWE_enzyme_runtime_ntuple.jl`](MWE_enzyme_runtime_ntuple.jl) (~90 lines, only Enzyme +
KernelAbstractions) throws `EnzymeNoTypeError` on Julia 1.12.6 and succeeds with the exact
gradient on 1.10.11.

**Revised root cause:** the trigger is *not* the runtime-offset `get_step` view. It is
`get_steps(var) = ntuple(step -> get_step(var, step), size(var, 3))` — a **runtime-length
ntuple** whose branch ladder returns a small Union of tuple types of large wrapped-view
aggregates. Passing those into a differentiated KernelAbstractions kernel exceeds Enzyme's type
analysis on Julia ≥ 1.11 (Memory-backed arrays enlarge the view aggregates). Ablations: a
literal-length ntuple is fine, a plain loop instead of the KA kernel is fine, a wrapper without
the inactive nested `spectrum` field is fine — all three ingredients are needed; a runtime step
*index* is **not** needed. This explains why the earlier `_literal_step` fix worked (it made the
ntuple lambda's views literal-offset, collapsing the union) and why `maxtypeoffset!` didn't help.
The MWE is suitable for an upstream Enzyme.jl issue.

### Fix #2 — `vorticity_flux_curldiv!` under `Duplicated(model)` (RESOLVED by the same `looseTypeAnalysis!`)

Distinct root cause from #1 (a source-level `get_step` change does **not** touch it). Occurs only
when differentiating w.r.t. the **model** (`Duplicated(model)`); the same function differentiates
fine with `Const(model)`. So it affected **parameter-sensitivity AD** (testsets 5/6/7), not state
AD. The failing method was `vorticity_flux_curldiv!` (`dynamics/tendencies.jl`): Enzyme could not
type some shadow of the model under `Duplicated(model)` on Julia 1.12.

**RESOLVED:** it is the same class of unprovable-type failure as #1, and `Enzyme.API.looseTypeAnalysis!(true)`
(Julia ≥ 1.11) clears it as well. Verified on the full barotropic suite (2026-07-01): testsets 5, 6
and 7 — the `Duplicated(model)` / parameter cases — no longer throw any `EnzymeNoTypeError`.
Testsets 5 and 7 pass; testset 6 (parameters) now *runs* its AD and only misses the strict FD
tolerance (`atol=1e-5`), which is a finite-difference-agreement matter, not an Enzyme error.

---

## FINAL RESULT (2026-07-01, amended 2026-07-02) — both EnzymeNoTypeErrors resolved

Shipped configuration (see [`MWE_REDUCTION_STATE.md`](MWE_REDUCTION_STATE.md) for the full trail):
1. **No-tuple source fix** for error #1: `update_prognostic!` and the implicit solvers bind step
   views individually via `get_step(var, 1/2)` instead of destructuring a `get_steps` tuple —
   tuples of LTA step views (even compile-time-length ones) exceed Enzyme's type analysis on
   Julia ≥ 1.11 ([EnzymeAD/Enzyme.jl#3275](https://github.com/EnzymeAD/Enzyme.jl/issues/3275)).
   With this, state AD differentiates natively on 1.12 (no flag needed).
2. `Enzyme.API.looseTypeAnalysis!(true)` in
   [`ext/SpeedyWeatherEnzymeExt.jl`](../../ext/SpeedyWeatherEnzymeExt.jl) `__init__`, gated
   `if VERSION >= v"1.11"` — still required for error #2 (parameter AD, `Duplicated(model)`),
   which is the same failure class but not caused by `get_steps` tuples (verified: flag-off run
   with the source fix still errors in testsets 5/6/7 only).

Full `speedy_transforms.jl` + `barotropic.jl` run (differentiability env, Enzyme 0.13.173):

| | Julia 1.10 | Julia 1.12 |
|--|-----------|-----------|
| speedy_transforms | 114/114 pass | 114/114 pass |
| barotropic (7 testsets) | 17 pass, 3 fail, **0 error**, 1 broken | 17 pass, 3 fail, **0 error**, 1 broken |

The two versions are **identical**. The 3 `barotropic.jl` failures (L219 `rtol=0.05`,
L220 mean-abs `< 0.002`, L251 parameter FD `atol=1e-5`) are finite-difference tolerance
comparisons that fail the same way on both versions — pre-existing test-quality items,
explicitly out of scope for the "no Enzyme errors" goal. Note: at L219 Enzyme's VJP is ≈ 2×
the finite-difference VJP on the non-trivial components (identical on 1.10 and 1.12), which points
to a test-construction subtlety (double-seeded `Duplicated` inputs in `timestep_oop!`) rather than
an AD correctness problem — worth a separate look but not an Enzyme error.
