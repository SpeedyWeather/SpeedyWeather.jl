# Recompute the semi-implicit reference temperature profile periodically

> Status: **in progress**. Implemented with a `ReferenceTemperatureSchedule`, uncommitted on `mg/instability`.
> T256L8 on the GPU is stable in all 10 runs with it: 3 × 2 years each with `Day(1)`, `Day(7)` and `Day(14)`,
> and 10 years with `Day(1)`. The default blew up in 8 of 8. GPU and differentiability tests are still to run.

Date of initial draft: 2026-10-19

Base revision: `e410eed0` (branch `mg/instability`, identical to `main`)

## Originating prompt

> continue. While we wait for the results, draft a possible fix that recomputes the reference temperature

## Revision log

- 2026-10-19: initial draft, with the code drafted alongside (`reference_update_period` option, `DateTime`
  of the last update compared on every step).

> We need to reimplement the reinitializing and recomputation of the reference_temperature. Implement it with a
> `schedule` that is set in a `Date` period when constructing the model. This schedule should precompute a
> `recompute_reference_temperature_every_n_steps`. We should already leave a possibility for Reactant to just
> do it at every time step in case that's more efficient than an explicit `@trace` conditions. But don't test
> this here yet, just make the implementation flexible enough for that. Once the implementation is done, test
> if the T256L8 GPU model is stable, start with 2y but also test 10y

- 2026-10-20: reimplemented with a `ReferenceTemperatureSchedule` (period `every`, precomputed
  number of steps). The decision of when to recompute is a separate function dispatching on the
  architecture, the hook for Reactant to recompute on every time step.

> why is the initialize!(schedule, time_stepping) called in EVERY reinitialize!, it should be only called in the
> beginning once, in the initialize! of the implicit correction

- 2026-10-20: the schedule is initialized once, in a new `initialize!(implicit, model)` called by
  `initialize!(model)` right after the time stepper, instead of in `reinitialize!` at every time-step change.

> before you start any validation, let me review the code and give a okay

> okay, go ahead with testing for stability. Once you've shown that it yields stable runs, also experiment with
> inceasing the schedule to Day(7) and Day(14) and see if that's stable as well.

- 2026-10-20: validation runs, see *Testing and verification*.

> change recompute_reference_temperature_every_n_steps to just `every_n_steps`

- 2026-10-20: field renamed to `ReferenceTemperatureSchedule.every_n_steps`.

> the unit test is too extensive. Dont actually run the simulations, this takes to long. Just make a test that it
> reinitialize! function does what it should based on different inputs of variables

- 2026-10-20: the unit test now calls `reinitialize!` directly with set step counters and profiles, without
  time stepping (~1 min instead of ~3 min).

> I don't understand the test, just do a test where you manaually set a new temperature, then reinitialize the
> implicit and check whethers it's really the averege and then do that twice

- 2026-10-20: simplified to that single scenario (see item 6 of *Summary of changes*).

## Problem description

The default `PrimitiveWetModel` at T256L8 blows up after 6 to 12 months. The investigation is in
`docs/dev/2026-10/2026-10-19-t256-instability.md` (untracked, on `mg/instability`). Its lead finding is that
the semi-implicit solver's reference temperature profile Tₖ is frozen at the initial conditions. The model
climate drifts far from it within weeks (k1: +32 K, k3–k8: −8 to −15 K). Top-layer cold spots then develop
over the Antarctic plateau and the run eventually produces NaNs.

## Background

`ImplicitPrimitiveEquation` linearizes the primitive equations about a resting atmosphere with the vertical
profile Tₖ (Hoskins and Simmons, 1975). `reinitialize!` builds all its operators (R, U, L, W and the inverted
S⁻¹ per degree) from `vars.dynamics.average_temperature_profile`. `transform!` updates that profile every time
step as the l = m = 0 coefficient of temperature, i.e. the global mean per layer. But `reinitialize!` only
rebuilds when Δt changes, which happens in the leapfrog start (Δt/2, Δt, 2Δt) and never afterwards. So
within one `run!`, Tₖ stays the initial-condition profile (200, 226, 243, 255, 265, 273, 279, 285 K
for the default initial conditions).

Every use of Tₖ in the dynamical core reads `implicit.temp_profile`: the anomaly T′ = T − Tₖ in
`dynamics_tendencies!`, the linear term R_d Tₖ ln pₛ in `linear_pressure_gradient!`, the pressure gradient in
`vordiv_grid_tendencies!`, the adiabatic term in `_temperature_tendency_kernel!`, and the implicit operators.
So the linear/nonlinear split is consistent for any Tₖ, and replacing `implicit.temp_profile` (with the
operators) between two time steps changes only how the terms are split, not what they sum to. The one
exception is `linear_virtual_temperature!`, which uses the current `average_temperature_profile` for the
small Tₖμq humidity correction (≈ 2 K). It is unaffected.

Evidence from the investigation (T256L8, driver-level refresh of the reference, one realization per row
unless noted):

| reference profile | outcome | lowest T_min | max. points < 170 K | max \|u\| |
|---|---|---|---|---|
| frozen at t = 0 (default) | NaN in 8/8 runs, day 181–358 | 126–166 K | 44–1949 | 105–126 m/s |
| refreshed daily | stable 400 d | 165.5 K | 18 | 108 m/s |
| refreshed every 30 days | stable 365 d, and 730 d in 2/2 seeds | 161–167 K | 15–78 | 94–118 m/s |
| isothermal 300 K | stable 365 d, and 730 d in 2/2 seeds | 180–186 K | 0 | 125–129 m/s |
| frozen + 30 K (same shape) | NaN day 286 | 148 K | 1138 | 111 m/s |

Refreshing the profile cures the instability and keeps the climate close to the default (jets). The
isothermal 300 K reference (the classical choice of Simmons, Hoskins and Burridge, 1978) also cures it, but
strengthens the subtropical jets by ~15–25%. The mechanism is not fully understood: a linear analysis of the
scheme finds the stale reference neither unstable nor differently damped. The cure is nevertheless robust
across realizations, and it is also what the model does at every new `run!`.

## Summary of changes

1. New `ReferenceTemperatureSchedule <: AbstractSchedule` (exported, mutable, in
   `SpeedyWeather/src/time_stepping/implicit/implicit_primitive_equations.jl`):
   - `every::Second = Day(1)`, the period set when the model is constructed (`Second(0)` = never, the old
     behaviour), also as `ReferenceTemperatureSchedule(Day(1))`;
   - `every_n_steps::Int`, precomputed from `every` and the time step by
     `initialize!(schedule, time_stepping)` (rounded, at least 1, 0 = never);
   - `counter::Int`, the number of recomputations, for diagnostics and tests.
2. `recompute_reference_temperature(schedule, clock, architecture)` decides whether to recompute on this time
   step: `step_counter > 0 && step_counter % n == 0`. It dispatches on the architecture, so that an
   architecture which compiles the time step (Reactant) can add a method that recomputes on every time step
   instead of tracing the condition. No such method is added yet.
3. `ImplicitPrimitiveEquation` gets the option
   `reference_temperature_schedule::ReferenceTemperatureSchedule = ReferenceTemperatureSchedule()` and a new
   `initialize!(implicit, model)`, which initializes the schedule once with the model's time step. It is called
   in `initialize!(model)` of `PrimitiveWet` and `PrimitiveDry` right after the time stepper; like horizontal
   diffusion, a time step changed with `set!` needs a new `initialize!(model)`. The operators are still built
   in the first time steps by `reinitialize!`.
4. `reinitialize!(implicit, model, vars)` rebuilds the operators when Δt has changed (leapfrog start, as
   before) **or** when the schedule says so, from the current `average_temperature_profile`.
5. `AbstractSchedule` moved from `output/schedule.jl` to `time_stepping/abstract_types.jl`, as the implicit
   solver is included before the output.
6. Unit test `"Implicit reference temperature profile is the global-mean temperature"` in
   `SpeedyWeather/test/dynamics/time_stepping.jl` (T22L4 PrimitiveWetModel, no time stepping, < 1 min).
   It does the following twice, with ΔT = 0 and ΔT = 20 K:
   - set the temperature with `set!` to 200 + 100σ + ΔT plus a wave 10 cos λ cos φ, which has zero global
     mean;
   - update the global mean with `transform!`;
   - call `reinitialize!` on a step where the schedule recomputes (n, then 2n);
   - check that the reference profile is 200 + 100σ + ΔT, the global mean of the temperature, to 0.01 K.
7. Docs: `docs/src/primitiveequation.md` (semi-implicit pressure gradient) states how Tₖ is chosen and updated.
8. `CHANGELOG.md` entry under Unreleased (PR number to be filled in). No version bump: `SpeedyWeather` is
   already at `0.23.0+DEV`, and no other package is touched.

Choice of the default period. Both 1 and 30 days cured the instability. Daily keeps the reference closest to
the climate: at most 18 points below 170 K, against 15–78 with 30 days. The cost is negligible either way.
One rebuild is 256 LU factorizations of 8×8 matrices on the CPU plus a few small host–device copies, once
every 288 steps at T256. Between rebuilds the per-step overhead is an integer modulo on the host,
with no device synchronization.

Alternatives considered:

- **Fixed isothermal reference** (e.g. 300 K): stable and cold-spot free, but it changes the jets by ~20%,
  i.e. the model climate, through the time-truncation error of the implicit gravity-wave treatment.
- **A climatological reference profile** (e.g. a fixed "model climate" profile) instead of the initial
  conditions. It is model-configuration dependent, and it still goes stale when users change physics,
  orbit, greenhouse gases, etc.
- **Exponential relaxation** of Tₖ towards the current mean. Smoother, but needs a rebuild at every change
  anyway, so a periodic rebuild is simpler for the same effect.

## Testing and verification

`ReferenceTemperatureSchedule` implementation (2026-10-20):

- [x] `test/dynamics/time_stepping.jl` (`--check-bounds=yes`): all pass, the new test
  "Implicit reference temperature profile is the global-mean temperature" 2/2 in 53 s. Runic: clean. `test/dynamics/dispatch.jl`
  (JET): `PrimitiveWet time_step! free of runtime dispatch` passes, so the per-step check adds no dynamic
  dispatch.
- [x] T256L8 default model on the GPU (H100), 730 days, unperturbed and 2 seeds (`PERTURB_SEED` 1, 2), with
  `Day(1)` (the default), `Day(7)` and `Day(14)` (driver option `REF_SCHEDULE_DAYS`). **All 9 runs are stable
  for 730 days**; the default blew up in 8 of 8. The logs confirm the precomputed periods (288, 2016 and
  4032 steps) and 730, 104 and 52 recomputations. Statistics over days 150–730 (min/max over the 3 runs):

  | schedule | lowest T_min | max. points < 170 K | days with points < 170 K | max \|u\| | year-2 T_top, T_sfc |
  |---|---|---|---|---|---|
  | `Day(1)` | 155.9–163.6 K | 18–129 | 9–15 | 108–123 m/s | 231.9–232.2, 276.8–276.9 K |
  | `Day(7)` | 165.0–168.0 K | 11–25 | 3–9 | 104–107 m/s | 232.0–232.1, 276.8 K |
  | `Day(14)` | 155.1–168.3 K | 1–202 | 2–17 | 103–139 m/s | 232.0–232.1, 276.9–277.0 K |
  | default (frozen) | 125.8–166.2 K | 44–1949 | 3–63 | 105–109 m/s | NaN within day 181–280 |

  With 3 runs each, `Day(1)` and `Day(7)` are indistinguishable. `Day(14)` leaves the largest residual
  excursions (one run with 202 points < 170 K). Single-day wind spikes of 121–154 m/s occur rarely: 1–2 days
  in a 2-year run, in layer k2 at 52–73°S in southern spring, with a cold point at the same time. They
  decay within a day and occur with `Day(1)` and `Day(14)`, but in none of the first-draft or `update30`
  runs, so they are sampling variability, not a trend.
- [x] T256L8 with `Day(1)`, 10 years (Slurm 2614622, `gpumedium`, 187 min, 77 SYPD): **stable for 3650 days**,
  3650 recomputations. Over days 150–3650: lowest T_min 162.5 K (day 997, k3), at most 85 points < 170 K,
  on 61 of 3501 days, and 18 single days with max |u| > 120 m/s (highest 154 m/s on day 1370, back to
  126 m/s the next day), with no trend over the 10 years. Mean of years 2–10: T_top 231.9 K, T_sfc 276.7 K.

First draft (`reference_update_period`, 2026-10-19):

- [x] The new unit test, and the rest of `test/dynamics/time_stepping.jl`, including the #1297 test
  "Single time steps match run! from a reset clock" (`--check-bounds=yes`): all pass (new test 16/16,
  #1297 test 10/10). Runic: clean.
- [x] T256L8 default model with the draft, 730 days, unperturbed and 2 seeds (`PERTURB_SEED` 1, 2), GPU, with
  the unchanged long-run driver, so the fix comes from the library default only (Slurm 2588109–2588111).
  **All three are stable for 730 days**, through two southern winters. Over days 150–730: lowest T_min
  166.6–167.3 K (k3/k4, not k1), at most 3–27 points < 170 K on 3–5 days, max |u| 104–109 m/s. Year-2
  global means: T_top 232.0–232.1 K, T_sfc 276.9–277.0 K, the same as the `update30` ensemble (232.1 K,
  276.8–277.0 K).
  First southern winter (days 150–288), compared with the default runs from the same initial conditions:

  | run | with the draft | default (frozen reference) |
  |---|---|---|
  | unperturbed | past day 288; lowest T_min 169.7 K (k3); at most 2 points < 170 K | NaN day 280; 125.8 K (k1); 1949 points |
  | seed 1 | past day 287; 169.5 K (k3); 2 points | NaN day 257; 141.9 K (k1); 1305 points |
  | seed 2 | past day 306; 171.9 K (k1); 0 points | NaN day 181; 166.2 K (k1); 44 points |

  Global means are unchanged (T_top 231.6–231.8 vs 231.3–231.5 K, T_sfc 277.4–277.6 vs 277.8–277.9 K), and
  max |u| is 91–101 m/s against 105–109 m/s.
- [x] T128L8 default with the draft, 730 days (Slurm 2588112): stable. The cold spots of the T128 default
  mostly disappear: lowest T_min 162.3 K (k3) vs 142.4 K (k1), and at most 6 points < 170 K on 8 days vs
  up to 375 points on 116 days (400-day default run). Global means unchanged.
- [x] Full `SpeedyWeather` test suite (`runtests.jl`, as `Pkg.test`, `--check-bounds=yes`, Slurm 2590250):
  118,252 passed, **0 failed**, 33 errors, all from the compute node's environment. 31 come from the asset
  `albedo.nc`, which could not be downloaded without internet; 2 from the ecRad data artifact of the
  optional NumericalRadiation extension, which is not installed. With the assets available (login node),
  `parameterizations/albedo.jl` and `parameterizations/ocean_sea_ice.jl` pass 320/320. The default primitive
  dry and wet long integrations pass, so no test needed updated reference values.
- [ ] GPU tests (`test/GPU`) and the extended differentiability tests (Enzyme). `reinitialize!` is
  `EnzymeRules.inactive`, so a rebuild inside a differentiated loop runs on the primal only, as for the
  leapfrog start.
- [ ] Benchmark: no measurable change expected (`SpeedyWeather/benchmark`).

## Documentation changes

- `ImplicitPrimitiveEquation` docstring (new option, via `TYPEDFIELDS`).
- `docs/src/primitiveequation.md`: one paragraph on how the reference profile Tₖ is chosen and updated.

## Known limitations

- **Reactant**: `time_stepping!` for a `ReactantSimulation` compiles `time_step!` once and calls it in a loop.
  The host-side decision in `reinitialize!` is evaluated at trace time only, so compiled Reactant runs keep
  the reference from tracing, exactly like the Δt-based rebuild today. The hook for Reactant is a method of
  `recompute_reference_temperature(schedule, clock, ::ReactantDevice)` that returns `true` (recompute on every
  time step, no traced condition). The recomputation itself (`initialize!(implicit, Δt, Tₖ, …)`, an LU per
  degree on the CPU) would then also need a traceable implementation. Neither is done or tested here.
- **Restarts**: a restart rebuilds the reference from the restart state (as before) and then updates it
  periodically. Bit-identical restarts were not possible before either.
- **Reproducibility**: results of primitive equation runs longer than the schedule's period change
  compared with v0.23.0. `ReferenceTemperatureSchedule(Second(0))` restores the old behaviour.
- **Mechanism**: why the stale reference produces the Antarctic cold spots is still open (see the
  investigation). The fix removes the trigger without a full explanation.

## Future work

- Reactant support for the periodic rebuild (rebuild between compiled steps).
- A temperature tendency budget at the cold spots with the stale vs the refreshed reference, to explain the
  mechanism.
