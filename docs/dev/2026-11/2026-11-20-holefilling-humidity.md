# Hole filling of the PrimitiveWet humidity: integrate the clipped humidity

> Status: **in progress**. Option 1 (the fill as a tendency, from the previous leapfrog step) is implemented and
> passes the unit, JET, Enzyme, Reactant (CPU) and 10-year CPU checks; the 100-year run and the project's AD check
> are pending.

Date of initial draft: 2026-11-20

Base revision: fa90cc8129a9423286731721b0452264aa1bb8c0

## Originating prompt

> Work on this plan to adjust the hole filling

The plan started as a handoff from the project side (FokkerPlanckSensitivities.jl): problem, cause, fix options,
constraints and acceptance tests below, committed in fa90cc81.

## Revision log

* 2026-11-20: Handoff written (sections Problem description and Background).
* 2026-11-20: Option 1 taken and refined: the fill is stored per leapfrog step, and the fill of the **previous** step
  (the step the leapfrog steps forward from, as the parameterizations) is added, not that of the current step.
  Added to the grid humidity tendency after the parameterizations (so that SPPT does not perturb it) and
  before the dynamics. Implemented on `mg/holefilling-revised` (on main), not on `mg/gfdt-upstream` as the handoff
  asked.
* 2026-11-20: The unit test found that the Robert–Asselin–Williams filter sees the unfilled previous step;
  accepted as a known limitation (≈ 2.5 % of the fill, refilled two steps later).
* 2026-11-20: Asked whether the problem is specific to `StochasticLeapfrog` (and the fix should then be limited
  to it): audited deterministic and SPPT runs, it is general (see Audit below), so the fix stays general.
* 2026-11-20: Enzyme reverse mode failed on the tendency written as a broadcast over `Field`s
  (`IllegalTypeAnalysisException`, not covered by the existing differentiability tests); both new broadcasts
  now act on `.data`, and a differentiability test was added.

## Problem description

**Goal.** A hole filling in SpeedyWeather's `PrimitiveWetModel` that keeps the spectral prognostic humidity
consistent with the grid humidity the model actually uses, so that the upper-level spectral humidity stops
drifting. This blocks the project's next control run (`docs/dev/plan-v025-recalibration-2026-11-20.md`,
stage 0.2): the exported state is the spectral coefficients, and a coordinate that drifts linearly over 1500 years is
a clock for the score network.

1500-yr control run (T32 L8 `PrimitiveWetModel`, `StochasticLeapfrog` + `AdditiveWienerNoise` with the default
amplitudes, default `SpeedLimitDrag`, Δt 40 min, SpeedyWeather `c19adc3d9`; `/p/tmp/maxgelbr/scratch_v025/control_t32l8_white`):

* the layer-2 spectral humidity's global mean (l = 0 coefficient / √4π) falls **linearly from ≈ 0 to −0.53 kg/kg**
  (R² 1.000); layer 1 goes slightly negative (run mean −1e-4 kg/kg); layers 3–8 and every temperature are flat;
* the grid humidity the physics uses is **never negative**: layer 2 between 0 and 2.2e-3 kg/kg, no point below zero;
* the old model (SpeedyWeather-main `496f1c6c`, 0.22.1+DEV) has the same `ClipNegatives` and the same noise
  amplitudes but no drift: layer 2 +3.4e-4 kg/kg and stationary. The new model's layers 1 and 2 are 18 K and 11 K
  colder (222.5 / 221.7 K against 240.5 / 233.1 K), so they hold almost no moisture, and noise and transport push that
  near-zero field below zero far more often.

## Background

### Cause

* `PrimitiveWetModel` defaults to `hole_filling = ClipNegatives(spectral_grid)` (`src/models/primitive_wet.jl:107`).
* `transform!` (`src/time_stepping/transform.jl:158–161`) transforms the prognostic step to grid space and then calls
  `hole_filling!(humid_grid, model.hole_filling, model)`, which is `q .= max.(q, 0)` (`src/dynamics/hole_filling.jl`).
  It changes **only the grid field**. The spectral prognostic humidity keeps its negative values, and nothing
  transforms the clipped field back.
* Everything downstream reads the clipped grid humidity, and their tendencies are added to the spectral humidity:
  * the humidity advection — `q·div` and the `(uq, vq)` fluxes in `humidity_grid_tendency!`
    (`src/dynamics/tendencies.jl:773`), and vertical advection;
  * the virtual temperature in the geopotential (`src/dynamics/geopotential.jl:78`);
  * large-scale condensation, convection, vertical diffusion, the surface fluxes, the clouds and both radiation
    schemes.
* So the field the model integrates (spectral) and the field its tendencies see (clipped grid) differ by the holes,
  `max(q, 0) − q ≥ 0`. Moisture carried into a hole region is absorbed by the negative spectral value and stays
  invisible, because the clipped field there is still zero. Sinks act on the moister clipped field. The holes deepen
  and spread, and the drift rate grows with their depth.
* Once established, the drift runs without any noise (table below). The noise seeds it: from rest, the same model
  without noise stays positive for at least 10 years.

### Evidence

10-yr runs (`speedy/checks/humidity_drift_test.jl`, jobs 3349322 and 3349373; daily TSVs in
`/p/tmp/maxgelbr/scratch_v025/check/humidity_drift/`). Layer-2 spectral global mean, year-1 mean → year-10 mean:

| start | noise | layer 2 [kg/kg] |
|---|---|---|
| end of the control run | default | −0.5345 → −0.5378 |
| end of the control run | humidity noise off, rest default | −0.5345 → −0.5377 |
| end of the control run | none (amplitudes × 1e-9; a zero-amplitude model has one array fewer and cannot load the saved state) | −0.5344 → −0.5370 |
| rest (initial conditions) | default | −2.6e-5 → −1.6e-3 |
| rest | none | +2.5e-5 → +7.6e-6 |

MWE (job 3349538: 3 years from rest, GPU, 6 min). The unclipped field is `transform(spectral humidity)`; the clipped
field is what the model uses:

| case | year | layer 1: spectral mean, unclipped < 0 | layer 2: spectral mean, unclipped min, unclipped < 0 | clipped min |
|---|---|---|---|---|
| noise | 1 | −2.3e-5, 64 % | −8.4e-5, −1.0e-2, 29 % | 0 |
| noise | 3 | −6.9e-5, 73 % | −4.2e-4, −2.1e-2, 48 % | 0 |
| no noise | 1 | +1.1e-5, 27 % | +5.8e-5, −6.2e-3, 8 % | 0 |
| no noise | 3 | +3.8e-5, 8 % | +8.6e-5, −9.8e-4, 9 % | 0 |

So holes exist without noise too (8–9 % of layer-2 points); with noise they cover half the layer and the mean runs off.

### Reproduce

```bash
# login node: environments speedy-upstream(-gpu) are instantiated; compute nodes have no internet
sbatch --account=flai --partition=gpu --qos=gpushort --gres=gpu:1 --cpus-per-task=8 --mem=32G --time=01:00:00 \
  --wrap "cd /p/tmp/maxgelbr/code/FokkerPlanckSensitivities.jl && module load julia/1.11.7 && \
          julia --project=speedy-upstream-gpu docs/upstream/hole_filling_drift_mwe.jl gpu 3"
```

On CPU, drop `gpu` and use `--project=speedy-upstream`; it is much slower. Longer and per-day diagnostics come from the project
script, with `HQ_VARIANTS` from `default, humoff, det` (end state) and `default_ic, det_ic` (from rest):

```bash
P4C_ARCH=gpu HQ_YEARS=10 HQ_VARIANTS=default_ic,det_ic julia --project=speedy-upstream-gpu speedy/checks/humidity_drift_test.jl
```

A useful extra diagnostic while developing a fix is the water added by the filling per year, and the global mean of
the clipped grid field next to the spectral one.

### Fix options

1. **The fill as a tendency (suggested first).** Keep the clip for the grid field. Write the fill
   `(max(q, 0) − q)` into a scratch field, and add it divided by the step's effective time step (leapfrog: mind the
   2Δt) to the grid humidity tendency before that tendency's grid→spectral transform. That transform happens anyway,
   so there is no extra one. The spectral state then receives the moisture the tendencies assumed. It does not
   conserve water: it adds the hole volume, which stays small once the holes stop deepening.
2. **Transform the clipped field back.** After the clip, write `transform(max(q, 0))` into the current spectral level.
   Grid and spectral then agree at once. It costs one extra 8-layer grid→spectral transform per step. Truncation
   reintroduces small negatives (Gibbs), so the clip–truncate cycle is a small systematic moisture source; that
   source is bounded. It also edits a leapfrog level directly, so check the interplay with the
   Robert–Asselin–Williams filter.
3. **A conservative filler.** Fill holes from the same column's moister layers, mass-weighted by Δσ, then a global
   rescaling, in the spirit of Rasch and Williamson (1990). Applied through route 1. Conserves water, but is more
   code.

Not a fix, but complementary on the project side: the additive humidity noise (3e-8 at every layer) seeds the holes
in the dry upper layers. Scaling it down there would slow the seeding; it does not remove the inconsistency.

**Constraints for any fix:**
* GPU via KernelAbstractions; CPU = GPU to rounding.
* Differentiable:
  * plain Enzyme forward and reverse on Julia 1.11 with a `Duplicated` model (plain `Reverse`, see
    `docs/dev/handoff-enzyme-reverse-pwet-2026-11-18.md`);
  * Reactant (`raise_first`, no scalar indexing, no data-dependent control flow);
  * `max` is fine; avoid runtime-`Symbol` dispatch inside kernels.
* No change for `PrimitiveDry`, the barotropic and the shallow-water models.
* Should be fast and not slowdown neither CPU nor GPU model significantly

### Acceptance tests

1. The MWE with noise: the layer-1 and layer-2 spectral means stay ≥ 0 and show no trend over 3 years. The unclipped
   minimum does not deepen from year to year, and the hole fraction falls well below today's no-noise 8–9 %.
2. 100 years from rest with noise (`HQ_YEARS=100 HQ_VARIANTS=default_ic`, ~0.5 GPU-h): no trend in any layer's
   spectral humidity mean (R² of a linear fit < 0.01; the old model had 0.006).
3. Climate impact, reported: 10-yr global means per layer of temperature and of the grid humidity, before and after
   the fix. The upper layers will moisten. The default longwave (`ByrneOGormanLongwaveTransmissivity`) reads
   `max(q, 0)`, so expect a small water-vapour effect; also report the filling's water source against evaporation.
4. The project's AD check `speedy/checks/ad_backends_primitivewet.jl` passes. As of job 3307167: plain-Enzyme VJP vs
   forward 7.5e-6, Reactant vs plain Enzyme 2.2e-4. Upstream's differentiability and unit tests pass.

**Report back:** the commit hash on `mg/gfdt-upstream`, the option taken, and the outputs of tests 1–4.

## Summary of changes

Option 1, with the fill of the step the time stepper steps forward from.

* **`ClipNegatives` stores its fill** (`src/dynamics/hole_filling.jl`). It defines the variable
  `vars.dynamics.humidity_hole_fill` (grid, all layers, as many steps as the grid humidity: 2 for `Leapfrog`).
  `hole_filling!(vars, ::ClipNegatives, model)`, called in `transform!` instead of the grid-only clip, writes
  `f = max(−q, 0)` of the current step and clips `q = max(q, 0)` (two broadcasts over `.data`, no branches; a
  broadcast over the `Field`s fails Enzyme's type analysis in reverse mode).
* **The fill as tendency.** `hole_filling_tendency!(vars, model)` is called in `time_step!` of the primitive
  equation models after the parameterizations (stochastic physics should not perturb the fill), ocean, sea ice and
  land, and before the dynamics, which scale all tendencies for the dynamical core. It adds `f / Δt` to the grid
  humidity tendency, with `Δt = time_step(time_stepping, clock)` the time step of this step (leapfrog: Δt/2, Δt in
  the leapfrog start, then 2Δt). It runs whatever `dynamics_only` and `dynamics` are, because the clip does too.
  No-op for `PrimitiveDry` (dispatch) and for `hole_filling = nothing`.
* **Which step.** With `Leapfrog`, `q⁽ⁿ⁺¹⁾ = q⁽ⁿ⁻¹⁾ + 2Δt T`, and the parameterizations read the previous step too.
  The fill of the previous step is added (`which_prognostic_step(var, ::AbstractLeapfrog, ::AbstractHoleFilling) = 1`),
  so that in grid space `q⁽ⁿ⁺¹⁾ = max(q⁽ⁿ⁻¹⁾, 0) + 2Δt T(q⁽ⁿ⁾₊, q⁽ⁿ⁻¹⁾₊)`: every input of the step is clipped. With
  the implicit horizontal diffusion this is identical to filling the base and then diffusing. The fill of the
  current step would instead add one leapfrog sequence's hole to the other sequence, and in the leapfrog start, where
  the previous step is not advanced, it would leave the second step unfilled.
  `move_prognostic_grid_variables_back!` moves the fill back with the grid variables (same `move` flag), so the
  fill's step 1 always belongs to the grid humidity's step 1, including the leapfrog start and the initial
  transform. A time stepper with a single step adds the fill of that step, the step it steps forward from.
* **API.** `hole_filling!(::Variables, ::Nothing, ::AbstractModel)` is the no-op. A custom `AbstractHoleFilling`
  that only defines `hole_filling!(humid_grid::AbstractField, H, model)` still works, with the old grid-only behaviour
  (fallback for `hole_filling!(::Variables, ::AbstractHoleFilling, model)`, `hole_filling_tendency!` no-op).
  `hole_filling!(field, ClipNegatives())` still clips a field.
* Barotropic, ShallowWater and PrimitiveDry are unchanged: they have no hole filling component, and the fill's
  move-back is skipped by a `haskey` on `vars.dynamics` that is resolved at compile time.

Files: `src/dynamics/hole_filling.jl`, `src/time_stepping/transform.jl`, `src/time_stepping/time_integration.jl`,
`src/time_stepping/steppers/leapfrog.jl`, `test/dynamics/hole_filling.jl` (new),
`test/differentiability/primitivewet.jl`, `docs/src/primitiveequation.md`, `CHANGELOG.md`. SpeedyWeather is
already `0.26.0-DEV`, no version bump.

## Testing and verification

Local, CPU (Apple M3), Julia 1.12.7, Float32, T32 L8 (`truncation = 32`; the MWE uses 33), `StochasticLeapfrog`
with the default noise amplitudes or none, from rest, 10 years.
Script: `docs/dev/2026-11/hole_filling_drift_check.jl`; `nofix` overrides the new tendency with a no-op, which
gives exactly the old trajectory. The old behaviour reproduces the MWE table above closely (noise, year 1: layer 1
−2.3e-5 with 66 % of points < 0; layer 2 −8.9e-5, minimum −7.9e-3, 33 %).

**Unit tests** (`test/dynamics/hole_filling.jl`, pass): the clip and its fill are exact
(`clipped − fill == unclipped`); a constant negative humidity of −1e-5 without physics is refilled by the tendency
(spectral mean within 5 % of the hole after one day, without the fix it stays at −1e-5); after one day of the default
model the grid humidity is exactly `max(transform(spectral humidity), 0)` and the fill exactly `max(−transform(…), 0)`.
The JET dispatch tests (`test/dynamics/dispatch.jl`, PrimitiveWet and StochasticLeapfrog `time_step!` free of
runtime dispatch) pass.

**1. MWE with noise.** Layers 1 and 2 stay positive. The unclipped minimum does not deepen and holes cover
0.1–0.7 % of the points (old: 33–67 %, old without noise 3–31 %).

| version | year | layer 1: spectral mean, unclipped min, < 0 | layer 2: spectral mean, unclipped min, < 0 |
|---|---|---|---|
| fix | 1 | +1.7e-5, −4.5e-6, 0.1 % | +6.4e-5, −2.2e-5, 0.7 % |
| fix | 2 | +1.9e-5, −3.5e-6, 0.2 % | +8.1e-5, −4.6e-6, 0.3 % |
| fix | 3 | +2.1e-5, −2.4e-6, 0.1 % | +9.7e-5, −1.9e-5, 0.4 % |
| fix | 10 | +2.9e-5, −4.0e-6, 0.3 % | +1.1e-4, −2.0e-6, 0.1 % |
| old | 1 | −2.3e-5, −6.3e-4, 66 % | −8.9e-5, −7.9e-3, 33 % |
| old | 3 | −6.8e-5, −2.4e-3, 67 % | −4.2e-4, −2.0e-2, 53 % |
| old | 10 | −2.2e-4, −6.7e-3, 63 % | −1.7e-3, −7.7e-2, 52 % |

The means still rise in the first years. That is spin-up from rest, not the filling: the runs without noise moisten
the upper layers as well, with and without the fix. Linear trend over years 4–10:

| run | layer 1 [kg/kg/yr], R² | layer 2 [kg/kg/yr], R² |
|---|---|---|
| fix, noise | +1.5e-6, 0.45 | +4.8e-6, 0.60 |
| old, noise | −2.1e-5, 1.00 | −1.8e-4, 1.00 |
| fix, no noise | +4.6e-6, 0.78 | +2.0e-6, 0.24 |
| old, no noise | +5.0e-6, 0.94 | +2.9e-6, 0.39 |

**Audit: not specific to `StochasticLeapfrog`.** The clip runs in `transform!` for every PrimitiveWet model,
whatever the time stepper. Holes come from the spectral truncation (Gibbs), the non-monotone spectral
advection, the leapfrog computational mode and physics tendencies over 2Δt; stochastic forcing (the additive
noise of `StochasticLeapfrog`, but also SPPT) adds to them. Runs without `StochasticLeapfrog`, old code vs
fix, from rest (plain `Leapfrog`; SPPT with a `SpectralAR1Process`; the hole volume is the global mean column
integral of `max(−q, 0)` of the unclipped humidity, column water vapour is ≈ 25 mm):

| run | old: points < 0, layers 1/2 | old: deepest | old: hole volume | fix: points < 0 | fix: hole volume | fix: water added |
|---|---|---|---|---|---|---|
| deterministic T32 L8, 10 yr | 31/11 % (yr 1) → 0/3 % (yr 10) | −6.4e-3 | 0.011–0.031 mm | ≤ 0.6 % | 0.0004 mm | 2.6–3.1 mm/yr |
| deterministic T32 L16, 3 yr | 8–16 % in layers 1–4 (yr 1) | −9.9e-4 | 0.014 → 0.021 mm | ≤ 0.8 % | 0.0005 mm | 3.0–4.2 mm/yr |
| SPPT T32 L8, 3 yr | 20/8 % → 5/4 % | −6.4e-3 | 0.012–0.028 mm | ≤ 0.3 % | 0.0004 mm | 2.6 mm/yr |
| deterministic T64 L8, 1 yr | 14/8 % | −1.9e-2 | 0.055 mm | 0.2 % | 0.0002 mm | 2.8 mm/yr |
| StochasticLeapfrog T32 L8, 10 yr | 66/33 % → 63/52 % | −7.7e-2 | 0.07 → 1.25 mm | ≤ 0.7 % | 0.0006 mm | 3.8–4.4 mm/yr |

In the deterministic T64 L8 run the layer-2 spectral mean is already negative after one year (−3.0e-5 kg/kg,
while the clipped grid humidity averages +4.0e-5), the start of the drift without any noise. The additive noise
makes the holes grow without bound, the other runs keep a hole volume of 0.01–0.06 mm and an upper-layer
humidity that differs from the field the model uses (10-yr deterministic T32 L8: layer 1 −40 %, layer 2 −30 %
against the fix). With the fix the hole volume is 25–260× smaller (up to 2000× with the additive noise),
and every run adds 2.6–4.4 mm/yr of water.

**2. 100 years from rest with noise.** Pending: about 2 h on this CPU, ~0.5 GPU-h on the cluster. 10 years from rest
are too short to separate the spin-up from a trend.

**3. Climate impact.** 10-yr means (years 1–10, incl. spin-up), fix against old:

| layer | T fix [K] | T old [K] | ΔT [K] | grid q fix [kg/kg] | grid q old [kg/kg] | ratio |
|---|---|---|---|---|---|---|
| noise 1 | 222.87 | 222.79 | +0.08 | 2.55e-5 | 2.30e-6 | 11.1 |
| noise 2 | 221.96 | 221.57 | +0.40 | 9.40e-5 | 2.61e-5 | 3.6 |
| noise 3 | 231.35 | 230.89 | +0.46 | 3.23e-4 | 2.80e-4 | 1.15 |
| noise 4 | 249.98 | 249.54 | +0.44 | 1.45e-3 | 1.38e-3 | 1.05 |
| noise 5 | 266.00 | 265.61 | +0.39 | 3.40e-3 | 3.29e-3 | 1.03 |
| noise 6 | 275.75 | 275.38 | +0.36 | 5.40e-3 | 5.27e-3 | 1.02 |
| noise 7 | 281.58 | 281.25 | +0.33 | 6.91e-3 | 6.76e-3 | 1.02 |
| noise 8 | 285.20 | 284.86 | +0.34 | 9.15e-3 | 8.96e-3 | 1.02 |
| no noise 1 | 222.29 | 222.38 | −0.08 | 9.71e-5 | 5.73e-5 | 1.70 |
| no noise 2 | 221.54 | 221.35 | +0.19 | 1.32e-4 | 9.45e-5 | 1.40 |
| no noise 3–8 | | | −0.04 … +0.12 | | | 0.99–1.04 |

With noise the old model dried the upper layers through the holes: the fix moistens layer 1 by 11× and layer 2 by
3.6×, and the column is 0.3–0.5 K warmer. Without noise the changes are small.

Water added by the filling: 3.8–4.4 mm/yr with noise, 2.6–3.1 mm/yr without, stationary over the 10 years, i.e.
0.4–0.5 % of the evaporation (720–860 mm/yr). In the old model the water hidden by the clip (the hole depth) grows
from 460 to 8200 mm/yr with noise (0.6× to 9.8× the evaporation), and is 75–210 mm/yr without noise.

**4. AD.** The project's check (`ad_backends_primitivewet.jl`, Enzyme vs Reactant) is pending, project side.
Upstream, in a scratch environment (Enzyme 0.13, Julia 1.12, CPU), all pass:
* the existing adjoint check of `transform!` (forward and reverse vs finite differences, `test/differentiability/primitivewet.jl`), 7/7;
* the new test of `hole_filling_tendency!` in reverse mode, which gives the adjoint 1/Δt on the previous step's fill and 0 on the current one, 2/2;
* reverse mode through `hole_filling!(vars, ::ClipNegatives, model)`: adjoint −1 where q < 0, 0 where q > 0 (not in the test suite).

**Reactant.** The PrimitiveWet CPU-vs-Reactant correctness test of `test/reactant/test_correctness.jl`
(MatrixSpectralTransform, tendencies after one step and 100 steps of prognostic and grid variables) passes, 23/23,
with Reactant on the CPU of this Mac. The Reactant state is also internally consistent: its stored fill equals
`max(−transform(q), 0)` of its own spectral humidity to 4.9e-9 and its grid humidity equals `max(q, 0)` to 2.8e-8
(Float32 rounding). Reactant on GPU not tested; the new code is broadcasts only, no scalar indexing.

**Performance.** CPU, T32 L8, per time step: +36 μs of 5.9 ms (+0.6 %): the fill broadcast (+8 μs), the tendency
(3 μs) and moving the fill back (24 μs). GPU not measured; it adds three small kernels per step.

## Documentation changes

* `docs/src/primitiveequation.md`: new subsection "Hole filling" under "Humidity equation" (the clip, the fill as
  tendency from the previous step, not water conserving, the filter interplay), and the hole filling in steps 1 and
  19 of the algorithm.
* `CHANGELOG.md`: entry under Unreleased (PR number to be filled in).

## Known limitations

* **Robert–Asselin–Williams filter.** The filter sees the unfilled previous step as part of the curvature
  `q⁽ⁿ⁻¹⁾ − 2q⁽ⁿ⁾ + q⁽ⁿ⁺¹⁾`. It moves ≈ ν α/2 ≈ 2.7 % of the fill back into the previous step as a hole and adds
  ≈ 2.3 % to the new step. The moved part is refilled two steps later; the computational mode this excites is
  damped by the filter. In the unit test (constant hole of −1e-5) about 1 % of the hole remains as extra water. An
  exact treatment needs the fill in spectral space inside the leapfrog kernel (see Future work).
* **Not water conserving.** The filling adds ≈ 4 mm/yr (0.5 % of the evaporation) at T32 L8 with noise.
* **Truncation.** The fill reaches the spectral humidity through the spectral truncation of the tendency, so small
  negatives remain on the grid (0.1–0.7 % of the points in layers 1–2, minima of a few 1e-6 to 2e-5 kg/kg).
* **`NCycleLorenz`.** Its weighted tendency memory spreads the fill over the cycle, so the spectral humidity is
  stepped from the clipped state only approximately.
* **Variables layout.** `vars.dynamics.humidity_hole_fill` is a new variable: states saved as a whole `Variables`
  from the old model may not load into the new one. A continued `run!` recomputes the fill in the initial transform.

## Future work

Other approaches, beyond the three fix options above:
* no clipping (`hole_filling = nothing`) and parameterizations robust to q < 0: consistent, but the holes are
  unbounded and the physics sees negative humidity;
* a global multiplicative water fixer after the clip (conserving, can be combined with option 1);
* a positive-definite prognostic variable, e.g. ln q or √q in spectral space: positive by construction, but
  changes the advection equation and loses conservation;
* humidity as a grid-point variable with a monotone or positive-definite advection scheme (semi-Lagrangian or
  flux-limited), as operational spectral models such as ECMWF's IFS do: removes the holes at their source, but
  a redesign of the humidity advection.

Next steps:
* A conservative filler (option 3) through the same `hole_filling_tendency!` route.
* An exact treatment of the filter: transform the fill alongside the tendencies into its own spectral slot and add
  it to the previous step inside the leapfrog kernel (one more 8-layer slot in the batched transform).
* Project side: taper the humidity noise with height, the noise seeds the holes in the dry upper layers.

## After the fix (project side)

* Update the pinned environments (speedy-upstream, speedy-upstream-gpu, speedy-ad) on the login node; Pkg cannot
  reach the internet from compute nodes.
* Run the five-day bitwise check of the dry and barotropic models against `c19adc3d9` (they should not change).
* The new control run (plan stage 0.3) and its drift check (`slurm/drift_check.sh`, which flags any coordinate with a
  linear-trend R² > 0.5 and prints block means of the slow PCs).

## Versions

Handoff: SpeedyWeather `0.25.0-DEV`, branch `mg/gfdt-upstream` at `c19adc3d9` (git tree `55c4d100`; `fb912357a` was
checked bitwise identical for this model); Julia 1.11.7; NVIDIA GPU via CUDA.jl's own artifacts (no system CUDA module).
Fix: SpeedyWeather `0.26.0-DEV`, branch `mg/holefilling-revised` on fa90cc81; checks with Julia 1.12.7 on CPU.

## Files

* `docs/dev/2026-11/hole_filling_mwe.jl` — standalone reproducer (SpeedyWeather only).
* `docs/dev/2026-11/hole_filling_drift_check.jl` — the yearly diagnostics above, with and without the fix.
