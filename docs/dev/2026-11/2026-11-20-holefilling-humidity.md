# Hole filling of the PrimitiveWet humidity: integrate the clipped humidity

> Status: **in progress**. Option 1 (the fill as a tendency, from the previous leapfrog step) is chosen, opt-in via
> `ClipNegatives(spectral_grid; fill_tendency = true)` and on by default only with `StochasticLeapfrog`. It passes the
> unit, JET, Enzyme, Reactant (CPU) and 10-year CPU checks, and the GPU checks on an H100 (GPU = CPU to rounding, all
> drift variants, the 100-year run, timing, Reactant in strict Float32); the project's AD check is pending.

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
* 2026-11-20: All options written up with implementation sketches and a GPU measurement protocol, for a GPU node
  to estimate their cost before an option is chosen ("Options and their GPU cost").
* 2026-11-20: Enzyme reverse mode failed on the tendency written as a broadcast over `Field`s
  (`IllegalTypeAnalysisException`, not covered by the existing differentiability tests); both new broadcasts
  now act on `.data`, and a differentiability test was added.
* 2026-11-20: GPU check on an NVIDIA H100 (see GPU below): GPU versions of the unit checks, CPU vs GPU from the
  same state, the drift variants of the audit (fix against the old clip), the 100-year run (acceptance test 2) and
  the timing. The GPU reproduces the CPU results; the implementation is unchanged. The `test/GPU/CUDA` environment
  pins Reactant 0.2.291 (compat is ≥ 0.2.305), on which `transform!` does not compile on the GPU with the old clip
  either, so Reactant was checked in a freshly resolved environment (0.2.308). There Reactant on the GPU matches the
  CPU only in strict Float32 (with TF32 it fails for the old clip too).
* 2026-11-20: Option 1 chosen, without measuring the other options' GPU cost (option 1 costs ≈ 22 μs per step on an
  H100). Decided: option 1 is **opt-in, not the default**, and is used automatically with `StochasticLeapfrog` only.
  This overrides the earlier "the fix stays general": deterministic and SPPT runs keep the old grid-only clip unless
  they opt in. `ClipNegatives{F}` has a type parameter `fill_tendency` (`true`, `false`, or `nothing` = automatic),
  and `adds_fill_tendency(hole_filling, time_stepping)` resolves the automatic choice by dispatch on the
  `StochasticLeapfrog` type alias at compile time (no runtime branch, Reactant-traceable). The fill variable is only
  allocated with the fill tendency.

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

Option 1, with the fill of the step the time stepper steps forward from, opt-in and automatic with
`StochasticLeapfrog` only.

* **Opt-in.** `ClipNegatives(spectral_grid; fill_tendency = nothing)` returns a `ClipNegatives{fill_tendency}`.
  `adds_fill_tendency(::ClipNegatives{F}, time_stepping)` is `F` for `true`/`false`. For `nothing`, the default and
  also the `PrimitiveWetModel` default, it is `true` only for a `StochasticLeapfrog` (a `Leapfrog` whose `noise` is an
  `AbstractWienerNoise`; method in `stochastic_leapfrog.jl`) and `false` otherwise. Without the fill tendency,
  `ClipNegatives` behaves exactly as before the fix: it clips the grid humidity, defines no variable and adds no
  tendency. The bullets below describe the behaviour with the fill tendency.
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
`src/time_stepping/steppers/leapfrog.jl`, `src/time_stepping/steppers/stochastic_leapfrog.jl`,
`test/dynamics/hole_filling.jl` (new),
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

**2. 100 years from rest with noise** (GPU, T32 L8, `StochasticLeapfrog` with the default amplitudes, seed 1, H100,
32 min). No drift of the filled layers. Linear fits over years 11–100 (the first decade is the spin-up from rest):

| layer | mean [kg/kg] | slope [kg/kg/yr] | R², daily values | R², yearly means | temperature: slope [K/yr], R² (yearly) |
|---|---|---|---|---|---|
| 1 | 4.11e-5 | +4.0e-8 | 0.011 | 0.015 | +0.0003, 0.012 |
| 2 | 1.22e-4 | +3.8e-8 | 0.006 | 0.008 | +0.0003, 0.008 |
| 3 | 3.64e-4 | +5.3e-8 | 0.005 | 0.016 | +0.0007, 0.010 |
| 4 | 1.61e-3 | +1.8e-7 | 0.008 | 0.022 | +0.0010, 0.023 |
| 5 | 3.65e-3 | +2.6e-7 | 0.004 | 0.020 | +0.0008, 0.021 |
| 6 | 5.69e-3 | +3.1e-7 | 0.002 | 0.023 | +0.0007, 0.021 |
| 7 | 7.25e-3 | +3.6e-7 | 0.001 | 0.023 | +0.0004, 0.006 |
| 8 | 9.57e-3 | +4.6e-7 | 0.001 | 0.023 | +0.0007, 0.023 |

With daily values R² < 0.01 in layers 2–8 and 0.011 in layer 1 (the old model: 0.006, metric not recorded), with yearly
means 0.008–0.023. The small residual trend is not the hole filling: temperature trends with the same R² (+0.03 to
+0.1 K per century), and the layers without holes (4–8) trend as much as layers 1–2, a slow warming and moistening of
the whole model. Before the fix the layer-2 trend was −1.5e-4 kg/kg/yr (R² 1.00), now +3.8e-8 kg/kg/yr. Over the 100
years the hole volume stays at 5.1e-4 → 6.5e-4 mm, 0.6/0.8 % (year 1) → 0.1/0.3 % (year 100) of the points of
layers 1/2 are below zero (yearly means), the deepest unclipped value in any layer is −1.2e-3 kg/kg, and the filling adds 3.4–4.3 mm/yr (0.5 % of the evaporation, 780–860
mm/yr). Year 10 equals that of the separate 10-year run exactly (same seed: the GPU runs are reproducible).

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
(Float32 rounding). Reactant on the GPU: see GPU below.

**Performance.** CPU, T32 L8, per time step: +36 μs of 5.9 ms (+0.6 %): the fill broadcast (+8 μs), the tendency
(3 μs) and moving the fill back (24 μs). GPU: about +22 μs per step, see GPU below.

**GPU.** One NVIDIA H100 80GB per SLURM job, Julia 1.12.2, CUDA.jl 6.2 (environment `SpeedyWeather/test/GPU/CUDA`),
Float32, from rest, with the GPU's default transform (`MatrixSpectralTransform` up to T64). Script
`docs/dev/2026-11/hole_filling_gpu_check.jl`, modes `consistency`, `drift`, `perf` and `summary`. The old behaviour
is a hole filling `OldClip <: AbstractHoleFilling` that only defines the grid-field clip. It takes the fallback,
without the fill variable or the tendency, so fix and old run in one process. Logs and daily TSVs are in
`/p/tmp/maxgelbr/code/claude_holefilling_gpu_2026-11-20/`. These runs predate the opt-in. At the time
`ClipNegatives` always added the fill tendency, which is now `ClipNegatives(spectral_grid; fill_tendency = true)`,
and the scripts now request it explicitly.

*Consistency* (`--check-bounds=yes`, all pass). `clip_negatives!` on the GPU equals the CPU bitwise, and
`hole_filling_tendency!` adds `fill(step 1) / Δt`. Without physics, the constant hole of −1e-5 is refilled to
+1.3e-7 after one day: as on the CPU, the filter leaves ≈ 1 % of the hole as extra water. After one day of the default
model, the current step's grid humidity is exactly `max(transform(q), 0)` and its fill exactly
`max(−transform(q), 0)`, with both `MatrixSpectralTransform` and `SpectralTransform`. This holds for the current
step only. On the CPU too, the previous step differs by up to 3.5e-5, because the filter changes its spectral
humidity after its grid humidity was computed (see Known limitations).

*CPU vs GPU from the same state.* The GPU state after 10 days is copied to a CPU model, `SpectralTransform` on both;
maximum relative differences:

| hole filling | 1 step: spectral q, grid q, T | 36 steps: spectral q, grid q, T |
|---|---|---|
| fix | 8.5e-9, 3.6e-7, 5.0e-7 | 1.6e-4, 2.6e-2, 2.1e-3 |
| old | 8.6e-9, 3.6e-7, 5.0e-7 | 1.5e-4, 1.7e-2, 5.3e-4 |

After one step, the fills differ by at most 7e-10 (largest fill 2.5e-4) and are nonzero at the same 791 points. After a
day the trajectories diverge chaotically, the grid humidity through condensation switches, with and without the fix
alike.

*Drift variants*, fix | old (compare the Audit table; T32 L8 unless stated):

| run | points < 0, layers 1/2, yr 1 → last | deepest unclipped | hole volume [mm] | water added \| hidden [mm/yr] | evaporation [mm/yr] |
|---|---|---|---|---|---|
| StochasticLeapfrog, 10 yr | 0.6/0.8 → 0.1/0.3 % \| 50/26 → 67/56 % | −7.1e-4 \| −7.2e-2 | 5.1e-4 → 6.1e-4 \| 0.066 → 1.12 | 3.4 → 4.0 \| 430 → 7400 | 780–850 |
| deterministic, 10 yr | 1.0/0.6 → 0.0/0.2 % \| 23/15 → 0.1/9.7 % | −6.4e-4 \| −5.7e-3 | 4.2e-4 → 4.6e-4 \| 0.027 → 0.032 | 2.8 → 3.0 \| 180 → 210 | 710–770 |
| deterministic T32 L16, 3 yr | 1.0/1.7 → 0.0/0.1 % \| 20/28 → 3.0/4.7 % | −2.1e-3 \| −8.4e-3 | 4.1e-4 → 6.2e-4 \| 0.015 → 0.020 | 2.7 → 4.2 \| 97 → 134 | 670–840 |
| SPPT (SpectralAR1Process), 3 yr | 0.9/0.6 → 0.1/0.1 % \| 25/13 → 8.9/4.8 % | −6.7e-4 \| −6.5e-3 | 3.7e-4 → 4.2e-4 \| 0.027 → 0.012 | 2.5 → 2.7 \| 174 → 77 | 710–740 |
| deterministic T64 L8, 1 yr | 0.6/0.4 % \| 14.5/12.6 % | −8.1e-4 \| −2.1e-2 | 2.0e-4 \| 0.052 | 2.6 \| 685 | 740 |

For the fix, "water added" is the fill accumulated every time step (column integral of the fill / 2 per leapfrog
step). For the old clip, "hidden" is the hole volume × time steps per year / 2, as in the Audit. For the fix both
measures agree (4.02 against 4.0 mm/yr). With noise, year 10, layers 1/2:

* spectral means: fix +3.2e-5 / +1.2e-4, old −2.0e-4 / −1.5e-3 kg/kg (CPU: +2.9e-5 / +1.1e-4 and
  −2.2e-4 / −1.7e-3);
* trend over years 4–10: fix +2.8e-6 (R² 0.65) / +4.4e-6 (0.58), old −2.3e-5 (0.99) / −1.5e-4 (1.00) kg/kg/yr;
* climate impact (10-yr means): +0.10 to +0.50 K, layer 1 9.3× and layer 2 3.7× moister (CPU 11× and 3.6×).

In the old deterministic T64 run the layer-2 spectral mean is again negative after one year: −3.3e-5 against a
clipped mean of +6.7e-5. The GPU reproduces the CPU results within the spread of chaotic trajectories, which here
differ in hardware and transform.

*Timing* (`perf`): default model with `Leapfrog`, median of 5 timings of 1000, 600 and 300 steps, alternating fix
and old:

| resolution, transform | ms per step, fix | old | difference |
|---|---|---|---|
| T32 L8, MatrixSpectralTransform | 1.111 | 1.084 | +0.027 (+2.5 %) |
| T64 L8, MatrixSpectralTransform | 1.140 | 1.086 | +0.055 (+5.0 %) |
| T128 L8, SpectralTransform | 3.802 | 3.791 | +0.011 (+0.3 %) |

At every resolution the calls alone take 11–12 μs (`hole_filling!`, two broadcasts), 7 μs
(`hole_filling_tendency!`) and 3 μs (moving the fill back), ≈ 22 μs per step. The T64 difference is inflated by the
spread of the fix's timings (1.11–1.77 ms). Up to T64 the GPU time step is bound by kernel launches, so the fix costs
≈ 2 % there and nothing measurable from T128 on.

*Reactant on the GPU* (`docs/dev/2026-11/hole_filling_reactant_gpu.jl`: the PrimitiveWet correctness test of
`test/reactant/test_correctness.jl` for both transforms, the fill of CPU and Reactant after its 100 steps, and the
consistency of the Reactant state):
* Reactant 0.2.291 (the stale `test/GPU/CUDA` manifest): `transform!` fails to compile, "cannot raise op to stablehlo"
  on a complex zero constant (`llvm.mlir.constant` `array<2 x f32>`), with the old clip as well. Not the fix.
* Reactant 0.2.308 (freshly resolved environment), Reactant on the CPU (`CUDA_VISIBLE_DEVICES=`): the correctness
  test passes; after 100 steps CPU and Reactant hold the same fill (total 4.0363e-2 at 918 points); the current step
  of the Reactant state is consistent to 1.1e-8 (grid) and 5.6e-9 (fill).
* Reactant 0.2.308 on the GPU, default TF32 matrix products: the correctness test fails for the fix and the old clip
  alike (maximum relative difference 1.2–2.0 in all prognostic variables after 100 steps), the TF32 deviations
  described in `docs/src/architectures_gpu.md`.
* Reactant 0.2.308 on the GPU, strict Float32 (`PrecisionConfig.HIGHEST`, `XLA_FLAGS=--xla_gpu_enable_triton_gemm=false`,
  `NVIDIA_TF32_OVERRIDE=0`): all pass, 48/48 with the fix and 42/42 with the old clip, both transforms. After 100
  steps CPU and Reactant hold the same fill (`MatrixSpectralTransform`: 4.0446e-2 at 924 points; `SpectralTransform`:
  4.0681e-2 against 4.0680e-2 at 935 points), and the current step of the Reactant state is consistent to 2.8e-8 and
  1.3e-8 (grid humidity) and 4.0e-9 and 1.9e-9 (fill).

## Documentation changes

* `docs/src/primitiveequation.md`: new subsection "Hole filling" under "Humidity equation" (the clip, the fill as
  tendency from the previous step, not water conserving, the filter interplay), and the hole filling in steps 1 and
  19 of the algorithm.
* `CHANGELOG.md`: entry under Unreleased (PR number to be filled in).
* The `ClipNegatives` docstring and the docs subsection describe `fill_tendency` and its automatic default.

## Known limitations

* **Off by default except with `StochasticLeapfrog`.** Deterministic and SPPT runs keep the old inconsistency
  unless they opt in with `fill_tendency = true`. On the GPU at T32 L8 the old clip hides 77–210 mm/yr in holes
  (hole volume 0.01–0.03 mm), and at T64 the layer-2 spectral mean is already negative after one year (Audit, GPU).
* **Robert–Asselin–Williams filter.** The filter sees the unfilled previous step as part of the curvature
  `q⁽ⁿ⁻¹⁾ − 2q⁽ⁿ⁾ + q⁽ⁿ⁺¹⁾`. It moves ≈ ν α/2 ≈ 2.7 % of the fill back into the previous step as a hole and adds
  ≈ 2.3 % to the new step. The moved part is refilled two steps later; the computational mode this excites is
  damped by the filter. In the unit test (constant hole of −1e-5) about 1 % of the hole remains as extra water. An
  exact treatment needs the fill in spectral space, added to the leapfrog base (option 3 below).
* **Not water conserving.** The filling adds ≈ 4 mm/yr (0.5 % of the evaporation) at T32 L8 with noise.
* **GPU cost at low resolution.** Three small kernels and a copy per step, ≈ 22 μs on an H100: ≈ 2 % of a T32/T64
  time step, which is bound by kernel launches there. Fusing the two broadcasts of the clip into one kernel would save
  one launch.
* **Truncation.** The fill reaches the spectral humidity through the spectral truncation of the tendency, so small
  negatives remain on the grid (0.1–0.7 % of the points in layers 1–2, minima of a few 1e-6 to 2e-5 kg/kg).
* **`NCycleLorenz`.** Its weighted tendency memory spreads the fill over the cycle, so the spectral humidity is
  stepped from the clipped state only approximately.
* **Variables layout.** `vars.dynamics.humidity_hole_fill` is a new variable: states saved as a whole `Variables`
  from the old model may not load into the new one. A continued `run!` recomputes the fill in the initial transform.

## Options and their GPU cost (handover to a GPU node)

Option 1 is implemented (Summary of changes). Before settling on it, the GPU cost of every option below is to
be estimated on a GPU node. The options differ in cost, water conservation, accuracy and implementation effort.
All but option 6 make the spectral humidity integrate the humidity that the tendencies see.
Numbering here replaces the three fix options of the handoff (Background): handoff option 1 → 1,
handoff 2 → 2, handoff 3 → 4.

**Reference work per time step** of the PrimitiveWet model with L layers: three batched spectral transforms,
spectral→grid of the prognostic state (4L+1 layers: vorticity, divergence, temperature, humidity, ln pₛ),
spectral→grid of U, V and ∇ln pₛ (2L+2), and grid→spectral of the tendencies (9L+1: u, v, temperature,
humidity, ln pₛ, uT′, vT′, kinetic energy, uq, vq), i.e. ≈ 124 layer-transforms at L8. Humidity accounts for
4L of them (32 at L8). Then the column parameterizations, which are one fused kernel on GPU.

### How to measure

* **Code.** Option 1 is on branch `mg/holefilling-revised`. Option 0 (old behaviour, the baseline) is the base
  revision `fa90cc81`. Prototypes of other options go on their own commits or checkouts.
* **End to end** (options with an implementation: 0, 1 and every prototype). The debug sweep of the benchmark
  suite (PrimitiveWet, T32–T128 at L8, T32 and T128 at L16, `--length 2`) per revision or checkout, then a
  comparison table:

  ```bash
  julia SpeedyWeather/benchmark/regression/regression.jl benchmark fa90cc81 option0.json --arch gpu
  julia SpeedyWeather/benchmark/regression/regression.jl benchmark mg/holefilling-revised option1.json --arch gpu
  julia SpeedyWeather/benchmark/regression/regression.jl benchmark <checkout dir> option3.json --arch gpu
  julia SpeedyWeather/benchmark/regression/regression.jl table --result option0=option0.json \
      --result option1=option1.json --candidate option1
  ```

  ±20 % between processes is normal, so run the options to compare in the same session and repeat once.
* **Building blocks** (estimates without a full implementation). Time each block on the GPU after a warm-up:
  `CUDA.@elapsed` or `CUDA.@sync`, median of ≥ 20 calls. Do this at T32, T64 and T128 with L8, and T128 with L16
  (launch overhead dominates at low resolution, the transforms at high resolution). Compare against one full
  `SpeedyWeather.time_step!(simulation)`, after `initialize!(simulation; steps = n)` on a spun-up simulation.
  * **B1** (option 1): `SpeedyWeather.hole_filling!(vars, model.hole_filling, model)`,
    `SpeedyWeather.hole_filling_tendency!(vars, model)`,
    `SpeedyWeather.copy_step_back!(vars.dynamics.humidity_hole_fill, true)`, and the old clip
    `SpeedyWeather.hole_filling!(humid_grid, ClipNegatives())` for the difference.
  * **B2** (option 2): one unbatched L-layer grid→spectral transform,
    `transform!(spectral_L, grid_L, vars.scratch.transform_memory, model.spectral_transform)`.
  * **B3** (options 3, 8): the batched grid→spectral transform with 9L+1 layers against 10L+1 (option 3) and 5L+1
    (option 8, without q tendency, uq, vq); and the spectral→grid batch with 4L+1 against 3L+1 (option 8).
    Use standalone arrays and a `SpectralTransform` whose scratch memory covers the largest layer count.
  * **B4** (options 4, 5): a column kernel over the grid points with a loop over the L layers (the prototype of
    option 4), and a global mass-weighted sum of an L-layer field (`sum`/`mapreduce` with the area weights from
    `RingGrids.get_quadrature_weights` per ring), including the device synchronisation it implies.
  * **B5** (option 8): departure points for all grid points from u, v (two iterations) and the interpolation of
    an L-layer field to them with RingGrids' `AnvilInterpolator` (`update_locator!` + `interpolate!`).
* **Report** per option and resolution: the added time per time step in ms and in % of `time_step!`, and the
  SYPD of the debug sweep for the implemented options. Fill in the table at the end of this section.

### Option 0: clip the grid humidity only (old behaviour, baseline)

* **What.** `q = max(q, 0)` on the grid after the transform, the spectral humidity is left unchanged.
* **Cost.** One element-wise kernel over L×npoints. Has the inconsistency and drift described above.

### Option 1: the fill as tendency, from the previous step (implemented)

* **What.** Store `f = max(−q, 0)` per leapfrog step, add the previous step's `f/Δt` to the grid humidity
  tendency, which is transformed with the other tendencies (Summary of changes).
* **Where.** `src/dynamics/hole_filling.jl`, called from `transform!` and `time_step!`, fill moved back in
  `move_prognostic_grid_variables_back!` (`leapfrog.jl`).
* **Extra work per step.** No transform. Three element-wise kernels over L×npoints: writing the fill next to the
  clip, adding the tendency, moving the 2-step fill back. Memory: 2L×npoints. CPU T32 L8: +36 μs of 5.9 ms (+0.6 %).
  On GPU mostly launch overhead.
* **Measure.** B1 and end to end (option 1 against option 0).
* **Variants if the GPU cost matters.** (a) Fuse fill and clip into one kernel. (b) Use the current step's fill:
  one step of storage, no move-back, but it fills the other leapfrog sequence and skips a step in the leapfrog
  start.
* **Accuracy.** Water source 2.6–4.4 mm/yr. The time filter sees ≈ 2.5 % of the fill (Known limitations).

### Option 2: transform the clipped field back to spectral

* **What.** After the clip, overwrite the current spectral humidity with the transform of the clipped grid
  humidity, so grid and spectral agree at once. No fill storage, no tendency.
* **Where.** `transform!` (`src/time_stepping/transform.jl`): after the clip,
  `transform!(get_step(vars.prognostic.humidity, 2), humid_grid, scratch_memory, S)`. In the initial transform
  (`initialize = true`) also copy the result to the previous step (1), which the leapfrog start steps from.
* **Extra work per step.** One unbatched L-layer grid→spectral transform (≈ L/124 more transform work at L8, plus
  the launch overhead of a separate transform).
* **Measure.** B2, then end to end with a prototype (≈ 5 lines).
* **Accuracy.** The truncation of `max(q, 0)` makes new small negatives every step, so the clip–truncate cycle is a
  systematic water source: measure it with `docs/dev/2026-11/hole_filling_drift_check.jl`. It edits a leapfrog
  level outside the time stepper; the filter then sees the clipped level, which is consistent.

### Option 3: like 1, but the spectral fill is added to the previous step (exact with the time filter)

* **What.** Transform the previous step's fill to spectral together with the tendencies, and add it to the
  spectral humidity of the previous step (the leapfrog base) before the implicit terms and the leapfrog. Then
  the horizontal diffusion, the leapfrog and the Robert–Asselin–Williams filter all see the filled base. This
  removes the filter limitation of option 1.
* **Where.**
  * Add a grid slot for the fill to the `:grid_tendencies` fuse and a matching spectral slot to the
    `:spectral_tendencies` fuse (the fuse alignment check requires the same members in both).
  * In `dynamics_tendencies!` (`src/dynamics/tendencies.jl`), copy the fill of step 1 into the grid slot after
    `scale_tendencies!`, so that it is not scaled. The batched transform then produces the spectral fill. Add it
    to `get_step(vars.prognostic.humidity, 1)` right after `spectral_tendencies!`.
  * Zero the grid fill of step 1 once applied. In the leapfrog start the previous step is not advanced, and its
  fill would otherwise be added twice. The `dynamics = false` path (`parameterization_tendencies_only!`) needs the
  same.
* **Extra work per step.** The batched grid→spectral transform grows from 9L+1 to 10L+1 layers (+11 % of that
  batch), plus one spectral element-wise kernel. Option 1's storage and move-back stay, its tendency kernel goes.
* **Measure.** B3 (9L+1 against 10L+1) plus B1, then end to end with a prototype.
* **Accuracy.** Like option 1 without the filter effect. Water source unchanged.

### Option 4: conservative filler (borrow from the column, then a global fixer)

* **What.** Fill each hole from the same column, mass-weighted by Δσ: going down the column, a negative q_k takes
  its deficit from the layer below, `q_{k+1} += q_k Δσ_k/Δσ_{k+1}`, `q_k = 0` (pₛ cancels within a column). A
  deficit left in the lowest layer is clipped and paid back by a global multiplicative fixer (option 5), in the
  spirit of Rasch and Williamson (1990). It conserves water. Bottom-up borrowing is the alternative to evaluate:
  in the dry upper layers the layer below is often dry as well.
* **Where.** A new `AbstractHoleFilling` type in `src/dynamics/hole_filling.jl`. Its `hole_filling!` is a column
  kernel (`LinearWorkOrder` over npoints, loop over k, branches inside the kernel are fine). The change
  `q_filled − q_unclipped` (positive in holes, negative where water was borrowed) is stored as the fill and
  applied through option 1's tendency (or option 3's route). Guard the fixer's division for differentiability:
  `ifelse(c, x, 0) / ifelse(c, y, 1)`.
* **Extra work per step.** One column kernel (sequential in k), one or two global reductions with a device
  synchronisation, plus option 1's kernels.
* **Measure.** B4 plus B1, then end to end with a prototype.
* **Accuracy.** Conserves water (no 2.6–4.4 mm/yr source); moves water vertically, mostly within the dry upper
  layers.

### Option 5: clip plus a global water fixer

* **What.** Clip as in option 1, then rescale the clipped humidity so that its global mass-weighted water equals
  that of the unclipped humidity: `q ← q·W(q_unclipped)/W(q_clipped)`, globally or per layer. Store
  `q_new − q_unclipped` (zero global integral) as the fill and apply it through option 1's tendency.
* **Where.** `src/dynamics/hole_filling.jl`, a new type or a type parameter of `ClipNegatives` (not a Bool field,
  to avoid a runtime branch).
* **Extra work per step.** One global reduction (two sums) with a synchronisation, one element-wise kernel, plus
  option 1's kernels.
* **Measure.** B4 (the reduction) plus B1.
* **Accuracy.** Conserves water globally but not locally: dries moist regions in proportion to pay for the
  holes in dry regions.

### Option 6: no clipping

* **What.** `hole_filling = nothing`: grid and spectral humidity agree, and the physics sees negative humidity.
* **Where.** Model setting only. Every parameterization reading humidity needs an audit for q < 0: large-scale
  condensation, convection, clouds (relative humidity), surface fluxes, both radiation schemes (the default
  longwave already uses `max(q, 0)`), and virtual temperature.
* **Extra work per step.** None, one kernel less than option 0.
* **Measure.** No GPU cost to measure. The question is whether the holes stay bounded: run
  `hole_filling_drift_check.jl` with `hole_filling = nothing` and the default noise. Nothing restores negative
  values, so expect unbounded holes.

### Option 7: a positive-definite prognostic variable

* **What.** Prognose s = √q (or ln q) in spectral space, q = s² on the grid, positive by construction.
  Tendencies convert as ∂s/∂t = (∂q/∂t)/(2s), which is singular for s → 0 and needs a floor. The advection
  becomes the advective form u·∇s, so the flux form of q and the uq, vq slots go. Not water conserving.
* **Where.** The humidity variable and its advection (`humidity_grid_tendency!`, `humidity_spectral_tendency!`,
  vertical advection), the transform (q = s² instead of the clip), and the conversion of all parameterization
  tendencies.
* **Extra work per step.** About the same transform work: the advective form needs ∇s on the grid (2L
  spectral→grid), which replaces uq and vq (2L grid→spectral). Plus two element-wise kernels.
* **Measure.** Element-wise kernels only (B1 as a proxy). The cost is implementation effort and accuracy, not GPU
  time: prototype only if options 1–5 fail.

### Option 8: humidity as a grid-point variable with positive-definite advection

* **What.** Remove humidity from the spectral state. Step it in grid space with a semi-Lagrangian (or
  flux-limited finite-volume) horizontal advection with a monotone interpolation, as operational spectral models
  such as ECMWF's IFS do. Negative values then do not arise from the advection. Semi-Lagrangian advection is not
  conservative, so it needs a global fixer (option 5).
* **Where.** A redesign: humidity becomes a grid prognostic like the ocean and land variables, with a
  departure-point and interpolation step (`RingGrids.AnvilInterpolator`; `src/dynamics/particle_advection.jl`
  moves points with u, v and is a starting point). The vertical advection stays.
* **Extra work per step.** Saves 4L transform layers: q in the spectral→grid batch, and the q tendency, uq and vq
  in the grid→spectral batch, about a quarter of the transform layers at L8. Adds departure points, an L-layer
  interpolation and a reduction for the fixer. It may be cheaper than today on GPU.
* **Measure.** B3 (batches without humidity) against B5 + B4.
* **Accuracy.** Removes the holes at their source, but changes the humidity advection and its numerics
  everywhere. Largest effort.

### Results to fill in

| option | T32 L8 | T64 L8 | T128 L8 | T128 L16 | conserves water | removes the drift | effort |
|---|---|---|---|---|---|---|---|
| 0 clip only | baseline | baseline | baseline | baseline | no (hides holes) | no | — |
| 1 fill as tendency (prev. step) | +0.027 ms (+2.5 %) | +0.055 ms (+5.0 %, noisy) | +0.011 ms (+0.3 %) | | no (+2.6–4.4 mm/yr) | yes (10 yr CPU, GPU; 100 yr GPU) | done, chosen |
| 2 transform clipped back | | | | | no (truncation source) | expected | small |
| 3 spectral fill into the base | | | | | no (as 1) | expected | medium |
| 4 column borrowing + fixer | | | | | yes | expected | medium |
| 5 clip + global fixer | | | | | yes (globally) | expected | small |
| 6 no clipping | none | none | none | none | yes | no (holes unbounded) | audit |
| 7 √q or ln q | | | | | no | expected | large |
| 8 grid-point semi-Lagrangian | | | | | with fixer | expected | very large |

Entries: added time per step in ms and in % of `time_step!` (building blocks or end to end, say which).
Option 1: end to end on one H100, the `perf` mode of `docs/dev/2026-11/hole_filling_gpu_check.jl` (median ms per step of
`run!` over 1000/600/300 steps, against the old clip in the same process, not the benchmark suite), see GPU. B1
alone is ≈ 22 μs per step at every resolution (11–12 + 7 + 3 μs). Options 2–8 were not measured: option 1 was chosen
on 2026-11-20.

## Future work

* Option 1 is chosen (opt-in). If water conservation or the filter limitation matters later, options 4 or 5
  (applied through route 1 or 3) or option 3 remain.
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
GPU checks at `4681db88f`: Julia 1.12.2, CUDA.jl 6.2.2, KernelAbstractions 0.9.44 (`SpeedyWeather/test/GPU/CUDA`),
one H100 80GB per job (SLURM jobs 3369964–3372655); Reactant 0.2.308 with CUDA.jl 6.4.1 in a freshly resolved copy
of that environment (`/p/tmp/maxgelbr/code/claude_holefilling_gpu_2026-11-20/env_fresh`).

## Files

* `docs/dev/2026-11/hole_filling_mwe.jl` — standalone reproducer (SpeedyWeather only; not in this checkout).
* `docs/dev/2026-11/hole_filling_drift_check.jl` — the yearly diagnostics above, with and without the fix (not in
  this checkout).
* `docs/dev/2026-11/hole_filling_gpu_check.jl` — the GPU checks: consistency, the drift variants (fix against the old
  grid-only clip), timing, and the tables from the daily TSVs.
* `docs/dev/2026-11/hole_filling_reactant_gpu.jl` — Reactant against CPU with the fix (or `HF_VARIANT=old`), on the GPU
  or the CPU, optionally strict Float32 (`HF_STRICT=1`).
