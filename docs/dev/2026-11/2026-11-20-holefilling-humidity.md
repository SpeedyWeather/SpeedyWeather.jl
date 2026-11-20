# Handoff: the PrimitiveWet hole filling lets the spectral humidity run away (2026-11-20)

**Goal.** A hole filling in SpeedyWeather's `PrimitiveWetModel` (branch `mg/gfdt-upstream`) that keeps the spectral
prognostic humidity consistent with the grid humidity the model actually uses, so that the upper-level spectral
humidity stops drifting. This blocks the project's next control run (`docs/dev/plan-v025-recalibration-2026-11-20.md`,
stage 0.2): the exported state is the spectral coefficients, and a coordinate that drifts linearly over 1500 years is
a clock for the score network.

**Status.** Cause found and reproduced with SpeedyWeather alone (`docs/upstream/hole_filling_drift_mwe.jl`); no fix yet.

## Symptom

1500-yr control run (T32 L8 `PrimitiveWetModel`, `StochasticLeapfrog` + `AdditiveWienerNoise` with the default
amplitudes, default `SpeedLimitDrag`, Δt 40 min, SpeedyWeather `c19adc3d9`; `/p/tmp/maxgelbr/scratch_v025/control_t32l8_white`):

* the layer-2 spectral humidity's global mean (l = 0 coefficient / √4π) falls **linearly from ≈ 0 to −0.53 kg/kg**
  (R² 1.000); layer 1 goes slightly negative (run mean −1e-4 kg/kg); layers 3–8 and every temperature are flat;
* the grid humidity the physics uses is **never negative**: layer 2 between 0 and 2.2e-3 kg/kg, no point below zero;
* the old model (SpeedyWeather-main `496f1c6c`, 0.22.1+DEV) has the same `ClipNegatives` and the same noise
  amplitudes but no drift: layer 2 +3.4e-4 kg/kg and stationary. The new model's layers 1 and 2 are 18 K and 11 K
  colder (222.5 / 221.7 K against 240.5 / 233.1 K), so they hold almost no moisture, and noise and transport push that
  near-zero field below zero far more often.

## Cause

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

## Evidence

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

## Reproduce

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

## Fix options

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

## Acceptance tests

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

## After the fix (project side)

* Update the pinned environments (speedy-upstream, speedy-upstream-gpu, speedy-ad) on the login node; Pkg cannot
  reach the internet from compute nodes.
* Run the five-day bitwise check of the dry and barotropic models against `c19adc3d9` (they should not change).
* The new control run (plan stage 0.3) and its drift check (`slurm/drift_check.sh`, which flags any coordinate with a
  linear-trend R² > 0.5 and prints block means of the slow PCs).

## Versions

SpeedyWeather `0.25.0-DEV`, branch `mg/gfdt-upstream` at `c19adc3d9` (git tree `55c4d100`; `fb912357a` was checked
bitwise identical for this model); Julia 1.11.7; NVIDIA GPU via CUDA.jl's own artifacts (no system CUDA module).

## Files

* `docs/upstream/hole_filling_drift_mwe.jl` — standalone reproducer (SpeedyWeather only).
* `speedy/checks/humidity_drift_test.jl` — per-day, per-layer diagnostics; variants from the end state and from rest.
* `speedy/pwet_v025.jl` — the project's model builder (`hum_frac` scales the humidity noise).
* `python/checks/drift_check.py`, `docs/figures/drift-v025-control.png` — the drift check that found it.
* `docs/dev/PROGRESS.md`, entries of 2026-11-20 — the full record.
