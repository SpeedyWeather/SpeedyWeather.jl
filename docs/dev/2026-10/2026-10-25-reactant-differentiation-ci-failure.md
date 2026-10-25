# Reactant differentiation test fails after merging main into `mg/new-reactant`

> Status: **in progress**. Reproduced locally; the forward difference comes from combining the new
> rational surface-flux time integration with the new Frierson sigma coordinates, and the exact
> spot inside the time step is still being located. The cause of the non-finite gradient is not
> yet investigated.

Date of initial draft: 2026-10-25

Base revision: 19247bfc2c7de9cbed5a609159153df1b32a9498 (branch `mg/new-reactant`, merge of `origin/main`)

## Originating prompt

> The Reactant correctness tests are failing after our latest merge of the corrected implicit stepping
> https://github.com/SpeedyWeather/SpeedyWeather.jl/actions/runs/40102090020/job/121025046501?pr=1304

## Revision log

> write a markdown file for the investigation that you continue to update

This document was started after the first two rounds of bisection. It is updated as the
investigation goes on.

## Problem description

In CI job 121025046501 (Reactant Correctness, Julia 1.12.7, ubuntu-latest, Reactant 0.2.299), one
test set fails: "Differentiability of 6 time steps with Reactant: PrimitiveWetModel"
(`SpeedyWeather/test/reactant/differentation.jl`).
- **Forward loss** (line 139): after 6 time steps the loss is 443188.75 on Reactant against 444370.88 on
  the CPU. That is 2.7e-3 relative, with a tolerance of `rtol = 1e-4`.
- **Gradients** (line 152): the Reactant gradient of the 6 steps contains non-finite values, both with
  `(;)` (every step stored) and with `checkpointing = Reactant.Binomial(2)`.

All other Reactant tests pass, including the CPU-vs-Reactant comparisons of the PrimitiveWet model and
`time_steps!` vs single compiled time steps. Those compare with `rtol = 0.05, atol = 0.002`, much looser
than the differentiation test's 1e-4.

The test passed on `ff459bd9` (Reactant Correctness, 2026-10-22) and fails on the merge `19247bfc`.

## Background

**What the merge brought in from main:**

| PR | Change relevant inside the time step |
| --- | --- |
| #1332 | The implicit solver also recomputes its reference temperature profile on a schedule, `PeriodicSchedule(every = Day(1))`, not only when the time step changes |
| #1336 | Surface fluxes are damped for implicit-like time integration, default `RationalFluxes`: factor `(6 + z)/(6 + 4z + z²)` with `z = Δt g C / Δp` of the lowermost layer, `Δt = default_time_step` (2Δt for leapfrog) |
| #1346 | Default vertical coordinates `FriersonSigmaCoordinates` (thinner lowermost layer); vertical diffusion off by default |
| #1184 | Boundary layer defaults: `NeutralWindSpeed` (was `nothing`) and `LearnedSurfaceRoughness` (was `ConstantSurfaceRoughness`) |

**What the test does.**
- `differentiability_simulations` builds a CPU and a Reactant PrimitiveWet model: T32L8,
  `SpectralTransform`, `convection = nothing`, constant longwave transmissivity (`setup.jl`).
- Both spin up 5 steps. Then the variables are synced Reactant → CPU, both are re-initialized for one
  step, and the variables are synced again.
- `loss_n_steps` runs 6 steps (`time_steps!`) and returns `sum(abs2, vorticity)`. On Reactant it is
  compiled with `reactant_compile_options(; raise_first = true)`.
- Model components (e.g. the implicit solver's state) are not synced, only the variables.

**A known gap, not the cause here.** The Reactant extension overrides `reinitialize!` for the implicit
solver (`SpeedyWeatherReactantExt.jl`): it updates the reference temperature profile only when the time
step changed, and doesn't know about #1332's schedule. With Δt = 2400 s the daily schedule fires every
36 steps, so not within this 6-step test. It still needs fixing for longer runs.

## Setup

- **Environment.** `SpeedyWeather/test/reactant` rebuilt as CI does: Julia 1.12.7, `Pkg.develop` of the
  five monorepo packages, Reactant 0.2.299, Reactant_jll 0.0.428, Enzyme 0.13.216. The old manifest
  (Julia 1.11.9, Reactant 0.2.292) is backed up in the session scratchpad.
  - `Pkg.develop` adds the monorepo packages to `test/reactant/Project.toml`; restore it with
    `git checkout SpeedyWeather/test/reactant/Project.toml` when done.
- **Scripts.** In the session scratchpad, `ci_repro/`, not committed:
  - `repro_setup.jl`: the test's model setup plus variants that swap the merged components back to
    their old versions.
  - `repro.jl <variants> [gradient]`: CPU vs Reactant loss after 6 steps per variant, optionally the
    Reactant gradient.
  - `one_step.jl <variant> [nsteps]`: from the synced state, the largest CPU-vs-Reactant differences of
    all variables after `nsteps` steps.

## Protocol

### 1. Reproduction and bisection of the forward difference

Relative difference of the 6-step loss, CPU vs Reactant, with one or several merged changes reverted
through model components:

| Variant | Reverted | Relative loss difference |
| --- | --- | --- |
| `default` | nothing (as CI) | **2.66e-3** (CPU 4.443954e5, Reactant 4.432127e5; CI: 2.7e-3) |
| `noschedule` | implicit schedule (`every = Second(0)`) | 2.66e-3, identical to `default` |
| `constroughness` | `ConstantSurfaceRoughness` | 2.52e-3 |
| `sigma` | `SigmaCoordinates` | **1.7e-5** |
| `explicitfluxes` | `ExplicitFluxes` for momentum, heat and humidity | **1.3e-6** |
| `all-old`, first version | schedule, roughness and drag constant, `SigmaCoordinates`, `ExplicitFluxes`; neutral wind speed still on | 3.3e-6 |

Findings:
- The difference reproduces locally with the same numbers as CI.
- The implicit schedule (#1332) and the learned roughness (#1184) are not the cause of the forward
  difference. Neither is `NeutralWindSpeed`: it was still on in `all-old`.
- Reverting either the flux integration (#1336) or the vertical coordinates (#1346) alone brings the
  difference below the tolerance. The cause is therefore their combination: rational damping of the
  surface fluxes into a thin lowermost layer.

### 2. Where inside the time step

`one_step.jl default 1` compares all variables (prognostic, grid, tendencies, dynamics,
parameterizations) after one step from the synced state. Largest differences, as max|Δ| / max|x|:

| Variable | CPU vs Reactant after 1 step |
| --- | --- |
| `tendencies.grid.temperature` | 4.2e-1 |
| `parameterizations.land.sensible_heat_flux` (and the total) | 3.7e-1 |
| `parameterizations.land.surface_humidity_flux` (and the total) | 3.7e-1 |
| `tendencies.land.soil_temperature`, `tendencies.land.soil_moisture` | 3.1e-1, 3.0e-1 |
| `parameterizations.ocean.sensible_heat_flux`, `ocean.surface_humidity_flux` | 3.4e-2, 3.3e-2 |
| `grid.vorticity`, `tendencies.vorticity` | 2.5e-2, 1.6e-2 |

All values are finite. The **land** surface fluxes differ by 37 % after a single step, and the ocean
fluxes by 3 %. That is not a rounding difference: the Reactant path computes the (damped) land surface
fluxes differently.

### 3. Land sensible heat flux: inputs or kernel?

`land_flux.jl` saves the synced start state (soil temperature, snow depth), takes one step on both
architectures, and recomputes the flux of `SurfaceLandHeatFlux` on the host from each architecture's
own inputs. The inputs are surface pressure, air density, wind speed, air temperature and drag.

| | max\|Δ\| / max\|x\| (land points) |
| --- | --- |
| Inputs, CPU vs Reactant | ≤ 4.2e-6 (drag), the others ≤ 6.5e-7 |
| CPU: stored flux vs flux recomputed from its own inputs | 5.0e-8 |
| Reactant: stored flux vs flux recomputed from its own inputs | **3.7e-1** |

The inputs are the same, but **Reactant's kernel computes a different flux from them**. The damping
factor is the prime suspect: z = Δt/τ reaches 1.66 at land points (Δσ of the lowermost Frierson layer
0.040, Δt = 2Δt = 4800 s), where the rational factor is about 0.5. Next: the damping factor each
architecture actually applied (stored / undamped flux) against z.

## Open questions and next steps

1. **Forward difference.** Find the quantity where CPU and Reactant first diverge (section 2), then
   whether it's a bug in the Reactant path or a difference amplified by the stiff coupling.
2. **Non-finite gradient.** Not investigated yet. Hypothesis: inside kernels raised to StableHLO an early
   `return` (land or ocean points) becomes a select that evaluates both sides, and the side not taken
   still gets differentiated. The new boundary-layer code has terms with infinite derivatives at zero
   wind speed, e.g. `exp(h_c4 / uₙ)` and `uₙ^m_c1` in `LearnedOceanRoughness`, evaluated there with
   uₙ = 0 at land points. 0 × ∞ gives NaN. Enzyme on the CPU follows the actual branch and never sees
   them. To test: `repro.jl ... gradient` with `oldbl` (old boundary layer), `noneutral` and
   `constroughness`.
3. **Implicit schedule on Reactant: addressed separately** (2026-10-25, independent of this failure).
   - The Reactant-specific override of `reinitialize!` is gone. The core method now decides, without
     branching on the clock, `update = (implicit.Δt[] != Δt) | isscheduled(implicit.schedule, clock)`,
     and dispatches on the type of `update`:
     - `update::Bool` (CPU, GPU) rebuilds the operators only when needed, as before.
     - `update::TracedRNumber{Bool}` (SpeedyWeatherReactantExt) updates the reference profile with
       `ifelse` and rebuilds the operators every step.
   - `isscheduled(every_n_steps, clock)` uses `&` on the step counter.
   - `PeriodicSchedule` counts executions through `count_execution!`, a no-op for a traced decision.
   - Tested in `test/reactant/test_implicit.jl` against the CPU: `StepSchedule`, `PeriodicSchedule`,
     and a schedule that never executes. All 12 checks pass locally (committed as `e6044288`).
   - CPU: the schedule and implicit-solver tests pass, and a 40-step run stays finite with one
     scheduled reinitialization (step 36).
   - **Pitfall:** run with `--check-bounds=yes`, the Reactant test fails with "MLIR pass pipeline all
     failed". Forced bounds checks keep the throwing branches in the KernelAbstractions kernels (here
     the Gauss-Jordan inversion), and those can't be raised. CI runs the Reactant tests without the
     flag, and so must local runs.
