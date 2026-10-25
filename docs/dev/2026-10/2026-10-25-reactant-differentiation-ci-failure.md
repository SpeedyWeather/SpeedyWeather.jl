# Reactant differentiation test fails after merging main into `mg/new-reactant`

> Status: **in progress, paused 2026-10-25**. The forward difference is found and fixed (section 4).
> The non-finite gradient comes from the physics (`dynamics_only` makes it finite), but the
> parameterization isn't identified yet. See "Handoff" for where to restart.

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
0.040, Δt = 2Δt = 4800 s), where the rational factor is about 0.5.

The damping factor each architecture actually applied (stored flux / undamped flux): the CPU applies
the expected `flux_damping(z)`. **Reactant applies exactly `flux_damping(z/2)`**, at every land point.

### 4. Cause of the forward difference: `default_time_step` inside kernels (fixed)

- Inside kernels the `Leapfrog` is adapted to a `LeapfrogCore` (`leapfrog.jl`,
  `Adapt.adapt_structure`).
- `default_time_step` was only defined as 2Δt for `Leapfrog`. `LeapfrogCore` fell back to
  `default_time_step(::AbstractTimeStepper) = L.Δt`.
- So every kernel calling `default_time_step(model.time_stepping)` got Δt instead of 2Δt. The only
  such kernel is the surface-flux damping of #1336.
- The CPU runs the column parameterizations as a plain loop with the unadapted model, so only
  kernel-launched paths were affected: **Reactant and the regular GPU path**. This is a bug from
  #1336 on main too.
- **Fix** (commit `0e9e523d`): `default_time_step(L::AbstractLeapfrog) = 2 * L.Δt`. A test in
  `test/dynamics/time_stepping.jl` checks that the adapted time steppers return the same default time
  step.
- **Verified:** CPU vs Reactant loss after 6 steps: **1.2e-6** (before 2.7e-3, tolerance 1e-4).
- Latent, unchanged: `LeapfrogCore` also drops `Leapfrog`'s separate `ocean`/`land` time steppers that
  `time_stepper` returns. With the defaults (`nothing`, the main stepper is used) this makes no
  difference.

### 5. Non-finite gradient (in progress)

With the forward fix the Reactant gradient of the 6 steps is still non-finite: 4480 of 8960 vorticity
entries. According to the user it is also non-finite with Enzyme on the CPU, so the cause is in code
that is actually executed, not only in the skipped side of a branch under Reactant.

**Scalar derivatives with Enzyme on the CPU** (`edge_derivatives.jl`):

| Formula | Derivative at the edge | Just above it |
| --- | --- | --- |
| `ocean_heat_roughness`, d/duₙ at uₙ = 0 (`surface_roughness.jl:541`, `exp(h_c4 / uₙ)`) | NaN | −3.0e-3 at uₙ = 1e-6 |
| `ocean_momentum_roughness`, d/duₙ at uₙ = 0 (`surface_roughness.jl:529`, `uₙ^m_c1`, m_c1 = 0.488) | Inf | 1.6e-28 at uₙ = 1e-6 |
| drag `max(drag_min, κ²/log(z/z₀)²)`, d/dz₀ at z₀ = 0 (`drag.jl:104`) | NaN | 1.6 at z₀ = 1e-4 |

**These edge inputs are common in the model** (CPU, default PrimitiveWet, after 10 steps):
- uₙ = 0 exactly at **40 % of ocean points** (1558 of 3933). There the clamp
  `max(…, 0)` at `neutral_wind.jl:70` is active, so the ocean roughness formulas are differentiated at
  uₙ = 0 in the branch that is actually taken.
- The drag's `log(z/z₀)` sees z₀ = 0 at 2758 ocean-only points (the land roughness) and 747 land-only
  points (the ocean roughness). The drag computes both coefficients everywhere.

**Mechanism.** `max`, `min`, `clamp` and `ifelse` are selects. In the reverse pass the side that isn't
selected gets a zero adjoint, which is then multiplied by that side's infinite local derivative:
0 × ∞ = NaN. The forward pass is finite throughout.

**Proposed fix**, to be checked first in a scratch script with method overrides: sanitize the input
before the singular formula, not only the result.
- Ocean roughness: evaluate at `max(uₙ, u_min)`.
- Drag: `ifelse(z₀ > 0, z₀, one(z₀))` inside the `log`.
- The threshold `u_min` slightly changes `uₙ^m_c1` at calm points (a factor `exp(u_min^0.488)`), so
  it is a modelling choice for the scheme's authors. The drag fix leaves forward values unchanged:
  those coefficients are weighted by a surface fraction of 0.

**Tested with both fixes applied** (`patched_gradient.jl both`, method overrides, u_min = 1e-6):
**still non-finite**, 4480 of 8960, unchanged.
- The structure is telling. 8960 = 560 spectral coefficients × 8 layers × 2 leapfrog time levels.
- The gradient w.r.t. one time level is non-finite in **every** coefficient. The other time level's
  gradient is all zero (the largest finite entry is 0).
- A NaN at even one grid point is spread over all spectral coefficients by the transform, so this
  doesn't localize the source.
- These edges are real non-differentiable spots, but not the (only) cause.

**Gradient bisection over the merged changes** (6 steps, Reactant, with the forward fix; each
variant in its own process):

| Variant | Reverted | Loss rel diff | Gradient |
| --- | --- | --- | --- |
| `default` | nothing | 1.2e-6 | non-finite (4480 of 8960) |
| `oldbl` | boundary layer as before the merge (no neutral wind speed, constant roughness) | 2.5e-6 | non-finite (4480 of 8960) |
| `noneutral` | `neutral_wind_speed = nothing` | 3.1e-6 | non-finite (4480 of 8960) |
| `explicitfluxes` | `ExplicitFluxes` for momentum, heat and humidity | 1.3e-6 | non-finite (4480 of 8960) |
| `all-old` | all of the above plus `SigmaCoordinates` and no schedule | — | compile didn't finish in 64 min (stuck in LLVM GVN), aborted |

So the new boundary layer (#1184) and the flux time integration (#1336) are **not** the cause of the
non-finite gradient under Reactant. The runs `sigma`, `vdiff` and `noschedule` were stopped in favour of
the following lead.

**Not covered by any variant: the new default surface condition** (also from #1184).
- Before the merge, the surface wind speed was `sqrt(uₛ² + vₛ² + 1²)` with a constant gust of 1 m/s.
- Now the gust is `BeljaarsWindGust` (`surface_condition.jl:294-298`):

  ```julia
  ΔTᵥ = max(Tᵥ_skin - Tᵥ, zero(Tᵥ))      # 0 in every stable column
  gust = sqrt(β³ * C_H * z_i * (gravity / max(Tᵥ, 200)) * ΔTᵥ)
  gust = max(gust, min_gust)
  ```

- `sqrt(0)` in stable columns has an infinite derivative. The adjoint of the `sqrt` argument is
  0 × ∞ = NaN. Through the product it reaches `gravity / max(Tᵥ, 200)` (NaN × ΔTᵥ = NaN × 0 = NaN), and
  through `max(Tᵥ, 200)` the lowermost-layer temperature, i.e. the model state. The spectral transform
  then spreads it over every coefficient.
- Scalar check: Enzyme on the CPU gives a derivative of 0 here, not NaN; it special-cases `sqrt` at 0.
  Enzyme-MLIR's reverse rule for `stablehlo.sqrt`, d/(2·sqrt(x)), gives 0/0. So this is a candidate
  for the **Reactant** NaN.
- The CPU NaN would then come from the spots that do give NaN with CPU Enzyme: the ocean roughness at
  uₙ = 0 and the drag at z₀ = 0, see above.
- **Exact fix:** take the `max` before the `sqrt`,
  `gust = sqrt(max(β³ C_H z_i g/Tᵥ · ΔTᵥ, min_gust²))`. The forward value is identical (`sqrt` is
  monotonic), and `sqrt` is never evaluated at 0.
- **Tested, not the (only) cause:** with `patched_gradient.jl gust` and with
  `gust,roughness,drag` the gradient is **still non-finite** (4480 of 8960, finite part all zero).
  The loss changes slightly, so the patches are applied.
- These edges are real non-differentiable spots worth fixing, but the main source is elsewhere.

**Systematic bisection** (gradient over 6 steps, Reactant):

| Variant | Change | Gradient |
| --- | --- | --- |
| `nophysics` | `dynamics_only = true`: no parameterizations, ocean, sea-ice or land time steps | **finite** |
| `sigma` | `SigmaCoordinates` (before #1346) | non-finite |
| `vdiff` | `BulkRichardsonDiffusion` on (before #1346) | non-finite |
| `noschedule` | implicit schedule never (#1332) | non-finite |

**The NaN comes from the physics.** Of the merged changes, the coordinates, the vertical diffusion,
the schedule, the flux integration and the boundary layer apart from its surface condition are each
ruled out on their own. Next round, splitting the physics: `constsc` (surface condition as before:
constant slowdown and gust), `nosurface` (boundary layer and surface fluxes off), `noveg` (no
vegetation), `nolsc` (no large-scale condensation).

## Handoff (2026-10-25, session closed): where we are and how to restart

### State

- **Forward difference: fixed and committed** (`0e9e523d`, `default_time_step(::AbstractLeapfrog)`,
  section 4). CPU vs Reactant loss after 6 steps: 1.2e-6.
- **Implicit `reinitialize!` with #1332's schedule: reworked, committed and tested** (`e6044288`,
  `test/reactant/test_implicit.jl`; see "Open questions" item 3).
- **Non-finite gradient: open.**
  - Known: the source is in the **physics**. With `dynamics_only = true` the gradient is finite.
  - Ruled out, each reverted on its own: the boundary layer (neutral wind speed, roughness), the flux
    time integration, the Frierson coordinates, vertical diffusion off, and the implicit schedule.
  - Fixing the three non-differentiable edges found (Beljaars gust `sqrt(0)`, ocean roughness at
    uₙ = 0, drag `log(z/0)`) does not make it finite either.
  - **Hint:** the gradient w.r.t. one leapfrog time level is non-finite in every coefficient, the other
    is all zero. So the NaN reaches the state and the transform spreads it everywhere.
- **Interrupted at closing** (results lost, rerun): the physics round `constsc`, `noveg`, `nolsc`.
  `nosurface` failed by construction: other components need the surface-flux variables (`FieldError:
  … no field sensible_heat_flux`). Turning the surface fluxes off needs a different variant, e.g. zero
  drag or prescribed fluxes.
- **The user reports non-finite gradients with Enzyme on the CPU too.** Not reproduced here yet:
  plain Enzyme on the CPU needs Julia 1.10 or 1.11 (see `differentiation_primitive_wet.jl`).
- **Repository:** all code changes are committed (the Runic whitespace cleanup of `leapfrog.jl` as
  `2e303303`). `SpeedyWeather/test/reactant/Project.toml` was restored after rebuilding the
  environment. Only this document has changes since the last commit.

### Restart

1. **Environment.** Julia 1.12, as in CI. The manifest is gitignored; if it's stale (Reactant
   < 0.2.299), regenerate it:

   ```bash
   rm SpeedyWeather/test/reactant/Manifest.toml
   julia +1.12 --project=SpeedyWeather/test/reactant -e 'using Pkg; Pkg.develop([PackageSpec(path=p) for p in ("SpeedyWeatherInternals", "LowerTriangularArrays", "RingGrids", "SpeedyTransforms", "SpeedyWeather")]); Pkg.instantiate()'
   git checkout SpeedyWeather/test/reactant/Project.toml   # Pkg.develop adds the monorepo packages to [deps]
   ```

   Run the Reactant scripts **without `--check-bounds=yes`**. Forced bounds checks make the kernels
   impossible to raise ("MLIR pass pipeline all failed").
2. **Scripts.** Recreate `repro_setup.jl` and `repro.jl` from the appendix in a scratch directory. One
   variant per process:

   ```bash
   julia +1.12 --project=SpeedyWeather/test/reactant repro.jl <variant> gradient
   ```

   - It prints the CPU vs Reactant loss after 6 steps, then whether the Reactant gradient of the 6 steps
     is finite.
   - It takes about 10 to 15 minutes per variant, mostly the gradient compilation. Up to 4 processes ran
     in parallel on a 16 GB M3.
   - Some combinations compile pathologically: `all-old` was stuck in LLVM GVN for over an hour. Cap
     runs, e.g. with `perl -e 'alarm shift; exec @ARGV' 1800 julia …`.
3. **Next bisection round inside the physics:**
   - Run `constsc` (surface condition as before the merge: `ConstantWindSlowdown`, `ConstantWindGust`),
     `noveg` and `nolsc`.
   - Then the remaining parameterizations one at a time (component `= nothing` where the model allows
     it): `radiation` (with `albedo`), `solar_zenith`, the ocean, sea-ice and land time steps
     (`ocean`, `sea_ice`, `land`), `greenhouse_gases`.
   - Decision tree:
     - If `constsc` is finite, the cause is in `LearnedWindSlowdown`: `U = sqrt(u² + v²)`,
       `log(z₀)` with the roughness of the previous step, `clamp`. The Beljaars gust alone is already
       ruled out by the patch test.
     - If a group is finite, bisect inside it.
     - Once the parameterization is known, locate the term by scalar derivatives at the inputs that
       occur in the model, as in `edge_derivatives.jl` (section 5).
4. **Fixes to make regardless** (all with identical forward values or tiny changes, see section 5):
   - Beljaars gust: `sqrt(max(arg, min_gust²))`.
   - Ocean roughness: evaluate at `max(uₙ, u_min)`; the threshold is for the scheme's authors.
   - Drag: sanitize z₀ inside `log`.

   These also concern Enzyme on the CPU, which does give NaN for the roughness and drag edges.

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

## Appendix: reproduction scripts

Not committed as code. Recreate them in a scratch directory to restart (see Handoff). `REPO` is the
local repository path.

### `repro_setup.jl`

```julia
#=
Reproduce the CI failure of "Differentiability of 6 time steps with Reactant: PrimitiveWetModel" after
merging main (#1332 implicit reference temperature schedule, #1336 implicit surface fluxes, #1346
FriersonSigmaCoordinates, #1184 learned surface roughness) and isolate the cause by swapping components.

    julia +1.12 --project=SpeedyWeather/test/reactant repro.jl <variant>[,<variant>...] [gradient]

Variants: default, noschedule, constroughness, constdrag, sigma, explicitfluxes, all-old
=#
using SpeedyWeather, Reactant, CUDA, Enzyme, Test, Printf

const TRUNCATION = 32
const REPO = "/Users/max/Nextcloud/SpeedyWeather/alt-version/SpeedyWeather.jl"
include(joinpath(REPO, "SpeedyWeather/test/reactant/setup.jl"))
const DIFFERENTIATION_FUNCTIONS_ONLY = true
include(joinpath(REPO, "SpeedyWeather/test/reactant/differentation.jl"))

# components replacing the merged defaults, by variant, as functions of the spectral grid
explicit_fluxes(sg) = (;
    surface_momentum_flux = SurfaceMomentumFlux(sg; time_integration = ExplicitFluxes()),
    surface_heat_flux = SurfaceHeatFlux(
        sg; ocean = SurfaceOceanHeatFlux(sg; time_integration = ExplicitFluxes()),
        land = SurfaceLandHeatFlux(sg; time_integration = ExplicitFluxes()),
    ),
    surface_humidity_flux = SurfaceHumidityFlux(
        sg; ocean = SurfaceOceanHumidityFlux(sg; time_integration = ExplicitFluxes()),
        land = SurfaceLandHumidityFlux(sg; time_integration = ExplicitFluxes()),
    ),
)
const VARIANTS = Dict(
    "default" => sg -> (;),
    "noschedule" => sg -> (; implicit = ImplicitPrimitiveEquation(sg; schedule = PeriodicSchedule(every = Second(0)))),
    "constroughness" => sg -> (; boundary_layer = BoundaryLayer(sg; surface_roughness = ConstantSurfaceRoughness(sg))),
    "constdrag" => sg -> (; boundary_layer = BoundaryLayer(sg; drag = ConstantDrag(sg))),
    "noneutral" => sg -> (; boundary_layer = BoundaryLayer(sg; neutral_wind_speed = nothing)),
    # the boundary layer before the merge: no neutral wind speed, constant roughness
    "oldbl" => sg -> (; boundary_layer = BoundaryLayer(sg; neutral_wind_speed = nothing, surface_roughness = ConstantSurfaceRoughness(sg))),
    "sigma" => sg -> (; geometry = Geometry(sg; vertical_coordinates = SigmaCoordinates(sg))),
    "explicitfluxes" => explicit_fluxes,
    # all parameterizations off (time_integration.jl: `if ~model.dynamics_only`)
    "nophysics" => sg -> (; dynamics_only = true),
    # surface condition as before the merge: constant wind slowdown and gust
    "constsc" => sg -> (; boundary_layer = BoundaryLayer(sg; surface_condition = SurfaceCondition(sg; wind_slowdown = ConstantWindSlowdown(sg), wind_gust = ConstantWindGust(sg)))),
    # no surface fluxes at all
    "nosurface" => sg -> (; boundary_layer = nothing, surface_momentum_flux = nothing, surface_heat_flux = nothing, surface_humidity_flux = nothing),
    "noveg" => sg -> (; land = LandModel(sg; vegetation = NoVegetation(sg))),
    "nolsc" => sg -> (; large_scale_condensation = nothing),
    # vertical diffusion on, the default before #1346
    "vdiff" => sg -> (; vertical_diffusion = BulkRichardsonDiffusion(sg)),
    "all-old" => sg -> (;
        implicit = ImplicitPrimitiveEquation(sg; schedule = PeriodicSchedule(every = Second(0))),
        boundary_layer = BoundaryLayer(sg; neutral_wind_speed = nothing, surface_roughness = ConstantSurfaceRoughness(sg)),
        geometry = Geometry(sg; vertical_coordinates = SigmaCoordinates(sg)),
        explicit_fluxes(sg)...,
    ),
)

# as create_model in setup.jl, plus the components of `variant`
function make_model(ModelType, arch, variant)
    sg = SpectralGrid(; architecture = arch, nlayers = nlayers_for_model(ModelType), truncation = TRUNCATION)
    return ModelType(
        sg; spectral_transform = SpectralTransform(sg), feedback = nothing,
        initial_conditions = InitialConditions(sg, ModelType),
        model_kwargs(ModelType, sg)..., VARIANTS[variant](sg)...,
    )
end

# as differentiability_simulations in differentation.jl
function simulations(variant; spinup_steps = 5)
    simulation_cpu = initialize!(make_model(PrimitiveWetModel, SpeedyWeather.CPU(), variant))
    simulation_reactant = initialize!(make_model(PrimitiveWetModel, SpeedyWeather.ReactantDevice(), variant))
    run!(simulation_cpu, steps = spinup_steps)
    run!(simulation_reactant, steps = spinup_steps)
    sync_variables!(simulation_cpu, simulation_reactant)
    initialize!(simulation_cpu, steps = 1)
    initialize!(simulation_reactant, steps = 1)
    sync_variables!(simulation_cpu, simulation_reactant)
    return simulation_cpu, simulation_reactant
end

```

### `repro.jl`

```julia
#=
Reproduce the CI failure of "Differentiability of 6 time steps with Reactant: PrimitiveWetModel" after
merging main (#1332 implicit reference temperature schedule, #1336 implicit surface fluxes, #1346
FriersonSigmaCoordinates, #1184 learned surface roughness) and isolate the cause by swapping components.

    julia +1.12 --project=SpeedyWeather/test/reactant repro.jl <variant>[,<variant>...] [gradient]

Variants: default, noschedule, constroughness, constdrag, sigma, explicitfluxes, all-old
=#
include(joinpath(@__DIR__, "repro_setup.jl"))

variants = split(ARGS[1], ",")
with_gradient = length(ARGS) >= 2 && ARGS[2] == "gradient"
n = 6
for variant in variants
    try
        simulation_cpu, simulation_reactant = simulations(variant)
        vars_cpu, model_cpu = SpeedyWeather.unpack(simulation_cpu)
        vars_reactant, model_reactant = SpeedyWeather.unpack(simulation_reactant)
        loss_cpu = loss_n_steps(deepcopy(vars_cpu), model_cpu, n)
        loss_reactant = Float64(@jit compile_options = REACTANT_AD_COMPILE_OPTIONS loss_n_steps(deepcopy(vars_reactant), model_reactant, n))
        s = @sprintf("RESULT %-15s loss CPU %.6e Reactant %.6e rel diff %.2e", variant, loss_cpu, loss_reactant, abs(loss_reactant - loss_cpu) / loss_cpu)
        println(s); flush(stdout)
        if with_gradient
            gradient = gradient_vorticity_reactant(deepcopy(vars_reactant), model_reactant, n, (;))
            println(@sprintf("RESULT %-15s gradient finite %s (%d non-finite of %d)", variant, all(isfinite, gradient), count(!isfinite, gradient), length(gradient)))
        end
    catch err
        println("RESULT $variant FAILED: ", first(sprint(showerror, err), 1500))
    end
    flush(stdout)
end
```
