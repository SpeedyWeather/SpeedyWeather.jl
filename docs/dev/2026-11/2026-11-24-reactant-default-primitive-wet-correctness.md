# Reactant correctness of the default PrimitiveWetModel

> Status: **in progress**. Default PrimitiveWetModel on Reactant reproduces the CPU to rounding; the largest
> deviations (divergence) came from the Float32 transform of the full geopotential in the Bernoulli
> potential, fixed (§12), tolerances tightened; the Bernoulli fix and new tolerances uncommitted (on top of 70cfdf95).

Date of initial draft: 2026-11-24

Base revision: 75bf89d9aeac994eb8857a06d136a5ceafec0099

## Originating prompt

> I want to achieve full Reactant correctness of the default PrimitiveWetModel, not just one without
> convection and with constant longwave trasnmissitivty

## Revision log

> the test enviroment should also use the most recent Reactant environment

`SpeedyWeather/test/reactant/Project.toml`: compat `Reactant = "0.2.309"` (was `0.2.305`). 0.2.309
is the latest registered version (2026-11-24); CI's last run on `main` already resolved to it. The
packages' own Reactant compat (`0.2.305` in all five) is unchanged.

> keep a markdown log of this investigation.

This document.

> okay keep going, but before that make the adjustment in the `mg/reactant-restart` branch that we
> test the default model there and don't pin Reactant anymore to a specific version

Commit 6656c88e on `mg/reactant-restart` (#1437, not pushed): `model_kwargs(::Type{PrimitiveWetModel}, _) = (;)`
in `test/reactant/setup.jl` (the PrimitiveDryModel keeps `convection = nothing`), the `[compat]` section
removed from `test/reactant/Project.toml` (the test environment then resolves to the latest Reactant the
packages' own compat allows, ≥ 0.2.305), CHANGELOG line. The same `setup.jl` and `Project.toml` changes
on this branch (replacing the 0.2.309 compat above), so both merge cleanly. Note: the CI differentiation
test (`test_differentiability_time_steps(PrimitiveWetModel; ...)` in `differentation.jl`) also builds its
model via `model_kwargs`, so it now differentiates the default model too; not yet checked locally.

> apply the same fix for the climatolgies, and then try to introduce somewhat tighter tolerences for
> the existing tests to pass

§10 (climatologies, all five file-based climatologies converted to the model's number format before
interpolating) and §11 (tolerances).

> okay your task is now to find the culprit for the biggest deviations between both models, this seems to
> be in the divergence. Try to work on use decreasing these tolerences in the end

§12.

## Problem description

`SpeedyWeather/test/reactant/setup.jl` compares CPU and Reactant only for a modified
PrimitiveWetModel: `ConstantLongwaveTransmissivity` and `convection = nothing`, with the comment
"Frierson's surface-concentrated absorption without convection amplifies ULP-level CPU-vs-Reactant
differences exponentially (chaotic, not a Reactant bug)". The PrimitiveDryModel is tested with
`convection = nothing`. The goal is CPU-vs-Reactant correctness of the default PrimitiveWetModel.

## Background

- The overrides come from #1304 (2026-10-27, `docs/dev/2026-10/2026-10-12-reactant-full-compatibility.md`).
  That doc gives no evidence for the "chaotic amplification" claim. Since then the default longwave
  transmissivity changed from Frierson to `ByrneOGormanLongwaveTransmissivity` (#1279, retuned in
  #1414, `docs/dev/2026-10/2026-10-30-longwave-byrne-ogorman.md`), so the comment is outdated.
- Relative to the test configuration, the default model adds only two column parameterizations:
  `BettsMillerConvection` (with `LinearEntrainment`) and `ByrneOGormanLongwaveTransmissivity` (which
  reads CO₂ from `vars.prognostic.greenhouse_gases.co2[]`, a `Ref` of a Reactant number, inside the
  column kernel). Shortwave, clouds, land, ocean, sea ice, surface fluxes and large-scale condensation
  are already in the tested configuration.
- `compare_arrays` uses `isapprox` on arrays, i.e. `norm(a - b) <= max(atol, rtol * max(norm(a), norm(b)))`.
  The `atol`s of `TOLERANCES` are in units of the prognostic variables (temperature 5.0, pressure 5.0,
  humidity 2e-3, geopotential 500), so for the tendencies (O(1e-5)) the temperature, pressure and
  humidity checks are vacuous.
- `compare_arrays` only compares top-level arrays of a group: the `ocean`, `land` (and other
  namespace) prognostic variables, the grid tendencies (`tendencies.grid`) and all of
  `parameterizations` (precipitation, radiation fluxes, cloud top) are never compared.

## Investigation

Environment: macOS arm64, Julia 1.11.9, Reactant CPU backend, `--project=SpeedyWeather/test/reactant`.
Scripts are in the session scratchpad (not kept); their essence is described here.

### 1. Stale test manifest (2026-11-24)

The local `SpeedyWeather/test/reactant/Manifest.toml` (gitignored, from 2026-10-16) pinned Reactant
0.2.294, violating the project's compat (0.2.305) and failing with `could not load symbol
"PjRtBufferGetReadyFuture"`. Regenerated (`rm Manifest.toml; Pkg.resolve(); Pkg.instantiate()`,
then `Pkg.Registry.update(); Pkg.update()`): Reactant 0.2.309, Reactant_jll 0.0.447, Enzyme 0.13.220,
CUDA 6.4.1 (Enzyme and CUDA held back by Reactant's compat).

### 2. The default model passes the existing test

`test_model(PrimitiveWetModel)` with `model_kwargs(::Type{PrimitiveWetModel}, _) = (;)`
(MatrixSpectralTransform, 10-day spin-up on Reactant, one-step tendencies, 100 steps): **24/24 pass**,
6m46s. Largest differences after 100 steps: grid temperature max abs 0.35 K (max rel 1.3e-3), grid
u/v max abs 0.2 m/s, grid humidity max abs 3.4e-4 kg/kg. Given the loose `atol`s (see Background)
this alone is weak evidence.

### 3. Error growth, Float32: CPU vs Reactant is indistinguishable from a 1-ulp perturbation

Default model, Float32, Reactant spun up 10 days, state copied to two CPU simulations, one of them
with the spectral temperature multiplied by `1 + eps(Float32)`. Norm-relative differences
`norm(a - b) / norm(a)` of the prognostic variables, and the number of grid points where convection
(`rain_rate_convection > 0`) or the cloud top differ:

| steps | pair | vorticity | divergence | temperature | humidity | pressure | conv. flips | cloud-top flips |
|---|---|---|---|---|---|---|---|---|
| 1 | CPU vs Reactant | 1.1e-7 | 6.1e-5 | 2.9e-8 | 6.5e-8 | 2.1e-9 | 0 | 0 |
| 1 | CPU vs CPU+1ulp | 8.5e-8 | 2.8e-5 | 1.3e-7 | 7.8e-8 | 1.1e-9 | 0 | 0 |
| 10 | CPU vs Reactant | 3.0e-5 | 3.3e-4 | 1.1e-6 | 8.5e-5 | 2.9e-8 | 2 | 1 |
| 10 | CPU vs CPU+1ulp | 1.2e-5 | 1.3e-4 | 2.6e-7 | 9.1e-7 | 6.9e-9 | 0 | 0 |
| 50 | CPU vs Reactant | 5.1e-4 | 3.5e-3 | 1.6e-5 | 1.0e-3 | 2.6e-7 | 7 | 7 |
| 50 | CPU vs CPU+1ulp | 3.3e-4 | 2.9e-3 | 7.9e-6 | 8.0e-4 | 2.0e-7 | 4 | 6 |
| 100 | CPU vs Reactant | 1.6e-3 | 9.2e-3 | 2.8e-5 | 1.2e-3 | 7.3e-7 | 20 | 20 |
| 100 | CPU vs CPU+1ulp | 9.6e-4 | 7.8e-3 | 2.4e-5 | 1.1e-3 | 6.8e-7 | 12 | 17 |

In Float32 the CPU-vs-Reactant differences follow the growth of a 1-ulp perturbation within a
factor of ~2 at every step. Convection switching on/off at a few points is the main amplifier
(threshold flips), which is the expected behaviour of the scheme, not a Reactant bug.

### 4. Error growth, Float64: a systematic difference in divergence

Same experiment with `NF = Float64` (10-day spin-up 453 s):

| steps | pair | vorticity | divergence | temperature | humidity | pressure |
|---|---|---|---|---|---|---|
| 1 | CPU vs Reactant | 2.3e-16 | **2.8e-8** | **1.5e-11** | 1.6e-16 | **4.3e-12** |
| 1 | CPU vs CPU+1ulp | 1.4e-16 | 5.0e-14 | 2.5e-16 | 1.5e-16 | 1.7e-18 |
| 5 | CPU vs Reactant | 3.1e-9 | 4.8e-8 | 5.3e-11 | 6.6e-11 | 2.3e-11 |
| 5 | CPU vs CPU+1ulp | 9.2e-15 | 1.6e-13 | 3.0e-16 | 7.1e-16 | 5.3e-18 |
| 100 | CPU vs Reactant | 5.5e-9 | 1.9e-8 | 5.9e-11 | 5.7e-10 | 2.3e-11 |
| 100 | CPU vs CPU+1ulp | 6.6e-14 | 6.3e-13 | 9.1e-16 | 4.0e-15 | 5.1e-17 |

The 1-ulp perturbation stays at 1e-16..1e-13 for 100 steps (no chaotic growth in Float64 on this
time scale, no convection flips), but CPU vs Reactant differs by 2.8e-8 in divergence **after one
step**, as a constant offset that does not grow. Divergence, temperature and pressure are affected,
vorticity and humidity are not (at step 1): the signature of the semi-implicit correction, which
only acts on divergence, temperature and surface pressure.

### 5. S⁻¹: Gauss-Jordan without pivoting is as accurate as LU (hypothesis rejected)

On Reactant S⁻¹ is computed every time step by `implicit_invert_S_kernel!` (Gauss-Jordan without
pivoting), on the CPU by `invert_S_lu!` (LU with partial pivoting). Both on the CPU for the operators
of a model after 3 steps, against a BigFloat inverse (relative ∞-norm error):

| NF | max cond(S) (l = 32) | worst error LU | worst error Gauss-Jordan | smallest pivot / max entry |
|---|---|---|---|---|
| Float64 | 20 | 7.2e-16 | 9.1e-16 | 0.18 |
| Float32 | 20 | 3.4e-7 | 2.9e-7 | 0.18 |

S is well conditioned and diagonally dominant enough, pivoting is not needed. The 2.8e-8 must come
from elsewhere.

### 6. The models differ before the first time step: Earth's orography

The time step ξ of the implicit solver is `Δt::NF` on both paths (no Float32 there). Instead, a
recursive comparison of every array and number in the CPU and the Reactant model right after
`initialize!(model)` (Float64, default PrimitiveWetModel, `output`, `feedback`, `callbacks`,
`spectral_grid`, `diagnostics` skipped) shows that only the orography differs:

| model field | norm-relative difference | entries that differ |
|---|---|---|
| `orography.orography` | 4.5e-9 | 3168/3168 |
| `orography.surface_geopotential` | 4.4e-9 | 465/560 |
| `spectral_transform.scratch_memory` | 3.7e-10 | (scratch, irrelevant) |

The orography is a model constant, `sync_variables!` copies only `Variables`, so CPU and Reactant
integrate over slightly different mountains. This explains the Float64 pattern of §4: the surface
geopotential enters the divergence tendency (∇²Φ) and the boundary layer, not the vorticity directly.

One time step from identical `Variables` (2-day CPU spin-up copied to Reactant, Float64): the largest
differences are in the boundary layer drag coefficients (2.6e-8, they depend on the height above the
surface), surface heat and humidity fluxes (1e-8), surface winds and tendencies (1e-9..4e-9), the
dynamics 1e-10. `model.implicit.S` (only the CPU's LU workspace) and `implicit.schedule.counter` (not
counted with Reactant) differ by design.

### 7. Which step of `initialize!(::EarthOrography)` differs

`initialize!` on Reactant: `get_asset` (the high-resolution orography file), `on_architecture`,
`interpolate!` (jitted, with the Reactant-only branchless `find_rings_kernel_branchless!` locator),
`transform!` to spectral, `spectral_smoothing!`, `transform!` back, `*gravity`, `spectral_truncation!`.
Each step on identical input, Float64: interpolation with a plain `@jit` bitwise identical, transforms
1.5e-15 and 9.6e-16, smoothing bitwise identical — but the full `initialize!` 4.5e-9. The difference:
`@maybe_jit` compiles with `reactant_compile_options()` (`raise = true`), the plain `@jit` doesn't.

The interpolation of the orography (`interpolate!` → `_interpolate_kernel!` → `anvil_average`) runs in
Float32: the file's high-resolution orography is a `Field{Float32}` (FullGaussianGrid, nlat_half = 512)
and `anvil_average` computes in the type of the data values, three `x + (y - x) * w` per point; only
the result is stored in the Float64 field. Raised to StableHLO, XLA evaluates these Float32 operations
differently (e.g. fused multiply-adds). Orography interpolated on the Reactant CPU backend vs the CPU:

| input | compilation | norm-relative | points that differ | max per point |
|---|---|---|---|---|
| Float32 (file) | plain `@jit` | 0 | 0/3168 | 0 |
| Float32 (file) | `reactant_compile_options()` (raised) | 1.2e-8 | 326/3168 | 1.9e-7 rel. (≤ 2 Float32 ulps), 0.24 mm |
| converted to Float64 | plain `@jit` | 0 | 0/3168 | 0 |
| converted to Float64 | `reactant_compile_options()` (raised) | 6.8e-22 | 1/3168 | 2.8e-17 m |

Not a Reactant bug: Float32 arithmetic compiled differently, the CPU result is not more accurate.
But the CPU and the Reactant model get different orographies, also in Float64.

**Fix** (`SpeedyWeather/src/dynamics/orography.jl`): convert the high-resolution orography to the
model's number format before interpolating (`Field(convert.(NF, field.data), field.grid)` on the CPU,
then `on_architecture`). Float32 models are unchanged (still interpolated in Float32), Float64 models
interpolate in Float64.

Same pattern, not fixed (not part of the default model's constants): `interpolate_ignoring_nan!` of
the Float32 climatologies for SST (`ocean.jl`), soil temperature (`land/soil_temperature.jl`) and
soil moisture (`land/soil_moisture.jl`) with `@maybe_jit`.

### 8. With the fix: Float64 CPU and Reactant agree to rounding

Same comparison as §6 (Float64, default model, 2-day CPU spin-up copied to Reactant), with the
orography fix:

- Model after `initialize!`: orography 1.9e-15, surface geopotential 1.6e-15 (transform rounding),
  nothing else differs (apart from transform scratch memory).
- One step from identical `Variables`: every array ≤ 1.1e-12 (spectral divergence tendency, which
  cancels heavily), all parameterizations ≤ 5e-14 (convective rain and snow, cloud cover, radiation,
  boundary layer, land and ocean tendencies), most 1e-15.
- 100 steps, norm-relative differences of the prognostic variables (§4 in brackets: before the fix;
  CPU vs CPU+1ulp in braces):

| steps | vorticity | divergence | temperature | humidity | pressure |
|---|---|---|---|---|---|
| 2 | 4.8e-16 | 1.7e-13 | 5.5e-17 | 1.3e-16 | 3.3e-18 |
| 10 | 1.5e-13 | 8.4e-13 | 9.1e-16 | 1.4e-15 | 2.8e-17 |
| 100 | 2.8e-13 [5.5e-9] {6.6e-14} | 2.2e-12 [1.9e-8] {6.3e-13} | 2.0e-15 [5.9e-11] {9.1e-16} | 6.2e-15 [5.7e-10] {4.0e-15} | 9.8e-17 [2.3e-11] {5.1e-17} |

CPU vs Reactant is now within a factor 2-4 of a 1-ulp perturbation in Float64 too: the default
PrimitiveWetModel on Reactant (CPU backend) reproduces the CPU to rounding, incl. convection and the
Byrne-O'Gorman longwave.

Regression test: `test/reactant/test_interpolation.jl`, "EarthOrography on Reactant in Float64",
orography and surface geopotential of a Float64 model, CPU vs Reactant, `rtol = 1e-12` (4.5e-9
before the fix).

### 9. Differentiation of the default model

The CI's `test_differentiability_time_steps(PrimitiveWetModel; compare_with_enzyme = false, options =
((;), (; checkpointing = Reactant.Binomial(2))))` with the default model (via `model_kwargs`): 5/5 pass,
13m12s locally (Reactant CPU backend, Julia 1.11). The same test with the previous configuration took
16m10s on CI (`main`, 2026-11-24, Linux runner), not comparable directly but no sign of a large
compile-time increase.

### 10. The climatologies, same fix

All file-based climatologies interpolated Float32 file data in Float32. Three of them with `@maybe_jit`
(raised with Reactant), so in Float64 models they differ between CPU and Reactant as the orography
did: `SeasonalOceanClimatology` (SST), `SeasonalLandTemperature`, `SeasonalSoilMoisture`, all via
`interpolate_ignoring_nan!`. SST file interpolated ignoring NaNs into Float64, CPU vs raised Reactant:
Float32 input 1.6e-8 norm-relative (229 of 31104 values differ), converted to Float64 first bitwise
identical. Two more without `@maybe_jit` (same on both architectures, but Float32 precision in Float64
models): `VegetationClimatology`, `AlbedoClimatology`. The land-sea mask already converted (`NF.(...)`)
and averages on the CPU.

Fix: all five convert the file data to the model's number format before interpolating
(`Field(convert.(NF, field.data), field.grid, field.dims)`, for the soil moisture the combined
two-layer field is allocated with the model's `NF`).

Test: `test/reactant/test_interpolation.jl`, "Climatologies on Reactant in Float64": the three monthly
climatologies of a Float64 model (land-sea mask initialized for the masking), CPU vs Reactant,
`rtol = 1e-12`, 3/3 pass. CPU unit tests of the components (`test/parameterizations/land.jl`,
`ocean_sea_ice.jl`, `albedo.jl`, `ocean_land_time_stepping.jl`, Julia 1.12, `--check-bounds=yes`) pass.

### 11. Tighter tolerances for the correctness tests

`compare_arrays` also returns the norm-relative difference `norm(cpu - reactant) / max(norm(cpu),
norm(reactant))` (what `isapprox` compares against `rtol`), printed for every variable
(`print_comparison`). Measured for all six `test_model` configurations of `test_correctness.jl` (default
models, all orography and climatology fixes, macOS arm64, Julia 1.11, Reactant 0.2.309; Baro = Barotropic,
SW = ShallowWater, Dry/Wet = PrimitiveDry/Wet, M = MatrixSpectralTransform, S = SpectralTransform):

| variable | Baro M | Baro S | SW M | Dry M | Wet M | Wet S |
|---|---|---|---|---|---|---|
| tend. vorticity | 2.9e-6 | 1.8e-6 | 5.6e-6 | 5.7e-6 | 7.1e-6 | 3.5e-6 |
| tend. divergence | | | 5.8e-5 | 5.6e-4 | 1.8e-3 | 1.1e-3 |
| tend. η / pressure | | | 3.2e-5 | 3.3e-5 | 1.4e-4 | 9.6e-5 |
| tend. temperature | | | | 1.1e-4 | 2.4e-4 | 1.5e-4 |
| tend. humidity | | | | | 8.1e-6 | 7.6e-6 |
| progn. vorticity | 2.4e-6 | 8.2e-7 | 3.0e-6 | 7.3e-5 | 1.5e-3 | 1.4e-3 |
| progn. divergence | | | 8.5e-5 | 6.6e-4 | 9.1e-3 | 9.6e-3 |
| progn. η / pressure | | | 2.6e-6 | 6.7e-8 | 7.6e-7 | 6.4e-7 |
| progn. temperature | | | | 1.1e-6 | 2.9e-5 | 2.6e-5 |
| progn. humidity | | | | | 1.6e-3 | 1.4e-3 |
| grid vorticity | 3.0e-6 | 1.0e-6 | 3.4e-6 | 7.8e-5 | 1.6e-3 | 1.5e-3 |
| grid divergence | | | 8.8e-5 | 6.0e-4 | 9.4e-3 | 9.9e-3 |
| grid u | 1.9e-6 | 9.4e-7 | 2.5e-6 | 4.1e-5 | 6.9e-4 | 7.5e-4 |
| grid v | 2.0e-6 | 8.0e-7 | 5.2e-6 | 1.7e-4 | 1.9e-3 | 1.9e-3 |
| grid η / pressure | | | 3.1e-6 | 2.7e-7 | 1.2e-6 | 9.2e-7 |
| grid temperature | | | | 1.4e-6 | 4.1e-5 | 3.7e-5 |
| grid humidity | | | | | 2.2e-3 | 2.0e-3 |

The one-layer models differ by ~1e-6, the PrimitiveWetModel by up to 1e-2 after 100 steps (§3:
amplified rounding, convection). One set of tolerances for all models would be set by the
PrimitiveWetModel, so two sets (`correctness_tolerances(ModelType)` in `runtests.jl`), each for the
tendencies (one step) and the time stepping (100 steps), norm-relative without `atol`, ~10x
(9x-17x) the largest measured difference:

| | one-layer: tendencies | one-layer: 100 steps | primitive: tendencies | primitive: 100 steps |
|---|---|---|---|---|
| default | 5e-5 | 5e-5 | 1e-4 | 2e-2 (vorticity, v, humidity) |
| divergence | 5e-4 | 1e-3 | 2e-2 | 1e-1 |
| η | 5e-4 | (default) | | |
| temperature | | | 3e-3 | 5e-4 |
| pressure | | | 2e-3 | 2e-5 |
| u | | | | 1e-2 |

Previously one set for everything, with `atol`s in units of the variables: e.g. temperature
`rtol = 2e-2, atol = 5`, pressure `1e-2, 5`, humidity `5e-2, 2e-3`, vorticity `5e-2`, u/v `5e-2`. Now the
tendency checks can fail, the one-layer models are checked ~1000x tighter, the primitive models'
temperature 40x, pressure 500x, u 5x, vorticity/v/humidity 2.5x tighter (divergence unchanged).
`RTOL`/`ATOL` remain for `test_maybe_jit.jl`.

Verification with the new tolerances, all six configurations (`test_model` as in `test_correctness.jl`):
Barotropic 8/8 (both transforms), ShallowWater 14/14, PrimitiveDry 17/17, PrimitiveWet 24/24 (both
transforms). The primitive models reproduce the measured differences bitwise (deterministic), the
Barotropic model's vary slightly between runs (random initial conditions, 2.2e-6 vs 2.4e-6 vorticity).

### 12. The divergence deviations

One time step of the default PrimitiveWetModel (Float32, MatrixSpectralTransform) split into its stages,
each stage from identical inputs (CPU directly, Reactant compiled with `reactant_compile_options()` per
stage, all variables compared, then Reactant resynced to the CPU; 10-day CPU spin-up). Largest
norm-relative difference per stage, in units of eps(Float32) = 1.2e-7:

| stage | largest difference |
|---|---|
| `reinitialize!` (implicit S⁻¹) | 0.6 eps |
| parameterizations (column kernel) | rain_rate_convection 212, tendencies.grid.temperature 153, cloud_cover 120, outgoing_shortwave 94, tendencies.grid.humidity 93, surface_shortwave_down 46, ocean albedo 45 eps |
| ocean, sea ice, land | ≤ 0.7 eps |
| pressure_gradient_flux!, vertical_integration!, vertical_velocity! | ≤ 0.7 eps |
| vertical_advection! | 6.3 eps (grid temperature tendency) |
| grid_tendencies! | ≤ 0.4 eps |
| transform grid → spectral | ≤ 4.1 eps |
| vordiv_spectral_tendencies! | vorticity 4.2, divergence 0.6 eps |
| temperature/humidity spectral tendency | 6.2 / 0.9 eps |
| bernoulli_spectral_potential! | divergence 3.1 eps |
| horizontal_diffusion!, implicit_correction! | ≤ 0.9 eps |
| update_prognostic!, transform spectral → grid | ≤ 1.3 eps (river runoff rate 38 eps) |

The dynamical core, the divergence included, agrees to a few eps per stage; the parameterizations are
the outlier by ~30x. Yet the test's one-step divergence tendency differs by 1.8e-3 (~15000 eps).

Cumulative over the stages of one step (no resync), the spectral tendencies after each stage:
divergence 4.9 eps after `vordiv_spectral_tendencies!`, then **1.46e-3 (12214 eps) after
`bernoulli_spectral_potential!`**, 1.83e-3 after diffusion, 1.64e-3 after the implicit correction, which
couples it into temperature (9.3e-6 → 1.85e-4) and pressure (1.4e-6 → 1.25e-4). The full compiled
`time_step!` from the same state gives the same (divergence 1.68e-3, temperature 2.3e-4, pressure
1.4e-4), so compiling stages separately doesn't matter.

**Culprit**: `bernoulli_spectral_potential!` subtracts ∇² of the Bernoulli potential
`½(u²+v²) + Φ (+ R Tₖ ln pₛ)`, which since #1381 ("Compute geopotential once in grid space", 2026-11-08)
contains the full grid-point geopotential Φ: up to ~9e5 m²/s² (#1381's plan), dominated by the layer
height and the surface geopotential, while only its horizontal variations (~1e3-1e4) matter for ∇².
In Float32 the grid → spectral transform of this large field has rounding errors of a few eps of its
norm (3.6e-7 of 1.3e6 for `dynamics.kinetic_energy`, CPU vs Reactant), spread over all wavenumbers;
∇² multiplies them by l(l+1) (up to ~1000 at T32), and the divergence tendency is ~17x smaller than the
Bernoulli potential. On identical inputs the stage itself is exact to 3 eps, but its input differs by
the transform's rounding. Not a Reactant bug: the CPU's Float32 divergence tendency has errors of the
same size, CPU and Reactant only round differently. Before #1381 the geopotential for ∇² was computed in
spectral space from the spectral temperature (each wavenumber with its own, small rounding error).

Experiment E1: transform only the horizontally varying part, `½(u²+v²) + (Φ - g·h - mean_k(Φ - g·h))`
on the grid (g·h the surface geopotential, mean over the grid points of each layer), and add the
spectral surface geopotential `orography.surface_geopotential` (= g × transform of the grid
orography, for all orographies and `set!`) to the spectral Bernoulli potential with R Tₖ ln pₛ in
`linear_pressure_gradient!`. Mathematically the same (∇² of a per-layer constant is 0 on the
Gaussian grids, Φₛ spectrally as before #1381). PrimitiveWet, MatrixSpectralTransform, CPU vs
Reactant (norm-relative):

| | before | E1 | factor |
|---|---|---|---|
| tendency divergence | 1.81e-3 | 1.72e-4 | 10.5 |
| tendency temperature | 2.37e-4 | 2.80e-5 | 8.4 |
| tendency pressure | 1.42e-4 | 4.06e-5 | 3.5 |
| 100 steps: divergence, vorticity | 9.4e-3, 1.6e-3 | 9.9e-3, 1.6e-3 | ~1 |

The 100-step differences don't change: a 1-ulp perturbation of the CPU model grows to the same size
(§3), they are the model's sensitivity (convection switching on/off), not the per-step rounding.

Accuracy on the CPU alone (not just CPU vs Reactant): ∇²B in Float32 vs a Float64 evaluation of the same
Float32 grid values (Φₛ from its spectral definition), 10-day spin-up, T32L8, relative error:

| | SpectralTransform | MatrixSpectralTransform |
|---|---|---|
| current: ½(u²+v²) + Φ | 1.1e-5 (90 eps) | 2.7e-5 (224 eps) |
| E1: (Φ - g·h) - mean | 2.7e-6 (22 eps) | 2.7e-6 (23 eps) |
| V2: (Φ - mean Φ) - g·h | 3.6e-7 (3 eps) | 6.2e-7 (5 eps) |
| anomaly rounded only once (floor) | 2.9e-7 (2 eps) | 5.8e-7 (5 eps) |

rms of the transformed grid field 1.3e5 → 3.5e3. E1 rounds Φ - g·h at the size of Φ (~2e5 at the
upper layers), V2 subtracts the layer mean first, exact there (Sterbenz lemma: Φ and its mean within a
factor 2), then g·h: at the floor, 30-43x more accurate than the current version.

**Fix** (V2, `src/dynamics/tendencies.jl`): `bernoulli_grid_potential!(vars, ::PrimitiveEquation, ...)`
writes `½(u²+v²) + (Φ - Φ̄ₖ) - g·orography` with Φ̄ₖ = `sum(Φ, dims = 1) / npoints` (a 112-byte allocation on
the CPU, the CPU time step allocates ~2.5 MB anyway; GPU graphs only capture the transforms), the new
`bernoulli_surface_geopotential!` adds `orography.surface_geopotential` to the spectral Bernoulli potential
after `linear_pressure_gradient!` in `spectral_tendencies!`. The `AbstractSpectralTransform` method (for the
generic `bernoulli_potential!(vars, S, TS)`) keeps the full Φ. Tests (`test/dynamics/geopotential.jl`):
the grid formula, the same ∇² as ½(u²+v²) + Φ with the spectral Φₛ (Float64 rtol 1e-10, Float32 1e-4), and
"Bernoulli potential: Float32 accuracy of its ∇²" (initial conditions: 3.3e-7 vs 1.1e-5 with the full Φ,
tested < 1e-6 and < 1/5 of the full-Φ error). All `geopotential.jl` tests pass (Julia 1.12).

Correctness tests with the fix, all six configurations pass with the previous tolerances; norm-relative
CPU-vs-Reactant differences before → after (largest over the configurations where it matters):

| | before | after | factor |
|---|---|---|---|
| tendencies: divergence (Dry, Wet M, Wet S) | 5.6e-4, 1.8e-3, 1.1e-3 | 2.8e-5, 8.8e-5, 4.2e-5 | 20-25 |
| tendencies: temperature | 1.1e-4, 2.4e-4, 1.5e-4 | 6.5e-6, 1.6e-5, 9.2e-6 | 15-17 |
| tendencies: pressure | 3.3e-5, 1.4e-4, 9.6e-5 | 1.0e-5, 3.9e-5, 1.5e-5 | 3-6.5 |
| PrimitiveDry 100 steps: divergence, vorticity, u, v, temperature | 6.6e-4, 7.8e-5, 4.1e-5, 1.7e-4, 1.4e-6 | 1.6e-4, 1.7e-5, 1.0e-5, 3.1e-5, 5.9e-7 | 2.4-5.5 |
| PrimitiveWet 100 steps (both transforms) | e.g. divergence 9.9e-3, v 1.9e-3 | 1.4e-2, 2.8e-3 | 0.7-1.2 |
| Barotropic, ShallowWater | | | 1 (unchanged code path) |

The dry model's 100-step differences were dominated by the Bernoulli noise too (no convection), the wet
model's remain those of its sensitivity to any perturbation.

New tolerances (`runtests.jl`), ~10x the largest measured difference (6.5x-20x):

| | before | after |
|---|---|---|
| primitive tendencies: divergence, temperature, pressure | 2e-2, 3e-3, 2e-3 | 1e-3, 2e-4, 5e-4 |
| primitive tendencies: default (vorticity, humidity) | 1e-4 | 1e-4 |
| PrimitiveDry 100 steps (own set, was shared with Wet) | 2e-2 (vorticity, v, humidity), u 1e-2, divergence 1e-1, temperature 5e-4, pressure 2e-5 | default (vorticity, u) 2e-4, v 5e-4, divergence 2e-3, temperature 1e-5, pressure 3e-6 |
| PrimitiveWet 100 steps | unchanged | unchanged |
| one-layer models | unchanged | unchanged |

Verification with these tolerances, all six configurations: Barotropic 8/8 (both transforms),
ShallowWater 14/14, PrimitiveDry 17/17, PrimitiveWet 24/24 (both transforms); the primitive models
reproduce the measured differences bitwise.

Further checks of the Bernoulli change: CPU tests (Julia 1.12, `--check-bounds=yes`) `dynamics/implicit.jl`,
`dynamics/vertical_advection.jl`, `dynamics/time_stepping.jl` (incl. bit reproducibility, single steps vs
`run!`), `dynamics/matrix_transform.jl`, `variables/fuse.jl` (128/128) pass; Enzyme reverse mode through
`bernoulli_grid_potential!` (reactant env, Julia 1.11): finite non-zero gradients, the gradient with respect
to Φ sums to zero on every layer (layer mean subtracted), 3/3. Not run: the differentiability test suite
(`test/differentiability`, separate CI), GPU tests.

## Summary of changes

- `SpeedyWeather/test/reactant/setup.jl`: `model_kwargs` is `(;)` for the PrimitiveWetModel, all
  Reactant tests use the default PrimitiveWetModel (the PrimitiveDryModel keeps `convection = nothing`).
- `SpeedyWeather/test/reactant/Project.toml`: no `[compat]` for Reactant, resolves to the latest
  version the packages allow.
- `SpeedyWeather/src/dynamics/orography.jl`: the high-resolution orography (Float32 file) is converted
  to the model's number format before it is interpolated, so that Float64 models interpolate in Float64
  and CPU and Reactant get the same orography.
- `SpeedyWeather/test/reactant/test_interpolation.jl`: test "EarthOrography on Reactant in Float64".
- `SpeedyWeather/test/reactant/runtests.jl`, `test_correctness.jl`: norm-relative tolerances per model
  (one-layer, PrimitiveDry, PrimitiveWet) for the tendencies and the time stepping, the norm-relative
  difference is printed for every variable (§11, §12).
- `SpeedyWeather/src/dynamics/tendencies.jl`: the grid Bernoulli potential of the primitive equation
  models without the layer mean of the geopotential and without the surface geopotential (added in
  spectral space by the new `bernoulli_surface_geopotential!`): 30-43x smaller Float32 rounding errors in
  ∇² of the Bernoulli potential, 20-25x smaller CPU-vs-Reactant differences in the divergence tendency
  (§12); tests in `test/dynamics/geopotential.jl`.
- `ocean.jl`, `land/soil_temperature.jl`, `land/soil_moisture.jl`, `land/vegetation.jl`,
  `radiation/albedo/general.jl`: the climatologies too are interpolated in the model's number format;
  test "Climatologies on Reactant in Float64".

## Testing and verification

- Default PrimitiveWetModel, `test_model(PrimitiveWetModel)` of `test_correctness.jl`: 24/24 (§2).
- Float32: CPU vs Reactant differences grow as a 1-ulp perturbation (§3).
- Float64: CPU vs Reactant within a factor 2-4 of a 1-ulp perturbation over 100 steps, all variables
  after one step ≤ 1.1e-12 (§8).
- `test_interpolation.jl` (incl. the new test) passes, CPU `test/dynamics/orography.jl` 4/4 (Julia 1.12),
  default-model differentiation 5/5 (§9).
- All six `test_model` configurations with the new tolerances pass (§11).
- Not run locally as a whole: the full Reactant test suite (`runtests.jl`), e.g. `test_output.jl`,
  `test_implicit.jl`, `test_time_steps.jl` with the default model.

## Known limitations

- The correctness tests compare only the top-level arrays of `prognostic`, `grid` and `tendencies`:
  nested namespaces (`ocean`, `land`), `tendencies.grid` and `parameterizations` are not compared.
- The tolerances are ~10x the differences measured on one machine (macOS arm64); CI (Linux x86, Julia
  1.12) may round differently, the 100-step differences of the PrimitiveWetModel depend on how many
  grid points switch convection on or off.
- Float32 comparisons over 100 steps are dominated by the amplification of rounding differences
  (mostly convection switching on/off at a few points), so they can only detect errors well above
  ~1e-3 (norm-relative) in vorticity/divergence/humidity.

## Future work

- The PrimitiveWetModel's 100-step comparison is limited by its sensitivity: any perturbation, also 1 ulp,
  grows to ~1e-3-1e-2 in 100 steps (§3). Tighter tolerances there need a shorter comparison (in §3 after
  20 steps vorticity 5e-5, divergence 4e-4, humidity 3e-6, i.e. 30-500x smaller than after 100) or a
  Float64 comparison (§8: ≤ 2e-12 after 100 steps).
- The parameterizations are now the largest per-stage CPU-vs-Reactant difference (~100-200 eps:
  convective rain, temperature tendency, cloud cover, shortwave; §12), most likely the same cancellation
  in the convection's T - T_ref of nearly equal temperatures.

- A Float64 CPU-vs-Reactant comparison of the default PrimitiveWetModel in CI (CPU spin-up copied to
  Reactant, one compiled time step, all variables incl. nested namespaces and parameterizations,
  tolerance ~1e-10): sensitive to errors that Float32 rounding hides, as the orography here.
- Comparing nested namespaces (`ocean`, `land`), `tendencies.grid` and `parameterizations`.
