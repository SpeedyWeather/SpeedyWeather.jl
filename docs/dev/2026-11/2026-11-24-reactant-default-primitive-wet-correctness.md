# Reactant correctness of the default PrimitiveWetModel

> Status: **in progress**. Investigation log; the default model already passes the existing
> correctness test, but a Float64 comparison shows a systematic ~1e-8 difference in divergence that is
> being tracked down.

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
