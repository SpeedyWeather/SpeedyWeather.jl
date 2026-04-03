# CLAUDE.md — SpeedyWeather.jl Codebase Summary

## Overview

SpeedyWeather.jl is a global atmospheric general circulation model (GCM) written in Julia.
It is organised as a **Julia monorepo** with 5 interdependent packages. The design philosophy
prioritises flexibility and interactivity over production forecasting.

## Package Structure

```
SpeedyWeather.jl/
├── LowerTriangularArrays/    # Efficient storage for spherical harmonic coefficients
├── RingGrids/                 # Iso-latitude ring grids, fields, interpolation
├── SpeedyTransforms/          # Spherical harmonic transforms (FFT + Legendre)
├── SpeedyWeatherInternals/    # Device abstraction (CPU/GPU), parameters, utilities
└── SpeedyWeather/             # Main atmospheric model (depends on all above)
```

Each package has its own `Project.toml`, `src/`, and `test/` directories.
Run tests with `julia --project=<PackageName> <PackageName>/test/runtests.jl`.

## SpeedyWeather Package Layout

```
SpeedyWeather/src/
├── dynamics/          # Core dynamics (tendencies, time stepping, diffusion, …)
├── parameterizations/ # Parameterisations (convection, radiation, land, ocean, …)
│   ├── land/          # 2-layer soil bucket model
│   ├── radiation/     # Shortwave & longwave
│   └── surface_fluxes/
├── models/            # Model definitions and simulation runner
├── variables/         # Unified Variables system (prognostic, grid, tendencies, …)
├── output/            # NetCDF, JLD2, callbacks, restart files
└── input/             # Asset loading
```

## Key Types

| Type | Purpose |
|------|---------|
| `SpectralGrid` | Central config: resolution (`trunc`), grid type, `nlayers`, `NF`, architecture |
| `Variables` | Unified container for all simulation variables (prognostic, grid, tendencies, dynamics, parameterizations, particles, scratch) |
| `Simulation{V,M}` | Top-level container holding model + variables |
| `Leapfrog` | Time stepper with Robert/Williams filters |
| `SpectralTransform` | FFT/Legendre transform machinery |

### Model Type Hierarchy

```
AbstractModel
├── Barotropic           # Single-layer vorticity only
├── ShallowWater         # Single-layer with divergence + pressure
└── PrimitiveEquation
    ├── PrimitiveDry     # Multi-layer without humidity
    └── PrimitiveWet     # Multi-layer with full physics
```

All model components (orography, ocean, convection, etc.) inherit from
`AbstractModelComponent` and implement `initialize!`, optionally `timestep!`
and `finalize!`.

### Variables Structure

`Variables` is a unified container with 7 groups, each a NamedTuple:

| Group | Purpose |
|-------|---------|
| `prognostic` | State variables subject to time stepping (vorticity, divergence, temperature, pressure, tracers, etc.) |
| `grid` | Grid-point copies of spectral prognostic variables (u, v, temperature, etc.) |
| `tendencies` | Time derivatives; spectral at top level, grid-space under `.grid` namespace |
| `dynamics` | Working arrays for the dynamical core (pressure gradients, vertical velocity, etc.) |
| `parameterizations` | Working arrays for parameterizations (fluxes, radiation, cloud properties, etc.) |
| `particles` | Particle advection variables |
| `scratch` | Temporary storage (write-before-read) |

Variables can be organized by namespace (e.g. `:ocean`, `:land`, `:tracers`):

```julia
vars.prognostic.vorticity                             # Spectral vorticity
vars.prognostic.ocean.sea_surface_temperature   # Ocean namespace
vars.tendencies.vorticity                             # Spectral vorticity tendency
vars.tendencies.grid.u                          # Grid-space u-wind tendency
vars.dynamics.w                                 # Vertical velocity
```

Model components define their variables via variable type categories:
`PrognosticVariable`, `GridVariable`, `TendencyVariable`, `DynamicsVariable`,
`ParameterizationVariable`, `ParticleVariable`, `ScratchVariable`.
These are collected and allocated automatically by the `Variables(model)` constructor.

## Typical Usage

```julia
arch = SpeedyWeather.CPU()
spectral_grid = SpectralGrid(trunc=31, nlayers=8, architecture=arch)
model = PrimitiveWetModel(spectral_grid)
simulation = initialize!(model)
run!(simulation, period=Day(10), output=true)
```

## Running on GPU

To run a model on GPU, load a GPU backend package (e.g. `CUDA`, `AMDGPU`, or `Metal`)
before using SpeedyWeather, then pass `GPU()` as the architecture:

```julia
using CUDA              # or AMDGPU, Metal
using SpeedyWeather

arch = SpeedyWeather.GPU()
spectral_grid = SpectralGrid(trunc=31, nlayers=8, architecture=arch)
model = PrimitiveWetModel(spectral_grid)
simulation = initialize!(model)
run!(simulation, period=Day(10))
```

The GPU backend must be loaded first so that SpeedyWeather's extension modules
register the correct array types (`CuArray`, `ROCArray`, `MtlArray`).

GPU tests live in `SpeedyWeather/test/GPU/` and can be run with:

```bash
julia --project=SpeedyWeather SpeedyWeather/test/GPU/runtests.jl
```

## Time-Step Information Flow

```
Spectral State (LowerTriangularArray)
    ↓  transform!()
Grid Variables (RingGrids Field)
    ↓  dynamics_tendencies!() + parameterization_tendencies!()
Grid Tendencies
    ↓  transform!() back to spectral
Spectral Tendencies
    ↓  timestep!() — Leapfrog + implicit correction
New Spectral State
```

Prognostic variables live in `vars.prognostic`, tendencies in `vars.tendencies` (spectral)
and `vars.tendencies.grid` (grid-space). Vorticity and divergence are scaled by `radius`
inside `run!` for spectral efficiency.

## Array Types 

Mainly two array types are used: `LowerTriangularArray` for spectral coefficients, and `Field` for gridded data. They can be initialized for testing as in the following: 

```julia 
arch = SpeedyWeather.CPU()
spectrum = Spectrum(trunc=10, architecture = arch)
nlayers = 5
coeffs_zero = zeros(ComplexF32, spectrum, nlayers)
coeffs_rand = rand(ComplexF32, spectrum, nlayers)

grid = HEALPixGrid(6, arch)
field_zero = zeros(Float32, grid, nlayers)
field_rand = rand(Float32, grid, nlayers)
```

In case there's already a `SpectralGrid` use it instead: 

```julia
arch = SpeedyWeather.CPU()
spectral_grid = SpectralGrid(trunc=10, architecture=arch)
coeffs = rand(ComplexF32, spectral_grid.spectrum)
field = rand(Float32, spectral_grid.grid)
```

## Key Source Files

| File | Notes |
|------|-------|
| `dynamics/tendencies.jl` (~1200 lines) | All dynamical tendency calculations |
| `dynamics/time_integration.jl` (~450 lines) | Main time loop (`time_stepping!`) |
| `dynamics/implicit.jl` | Semi-implicit gravity-wave treatment (allows CFL > 1) |
| `models/simulation.jl` | `run!`, `initialize!`, `finalize!` for `Simulation` |
| `parameterizations/convection.jl` | Betts-Miller convection |
| `parameterizations/large_scale_condensation.jl` | Implicit condensation |
| `output/netcdf_output.jl` (~540 lines) | NetCDF writing with interpolation |
| `output/callbacks.jl` | Hook system for runtime code injection |

## Numerical Methods

- **Spectral transform**: FFT in longitude + Associated Legendre polynomials in latitude
- **Time integration**: Leapfrog with Robert (1966) and Williams (Amezcua 2011) filters
- **Semi-implicit**: Gravity waves treated implicitly; allows larger timesteps
- **Grids**: Full and reduced (octahedral, HEALPix) Gaussian and Clenshaw-Curtis grids
- **Precision**: Float32 (default), Float64, experimental BFloat16 + stochastic rounding
- **GPU**: Device-agnostic via `KernelAbstractions` (CUDA, Metal, AMDGPU)

## Extending the Model

- **Custom component**: Inherit abstract type → implement `initialize!(c, model)` → pass as kwarg to model constructor
- **Custom parameterisation**: Use `CustomParameterization` slot in `PrimitiveWet/Dry`
- **Callbacks**: Implement `AbstractCallback`, attach via `model.callbacks`
- **Output variables**: Add/remove via `NetCDFOutput` variable list

## Testing

- Every new method that is implemented has to be tested in the unit tests of its respective submodule
- Keep the unit tests short and concise, they should finish quickly
- Call the test with the `--check-bounds=yes` flag activated

```bash
# Main model tests
julia --project=SpeedyWeather --check-bounds=yes -e 'using Pkg; Pkg.test("SpeedyWeather")'

# Extended tests
julia --project=SpeedyWeather SpeedyWeather/test/runtests.jl extended_tests

# Individual packages
julia --project=RingGrids --check-bounds=yes -e 'using Pkg; Pkg.test("RingGrids")'
```

Test subdirectories: `dynamics/`, `parameterizations/`, `output/`, `GPU/`, `differentiability/`.

## Documentation

- Stable docs: https://speedyweather.github.io/SpeedyWeatherDocumentation/stable/
- Dev docs:    https://speedyweather.github.io/SpeedyWeatherDocumentation/dev/
- Build: `julia --project=docs docs/make.jl`

## Kernel Launching

Kernel infrastructure lives in `SpeedyWeatherInternals/src/KernelLaunching/KernelLaunching.jl`
(a standalone submodule of `SpeedyWeatherInternals`) and `SpeedyWeatherInternals/src/Architectures/`.
Other packages import it via `using SpeedyWeatherInternals.KernelLaunching`.

### Architecture abstraction

```julia
CPU()        # wraps KernelAbstractions.CPU()
CPUStatic()  # wraps KernelAbstractions.CPU(; static = true)
GPU()        # wraps backend-specific GPU (CUDA default; Metal/AMDGPU via extensions)
```

Backend-specific array types (`CuArray`, `ROCArray`, `MtlArray`) are registered
in the extension modules under `SpeedyWeatherInternals/ext/`.
Use `on_architecture(device, data)` to transfer data between devices.

### Work order types

Six types encode what iteration space a kernel covers:

| Type | Iteration space |
|------|----------------|
| `SpectralWorkOrder` | all `lm` harmonics × vertical layers |
| `SpectralInnerWorkOrder` | as above, but skip `lm=1` |
| `DiagonalWorkOrder` | diagonal elements of a `LowerTriangularArray` |
| `RingGridWorkOrder` | all grid points `ij` × vertical layers |
| `Array3DWorkOrder` | regular 3D arrays |
| `LinearWorkOrder` | flattened 1D (eachindex) |

### Launching a kernel

```julia
launch!(arch, WorkOrderType, worksize, kernel!, args...)
```

Internally this calls `work_layout` to compute a sensible `workgroup` size
(capped at 256 threads, heuristic aspect ratios for 2D/3D), then dispatches
via the KernelAbstractions API:

```julia
loop! = kernel!(device(arch), workgroup, worksize)
loop!(args...)
```

### Defining a kernel

```julia
@kernel inbounds = true function my_kernel!(A, B, c)
    I = @index(Global, Linear)      # or NTuple / Cartesian for multi-dim
    A[I] = B[I] * c
end
```

- `inbounds = true` — standard on all performance kernels (28+ files)
- Use `@index(Global, Linear)` for 1D / flattened access,
  `@index(Global, NTuple)` for `ij, k` style,
  `@index(Global, Cartesian)` for `I[1], I[2]` style

### CPU vs GPU dispatch

Many functions branch on architecture type:

```julia
arch = architecture(field)
if arch isa GPU
    launch!(arch, LinearWorkOrder, (n,), my_kernel!, ...)
else
    my_cpu_loop!(...)   # plain Julia loops
end
```

Kernels are spread across ~36 files in `SpeedyWeather/src/dynamics/`,
`parameterizations/`, `SpeedyTransforms/src/`, `RingGrids/src/`, and
`LowerTriangularArrays/src/`.

## Benchmarks

Benchmarks live in `SpeedyWeather/benchmark/`. Run the full suite with:

```bash
julia --project=SpeedyWeather/benchmark SpeedyWeather/benchmark/manual_benchmarking.jl
```

Results are written to `SpeedyWeather/benchmark/README.md`. Previous benchmark results
in that file serve as the baseline. When running benchmarks, compare new results against
the existing README.md values: deviations of +/- 20% are normal and acceptable, but
larger regressions should be reported to the user.

## Code Style

The project uses [Runic.jl](https://github.com/fredrikekre/Runic.jl) for formatting.

## Pull Request Convention

Every PR **must** add a line to `CHANGELOG.md` under the `## Unreleased` section:

```
- Brief description of change [#NNN](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/NNN)
```

Add it as the first bullet under `## Unreleased`.
