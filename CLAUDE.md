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
├── physics/           # Parameterisations (convection, radiation, land, ocean, …)
│   ├── land/          # 2-layer soil bucket model
│   ├── radiation/     # Shortwave & longwave
│   └── surface_fluxes/
├── models/            # Model definitions and simulation runner
├── variables/         # Prognostic & diagnostic variable metadata
├── output/            # NetCDF, JLD2, callbacks, restart files
└── input/             # Asset loading
```

## Key Types

| Type | Purpose |
|------|---------|
| `SpectralGrid` | Central config: resolution (`trunc`), grid type, `nlayers`, `NF`, architecture |
| `PrognosticVariables` | All state variables (vor, div, temp, humid, pres, particles, tracers) |
| `DiagnosticVariables` | Tendencies and auxiliary grid-space fields |
| `Simulation{Model,Progn,Diagn}` | Top-level container holding model + variables |
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

## Typical Usage

```julia
spectral_grid = SpectralGrid(trunc=31, nlayers=8)
model = PrimitiveWetModel(spectral_grid)
simulation = initialize!(model)
run!(simulation, period=Day(10), output=true)
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

Prognostic variables are stored with 2 time levels (`lf=1,2`) for the leapfrog scheme.
Vorticity and divergence are scaled by `radius` inside `run!` for spectral efficiency.

## Key Source Files

| File | Notes |
|------|-------|
| `dynamics/tendencies.jl` (~1200 lines) | All dynamical tendency calculations |
| `dynamics/time_integration.jl` (~450 lines) | Main time loop (`time_stepping!`) |
| `dynamics/implicit.jl` | Semi-implicit gravity-wave treatment (allows CFL > 1) |
| `models/simulation.jl` | `run!`, `initialize!`, `finalize!` for `Simulation` |
| `physics/convection.jl` | Betts-Miller convection |
| `physics/large_scale_condensation.jl` | Implicit condensation |
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

```bash
# Main model tests
julia --project=SpeedyWeather SpeedyWeather/test/runtests.jl

# Extended tests
julia --project=SpeedyWeather SpeedyWeather/test/runtests.jl extended_tests

# Individual packages
julia --project=RingGrids RingGrids/test/runtests.jl
```

Test subdirectories: `dynamics/`, `physics/`, `output/`, `GPU/`, `differentiability/`.

## Documentation

- Stable docs: https://speedyweather.github.io/SpeedyWeatherDocumentation/stable/
- Dev docs:    https://speedyweather.github.io/SpeedyWeatherDocumentation/dev/
- Build: `julia --project=docs docs/make.jl`

## Kernel Launching

Kernel infrastructure lives in `SpeedyWeatherInternals/src/Utils/kernel_launching.jl`
and `SpeedyWeatherInternals/src/Architectures/`.

### Architecture abstraction

```julia
CPU()        # wraps KernelAbstractions.CPU()
CPUStatic()  # wraps KernelAbstractions.CPU(; static = true)
GPU()        # wraps backend-specific GPU (CUDA default; Metal/AMDGPU via extensions)
```

Backend-specific array types (`CuArray`, `ROCArray`, `MtlArray`) and
`convert_to_device` are registered in the extension modules under
`SpeedyWeatherInternals/ext/`.

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
@kernel inbounds = true function my_kernel!(A, B, @Const(c))
    I = @index(Global, Linear)      # or NTuple / Cartesian for multi-dim
    A[I] = B[I] * c
end
```

- `inbounds = true` — standard on all performance kernels (28+ files)
- `@Const(x)` — marks read-only arguments for compiler optimisation (23+ files)
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
`physics/`, `SpeedyTransforms/src/`, `RingGrids/src/`, and
`LowerTriangularArrays/src/`.

## Code Style

The project uses [Runic.jl](https://github.com/fredrikekre/Runic.jl) for formatting.

## Pull Request Convention

Every PR **must** add a line to `CHANGELOG.md` under the `## Unreleased` section:

```
- Brief description of change [#NNN](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/NNN)
```

Add it as the first bullet under `## Unreleased`.
