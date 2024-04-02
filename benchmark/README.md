# Benchmarks

created for SpeedyWeather.jl v0.9.0 on Tue, 02 Apr 2024 14:05:51. 

All simulations have been benchmarked over several seconds (wallclock time) without output. Benchmarking excludes initialization and is started just before the main time loop and finishes right after. The benchmarking results here are not very robust, timings that change with +-50% are not uncommon. Proper benchmarking for performance optimization uses the minimum or median of many executions, while we run a simulation for several time steps which effectively represents the mean, susceptible to outliers that slow down the simulation. However, this is what a user will experience in most situations anyway and the following therefore presents a rough idea of how fast a SpeedyWeather simulation will run, and how much memory it requires.

### Machine details

All benchmark simulation were single-threaded on a CPU:
```julia
julia> versioninfo()
Julia Version 1.10.2
Commit bd47eca2c8a (2024-03-01 10:14 UTC)
Build Info:
  Official https://julialang.org/ release
Platform Info:
  OS: macOS (x86_64-apple-darwin22.4.0)
  CPU: 8 × Intel(R) Core(TM) i5-1030NG7 CPU @ 1.10GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-15.0.7 (ORCJIT, icelake-client)
Threads: 1 default, 0 interactive, 1 GC (on 8 virtual cores)
Environment:
  LD_LIBRARY_PATH = /Users/milan/.julia/conda/3/lib:
```

### Explanation

Abbreviations in the tables below are as follows, omitted columns use defaults.
- NF: Number format, default: Float32
- T: Spectral resolution, maximum degree of spherical harmonics, default: T31
- L: Number of vertical layers, default: 8 (for 3D models)
- Grid: Horizontal grid, default: OctahedralGaussianGrid
- Rings: Grid-point resolution, number of latitude rings pole to pole
- Dynamics: With dynamics?, default: true
- Physics: With physical parameterizations?, default: true (for primitive equation models)
- Δt: time step [s].
- SYPD: Speed of simulation, simulated years per wallclock day.
- Memory: Memory footprint of simulation, variables and constants.

### Running the benchmarks

The benchmark suite here can be reproduced by executing:

```> julia manual_benchmarking.jl```

inside `the SpeedyWeather.jl/benchmark` folder. It will create this `README.md` which can be pushed to the repository for updates or comparison.
## Models, default setups

| Model | T | L | Physics | Δt | SYPD | Memory|
| - | - | - | - | - | - | - |
| BarotropicModel | 31 | 1 | false | 1800 | 19349 | 1.22 MB |
| ShallowWaterModel | 31 | 1 | false | 1800 | 10195 | 1.24 MB |
| PrimitiveDryModel | 31 | 8 | true | 1800 | 530 | 3.95 MB |
| PrimitiveWetModel | 31 | 8 | true | 1800 | 399 | 4.28 MB |

## Grids

| Model | T | L | Grid | Rings | Δt | SYPD | Memory|
| - | - | - | - | - | - | - | - |
| PrimitiveWetModel | 63 | 8 | FullGaussianGrid | 96 | 900 | 35 | 22.50 MB |
| PrimitiveWetModel | 63 | 8 | FullClenshawGrid | 95 | 900 | 34 | 22.29 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralGaussianGrid | 96 | 900 | 45 | 15.34 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralClenshawGrid | 95 | 900 | 54 | 15.12 MB |
| PrimitiveWetModel | 63 | 8 | HEALPixGrid | 95 | 900 | 77 | 11.46 MB |
| PrimitiveWetModel | 63 | 8 | OctaHEALPixGrid | 95 | 900 | 56 | 13.67 MB |

## Primitive wet model, resolution

| Model | T | L | Rings | Δt | SYPD | Memory|
| - | - | - | - | - | - | - |
| PrimitiveWetModel | 31 | 8 | 48 | 1800 | 378 | 4.28 MB |
| PrimitiveWetModel | 42 | 8 | 64 | 1350 | 161 | 7.22 MB |
| PrimitiveWetModel | 63 | 8 | 96 | 900 | 47 | 15.34 MB |
| PrimitiveWetModel | 85 | 8 | 128 | 675 | 20 | 26.73 MB |
| PrimitiveWetModel | 127 | 8 | 192 | 450 | 6 | 59.10 MB |
| PrimitiveWetModel | 170 | 8 | 256 | 338 | 2 | 105.37 MB |

## PrimitiveWetModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| - | - | - | - | - | - | - | - |
| PrimitiveWetModel | 31 | 8 | true | true | 1800 | 322 | 4.28 MB |
| PrimitiveWetModel | 31 | 8 | true | false | 1800 | 587 | 4.28 MB |
| PrimitiveWetModel | 31 | 8 | false | true | 1800 | 527 | 4.28 MB |

## Shallow water model, resolution

| Model | T | L | Rings | Δt | SYPD | Memory|
| - | - | - | - | - | - | - |
| ShallowWaterModel | 31 | 1 | 48 | 1800 | 9402 | 1.24 MB |
| ShallowWaterModel | 42 | 1 | 64 | 1350 | 4034 | 2.13 MB |
| ShallowWaterModel | 63 | 1 | 96 | 900 | 1132 | 4.71 MB |
| ShallowWaterModel | 85 | 1 | 128 | 675 | 415 | 8.45 MB |
| ShallowWaterModel | 127 | 1 | 192 | 450 | 105 | 19.60 MB |
| ShallowWaterModel | 170 | 1 | 256 | 338 | 39 | 36.41 MB |
| ShallowWaterModel | 255 | 1 | 384 | 225 | 9 | 89.44 MB |

## Primitive Equation, Float32 vs Float64

| Model | NF | T | L | Δt | SYPD | Memory|
| - | - | - | - | - | - | - |
| PrimitiveWetModel | Float32 | 31 | 8 | 1800 | 317 | 4.28 MB |
| PrimitiveWetModel | Float64 | 31 | 8 | 1800 | 343 | 8.03 MB |

## PrimitiveDryModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| - | - | - | - | - | - | - | - |
| PrimitiveDryModel | 31 | 8 | true | true | 1800 | 462 | 3.95 MB |
| PrimitiveDryModel | 31 | 8 | true | false | 1800 | 657 | 3.95 MB |
| PrimitiveDryModel | 31 | 8 | false | true | 1800 | 683 | 3.95 MB |

## Number of vertical layers

| Model | T | L | Δt | SYPD | Memory|
| - | - | - | - | - | - |
| PrimitiveWetModel | 31 | 4 | 1800 | 583 | 2.92 MB |
| PrimitiveWetModel | 31 | 8 | 1800 | 345 | 4.28 MB |
| PrimitiveWetModel | 31 | 12 | 1800 | 217 | 5.65 MB |
| PrimitiveWetModel | 31 | 16 | 1800 | 186 | 7.03 MB |
