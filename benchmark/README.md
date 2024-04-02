# Benchmarks

created for SpeedyWeather.jl v0.9.0 on Tue, 02 Apr 2024 13:47:54. 

All simulations have been benchmarked over several seconds (wallclock time) without output. Benchmarking excludes initialization and is started just before the main time loop and finishes right after. All simulations single-threaded on a CPU, more details:
<details><summary>Machine details</summary>
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
</details>
The benchmarking results here are not very robust, timings that change with +-50% are not uncommon. Proper benchmarking for performance optimization uses the minimum or median of many executions, while we run a simulation for several time steps which effectively represents the mean, susceptible to outliers that slow down the simulation. However, this is what a user will experience in most situations anyway and the following therefore presents a rough idea of how fast a SpeedyWeather simulation will run, and how much memory it requires.

Explanation
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

## Models, default setups

| Model | T | L | Physics | Δt | SYPD | Memory|
| - | - | - | - | - | - | - |
| BarotropicModel | 31 | 1 | false | 1800 | 26309 | 1.22 MB |
| ShallowWaterModel | 31 | 1 | false | 1800 | 14511 | 1.24 MB |
| PrimitiveDryModel | 31 | 8 | true | 1800 | 670 | 3.95 MB |
| PrimitiveWetModel | 31 | 8 | true | 1800 | 510 | 4.28 MB |

## Grids

| Model | T | L | Grid | Rings | Δt | SYPD | Memory|
| - | - | - | - | - | - | - | - |
| PrimitiveWetModel | 63 | 8 | FullGaussianGrid | 96 | 900 | 41 | 22.50 MB |
| PrimitiveWetModel | 63 | 8 | FullClenshawGrid | 95 | 900 | 41 | 22.29 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralGaussianGrid | 96 | 900 | 65 | 15.34 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralClenshawGrid | 95 | 900 | 65 | 15.12 MB |
| PrimitiveWetModel | 63 | 8 | HEALPixGrid | 95 | 900 | 74 | 11.46 MB |
| PrimitiveWetModel | 63 | 8 | OctaHEALPixGrid | 95 | 900 | 64 | 13.67 MB |

## Primitive wet model, resolution

| Model | T | L | Rings | Δt | SYPD | Memory|
| - | - | - | - | - | - | - |
| PrimitiveWetModel | 31 | 8 | 48 | 1800 | 483 | 4.28 MB |
| PrimitiveWetModel | 42 | 8 | 64 | 1350 | 214 | 7.22 MB |
| PrimitiveWetModel | 63 | 8 | 96 | 900 | 57 | 15.34 MB |
| PrimitiveWetModel | 85 | 8 | 128 | 675 | 27 | 26.73 MB |
| PrimitiveWetModel | 127 | 8 | 192 | 450 | 7 | 59.10 MB |
| PrimitiveWetModel | 170 | 8 | 256 | 338 | 3 | 105.37 MB |

## PrimitiveWetModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| - | - | - | - | - | - | - | - |
| PrimitiveWetModel | 31 | 8 | true | true | 1800 | 373 | 4.28 MB |
| PrimitiveWetModel | 31 | 8 | true | false | 1800 | 713 | 4.28 MB |
| PrimitiveWetModel | 31 | 8 | false | true | 1800 | 658 | 4.28 MB |

## Shallow water model, resolution

| Model | T | L | Rings | Δt | SYPD | Memory|
| - | - | - | - | - | - | - |
| ShallowWaterModel | 31 | 1 | 48 | 1800 | 12110 | 1.24 MB |
| ShallowWaterModel | 42 | 1 | 64 | 1350 | 5239 | 2.13 MB |
| ShallowWaterModel | 63 | 1 | 96 | 900 | 1418 | 4.71 MB |
| ShallowWaterModel | 85 | 1 | 128 | 675 | 534 | 8.45 MB |
| ShallowWaterModel | 127 | 1 | 192 | 450 | 129 | 19.60 MB |
| ShallowWaterModel | 170 | 1 | 256 | 338 | 47 | 36.41 MB |
| ShallowWaterModel | 255 | 1 | 384 | 225 | 11 | 89.44 MB |

## Primitive Equation, Float32 vs Float64

| Model | NF | T | L | Δt | SYPD | Memory|
| - | - | - | - | - | - | - |
| PrimitiveWetModel | Float32 | 31 | 8 | 1800 | 375 | 4.28 MB |
| PrimitiveWetModel | Float64 | 31 | 8 | 1800 | 399 | 8.03 MB |

## PrimitiveDryModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| - | - | - | - | - | - | - | - |
| PrimitiveDryModel | 31 | 8 | true | true | 1800 | 555 | 3.95 MB |
| PrimitiveDryModel | 31 | 8 | true | false | 1800 | 804 | 3.95 MB |
| PrimitiveDryModel | 31 | 8 | false | true | 1800 | 820 | 3.95 MB |

## Number of vertical layers

| Model | T | L | Δt | SYPD | Memory|
| - | - | - | - | - | - |
| PrimitiveWetModel | 31 | 4 | 1800 | 711 | 2.92 MB |
| PrimitiveWetModel | 31 | 8 | 1800 | 418 | 4.28 MB |
| PrimitiveWetModel | 31 | 12 | 1800 | 279 | 5.65 MB |
| PrimitiveWetModel | 31 | 16 | 1800 | 202 | 7.03 MB |
