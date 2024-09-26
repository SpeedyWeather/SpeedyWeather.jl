# Benchmarks

created for SpeedyWeather.jl v0.11.0 on Thu, 26 Sep 2024 13:52:34. 

All simulations have been benchmarked over several seconds (wallclock time) without output. Benchmarking excludes initialization and is started just before the main time loop and finishes right after. The benchmarking results here are not very robust, timings that change with +-50% are not uncommon. Proper benchmarking for performance optimization uses the minimum or median of many executions, while we run a simulation for several time steps which effectively represents the mean, susceptible to outliers that slow down the simulation. However, this is what a user will experience in most situations anyway and the following therefore presents a rough idea of how fast a SpeedyWeather simulation will run, and how much memory it requires.

### Machine details

All benchmark simulation were single-threaded on a CPU:
```julia
julia> versioninfo()
Julia Version 1.10.4
Commit 48d4fd48430 (2024-06-04 10:41 UTC)
Build Info:
  Official https://julialang.org/ release
Platform Info:
  OS: macOS (arm64-apple-darwin22.4.0)
  CPU: 8 × Apple M3
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-15.0.7 (ORCJIT, apple-m1)
Threads: 1 default, 0 interactive, 1 GC (on 4 virtual cores)
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
| BarotropicModel | 31 | 1 | false | 1800 | 62943 | 1.17 MB |
| ShallowWaterModel | 31 | 1 | false | 1800 | 33470 | 1.19 MB |
| PrimitiveDryModel | 31 | 8 | true | 1800 | 2007 | 4.11 MB |
| PrimitiveWetModel | 31 | 8 | true | 1800 | 1247 | 4.47 MB |

## Grids

| Model | T | L | Grid | Rings | Δt | SYPD | Memory|
| - | - | - | - | - | - | - | - |
| PrimitiveWetModel | 63 | 8 | FullGaussianGrid | 96 | 900 | 110 | 23.40 MB |
| PrimitiveWetModel | 63 | 8 | FullClenshawGrid | 95 | 900 | 139 | 23.19 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralGaussianGrid | 96 | 900 | 198 | 15.91 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralClenshawGrid | 95 | 900 | 184 | 15.68 MB |
| PrimitiveWetModel | 63 | 8 | HEALPixGrid | 95 | 900 | 201 | 11.86 MB |
| PrimitiveWetModel | 63 | 8 | OctaHEALPixGrid | 95 | 900 | 190 | 14.16 MB |

## Primitive wet model, resolution

| Model | T | L | Rings | Δt | SYPD | Memory|
| - | - | - | - | - | - | - |
| PrimitiveWetModel | 31 | 8 | 48 | 1800 | 1103 | 4.47 MB |
| PrimitiveWetModel | 42 | 8 | 64 | 1350 | 604 | 7.51 MB |
| PrimitiveWetModel | 63 | 8 | 96 | 900 | 199 | 15.91 MB |
| PrimitiveWetModel | 85 | 8 | 128 | 675 | 71 | 27.66 MB |
| PrimitiveWetModel | 127 | 8 | 192 | 450 | 19 | 60.99 MB |
| PrimitiveWetModel | 170 | 8 | 256 | 338 | 6 | 108.55 MB |

## PrimitiveWetModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| - | - | - | - | - | - | - | - |
| PrimitiveWetModel | 31 | 8 | true | true | 1800 | 1097 | 4.47 MB |
| PrimitiveWetModel | 31 | 8 | true | false | 1800 | 3345 | 4.47 MB |
| PrimitiveWetModel | 31 | 8 | false | true | 1800 | 2011 | 4.47 MB |

## Individual dynamics functions


### PrimitiveWetModel | Float32 | T31 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| - | - | - | - |
| pressure_gradient_flux! | 27.792 μs| 12.88 KiB| 260 |
| linear_virtual_temperature! | 952.870 ns| 0 bytes| 0 |
| temperature_anomaly! | 3.266 μs| 0 bytes| 0 |
| geopotential! | 2.995 μs| 0 bytes| 0 |
| vertical_integration! | 12.417 μs| 0 bytes| 0 |
| surface_pressure_tendency! | 10.750 μs| 6.44 KiB| 130 |
| vertical_velocity! | 4.250 μs| 0 bytes| 0 |
| linear_pressure_gradient! | 933.333 ns| 0 bytes| 0 |
| vertical_advection! | 86.417 μs| 2.81 KiB| 24 |
| vordiv_tendencies! | 178.500 μs| 98.38 KiB| 1940 |
| temperature_tendency! | 257.750 μs| 147.56 KiB| 2910 |
| humidity_tendency! | 254.125 μs| 147.56 KiB| 2910 |
| bernoulli_potential! | 86.917 μs| 49.19 KiB| 970 |

## Shallow water model, resolution

| Model | T | L | Rings | Δt | SYPD | Memory|
| - | - | - | - | - | - | - |
| ShallowWaterModel | 31 | 1 | 48 | 1800 | 36300 | 1.19 MB |
| ShallowWaterModel | 42 | 1 | 64 | 1350 | 14941 | 2.02 MB |
| ShallowWaterModel | 63 | 1 | 96 | 900 | 3816 | 4.40 MB |
| ShallowWaterModel | 85 | 1 | 128 | 675 | 1466 | 7.86 MB |
| ShallowWaterModel | 127 | 1 | 192 | 450 | 391 | 18.16 MB |
| ShallowWaterModel | 170 | 1 | 256 | 338 | 119 | 33.76 MB |
| ShallowWaterModel | 255 | 1 | 384 | 225 | 26 | 83.28 MB |

## Primitive Equation, Float32 vs Float64

| Model | NF | T | L | Δt | SYPD | Memory|
| - | - | - | - | - | - | - |
| PrimitiveWetModel | Float32 | 31 | 8 | 1800 | 1261 | 4.47 MB |
| PrimitiveWetModel | Float64 | 31 | 8 | 1800 | 1260 | 8.36 MB |

## PrimitiveDryModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| - | - | - | - | - | - | - | - |
| PrimitiveDryModel | 31 | 8 | true | true | 1800 | 2007 | 4.11 MB |
| PrimitiveDryModel | 31 | 8 | true | false | 1800 | 5105 | 4.11 MB |
| PrimitiveDryModel | 31 | 8 | false | true | 1800 | 2521 | 4.11 MB |

## Number of vertical layers

| Model | T | L | Δt | SYPD | Memory|
| - | - | - | - | - | - |
| PrimitiveWetModel | 31 | 4 | 1800 | 1833 | 3.00 MB |
| PrimitiveWetModel | 31 | 8 | 1800 | 1257 | 4.47 MB |
| PrimitiveWetModel | 31 | 12 | 1800 | 1087 | 5.94 MB |
| PrimitiveWetModel | 31 | 16 | 1800 | 819 | 7.42 MB |
