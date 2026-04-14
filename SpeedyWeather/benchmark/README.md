# Benchmarks

created for SpeedyWeather.jl v0.19.0 on Thu, 02 Apr 2026 12:43:18. 

All simulations have been benchmarked over several seconds (wallclock time) without output. Benchmarking excludes initialization and is started just before the main time loop and finishes right after. The benchmarking results here are not very robust, timings that change with +-50% are not uncommon. Proper benchmarking for performance optimization uses the minimum or median of many executions, while we run a simulation for several time steps which effectively represents the mean, susceptible to outliers that slow down the simulation. However, this is what a user will experience in most situations anyway and the following therefore presents a rough idea of how fast a SpeedyWeather simulation will run, and how much memory it requires.

### Machine details

All benchmark simulation were single-threaded on a CPU:
```julia
julia> versioninfo()
Julia Version 1.11.7
Commit f2b3dbda30a (2025-09-08 12:10 UTC)
Build Info:
  Official https://julialang.org/ release
Platform Info:
  OS: macOS (arm64-apple-darwin24.0.0)
  CPU: 8 × Apple M3
  WORD_SIZE: 64
  LLVM: libLLVM-16.0.6 (ORCJIT, apple-m2)
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
| BarotropicModel | 31 | 1 | false | 2400 | 43763 | 688.35 KB |
| ShallowWaterModel | 31 | 1 | false | 2400 | 26763 | 822.86 KB |
| PrimitiveDryModel | 31 | 8 | true | 2400 | 2240 | 4.16 MB |
| PrimitiveWetModel | 31 | 8 | true | 2400 | 1211 | 4.84 MB |

## Grids

| Model | T | L | Grid | Rings | Δt | SYPD | Memory|
| - | - | - | - | - | - | - | - |
| PrimitiveWetModel | 63 | 8 | FullGaussianGrid | 96 | 1200 | 89 | 25.36 MB |
| PrimitiveWetModel | 63 | 8 | FullClenshawGrid | 95 | 1200 | 104 | 25.14 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralGaussianGrid | 96 | 1200 | 188 | 17.42 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralClenshawGrid | 95 | 1200 | 148 | 17.18 MB |
| PrimitiveWetModel | 63 | 8 | HEALPixGrid | 95 | 1200 | 192 | 12.73 MB |
| PrimitiveWetModel | 63 | 8 | OctaHEALPixGrid | 95 | 1200 | 146 | 15.50 MB |

## Primitive wet model, resolution

| Model | T | L | Rings | Δt | SYPD | Memory|
| - | - | - | - | - | - | - |
| PrimitiveWetModel | 31 | 8 | 48 | 2400 | 1491 | 4.84 MB |
| PrimitiveWetModel | 42 | 8 | 64 | 1800 | 664 | 8.18 MB |
| PrimitiveWetModel | 63 | 8 | 96 | 1200 | 191 | 17.42 MB |
| PrimitiveWetModel | 85 | 8 | 128 | 900 | 65 | 30.35 MB |
| PrimitiveWetModel | 127 | 8 | 192 | 600 | 20 | 67.03 MB |
| PrimitiveWetModel | 170 | 8 | 256 | 450 | 7 | 119.25 MB |

## PrimitiveWetModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| - | - | - | - | - | - | - | - |
| PrimitiveWetModel | 31 | 8 | true | true | 2400 | 1493 | 4.84 MB |
| PrimitiveWetModel | 31 | 8 | true | false | 2400 | 3396 | 4.84 MB |
| PrimitiveWetModel | 31 | 8 | false | true | 2400 | 1903 | 4.84 MB |

## Individual dynamics functions


### PrimitiveWetModel | Float32 | T31 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| - | - | - | - |
| pressure_gradient_flux! | 39.292 μs| 100.28 KiB| 790 |
| linear_virtual_temperature! | 1.825 μs| 0 bytes| 0 |
| geopotential! | 6.258 μs| 5.61 KiB| 23 |
| vertical_integration! | 14.375 μs| 0 bytes| 0 |
| surface_pressure_tendency! | 14.666 μs| 24.66 KiB| 288 |
| vertical_velocity! | 48.667 μs| 384.31 KiB| 12 |
| linear_pressure_gradient! | 1.754 μs| 0 bytes| 0 |
| vertical_advection! | 110.625 μs| 8.62 KiB| 100 |
| vordiv_tendencies! | 289.209 μs| 259.66 KiB| 724 |
| temperature_tendency! | 314.833 μs| 381.95 KiB| 1027 |
| humidity_tendency! | 288.292 μs| 380.59 KiB| 1017 |
| bernoulli_potential! | 98.833 μs| 510.22 KiB| 345 |

## Shallow water model, resolution

| Model | T | L | Rings | Δt | SYPD | Memory|
| - | - | - | - | - | - | - |
| ShallowWaterModel | 31 | 1 | 48 | 2400 | 27667 | 822.86 KB |
| ShallowWaterModel | 42 | 1 | 64 | 1800 | 12256 | 1.44 MB |
| ShallowWaterModel | 63 | 1 | 96 | 1200 | 3691 | 3.25 MB |
| ShallowWaterModel | 85 | 1 | 128 | 900 | 1296 | 5.94 MB |
| ShallowWaterModel | 127 | 1 | 192 | 600 | 386 | 14.14 MB |
| ShallowWaterModel | 170 | 1 | 256 | 450 | 105 | 26.86 MB |
| ShallowWaterModel | 255 | 1 | 384 | 300 | 24 | 68.30 MB |

## Primitive Equation, Float32 vs Float64

| Model | NF | T | L | Δt | SYPD | Memory|
| - | - | - | - | - | - | - |
| PrimitiveWetModel | Float32 | 31 | 8 | 2400 | 1483 | 4.84 MB |
| PrimitiveWetModel | Float64 | 31 | 8 | 2400 | 1202 | 9.00 MB |

## PrimitiveDryModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| - | - | - | - | - | - | - | - |
| PrimitiveDryModel | 31 | 8 | true | true | 2400 | 2236 | 4.16 MB |
| PrimitiveDryModel | 31 | 8 | true | false | 2400 | 3305 | 4.16 MB |
| PrimitiveDryModel | 31 | 8 | false | true | 2400 | 3341 | 4.16 MB |

## Number of vertical layers

| Model | T | L | Δt | SYPD | Memory|
| - | - | - | - | - | - |
| PrimitiveWetModel | 31 | 4 | 2400 | 2289 | 3.11 MB |
| PrimitiveWetModel | 31 | 8 | 2400 | 1268 | 4.84 MB |
| PrimitiveWetModel | 31 | 12 | 2400 | 980 | 6.58 MB |
| PrimitiveWetModel | 31 | 16 | 2400 | 547 | 8.33 MB |
