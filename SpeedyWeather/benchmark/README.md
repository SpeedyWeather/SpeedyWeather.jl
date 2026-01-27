# Benchmarks

created for SpeedyWeather.jl v0.15.0 on Fri, 23 May 2025 12:48:31. 

All simulations have been benchmarked over several seconds (wallclock time) without output. Benchmarking excludes initialization and is started just before the main time loop and finishes right after. The benchmarking results here are not very robust, timings that change with +-50% are not uncommon. Proper benchmarking for performance optimization uses the minimum or median of many executions, while we run a simulation for several time steps which effectively represents the mean, susceptible to outliers that slow down the simulation. However, this is what a user will experience in most situations anyway and the following therefore presents a rough idea of how fast a SpeedyWeather simulation will run, and how much memory it requires.

### Machine details

All benchmark simulation were single-threaded on a CPU:
```julia
julia> versioninfo()
Julia Version 1.10.9
Commit 5595d20a287 (2025-03-10 12:51 UTC)
Build Info:
  Official https://julialang.org/ release
Platform Info:
  OS: macOS (arm64-apple-darwin24.0.0)
  CPU: 8 × Apple M3
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-15.0.7 (ORCJIT, apple-m1)
Threads: 1 default, 0 interactive, 1 GC (on 4 virtual cores)
Environment:
  JULIA_EDITOR = code
  JULIA_VSCODE_REPL = 1
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
| BarotropicModel | 31 | 1 | false | 2400 | 61406 | 1.66 MB |
| ShallowWaterModel | 31 | 1 | false | 2400 | 36373 | 1.68 MB |
| PrimitiveDryModel | 31 | 8 | true | 2400 | 1752 | 4.85 MB |
| PrimitiveWetModel | 31 | 8 | true | 2400 | 1540 | 5.18 MB |

## Grids

| Model | T | L | Grid | Rings | Δt | SYPD | Memory|
| - | - | - | - | - | - | - | - |
| PrimitiveWetModel | 63 | 8 | FullGaussianGrid | 96 | 1200 | 75 | 26.83 MB |
| PrimitiveWetModel | 63 | 8 | FullClenshawGrid | 95 | 1200 | 104 | 26.59 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralGaussianGrid | 96 | 1200 | 254 | 18.44 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralClenshawGrid | 95 | 1200 | 147 | 18.19 MB |
| PrimitiveWetModel | 63 | 8 | HEALPixGrid | 95 | 1200 | 127 | 13.50 MB |
| PrimitiveWetModel | 63 | 8 | OctaHEALPixGrid | 95 | 1200 | 173 | 16.41 MB |

## Primitive wet model, resolution

| Model | T | L | Rings | Δt | SYPD | Memory|
| - | - | - | - | - | - | - |
| PrimitiveWetModel | 31 | 8 | 48 | 2400 | 1893 | 5.18 MB |
| PrimitiveWetModel | 42 | 8 | 64 | 1800 | 787 | 8.70 MB |
| PrimitiveWetModel | 63 | 8 | 96 | 1200 | 259 | 18.44 MB |
| PrimitiveWetModel | 85 | 8 | 128 | 900 | 97 | 32.03 MB |
| PrimitiveWetModel | 127 | 8 | 192 | 600 | 26 | 70.50 MB |
| PrimitiveWetModel | 170 | 8 | 256 | 450 | 8 | 125.18 MB |

## PrimitiveWetModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| - | - | - | - | - | - | - | - |
| PrimitiveWetModel | 31 | 8 | true | true | 2400 | 1921 | 5.18 MB |
| PrimitiveWetModel | 31 | 8 | true | false | 2400 | 3400 | 5.18 MB |
| PrimitiveWetModel | 31 | 8 | false | true | 2400 | 2668 | 5.18 MB |

## Individual dynamics functions


### PrimitiveWetModel | Float32 | T31 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| - | - | - | - |
| pressure_gradient_flux! | 45.125 μs| 182.81 KiB| 2876 |
| linear_virtual_temperature! | 952.870 ns| 0 bytes| 0 |
| temperature_anomaly! | 3.255 μs| 0 bytes| 0 |
| geopotential! | 2.806 μs| 0 bytes| 0 |
| vertical_integration! | 12.542 μs| 0 bytes| 0 |
| surface_pressure_tendency! | 17.000 μs| 65.16 KiB| 1390 |
| vertical_velocity! | 4.315 μs| 0 bytes| 0 |
| linear_pressure_gradient! | 933.897 ns| 0 bytes| 0 |
| vertical_advection! | 86.542 μs| 3.45 KiB| 21 |
| vordiv_tendencies! | 161.416 μs| 15.75 KiB| 210 |
| temperature_tendency! | 233.333 μs| 19.88 KiB| 297 |
| humidity_tendency! | 225.917 μs| 19.88 KiB| 297 |
| bernoulli_potential! | 127.375 μs| 216.00 KiB| 13536 |

## Shallow water model, resolution

| Model | T | L | Rings | Δt | SYPD | Memory|
| - | - | - | - | - | - | - |
| ShallowWaterModel | 31 | 1 | 48 | 2400 | 34894 | 1.68 MB |
| ShallowWaterModel | 42 | 1 | 64 | 1800 | 17106 | 2.84 MB |
| ShallowWaterModel | 63 | 1 | 96 | 1200 | 4684 | 6.08 MB |
| ShallowWaterModel | 85 | 1 | 128 | 900 | 1974 | 10.72 MB |
| ShallowWaterModel | 127 | 1 | 192 | 600 | 440 | 24.32 MB |
| ShallowWaterModel | 170 | 1 | 256 | 450 | 158 | 44.45 MB |
| ShallowWaterModel | 255 | 1 | 384 | 300 | 23 | 106.77 MB |

## Primitive Equation, Float32 vs Float64

| Model | NF | T | L | Δt | SYPD | Memory|
| - | - | - | - | - | - | - |
| PrimitiveWetModel | Float32 | 31 | 8 | 2400 | 1647 | 5.18 MB |
| PrimitiveWetModel | Float64 | 31 | 8 | 2400 | 1457 | 9.60 MB |

## PrimitiveDryModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| - | - | - | - | - | - | - | - |
| PrimitiveDryModel | 31 | 8 | true | true | 2400 | 2617 | 4.85 MB |
| PrimitiveDryModel | 31 | 8 | true | false | 2400 | 4340 | 4.85 MB |
| PrimitiveDryModel | 31 | 8 | false | true | 2400 | 2741 | 4.85 MB |

## Number of vertical layers

| Model | T | L | Δt | SYPD | Memory|
| - | - | - | - | - | - |
| PrimitiveWetModel | 31 | 4 | 2400 | 2520 | 3.56 MB |
| PrimitiveWetModel | 31 | 8 | 2400 | 1922 | 5.18 MB |
| PrimitiveWetModel | 31 | 12 | 2400 | 1489 | 6.81 MB |
| PrimitiveWetModel | 31 | 16 | 2400 | 1612 | 8.44 MB |
