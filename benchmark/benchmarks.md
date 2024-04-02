# Benchmarks

created for SpeedyWeather.jl v0.9.0 on Tue, 02 Apr 2024 11:46:14. 

All simulations have been benchmarked over several seconds (wallclock time) without output. Benchmarking excludes initialization and is started just before the main time loop and finishes right after. All simulations single-threaded on a CPU. 

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
| BarotropicModel | 31 | 1 | false | 1800 | 22436 | 1.22 MB |
| ShallowWaterModel | 31 | 1 | false | 1800 | 13086 | 1.24 MB |
| PrimitiveDryModel | 31 | 8 | true | 1800 | 582 | 3.95 MB |
| PrimitiveWetModel | 31 | 8 | true | 1800 | 472 | 4.28 MB |

## Grids

| Model | T | L | Grid | Rings | Δt | SYPD | Memory|
| - | - | - | - | - | - | - | - |
| PrimitiveWetModel | 63 | 8 | FullGaussianGrid | 96 | 900 | 39 | 22.50 MB |
| PrimitiveWetModel | 63 | 8 | FullClenshawGrid | 95 | 900 | 39 | 22.29 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralGaussianGrid | 96 | 900 | 56 | 15.34 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralClenshawGrid | 95 | 900 | 57 | 15.12 MB |
| PrimitiveWetModel | 63 | 8 | HEALPixGrid | 95 | 900 | 78 | 11.46 MB |
| PrimitiveWetModel | 63 | 8 | OctaHEALPixGrid | 95 | 900 | 64 | 13.67 MB |

## Primitive wet model, resolution

| Model | T | L | Rings | Δt | SYPD | Memory|
| - | - | - | - | - | - | - |
| PrimitiveWetModel | 31 | 8 | 48 | 1800 | 420 | 4.28 MB |
| PrimitiveWetModel | 42 | 8 | 64 | 1350 | 189 | 7.22 MB |
| PrimitiveWetModel | 63 | 8 | 96 | 900 | 53 | 15.34 MB |
| PrimitiveWetModel | 85 | 8 | 128 | 675 | 23 | 26.73 MB |
| PrimitiveWetModel | 127 | 8 | 192 | 450 | 6 | 59.10 MB |
| PrimitiveWetModel | 170 | 8 | 256 | 338 | 2 | 105.37 MB |

## PrimitiveWetModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| - | - | - | - | - | - | - | - |
| PrimitiveWetModel | 31 | 8 | true | true | 1800 | 421 | 4.28 MB |
| PrimitiveWetModel | 31 | 8 | true | false | 1800 | 659 | 4.28 MB |
| PrimitiveWetModel | 31 | 8 | false | true | 1800 | 561 | 4.28 MB |

## Shallow water model, resolution

| Model | T | L | Rings | Δt | SYPD | Memory|
| - | - | - | - | - | - | - |
| ShallowWaterModel | 31 | 1 | 48 | 1800 | 10448 | 1.24 MB |
| ShallowWaterModel | 42 | 1 | 64 | 1350 | 4620 | 2.13 MB |
| ShallowWaterModel | 63 | 1 | 96 | 900 | 1329 | 4.71 MB |
| ShallowWaterModel | 85 | 1 | 128 | 675 | 502 | 8.45 MB |
| ShallowWaterModel | 127 | 1 | 192 | 450 | 128 | 19.60 MB |
| ShallowWaterModel | 170 | 1 | 256 | 338 | 45 | 36.41 MB |
| ShallowWaterModel | 255 | 1 | 384 | 225 | 10 | 89.44 MB |

## Primitive Equation, Float32 vs Float64

| Model | NF | T | L | Δt | SYPD | Memory|
| - | - | - | - | - | - | - |
| PrimitiveWetModel | Float32 | 31 | 8 | 1800 | 339 | 4.28 MB |
| PrimitiveWetModel | Float64 | 31 | 8 | 1800 | 393 | 8.03 MB |

## PrimitiveDryModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| - | - | - | - | - | - | - | - |
| PrimitiveDryModel | 31 | 8 | true | true | 1800 | 519 | 3.95 MB |
| PrimitiveDryModel | 31 | 8 | true | false | 1800 | 714 | 3.95 MB |
| PrimitiveDryModel | 31 | 8 | false | true | 1800 | 802 | 3.95 MB |

## Number of vertical layers

| Model | T | L | Δt | SYPD | Memory|
| - | - | - | - | - | - |
| PrimitiveWetModel | 31 | 4 | 1800 | 693 | 2.92 MB |
| PrimitiveWetModel | 31 | 8 | 1800 | 399 | 4.28 MB |
| PrimitiveWetModel | 31 | 12 | 1800 | 276 | 5.65 MB |
| PrimitiveWetModel | 31 | 16 | 1800 | 214 | 7.03 MB |
