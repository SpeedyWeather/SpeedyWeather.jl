# Benchmarks

Performance benchmarks for SpeedyWeather.jl, collected across multiple architectures. Each architecture's results live in its own section below; the overview table at the top compares the headline PrimitiveWet resolution sweep across all archs that have been benchmarked so far.

All simulations are benchmarked over several seconds (wallclock time) without output. Benchmarking excludes initialization and is started just before the main time loop and finishes right after. The benchmarking results here are not very robust; timings that change by ±50% are not uncommon. Proper benchmarking for performance optimization uses the minimum or median of many executions, while we run a simulation for several time steps which effectively represents the mean, susceptible to outliers. However, this is what a user will experience in most situations anyway and the following therefore presents a rough idea of how fast a SpeedyWeather simulation will run, and how much memory it requires.

### Explanation

Abbreviations in the tables below are as follows; omitted columns use defaults.
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

Reproduce the benchmark suite by running, from `SpeedyWeather/benchmark`:

```
julia --project=. manual_benchmarking.jl                # CPU (auto-labelled cpu-arm or cpu-x86)
julia --project=. manual_benchmarking.jl gpu            # CUDA GPU
julia --project=. manual_benchmarking.jl reactant-cpu   # Reactant on CPU
julia --project=. manual_benchmarking.jl reactant-gpu   # Reactant on CUDA GPU
```

Each run updates only its own architecture's section in this `README.md`; results for other architectures are preserved via `benchmark_results.json`.

## Overview: PrimitiveWet resolution across architectures

Simulated years per wallclock day (SYPD) for the `PrimitiveWetModel` resolution sweep, one column per architecture. Empty cells mean the architecture has not yet been benchmarked or that suite was skipped. Comparison figures across architectures are available on the documentation's `Benchmarks` page.

| T | L | cpu-arm | cpu-x86 |
| --- | --- | --- | --- |
| 31 | 8 | 1484 | 531 |
| 42 | 8 | 662 | 305 |
| 63 | 8 | 188 | 91 |
| 85 | 8 | 65 | 40 |
| 85 | 16 | 42 | 14 |
| 85 | 24 | 23 | — |
| 127 | 8 | 18 | 10 |
| 127 | 16 | 10.0 | 4.8 |
| 127 | 24 | 6.2 | — |
| 170 | 8 | 7.0 | 4.2 |
| 170 | 16 | 3.7 | 2.3 |
| 170 | 24 | 2.2 | — |
| 255 | 8 | 1.8 | 1.1 |
| 255 | 16 | 0.9 | 0.6 |
| 255 | 24 | 0.6 | — |

## Architecture: `cpu-arm`

Created for SpeedyWeather.jl v0.20.1+DEV on Thu, 28 May 2026 11:48:39.

### Machine details

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


### Models, default setups

| Model | T | L | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| BarotropicModel | 31 | 1 | false | 2400 | 47411 | 688.34 KB |
| ShallowWaterModel | 31 | 1 | false | 2400 | 26114 | 822.86 KB |
| PrimitiveDryModel | 31 | 8 | true | 2400 | 2677 | 4.16 MB |
| PrimitiveWetModel | 31 | 8 | true | 2400 | 1222 | 4.84 MB |

### Shallow water model, resolution

| Model | T | L | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 31 | 1 | 48 | 2400 | 27856 | 822.86 KB |
| ShallowWaterModel | 42 | 1 | 64 | 1800 | 13085 | 1.44 MB |
| ShallowWaterModel | 63 | 1 | 96 | 1200 | 3655 | 3.25 MB |
| ShallowWaterModel | 85 | 1 | 128 | 900 | 1441 | 5.94 MB |
| ShallowWaterModel | 127 | 1 | 192 | 600 | 316 | 14.14 MB |
| ShallowWaterModel | 170 | 1 | 256 | 450 | 132 | 26.86 MB |
| ShallowWaterModel | 255 | 1 | 384 | 300 | 33 | 68.30 MB |

### Primitive wet model, resolution

| Model | T | L | Rings | Transform | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 8 | 48 | default | 2400 | 1484 | 4.84 MB |
| PrimitiveWetModel | 42 | 8 | 64 | default | 1800 | 662 | 8.18 MB |
| PrimitiveWetModel | 63 | 8 | 96 | default | 1200 | 188 | 17.42 MB |
| PrimitiveWetModel | 85 | 8 | 128 | default | 900 | 65 | 30.35 MB |
| PrimitiveWetModel | 127 | 8 | 192 | default | 600 | 18 | 67.03 MB |
| PrimitiveWetModel | 170 | 8 | 256 | default | 450 | 7.0 | 119.25 MB |
| PrimitiveWetModel | 255 | 8 | 384 | default | 300 | 1.8 | 272.27 MB |
| PrimitiveWetModel | 85 | 16 | 128 | default | 900 | 42 | 51.65 MB |
| PrimitiveWetModel | 127 | 16 | 192 | default | 600 | 10.0 | 113.22 MB |
| PrimitiveWetModel | 170 | 16 | 256 | default | 450 | 3.7 | 200.00 MB |
| PrimitiveWetModel | 255 | 16 | 384 | default | 300 | 0.9 | 450.63 MB |
| PrimitiveWetModel | 85 | 24 | 128 | default | 900 | 23 | 73.00 MB |
| PrimitiveWetModel | 127 | 24 | 192 | default | 600 | 6.2 | 159.48 MB |
| PrimitiveWetModel | 170 | 24 | 256 | default | 450 | 2.2 | 280.85 MB |
| PrimitiveWetModel | 255 | 24 | 384 | default | 300 | 0.6 | 629.13 MB |
| PrimitiveWetModel | 31 | 8 | 48 | matrix | 2400 | 484 | 47.07 MB |
| PrimitiveWetModel | 42 | 8 | 64 | matrix | 1800 | 173 | 132.12 MB |
| PrimitiveWetModel | 63 | 8 | 96 | matrix | 1200 | 34 | 579.07 MB |
| PrimitiveWetModel | 85 | 8 | 128 | matrix | 900 | 8.7 | 1.74 GB |
| PrimitiveWetModel | 127 | 8 | 192 | matrix | 600 | 0.6 | 8.17 GB |
| PrimitiveWetModel | 85 | 16 | 128 | matrix | 900 | 7.1 | 1.76 GB |
| PrimitiveWetModel | 127 | 16 | 192 | matrix | 600 | 1.1 | 8.22 GB |
| PrimitiveWetModel | 85 | 24 | 128 | matrix | 900 | 5.2 | 1.78 GB |
| PrimitiveWetModel | 127 | 24 | 192 | matrix | 600 | 0.9 | 8.26 GB |

### Primitive Equation, Float32 vs Float64

| Model | NF | T | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | Float32 | 31 | 8 | 2400 | 1145 | 4.84 MB |
| PrimitiveWetModel | Float64 | 31 | 8 | 2400 | 1337 | 9.00 MB |

### Grids

| Model | T | L | Grid | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 63 | 8 | FullGaussianGrid | 96 | 1200 | 76 | 25.36 MB |
| PrimitiveWetModel | 63 | 8 | FullClenshawGrid | 95 | 1200 | 96 | 25.14 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralGaussianGrid | 96 | 1200 | 182 | 17.42 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralClenshawGrid | 95 | 1200 | 122 | 17.18 MB |
| PrimitiveWetModel | 63 | 8 | HEALPixGrid | 95 | 1200 | 257 | 12.73 MB |
| PrimitiveWetModel | 63 | 8 | OctaHEALPixGrid | 95 | 1200 | 190 | 15.50 MB |

### Number of vertical layers

| Model | T | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 4 | 2400 | 2184 | 3.11 MB |
| PrimitiveWetModel | 31 | 8 | 2400 | 1263 | 4.84 MB |
| PrimitiveWetModel | 31 | 12 | 2400 | 981 | 6.58 MB |
| PrimitiveWetModel | 31 | 16 | 2400 | 619 | 8.33 MB |

### PrimitiveDryModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveDryModel | 31 | 8 | true | true | 2400 | 2697 | 4.16 MB |
| PrimitiveDryModel | 31 | 8 | true | false | 2400 | 4494 | 4.16 MB |
| PrimitiveDryModel | 31 | 8 | false | true | 2400 | 3389 | 4.16 MB |

### PrimitiveWetModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 8 | true | true | 2400 | 1257 | 4.84 MB |
| PrimitiveWetModel | 31 | 8 | true | false | 2400 | 3354 | 4.84 MB |
| PrimitiveWetModel | 31 | 8 | false | true | 2400 | 1921 | 4.84 MB |

### Individual dynamics functions


#### PrimitiveWetModel | Float32 | T31 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| --- | --- | --- | --- |
| pressure_gradient_flux! | 40.041 μs| 100.28 KiB| 790 |
| linear_virtual_temperature! | 1.754 μs| 0 bytes| 0 |
| geopotential! | 6.308 μs| 5.61 KiB| 23 |
| vertical_integration! | 13.833 μs| 0 bytes| 0 |
| surface_pressure_tendency! | 14.791 μs| 24.66 KiB| 288 |
| vertical_velocity! | 48.125 μs| 384.31 KiB| 12 |
| linear_pressure_gradient! | 1.750 μs| 0 bytes| 0 |
| vertical_advection! | 98.083 μs| 8.62 KiB| 100 |
| vordiv_tendencies! | 287.333 μs| 259.66 KiB| 724 |
| temperature_tendency! | 313.875 μs| 381.95 KiB| 1027 |
| humidity_tendency! | 287.958 μs| 380.59 KiB| 1017 |
| bernoulli_potential! | 99.125 μs| 510.22 KiB| 345 |

## Architecture: `cpu-x86`

Created for SpeedyWeather.jl v0.20.1+DEV on Fri, 22 May 2026 15:41:35.

### Machine details

```julia
julia> versioninfo()
Julia Version 1.12.2
Commit ca9b6662be4 (2025-11-20 16:25 UTC)
Build Info:
  Official https://julialang.org release
Platform Info:
  OS: Linux (x86_64-linux-gnu)
  CPU: 128 × AMD EPYC 9354 32-Core Processor
  WORD_SIZE: 64
  LLVM: libLLVM-18.1.7 (ORCJIT, znver4)
  GC: Built with stock GC
Threads: 1 default, 1 interactive, 1 GC (on 128 virtual cores)
Environment:
  LD_LIBRARY_PATH = /usr/local/lib:/usr/local/lib:
```


### Models, default setups

| Model | T | L | Physics | Δt | SYPD | Memory|
| - | - | - | - | - | - | - |
| BarotropicModel | 31 | 1 | false | 2400 | 19080 | 688.32 KB |
| ShallowWaterModel | 31 | 1 | false | 2400 | 13005 | 822.84 KB |
| PrimitiveDryModel | 31 | 8 | true | 2400 | 722 | 4.16 MB |
| PrimitiveWetModel | 31 | 8 | true | 2400 | 703 | 4.84 MB |

### Shallow water model, resolution

| Model | T | L | Rings | Δt | SYPD | Memory|
| - | - | - | - | - | - | - |
| ShallowWaterModel | 31 | 1 | 48 | 2400 | 12672 | 822.84 KB |
| ShallowWaterModel | 42 | 1 | 64 | 1800 | 5832 | 1.44 MB |
| ShallowWaterModel | 63 | 1 | 96 | 1200 | 1766 | 3.25 MB |
| ShallowWaterModel | 85 | 1 | 128 | 900 | 723 | 5.94 MB |
| ShallowWaterModel | 127 | 1 | 192 | 600 | 198 | 14.14 MB |
| ShallowWaterModel | 170 | 1 | 256 | 450 | 72 | 26.86 MB |
| ShallowWaterModel | 255 | 1 | 384 | 300 | 17 | 68.30 MB |

### Primitive wet model, resolution

| Model | T | L | Rings | Δt | SYPD | Memory|
| - | - | - | - | - | - | - |
| PrimitiveWetModel | 31 | 8 | 48 | 2400 | 531 | 4.84 MB |
| PrimitiveWetModel | 42 | 8 | 64 | 1800 | 305 | 8.18 MB |
| PrimitiveWetModel | 63 | 8 | 96 | 1200 | 91 | 17.42 MB |
| PrimitiveWetModel | 85 | 8 | 128 | 900 | 40 | 30.35 MB |
| PrimitiveWetModel | 127 | 8 | 192 | 600 | 10 | 67.03 MB |
| PrimitiveWetModel | 170 | 8 | 256 | 450 | 4 | 119.25 MB |
| PrimitiveWetModel | 255 | 8 | 384 | 300 | 1 | 272.27 MB |
| PrimitiveWetModel | 85 | 16 | 128 | 900 | 14 | 51.65 MB |
| PrimitiveWetModel | 127 | 16 | 192 | 600 | 5 | 113.22 MB |
| PrimitiveWetModel | 170 | 16 | 256 | 450 | 2 | 200.00 MB |
| PrimitiveWetModel | 255 | 16 | 384 | 300 | 1 | 450.63 MB |

### Primitive Equation, Float32 vs Float64

| Model | NF | T | L | Δt | SYPD | Memory|
| - | - | - | - | - | - | - |
| PrimitiveWetModel | Float32 | 31 | 8 | 2400 | 736 | 4.84 MB |
| PrimitiveWetModel | Float64 | 31 | 8 | 2400 | 657 | 9.00 MB |

### Grids

| Model | T | L | Grid | Rings | Δt | SYPD | Memory|
| - | - | - | - | - | - | - | - |
| PrimitiveWetModel | 63 | 8 | FullGaussianGrid | 96 | 1200 | 71 | 25.36 MB |
| PrimitiveWetModel | 63 | 8 | FullClenshawGrid | 95 | 1200 | 66 | 25.14 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralGaussianGrid | 96 | 1200 | 93 | 17.42 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralClenshawGrid | 95 | 1200 | 88 | 17.18 MB |
| PrimitiveWetModel | 63 | 8 | HEALPixGrid | 95 | 1200 | 59 | 12.73 MB |
| PrimitiveWetModel | 63 | 8 | OctaHEALPixGrid | 95 | 1200 | 98 | 15.50 MB |

### Number of vertical layers

| Model | T | L | Δt | SYPD | Memory|
| - | - | - | - | - | - |
| PrimitiveWetModel | 31 | 4 | 2400 | 1162 | 3.11 MB |
| PrimitiveWetModel | 31 | 8 | 2400 | 731 | 4.84 MB |
| PrimitiveWetModel | 31 | 12 | 2400 | 510 | 6.58 MB |
| PrimitiveWetModel | 31 | 16 | 2400 | 203 | 8.33 MB |

### PrimitiveDryModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| - | - | - | - | - | - | - | - |
| PrimitiveDryModel | 31 | 8 | true | true | 2400 | 1214 | 4.16 MB |
| PrimitiveDryModel | 31 | 8 | true | false | 2400 | 1905 | 4.16 MB |
| PrimitiveDryModel | 31 | 8 | false | true | 2400 | 1625 | 4.16 MB |

### PrimitiveWetModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| - | - | - | - | - | - | - | - |
| PrimitiveWetModel | 31 | 8 | true | true | 2400 | 738 | 4.84 MB |
| PrimitiveWetModel | 31 | 8 | true | false | 2400 | 1338 | 4.84 MB |
| PrimitiveWetModel | 31 | 8 | false | true | 2400 | 1024 | 4.84 MB |

### Individual dynamics functions


#### PrimitiveWetModel | Float32 | T31 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| - | - | - | - |
| pressure_gradient_flux! | 84.206 μs| 100.30 KiB| 789 |
| linear_virtual_temperature! | 3.604 μs| 0 bytes| 0 |
| geopotential! | 15.162 μs| 5.51 KiB| 22 |
| vertical_integration! | 14.702 μs| 0 bytes| 0 |
| surface_pressure_tendency! | 27.591 μs| 24.66 KiB| 288 |
| vertical_velocity! | 90.576 μs| 346.84 KiB| 12 |
| linear_pressure_gradient! | 3.286 μs| 0 bytes| 0 |
| vertical_advection! | 295.173 μs| 8.69 KiB| 96 |
| vordiv_tendencies! | 488.233 μs| 246.36 KiB| 721 |
| temperature_tendency! | 620.212 μs| 361.98 KiB| 1024 |
| humidity_tendency! | 589.686 μs| 360.62 KiB| 1014 |
| bernoulli_potential! | 209.745 μs| 416.55 KiB| 340 |

