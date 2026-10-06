# Benchmarks

Performance benchmarks for SpeedyWeather.jl, collected across multiple architectures. Each architecture's results live in its own section below; the overview table at the top compares the headline PrimitiveWet resolution sweep across all archs that have been benchmarked so far.

All simulations are benchmarked over several seconds (wallclock time) without output. Benchmarking excludes initialization and is started just before the main time loop and finishes right after. The benchmarking results here are not very robust; timings that change by ±50% are not uncommon. Proper benchmarking for performance optimization uses the minimum or median of many executions, while we run a simulation for several time steps which effectively represents the mean, susceptible to outliers. However, this is what a user will experience in most situations anyway and the following therefore presents a rough idea of how fast a SpeedyWeather simulation will run, and how much memory it requires.

### Explanation

Abbreviations in the tables below are as follows; omitted columns use defaults.
- NF: Number format, default: Float32
- T: Spectral resolution, maximum degree of spherical harmonics (1-based), default: T32
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
julia --project=. manual_benchmarking.jl amdgpu         # AMDGPU (HIP graphs forced on)
julia --project=. manual_benchmarking.jl reactant-cpu   # Reactant on CPU
julia --project=. manual_benchmarking.jl reactant-gpu   # Reactant on CUDA GPU
```

Each run updates only its own architecture's section in this `README.md`; results for other architectures are preserved via `benchmark_results.json`.

## Overview: PrimitiveWet resolution across architectures

Simulated years per wallclock day (SYPD) for the `PrimitiveWetModel` resolution sweep, one column per architecture. Each (T, L) configuration is reported for both the standard Legendre transform and fast Fourier transform (LT+FFT) and the single matrix transform (MT). Empty cells mean the architecture has not yet been benchmarked or that suite was skipped. Comparison figures across architectures are available on the documentation's `Benchmarks` page.

| T | L | Transform | cpu-arm | cpu-x86 | gpu-nvidia | gpu-amd |
| --- | --- | --- | --- | --- | --- | --- |
| 32 | 8 | LT+FFT | 1308 | 830 | 5621 | 3763 |
| 32 | 8 | MT | 772 | 56 | 5538 | 3696 |
| 43 | 8 | LT+FFT | 559 | 360 | 3888 | 2149 |
| 43 | 8 | MT | 253 | 15 | 3899 | 2830 |
| 64 | 8 | LT+FFT | 152 | 104 | 1181 | 828 |
| 64 | 8 | MT | 52 | 2.2 | 1179 | 826 |
| 86 | 8 | LT+FFT | 58 | 39 | 653 | 246 |
| 86 | 8 | MT | 14 | 0.5 | 306 | 245 |
| 86 | 16 | LT+FFT | 29 | 20 | 573 | 202 |
| 86 | 16 | MT | 8.8 | 0.3 | 138 | 150 |
| 86 | 24 | LT+FFT | 20 | 14 | 537 | 196 |
| 86 | 24 | MT | 6.1 | 0.2 | 76 | 108 |
| 128 | 8 | LT+FFT | 14 | 10 | 261 | 86 |
| 128 | 8 | MT | 1.2 | 0.1 | 33 | 39 |
| 128 | 16 | LT+FFT | 7.6 | 5.3 | 237 | 80 |
| 128 | 16 | MT | 1.2 | 0.0 | 14 | 23 |
| 128 | 24 | LT+FFT | 5.5 | 3.6 | 221 | 76 |
| 128 | 24 | MT | 1.0 | 0.0 | 9.3 | 17 |
| 171 | 8 | LT+FFT | 5.7 | 3.8 | 136 | 44 |
| 171 | 16 | LT+FFT | 3.1 | 2.1 | 137 | 41 |
| 171 | 24 | LT+FFT | 1.9 | 1.4 | 110 | 38 |
| 256 | 8 | LT+FFT | 1.4 | 1.0 | 53 | 16 |
| 256 | 16 | LT+FFT | 0.8 | 0.5 | 45 | 13 |
| 256 | 24 | LT+FFT | 0.5 | 0.3 | 38 | 11 |

## Architecture: `cpu-arm`

Created for SpeedyWeather.jl v0.23.0 on Tue, 06 Oct 2026 12:25:03.

### Machine details

```julia
julia> versioninfo()
Julia Version 1.12.6
Commit 15346901f00 (2026-04-09 19:20 UTC)
Build Info:
  Official https://julialang.org release
Platform Info:
  OS: macOS (arm64-apple-darwin24.0.0)
  CPU: 8 × Apple M3
  WORD_SIZE: 64
  LLVM: libLLVM-18.1.7 (ORCJIT, apple-m3)
  GC: Built with stock GC
Threads: 1 default, 1 interactive, 1 GC (on 4 virtual cores)
```


### Models, default setups

| Model | truncation | L | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| BarotropicModel | 32 | 1 | false | 1800 | 33737 | 725.53 KB |
| ShallowWaterModel | 32 | 1 | false | 2400 | 30728 | 907.89 KB |
| PrimitiveDryModel | 32 | 8 | true | 2400 | 2063 | 4.78 MB |
| PrimitiveWetModel | 32 | 8 | true | 2400 | 1275 | 5.82 MB |

### Shallow water model, resolution

| Model | truncation | L | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 32 | 1 | 48 | 2400 | 35956 | 907.89 KB |
| ShallowWaterModel | 43 | 1 | 64 | 1800 | 14935 | 1.59 MB |
| ShallowWaterModel | 64 | 1 | 96 | 1200 | 4135 | 3.57 MB |
| ShallowWaterModel | 86 | 1 | 128 | 900 | 1595 | 6.49 MB |
| ShallowWaterModel | 128 | 1 | 192 | 600 | 376 | 15.36 MB |
| ShallowWaterModel | 171 | 1 | 256 | 450 | 127 | 29.00 MB |
| ShallowWaterModel | 256 | 1 | 384 | 300 | 27 | 73.09 MB |

### Primitive wet model, resolution

| Model | truncation | L | Rings | Transform | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | 48 | default | 2400 | 1308 | 5.82 MB |
| PrimitiveWetModel | 43 | 8 | 64 | default | 1800 | 559 | 9.82 MB |
| PrimitiveWetModel | 64 | 8 | 96 | default | 1200 | 152 | 20.87 MB |
| PrimitiveWetModel | 86 | 8 | 128 | default | 900 | 58 | 36.31 MB |
| PrimitiveWetModel | 128 | 8 | 192 | default | 600 | 14 | 79.95 MB |
| PrimitiveWetModel | 171 | 8 | 256 | default | 450 | 5.7 | 141.86 MB |
| PrimitiveWetModel | 256 | 8 | 384 | default | 300 | 1.4 | 322.22 MB |
| PrimitiveWetModel | 86 | 16 | 128 | default | 900 | 29 | 62.14 MB |
| PrimitiveWetModel | 128 | 16 | 192 | default | 600 | 7.6 | 135.95 MB |
| PrimitiveWetModel | 171 | 16 | 256 | default | 450 | 3.1 | 239.79 MB |
| PrimitiveWetModel | 256 | 16 | 384 | default | 300 | 0.8 | 538.55 MB |
| PrimitiveWetModel | 86 | 24 | 128 | default | 900 | 20 | 88.01 MB |
| PrimitiveWetModel | 128 | 24 | 192 | default | 600 | 5.5 | 192.01 MB |
| PrimitiveWetModel | 171 | 24 | 256 | default | 450 | 1.9 | 337.81 MB |
| PrimitiveWetModel | 256 | 24 | 384 | default | 300 | 0.5 | 755.00 MB |
| PrimitiveWetModel | 32 | 8 | 48 | matrix | 2400 | 772 | 34.26 MB |
| PrimitiveWetModel | 43 | 8 | 64 | matrix | 1800 | 253 | 92.95 MB |
| PrimitiveWetModel | 64 | 8 | 96 | matrix | 1200 | 52 | 396.37 MB |
| PrimitiveWetModel | 86 | 8 | 128 | matrix | 900 | 14 | 1.18 GB |
| PrimitiveWetModel | 128 | 8 | 192 | matrix | 600 | 1.2 | 5.49 GB |
| PrimitiveWetModel | 86 | 16 | 128 | matrix | 900 | 8.8 | 1.21 GB |
| PrimitiveWetModel | 128 | 16 | 192 | matrix | 600 | 1.2 | 5.55 GB |
| PrimitiveWetModel | 86 | 24 | 128 | matrix | 900 | 6.1 | 1.23 GB |
| PrimitiveWetModel | 128 | 24 | 192 | matrix | 600 | 1.0 | 5.60 GB |

### Primitive Equation, Float32 vs Float64

| Model | NF | truncation | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | Float32 | 32 | 8 | 2400 | 1306 | 5.82 MB |
| PrimitiveWetModel | Float64 | 32 | 8 | 2400 | 1206 | 11.04 MB |

### Grids

| Model | truncation | L | Grid | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 64 | 8 | FullGaussianGrid | 96 | 1200 | 86 | 30.59 MB |
| PrimitiveWetModel | 64 | 8 | FullClenshawGrid | 127 | 1200 | 61 | 51.09 MB |
| PrimitiveWetModel | 64 | 8 | OctahedralGaussianGrid | 96 | 1200 | 154 | 20.87 MB |
| PrimitiveWetModel | 64 | 8 | OctahedralClenshawGrid | 127 | 1200 | 92 | 32.77 MB |
| PrimitiveWetModel | 64 | 8 | HEALPixGrid | 127 | 1200 | 95 | 24.19 MB |
| PrimitiveWetModel | 64 | 8 | OctaHEALPixGrid | 127 | 1200 | 72 | 30.06 MB |

### Number of vertical layers

| Model | truncation | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 4 | 2400 | 2226 | 3.72 MB |
| PrimitiveWetModel | 32 | 8 | 2400 | 1388 | 5.82 MB |
| PrimitiveWetModel | 32 | 12 | 2400 | 972 | 7.93 MB |
| PrimitiveWetModel | 32 | 16 | 2400 | 754 | 10.04 MB |

### PrimitiveDryModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveDryModel | 32 | 8 | true | true | 2400 | 2363 | 4.78 MB |
| PrimitiveDryModel | 32 | 8 | true | false | 2400 | 3957 | 4.78 MB |
| PrimitiveDryModel | 32 | 8 | false | true | 2400 | 2631 | 4.78 MB |

### PrimitiveWetModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | true | true | 2400 | 1371 | 5.82 MB |
| PrimitiveWetModel | 32 | 8 | true | false | 2400 | 3048 | 5.82 MB |
| PrimitiveWetModel | 32 | 8 | false | true | 2400 | 1531 | 5.82 MB |

### Individual dynamics functions


#### PrimitiveWetModel | Float32 | T31 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| --- | --- | --- | --- |
| pressure_gradient_flux! | 41.458 μs| 31.98 KiB| 200 |
| linear_virtual_temperature! | 2.069 μs| 0 bytes| 0 |
| geopotential! | 7.521 μs| 384 bytes| 6 |
| vertical_integration! | 14.583 μs| 0 bytes| 0 |
| surface_pressure_tendency! | 11.625 μs| 15.66 KiB| 96 |
| vertical_velocity! | 23.333 μs| 0 bytes| 0 |
| linear_pressure_gradient! | 2.079 μs| 0 bytes| 0 |
| vertical_advection! | 110.167 μs| 2.44 KiB| 32 |
| vordiv_tendencies! | 224.584 μs| 231.83 KiB| 284 |
| temperature_tendency! | 291.333 μs| 344.77 KiB| 401 |
| humidity_tendency! | 280.000 μs| 344.09 KiB| 396 |
| bernoulli_potential! | 91.833 μs| 114.28 KiB| 129 |

## Architecture: `cpu-x86`

Created for SpeedyWeather.jl v0.23.0 on Fri, 02 Oct 2026 12:26:22.

### Machine details

```julia
julia> versioninfo()
Julia Version 1.12.2
Commit ca9b6662be4 (2025-11-20 16:25 UTC)
Build Info:
  Official https://julialang.org release
Platform Info:
  OS: Linux (x86_64-linux-gnu)
  CPU: 128 × AMD EPYC 9554 64-Core Processor
  WORD_SIZE: 64
  LLVM: libLLVM-18.1.7 (ORCJIT, znver4)
  GC: Built with stock GC
Threads: 1 default, 1 interactive, 1 GC (on 128 virtual cores)
Environment:
  LD_LIBRARY_PATH = /usr/local/lib:/usr/local/lib:
```


### Models, default setups

| Model | truncation | L | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| BarotropicModel | 32 | 1 | false | 1800 | 23839 | 725.64 KB |
| ShallowWaterModel | 32 | 1 | false | 2400 | 18452 | 908.00 KB |
| PrimitiveDryModel | 32 | 8 | true | 2400 | 1381 | 4.78 MB |
| PrimitiveWetModel | 32 | 8 | true | 2400 | 839 | 5.77 MB |

### Shallow water model, resolution

| Model | truncation | L | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 32 | 1 | 48 | 2400 | 17815 | 908.00 KB |
| ShallowWaterModel | 43 | 1 | 64 | 1800 | 8191 | 1.59 MB |
| ShallowWaterModel | 64 | 1 | 96 | 1200 | 2153 | 3.57 MB |
| ShallowWaterModel | 86 | 1 | 128 | 900 | 905 | 6.49 MB |
| ShallowWaterModel | 128 | 1 | 192 | 600 | 229 | 15.36 MB |
| ShallowWaterModel | 171 | 1 | 256 | 450 | 84 | 29.00 MB |
| ShallowWaterModel | 256 | 1 | 384 | 300 | 18 | 73.09 MB |

### Primitive wet model, resolution

| Model | truncation | L | Rings | Transform | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | 48 | default | 2400 | 830 | 5.77 MB |
| PrimitiveWetModel | 43 | 8 | 64 | default | 1800 | 360 | 9.74 MB |
| PrimitiveWetModel | 64 | 8 | 96 | default | 1200 | 104 | 20.70 MB |
| PrimitiveWetModel | 86 | 8 | 128 | default | 900 | 39 | 36.01 MB |
| PrimitiveWetModel | 128 | 8 | 192 | default | 600 | 10 | 79.31 MB |
| PrimitiveWetModel | 171 | 8 | 256 | default | 450 | 3.8 | 140.74 MB |
| PrimitiveWetModel | 256 | 8 | 384 | default | 300 | 1.0 | 319.75 MB |
| PrimitiveWetModel | 86 | 16 | 128 | default | 900 | 20 | 61.84 MB |
| PrimitiveWetModel | 128 | 16 | 192 | default | 600 | 5.3 | 135.30 MB |
| PrimitiveWetModel | 171 | 16 | 256 | default | 450 | 2.1 | 238.67 MB |
| PrimitiveWetModel | 256 | 16 | 384 | default | 300 | 0.5 | 536.08 MB |
| PrimitiveWetModel | 86 | 24 | 128 | default | 900 | 14 | 87.71 MB |
| PrimitiveWetModel | 128 | 24 | 192 | default | 600 | 3.6 | 191.37 MB |
| PrimitiveWetModel | 171 | 24 | 256 | default | 450 | 1.4 | 336.69 MB |
| PrimitiveWetModel | 256 | 24 | 384 | default | 300 | 0.3 | 752.53 MB |
| PrimitiveWetModel | 32 | 8 | 48 | matrix | 2400 | 56 | 48.15 MB |
| PrimitiveWetModel | 43 | 8 | 64 | matrix | 1800 | 15 | 133.94 MB |
| PrimitiveWetModel | 64 | 8 | 96 | matrix | 1200 | 2.2 | 582.94 MB |
| PrimitiveWetModel | 86 | 8 | 128 | matrix | 900 | 0.5 | 1.75 GB |
| PrimitiveWetModel | 128 | 8 | 192 | matrix | 600 | 0.1 | 8.19 GB |
| PrimitiveWetModel | 86 | 16 | 128 | matrix | 900 | 0.3 | 1.78 GB |
| PrimitiveWetModel | 128 | 16 | 192 | matrix | 600 | 0.0 | 8.24 GB |
| PrimitiveWetModel | 86 | 24 | 128 | matrix | 900 | 0.2 | 1.80 GB |
| PrimitiveWetModel | 128 | 24 | 192 | matrix | 600 | 0.0 | 8.30 GB |

### Primitive Equation, Float32 vs Float64

| Model | NF | truncation | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | Float32 | 32 | 8 | 2400 | 801 | 5.77 MB |
| PrimitiveWetModel | Float64 | 32 | 8 | 2400 | 748 | 10.93 MB |

### Grids

| Model | truncation | L | Grid | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 64 | 8 | FullGaussianGrid | 96 | 1200 | 73 | 30.30 MB |
| PrimitiveWetModel | 64 | 8 | FullClenshawGrid | 127 | 1200 | 41 | 50.57 MB |
| PrimitiveWetModel | 64 | 8 | OctahedralGaussianGrid | 96 | 1200 | 101 | 20.70 MB |
| PrimitiveWetModel | 64 | 8 | OctahedralClenshawGrid | 127 | 1200 | 58 | 32.48 MB |
| PrimitiveWetModel | 64 | 8 | HEALPixGrid | 127 | 1200 | 79 | 23.99 MB |
| PrimitiveWetModel | 64 | 8 | OctaHEALPixGrid | 127 | 1200 | 59 | 29.79 MB |

### Number of vertical layers

| Model | truncation | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 4 | 2400 | 1358 | 3.67 MB |
| PrimitiveWetModel | 32 | 8 | 2400 | 818 | 5.77 MB |
| PrimitiveWetModel | 32 | 12 | 2400 | 606 | 7.88 MB |
| PrimitiveWetModel | 32 | 16 | 2400 | 403 | 9.99 MB |

### PrimitiveDryModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveDryModel | 32 | 8 | true | true | 2400 | 1335 | 4.78 MB |
| PrimitiveDryModel | 32 | 8 | true | false | 2400 | 2177 | 4.78 MB |
| PrimitiveDryModel | 32 | 8 | false | true | 2400 | 1526 | 4.78 MB |

### PrimitiveWetModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | true | true | 2400 | 825 | 5.77 MB |
| PrimitiveWetModel | 32 | 8 | true | false | 2400 | 1624 | 5.77 MB |
| PrimitiveWetModel | 32 | 8 | false | true | 2400 | 901 | 5.77 MB |

### Individual dynamics functions


#### PrimitiveWetModel | Float32 | T31 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| --- | --- | --- | --- |
| pressure_gradient_flux! | 68.563 μs| 31.98 KiB| 200 |
| linear_virtual_temperature! | 3.549 μs| 0 bytes| 0 |
| geopotential! | 9.981 μs| 384 bytes| 6 |
| vertical_integration! | 13.550 μs| 0 bytes| 0 |
| surface_pressure_tendency! | 18.371 μs| 15.66 KiB| 96 |
| vertical_velocity! | 59.283 μs| 0 bytes| 0 |
| linear_pressure_gradient! | 3.085 μs| 0 bytes| 0 |
| vertical_advection! | 160.917 μs| 2.44 KiB| 32 |
| vordiv_tendencies! | 403.337 μs| 218.48 KiB| 284 |
| temperature_tendency! | 525.762 μs| 324.75 KiB| 401 |
| humidity_tendency! | 500.931 μs| 324.08 KiB| 396 |
| bernoulli_potential! | 167.057 μs| 107.61 KiB| 129 |

## Architecture: `gpu-nvidia`

Created for SpeedyWeather.jl v0.23.0 on Thu, 01 Oct 2026 17:21:36.

### Machine details

```julia
julia> versioninfo()
Julia Version 1.12.2
Commit ca9b6662be4 (2025-11-20 16:25 UTC)
Build Info:
  Official https://julialang.org release
Platform Info:
  OS: Linux (x86_64-linux-gnu)
  CPU: 128 × AMD EPYC 9554 64-Core Processor
  WORD_SIZE: 64
  LLVM: libLLVM-18.1.7 (ORCJIT, znver4)
  GC: Built with stock GC
Threads: 1 default, 1 interactive, 1 GC (on 128 virtual cores)
Environment:
  LD_LIBRARY_PATH = /usr/local/lib:/usr/local/lib:
```

```julia
julia> CUDA.versioninfo()
CUDA toolchain: 
- runtime 13.3.0, artifact installation
- driver 580.126.9 for 13.3
- compiler 13.3.33, artifact installation

CUDA libraries: 
- cuBLAS: 13.6.0
- cuSPARSE: 12.8.2
- cuSOLVER: 12.2.6
- cuFFT: 12.3.0
- cuRAND: 10.4.3
- CUPTI: 2026.2.1 (API 13.3.1)
- NVML: 13.0.0+580.126.9

Julia packages: 
- CUDACore: 6.2.2
- GPUArrays: 11.5.13
- GPUCompiler: 1.23.0
- KernelAbstractions: 0.9.42
- CUDA_Driver_jll: 13.3.1+0
- CUDA_Compiler_jll: 0.4.4+1
- CUDA_Runtime_jll: 0.23.0+1
- NVPTX_LLVM_Backend_jll: 22.1.7+1

Toolchain:
- Julia: 1.12.2
- LLVM: 18.1.7

1 device:
  0: NVIDIA H100 80GB HBM3 (sm_90, 77.878 GiB / 79.647 GiB available)
     compiles to sm_90a / PTX 9.3 (LLVM: sm_90a / PTX 9.0)
```


### Models, default setups

| Model | truncation | L | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| BarotropicModel | 32 | 1 | false | 1800 | 8668 | 429.08 KB |
| ShallowWaterModel | 32 | 1 | false | 2400 | 9067 | 432.66 KB |
| PrimitiveDryModel | 32 | 8 | true | 2400 | 6623 | 576.95 KB |
| PrimitiveWetModel | 32 | 8 | true | 2400 | 5628 | 582.49 KB |

### Shallow water model, resolution

| Model | truncation | L | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 32 | 1 | 48 | 2400 | 9097 | 432.66 KB |
| ShallowWaterModel | 43 | 1 | 64 | 1800 | 2962 | 745.24 KB |
| ShallowWaterModel | 64 | 1 | 96 | 1200 | 943 | 1.63 MB |
| ShallowWaterModel | 86 | 1 | 128 | 900 | 688 | 3.07 MB |
| ShallowWaterModel | 128 | 1 | 192 | 600 | 293 | 6.70 MB |
| ShallowWaterModel | 171 | 1 | 256 | 450 | 129 | 11.74 MB |
| ShallowWaterModel | 256 | 1 | 384 | 300 | 51 | 26.06 MB |

### Primitive wet model, resolution

| Model | truncation | L | Rings | Transform | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | 48 | default | 2400 | 5621 | 582.49 KB |
| PrimitiveWetModel | 43 | 8 | 64 | default | 1800 | 3888 | 995.42 KB |
| PrimitiveWetModel | 64 | 8 | 96 | default | 1200 | 1181 | 2.17 MB |
| PrimitiveWetModel | 86 | 8 | 128 | default | 900 | 653 | 4.07 MB |
| PrimitiveWetModel | 128 | 8 | 192 | default | 600 | 261 | 8.88 MB |
| PrimitiveWetModel | 171 | 8 | 256 | default | 450 | 136 | 15.56 MB |
| PrimitiveWetModel | 256 | 8 | 384 | default | 300 | 53 | 34.53 MB |
| PrimitiveWetModel | 86 | 16 | 128 | default | 900 | 573 | 5.12 MB |
| PrimitiveWetModel | 128 | 16 | 192 | default | 600 | 237 | 11.24 MB |
| PrimitiveWetModel | 171 | 16 | 256 | default | 450 | 137 | 19.76 MB |
| PrimitiveWetModel | 256 | 16 | 384 | default | 300 | 45 | 43.96 MB |
| PrimitiveWetModel | 86 | 24 | 128 | default | 900 | 537 | 6.17 MB |
| PrimitiveWetModel | 128 | 24 | 192 | default | 600 | 221 | 13.60 MB |
| PrimitiveWetModel | 171 | 24 | 256 | default | 450 | 110 | 23.95 MB |
| PrimitiveWetModel | 256 | 24 | 384 | default | 300 | 38 | 53.40 MB |
| PrimitiveWetModel | 32 | 8 | 48 | matrix | 2400 | 5538 | 582.49 KB |
| PrimitiveWetModel | 43 | 8 | 64 | matrix | 1800 | 3899 | 995.42 KB |
| PrimitiveWetModel | 64 | 8 | 96 | matrix | 1200 | 1179 | 2.17 MB |
| PrimitiveWetModel | 86 | 8 | 128 | matrix | 900 | 306 | 3.81 MB |
| PrimitiveWetModel | 128 | 8 | 192 | matrix | 600 | 33 | 8.50 MB |
| PrimitiveWetModel | 86 | 16 | 128 | matrix | 900 | 138 | 4.86 MB |
| PrimitiveWetModel | 128 | 16 | 192 | matrix | 600 | 14 | 10.86 MB |
| PrimitiveWetModel | 86 | 24 | 128 | matrix | 900 | 76 | 5.91 MB |
| PrimitiveWetModel | 128 | 24 | 192 | matrix | 600 | 9.3 | 13.22 MB |

### Primitive Equation, Float32 vs Float64

| Model | NF | truncation | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | Float32 | 32 | 8 | 2400 | 5609 | 582.49 KB |
| PrimitiveWetModel | Float64 | 32 | 8 | 2400 | 4578 | 583.08 KB |

### Grids

| Model | truncation | L | Grid | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 64 | 8 | FullGaussianGrid | 96 | 1200 | 759 | 2.26 MB |
| PrimitiveWetModel | 64 | 8 | FullClenshawGrid | 127 | 1200 | 481 | 3.95 MB |
| PrimitiveWetModel | 64 | 8 | OctahedralGaussianGrid | 96 | 1200 | 1184 | 2.17 MB |
| PrimitiveWetModel | 64 | 8 | OctahedralClenshawGrid | 127 | 1200 | 761 | 3.78 MB |
| PrimitiveWetModel | 64 | 8 | HEALPixGrid | 127 | 1200 | 1080 | 3.71 MB |
| PrimitiveWetModel | 64 | 8 | OctaHEALPixGrid | 127 | 1200 | 865 | 3.76 MB |

### Number of vertical layers

| Model | truncation | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 4 | 2400 | 5641 | 508.76 KB |
| PrimitiveWetModel | 32 | 8 | 2400 | 6013 | 582.49 KB |
| PrimitiveWetModel | 32 | 12 | 2400 | 5955 | 656.22 KB |
| PrimitiveWetModel | 32 | 16 | 2400 | 5790 | 729.95 KB |

### PrimitiveDryModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveDryModel | 32 | 8 | true | true | 2400 | 7079 | 576.95 KB |
| PrimitiveDryModel | 32 | 8 | true | false | 2400 | 8582 | 576.95 KB |
| PrimitiveDryModel | 32 | 8 | false | true | 2400 | 11047 | 576.95 KB |

### PrimitiveWetModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | true | true | 2400 | 6044 | 582.49 KB |
| PrimitiveWetModel | 32 | 8 | true | false | 2400 | 7835 | 582.49 KB |
| PrimitiveWetModel | 32 | 8 | false | true | 2400 | 9725 | 582.49 KB |

### Individual dynamics functions


#### PrimitiveWetModel | Float32 | T31 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| --- | --- | --- | --- |
| pressure_gradient_flux! | 75.900 μs| 14.42 KiB| 348 |
| linear_virtual_temperature! | 12.789 μs| 2.77 KiB| 52 |
| geopotential! | 19.100 μs| 3.73 KiB| 101 |
| vertical_integration! | 27.250 μs| 7.17 KiB| 141 |
| surface_pressure_tendency! | N/A| N/A| N/A |
| vertical_velocity! | 24.250 μs| 8.72 KiB| 192 |
| linear_pressure_gradient! | 12.540 μs| 2.34 KiB| 50 |
| vertical_advection! | 24.160 μs| 9.09 KiB| 142 |
| vordiv_tendencies! | N/A| N/A| N/A |
| temperature_tendency! | N/A| N/A| N/A |
| humidity_tendency! | N/A| N/A| N/A |
| bernoulli_potential! | N/A| N/A| N/A |

## Architecture: `gpu-amd`

Created for SpeedyWeather.jl v0.22.1 on Tue, 29 Sep 2026 16:46:03.

### Machine details

```julia
julia> versioninfo()
Julia Version 1.12.7
Commit 6d172b025e4 (2026-08-15 08:05 UTC)
Build Info:
  Official https://julialang.org release
Platform Info:
  OS: Linux (x86_64-linux-gnu)
  CPU: 128 × AMD EPYC 7A53 64-Core Processor
  WORD_SIZE: 64
  LLVM: libLLVM-18.1.7 (ORCJIT, znver3)
  GC: Built with stock GC
Threads: 1 default, 1 interactive, 1 GC (on 128 virtual cores)
Environment:
  LD_LIBRARY_PATH = /opt/cray/pe/papi/7.2.0.1/lib64:/opt/cray/libfabric/1.22.0/lib64
  JULIA_DEPOT_PATH = /projappl/project_462000008/decristoforo/.julia:
```

```julia
julia> AMDGPU.versioninfo()
AMDGPU versioninfo
(AMDGPU.versioninfo() failed: ErrorException("could not load symbol \"hiptensorGetVersion\":\n/opt/rocm-6.3.4/lib/libhiptensor.so: undefined symbol: hiptensorGetVersion"))
```


### Models, default setups

| Model | truncation | L | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| BarotropicModel | 32 | 1 | false | 1800 | 9483 | 436.34 KB |
| ShallowWaterModel | 32 | 1 | false | 2400 | 7389 | 440.48 KB |
| PrimitiveDryModel | 32 | 8 | true | 2400 | 4148 | 590.11 KB |
| PrimitiveWetModel | 32 | 8 | true | 2400 | 3790 | 596.78 KB |

### Shallow water model, resolution

| Model | truncation | L | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 32 | 1 | 48 | 2400 | 8214 | 440.48 KB |
| ShallowWaterModel | 43 | 1 | 64 | 1800 | 4029 | 753.06 KB |
| ShallowWaterModel | 64 | 1 | 96 | 1200 | 1038 | 1.64 MB |
| ShallowWaterModel | 86 | 1 | 128 | 900 | 84 | 3.26 MB |
| ShallowWaterModel | 128 | 1 | 192 | 600 | 37 | 6.98 MB |
| ShallowWaterModel | 171 | 1 | 256 | 450 | 20 | 12.12 MB |
| ShallowWaterModel | 256 | 1 | 384 | 300 | 9.0 | 26.62 MB |

### Primitive wet model, resolution

| Model | truncation | L | Rings | Transform | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | 48 | default | 2400 | 3763 | 596.78 KB |
| PrimitiveWetModel | 43 | 8 | 64 | default | 1800 | 2149 | 1.01 MB |
| PrimitiveWetModel | 64 | 8 | 96 | default | 1200 | 828 | 2.19 MB |
| PrimitiveWetModel | 86 | 8 | 128 | default | 900 | 246 | 4.35 MB |
| PrimitiveWetModel | 128 | 8 | 192 | default | 600 | 86 | 9.30 MB |
| PrimitiveWetModel | 171 | 8 | 256 | default | 450 | 44 | 16.12 MB |
| PrimitiveWetModel | 256 | 8 | 384 | default | 300 | 16 | 35.35 MB |
| PrimitiveWetModel | 86 | 16 | 128 | default | 900 | 202 | 5.40 MB |
| PrimitiveWetModel | 128 | 16 | 192 | default | 600 | 80 | 11.66 MB |
| PrimitiveWetModel | 171 | 16 | 256 | default | 450 | 41 | 20.31 MB |
| PrimitiveWetModel | 256 | 16 | 384 | default | 300 | 13 | 44.79 MB |
| PrimitiveWetModel | 86 | 24 | 128 | default | 900 | 196 | 6.45 MB |
| PrimitiveWetModel | 128 | 24 | 192 | default | 600 | 76 | 14.02 MB |
| PrimitiveWetModel | 171 | 24 | 256 | default | 450 | 38 | 24.50 MB |
| PrimitiveWetModel | 256 | 24 | 384 | default | 300 | 11 | 54.22 MB |
| PrimitiveWetModel | 32 | 8 | 48 | matrix | 2400 | 3696 | 596.78 KB |
| PrimitiveWetModel | 43 | 8 | 64 | matrix | 1800 | 2830 | 1.01 MB |
| PrimitiveWetModel | 64 | 8 | 96 | matrix | 1200 | 826 | 2.19 MB |
| PrimitiveWetModel | 86 | 8 | 128 | matrix | 900 | 245 | 3.83 MB |
| PrimitiveWetModel | 128 | 8 | 192 | matrix | 600 | 39 | 8.52 MB |
| PrimitiveWetModel | 86 | 16 | 128 | matrix | 900 | 150 | 4.88 MB |
| PrimitiveWetModel | 128 | 16 | 192 | matrix | 600 | 23 | 10.87 MB |
| PrimitiveWetModel | 86 | 24 | 128 | matrix | 900 | 108 | 5.93 MB |
| PrimitiveWetModel | 128 | 24 | 192 | matrix | 600 | 17 | 13.23 MB |

### Primitive Equation, Float32 vs Float64

| Model | NF | truncation | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | Float32 | 32 | 8 | 2400 | 3701 | 596.78 KB |
| PrimitiveWetModel | Float64 | 32 | 8 | 2400 | 3101 | 597.37 KB |

### Grids

| Model | truncation | L | Grid | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 64 | 8 | FullGaussianGrid | 96 | 1200 | 515 | 2.28 MB |
| PrimitiveWetModel | 64 | 8 | FullClenshawGrid | 127 | 1200 | 325 | 3.97 MB |
| PrimitiveWetModel | 64 | 8 | OctahedralGaussianGrid | 96 | 1200 | 826 | 2.19 MB |
| PrimitiveWetModel | 64 | 8 | OctahedralClenshawGrid | 127 | 1200 | 509 | 3.80 MB |
| PrimitiveWetModel | 64 | 8 | HEALPixGrid | 127 | 1200 | 742 | 3.72 MB |
| PrimitiveWetModel | 64 | 8 | OctaHEALPixGrid | 127 | 1200 | 566 | 3.77 MB |

### Number of vertical layers

| Model | truncation | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 4 | 2400 | 3890 | 523.05 KB |
| PrimitiveWetModel | 32 | 8 | 2400 | 3935 | 596.78 KB |
| PrimitiveWetModel | 32 | 12 | 2400 | 3801 | 670.51 KB |
| PrimitiveWetModel | 32 | 16 | 2400 | 3653 | 744.24 KB |

### PrimitiveDryModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveDryModel | 32 | 8 | true | true | 2400 | 4635 | 590.11 KB |
| PrimitiveDryModel | 32 | 8 | true | false | 2400 | 5238 | 590.11 KB |
| PrimitiveDryModel | 32 | 8 | false | true | 2400 | 8782 | 590.11 KB |

### PrimitiveWetModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | true | true | 2400 | 3830 | 596.78 KB |
| PrimitiveWetModel | 32 | 8 | true | false | 2400 | 4740 | 596.78 KB |
| PrimitiveWetModel | 32 | 8 | false | true | 2400 | 7114 | 596.78 KB |

### Individual dynamics functions


#### PrimitiveWetModel | Float32 | T31 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| --- | --- | --- | --- |
| pressure_gradient_flux! | 127.728 μs| 19.23 KiB| 421 |
| linear_virtual_temperature! | 26.051 μs| 3.64 KiB| 66 |
| geopotential! | 38.305 μs| 5.06 KiB| 129 |
| vertical_integration! | 55.427 μs| 9.19 KiB| 188 |
| surface_pressure_tendency! | N/A| N/A| N/A |
| vertical_velocity! | 48.564 μs| 10.22 KiB| 211 |
| linear_pressure_gradient! | 25.540 μs| 3.27 KiB| 66 |
| vertical_advection! | 52.752 μs| 13.06 KiB| 224 |
| vordiv_tendencies! | N/A| N/A| N/A |
| temperature_tendency! | N/A| N/A| N/A |
| humidity_tendency! | N/A| N/A| N/A |
| bernoulli_potential! | N/A| N/A| N/A |

