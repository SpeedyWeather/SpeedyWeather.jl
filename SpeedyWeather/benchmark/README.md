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

Simulated years per wallclock day (SYPD) for the `PrimitiveWetModel` resolution sweep, one column per architecture. Each (T, L) configuration is reported for both the standard Legendre transform and fast Fourier transform (LT+FFT) and the single matrix transform (MT). Empty cells mean the architecture has not yet been benchmarked or that suite was skipped. Comparison figures across architectures are available on the documentation's `Benchmarks` page.

| T | L | Transform | cpu-arm | cpu-x86 | gpu-nvidia |
| --- | --- | --- | --- | --- | --- |
| 31 | 8 | LT+FFT | 1326 | 864 | 2411 |
| 31 | 8 | MT | 678 | 107 | 4660 |
| 42 | 8 | LT+FFT | 556 | 372 | 1443 |
| 42 | 8 | MT | 157 | 28 | 2358 |
| 63 | 8 | LT+FFT | 146 | 107 | 665 |
| 63 | 8 | MT | 33 | 3.9 | 995 |
| 85 | 8 | LT+FFT | 63 | 42 | 397 |
| 85 | 8 | MT | 7.4 | 1.0 | 267 |
| 85 | 16 | LT+FFT | 31 | 50 | 1196 |
| 85 | 16 | MT | 5.1 | 1.1 | 576 |
| 85 | 24 | LT+FFT | 15 | 39 | 1210 |
| 85 | 24 | MT | 4.1 | 0.9 | 371 |
| 127 | 8 | LT+FFT | 14 | 11 | 157 |
| 127 | 8 | MT | 0.5 | 0.1 | 57 |
| 127 | 16 | LT+FFT | 7.9 | 15 | 381 |
| 127 | 16 | MT | 0.4 | 0.2 | 73 |
| 127 | 24 | LT+FFT | 4.8 | 11 | 374 |
| 127 | 24 | MT | 0.4 | 0.1 | 54 |
| 170 | 8 | LT+FFT | 5.3 | 4.0 | 72 |
| 170 | 16 | LT+FFT | 2.9 | 5.7 | 144 |
| 170 | 24 | LT+FFT | 1.7 | 4.5 | 140 |
| 255 | 8 | LT+FFT | 1.4 | 1.0 | 23 |
| 255 | 16 | LT+FFT | 0.7 | 1.3 | 40 |
| 255 | 24 | LT+FFT | 0.5 | 1.0 | 36 |

## Architecture: `cpu-arm`

Created for SpeedyWeather.jl v0.21.1+DEV on Tue, 30 Jun 2026 18:02:48.

### Machine details

```julia
julia> versioninfo()
Julia Version 1.11.9
Commit 53a02c0720c (2026-02-06 00:27 UTC)
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
| BarotropicModel | 31 | 1 | false | 1800 | 35318 | 746.37 KB |
| ShallowWaterModel | 31 | 1 | false | 2400 | 24399 | 880.97 KB |
| PrimitiveDryModel | 31 | 8 | true | 2400 | 2247 | 4.70 MB |
| PrimitiveWetModel | 31 | 8 | true | 2400 | 1212 | 5.38 MB |

### Shallow water model, resolution

| Model | T | L | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 31 | 1 | 48 | 2400 | 25538 | 880.97 KB |
| ShallowWaterModel | 42 | 1 | 64 | 1800 | 11692 | 1.54 MB |
| ShallowWaterModel | 63 | 1 | 96 | 1200 | 3663 | 3.47 MB |
| ShallowWaterModel | 85 | 1 | 128 | 900 | 1411 | 6.32 MB |
| ShallowWaterModel | 127 | 1 | 192 | 600 | 294 | 14.97 MB |
| ShallowWaterModel | 170 | 1 | 256 | 450 | 137 | 28.31 MB |
| ShallowWaterModel | 255 | 1 | 384 | 300 | 27 | 71.53 MB |

### Primitive wet model, resolution

| Model | T | L | Rings | Transform | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 8 | 48 | default | 2400 | 1326 | 5.38 MB |
| PrimitiveWetModel | 42 | 8 | 64 | default | 1800 | 556 | 9.10 MB |
| PrimitiveWetModel | 63 | 8 | 96 | default | 1200 | 146 | 19.35 MB |
| PrimitiveWetModel | 85 | 8 | 128 | default | 900 | 63 | 33.68 MB |
| PrimitiveWetModel | 127 | 8 | 192 | default | 600 | 14 | 74.25 MB |
| PrimitiveWetModel | 170 | 8 | 256 | default | 450 | 5.3 | 131.86 MB |
| PrimitiveWetModel | 255 | 8 | 384 | default | 300 | 1.4 | 300.09 MB |
| PrimitiveWetModel | 85 | 16 | 128 | default | 900 | 31 | 57.92 MB |
| PrimitiveWetModel | 127 | 16 | 192 | default | 600 | 7.9 | 126.80 MB |
| PrimitiveWetModel | 170 | 16 | 256 | default | 450 | 2.9 | 223.74 MB |
| PrimitiveWetModel | 255 | 16 | 384 | default | 300 | 0.7 | 502.99 MB |
| PrimitiveWetModel | 85 | 24 | 128 | default | 900 | 15 | 82.21 MB |
| PrimitiveWetModel | 127 | 24 | 192 | default | 600 | 4.8 | 179.42 MB |
| PrimitiveWetModel | 170 | 24 | 256 | default | 450 | 1.7 | 315.70 MB |
| PrimitiveWetModel | 255 | 24 | 384 | default | 300 | 0.5 | 706.03 MB |
| PrimitiveWetModel | 31 | 8 | 48 | matrix | 2400 | 678 | 47.28 MB |
| PrimitiveWetModel | 42 | 8 | 64 | matrix | 1800 | 157 | 132.47 MB |
| PrimitiveWetModel | 63 | 8 | 96 | matrix | 1200 | 33 | 579.81 MB |
| PrimitiveWetModel | 85 | 8 | 128 | matrix | 900 | 7.4 | 1.75 GB |
| PrimitiveWetModel | 127 | 8 | 192 | matrix | 600 | 0.5 | 8.18 GB |
| PrimitiveWetModel | 85 | 16 | 128 | matrix | 900 | 5.1 | 1.77 GB |
| PrimitiveWetModel | 127 | 16 | 192 | matrix | 600 | 0.4 | 8.22 GB |
| PrimitiveWetModel | 85 | 24 | 128 | matrix | 900 | 4.1 | 1.79 GB |
| PrimitiveWetModel | 127 | 24 | 192 | matrix | 600 | 0.4 | 8.26 GB |

### Primitive Equation, Float32 vs Float64

| Model | NF | T | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | Float32 | 31 | 8 | 2400 | 677 | 5.38 MB |
| PrimitiveWetModel | Float64 | 31 | 8 | 2400 | 1055 | 9.69 MB |

### Grids

| Model | T | L | Grid | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 63 | 8 | FullGaussianGrid | 96 | 1200 | 86 | 28.09 MB |
| PrimitiveWetModel | 63 | 8 | FullClenshawGrid | 95 | 1200 | 115 | 27.86 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralGaussianGrid | 96 | 1200 | 145 | 19.35 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralClenshawGrid | 95 | 1200 | 144 | 19.11 MB |
| PrimitiveWetModel | 63 | 8 | HEALPixGrid | 95 | 1200 | 266 | 14.14 MB |
| PrimitiveWetModel | 63 | 8 | OctaHEALPixGrid | 95 | 1200 | 181 | 17.21 MB |

### Number of vertical layers

| Model | T | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 4 | 2400 | 2184 | 3.41 MB |
| PrimitiveWetModel | 31 | 8 | 2400 | 1500 | 5.38 MB |
| PrimitiveWetModel | 31 | 12 | 2400 | 973 | 7.36 MB |
| PrimitiveWetModel | 31 | 16 | 2400 | 803 | 9.35 MB |

### PrimitiveDryModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveDryModel | 31 | 8 | true | true | 2400 | 2342 | 4.70 MB |
| PrimitiveDryModel | 31 | 8 | true | false | 2400 | 3058 | 4.70 MB |
| PrimitiveDryModel | 31 | 8 | false | true | 2400 | 3352 | 4.70 MB |

### PrimitiveWetModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 8 | true | true | 2400 | 1482 | 5.38 MB |
| PrimitiveWetModel | 31 | 8 | true | false | 2400 | 3360 | 5.38 MB |
| PrimitiveWetModel | 31 | 8 | false | true | 2400 | 1913 | 5.38 MB |

### Individual dynamics functions


#### PrimitiveWetModel | Float32 | T31 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| --- | --- | --- | --- |
| pressure_gradient_flux! | 38.750 μs| 49.06 KiB| 594 |
| linear_virtual_temperature! | 1.754 μs| 0 bytes| 0 |
| geopotential! | 5.868 μs| 800 bytes| 14 |
| vertical_integration! | 13.666 μs| 0 bytes| 0 |
| surface_pressure_tendency! | 14.375 μs| 25.41 KiB| 288 |
| vertical_velocity! | 21.541 μs| 0 bytes| 0 |
| linear_pressure_gradient! | 1.750 μs| 0 bytes| 0 |
| vertical_advection! | 115.834 μs| 3.88 KiB| 64 |
| vordiv_tendencies! | 229.208 μs| 260.75 KiB| 714 |
| temperature_tendency! | 300.292 μs| 382.69 KiB| 1015 |
| humidity_tendency! | 287.875 μs| 381.44 KiB| 1005 |
| bernoulli_potential! | 90.042 μs| 126.03 KiB| 328 |

## Architecture: `cpu-x86`

Created for SpeedyWeather.jl v0.21.1+DEV on Tue, 21 Jul 2026 12:52:19.

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

| Model | T | L | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| BarotropicModel | 31 | 1 | false | 1800 | 29472 | 780.69 KB |
| ShallowWaterModel | 31 | 1 | false | 2400 | 35148 | 962.97 KB |
| PrimitiveDryModel | 31 | 8 | true | 2400 | 1472 | 5.27 MB |
| PrimitiveWetModel | 31 | 8 | true | 2400 | 865 | 6.22 MB |

### Shallow water model, resolution

| Model | T | L | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 31 | 1 | 48 | 2400 | 39280 | 962.97 KB |
| ShallowWaterModel | 42 | 1 | 64 | 1800 | 9698 | 1.68 MB |
| ShallowWaterModel | 63 | 1 | 96 | 1200 | 2437 | 3.77 MB |
| ShallowWaterModel | 85 | 1 | 128 | 900 | 971 | 6.84 MB |
| ShallowWaterModel | 127 | 1 | 192 | 600 | 234 | 16.12 MB |
| ShallowWaterModel | 170 | 1 | 256 | 450 | 81 | 30.33 MB |
| ShallowWaterModel | 255 | 1 | 384 | 300 | 19 | 76.02 MB |

### Primitive wet model, resolution

| Model | T | L | Rings | Transform | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 8 | 48 | default | 2400 | 864 | 6.22 MB |
| PrimitiveWetModel | 42 | 8 | 64 | default | 1800 | 372 | 10.51 MB |
| PrimitiveWetModel | 63 | 8 | 96 | default | 1200 | 107 | 22.34 MB |
| PrimitiveWetModel | 85 | 8 | 128 | default | 900 | 42 | 38.87 MB |
| PrimitiveWetModel | 127 | 8 | 192 | default | 600 | 11 | 85.50 MB |
| PrimitiveWetModel | 170 | 8 | 256 | default | 450 | 4.0 | 151.57 MB |
| PrimitiveWetModel | 255 | 8 | 384 | default | 300 | 1.0 | 343.69 MB |
| PrimitiveWetModel | 85 | 16 | 128 | default | 900 | 50 | 67.81 MB |
| PrimitiveWetModel | 127 | 16 | 192 | default | 600 | 15 | 148.26 MB |
| PrimitiveWetModel | 170 | 16 | 256 | default | 450 | 5.7 | 261.34 MB |
| PrimitiveWetModel | 255 | 16 | 384 | default | 300 | 1.3 | 586.14 MB |
| PrimitiveWetModel | 85 | 24 | 128 | default | 900 | 39 | 96.80 MB |
| PrimitiveWetModel | 127 | 24 | 192 | default | 600 | 11 | 211.09 MB |
| PrimitiveWetModel | 170 | 24 | 256 | default | 450 | 4.5 | 371.20 MB |
| PrimitiveWetModel | 255 | 24 | 384 | default | 300 | 1.0 | 828.73 MB |
| PrimitiveWetModel | 31 | 8 | 48 | matrix | 2400 | 107 | 48.11 MB |
| PrimitiveWetModel | 42 | 8 | 64 | matrix | 1800 | 28 | 133.88 MB |
| PrimitiveWetModel | 63 | 8 | 96 | matrix | 1200 | 3.9 | 582.80 MB |
| PrimitiveWetModel | 85 | 8 | 128 | matrix | 900 | 1.0 | 1.75 GB |
| PrimitiveWetModel | 127 | 8 | 192 | matrix | 600 | 0.1 | 8.19 GB |
| PrimitiveWetModel | 85 | 16 | 128 | matrix | 900 | 1.1 | 1.78 GB |
| PrimitiveWetModel | 127 | 16 | 192 | matrix | 600 | 0.2 | 8.24 GB |
| PrimitiveWetModel | 85 | 24 | 128 | matrix | 900 | 0.9 | 1.80 GB |
| PrimitiveWetModel | 127 | 24 | 192 | matrix | 600 | 0.1 | 8.30 GB |

### Primitive Equation, Float32 vs Float64

| Model | NF | T | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | Float32 | 31 | 8 | 2400 | 864 | 6.22 MB |
| PrimitiveWetModel | Float64 | 31 | 8 | 2400 | 792 | 11.35 MB |

### Grids

| Model | T | L | Grid | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 63 | 8 | FullGaussianGrid | 96 | 1200 | 70 | 32.40 MB |
| PrimitiveWetModel | 63 | 8 | FullClenshawGrid | 95 | 1200 | 75 | 32.13 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralGaussianGrid | 96 | 1200 | 107 | 22.34 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralClenshawGrid | 95 | 1200 | 94 | 22.06 MB |
| PrimitiveWetModel | 63 | 8 | HEALPixGrid | 95 | 1200 | 165 | 16.42 MB |
| PrimitiveWetModel | 63 | 8 | OctaHEALPixGrid | 95 | 1200 | 117 | 19.89 MB |

### Number of vertical layers

| Model | T | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 4 | 2400 | 1419 | 3.87 MB |
| PrimitiveWetModel | 31 | 8 | 2400 | 872 | 6.22 MB |
| PrimitiveWetModel | 31 | 12 | 2400 | 613 | 8.58 MB |
| PrimitiveWetModel | 31 | 16 | 2400 | 480 | 10.94 MB |

### PrimitiveDryModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveDryModel | 31 | 8 | true | true | 2400 | 1266 | 5.27 MB |
| PrimitiveDryModel | 31 | 8 | true | false | 2400 | 2388 | 5.27 MB |
| PrimitiveDryModel | 31 | 8 | false | true | 2400 | 1654 | 5.27 MB |

### PrimitiveWetModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 8 | true | true | 2400 | 860 | 6.22 MB |
| PrimitiveWetModel | 31 | 8 | true | false | 2400 | 1824 | 6.22 MB |
| PrimitiveWetModel | 31 | 8 | false | true | 2400 | 957 | 6.22 MB |

### Individual dynamics functions


#### PrimitiveWetModel | Float32 | T31 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| --- | --- | --- | --- |
| pressure_gradient_flux! | 68.593 μs| 31.98 KiB| 200 |
| linear_virtual_temperature! | 3.280 μs| 0 bytes| 0 |
| geopotential! | 9.780 μs| 384 bytes| 6 |
| vertical_integration! | 14.120 μs| 0 bytes| 0 |
| surface_pressure_tendency! | 18.671 μs| 15.66 KiB| 96 |
| vertical_velocity! | 59.422 μs| 0 bytes| 0 |
| linear_pressure_gradient! | 3.113 μs| 0 bytes| 0 |
| vertical_advection! | 160.175 μs| 2.44 KiB| 32 |
| vordiv_tendencies! | 402.852 μs| 218.48 KiB| 284 |
| temperature_tendency! | 524.946 μs| 324.75 KiB| 401 |
| humidity_tendency! | 501.785 μs| 324.08 KiB| 396 |
| bernoulli_potential! | 172.966 μs| 107.62 KiB| 129 |

## Architecture: `gpu-nvidia`

Created for SpeedyWeather.jl v0.21.1+DEV on Tue, 21 Jul 2026 12:36:45.

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
- CUDACore: 6.2.1
- GPUArrays: 11.5.8
- GPUCompiler: 1.23.0
- KernelAbstractions: 0.9.42
- CUDA_Driver_jll: 13.3.0+1
- CUDA_Compiler_jll: 0.4.4+1
- CUDA_Runtime_jll: 0.23.0+1
- NVPTX_LLVM_Backend_jll: 22.1.7+1

Toolchain:
- Julia: 1.12.2
- LLVM: 18.1.7

1 device:
  0: NVIDIA H100 80GB HBM3 (sm_90, 66.036 GiB / 79.647 GiB available)
     compiles to sm_90a / PTX 9.3 (LLVM: sm_90a / PTX 9.0)
```


### Models, default setups

| Model | T | L | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| BarotropicModel | 31 | 1 | false | 1800 | 2439 | 489.24 KB |
| ShallowWaterModel | 31 | 1 | false | 2400 | 1930 | 492.74 KB |
| PrimitiveDryModel | 31 | 8 | true | 2400 | 337 | 648.86 KB |
| PrimitiveWetModel | 31 | 8 | true | 2400 | 2417 | 654.00 KB |

### Shallow water model, resolution

| Model | T | L | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 31 | 1 | 48 | 2400 | 1734 | 492.74 KB |
| ShallowWaterModel | 42 | 1 | 64 | 1800 | 1051 | 824.77 KB |
| ShallowWaterModel | 63 | 1 | 96 | 1200 | 429 | 1.75 MB |
| ShallowWaterModel | 85 | 1 | 128 | 900 | 245 | 3.03 MB |
| ShallowWaterModel | 127 | 1 | 192 | 600 | 110 | 6.65 MB |
| ShallowWaterModel | 170 | 1 | 256 | 450 | 60 | 11.68 MB |
| ShallowWaterModel | 255 | 1 | 384 | 300 | 28 | 25.96 MB |

### Primitive wet model, resolution

| Model | T | L | Rings | Transform | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 8 | 48 | default | 2400 | 2411 | 654.00 KB |
| PrimitiveWetModel | 42 | 8 | 64 | default | 1800 | 1443 | 1.09 MB |
| PrimitiveWetModel | 63 | 8 | 96 | default | 1200 | 665 | 2.31 MB |
| PrimitiveWetModel | 85 | 8 | 128 | default | 900 | 397 | 4.00 MB |
| PrimitiveWetModel | 127 | 8 | 192 | default | 600 | 157 | 8.78 MB |
| PrimitiveWetModel | 170 | 8 | 256 | default | 450 | 72 | 15.43 MB |
| PrimitiveWetModel | 255 | 8 | 384 | default | 300 | 23 | 34.33 MB |
| PrimitiveWetModel | 85 | 16 | 128 | default | 900 | 1196 | 5.05 MB |
| PrimitiveWetModel | 127 | 16 | 192 | default | 600 | 381 | 11.14 MB |
| PrimitiveWetModel | 170 | 16 | 256 | default | 450 | 144 | 19.63 MB |
| PrimitiveWetModel | 255 | 16 | 384 | default | 300 | 40 | 43.77 MB |
| PrimitiveWetModel | 85 | 24 | 128 | default | 900 | 1210 | 6.10 MB |
| PrimitiveWetModel | 127 | 24 | 192 | default | 600 | 374 | 13.50 MB |
| PrimitiveWetModel | 170 | 24 | 256 | default | 450 | 140 | 23.82 MB |
| PrimitiveWetModel | 255 | 24 | 384 | default | 300 | 36 | 53.21 MB |
| PrimitiveWetModel | 31 | 8 | 48 | matrix | 2400 | 4660 | 581.36 KB |
| PrimitiveWetModel | 42 | 8 | 64 | matrix | 1800 | 2358 | 994.16 KB |
| PrimitiveWetModel | 63 | 8 | 96 | matrix | 1200 | 995 | 2.17 MB |
| PrimitiveWetModel | 85 | 8 | 128 | matrix | 900 | 267 | 3.81 MB |
| PrimitiveWetModel | 127 | 8 | 192 | matrix | 600 | 57 | 8.50 MB |
| PrimitiveWetModel | 85 | 16 | 128 | matrix | 900 | 576 | 4.86 MB |
| PrimitiveWetModel | 127 | 16 | 192 | matrix | 600 | 73 | 10.86 MB |
| PrimitiveWetModel | 85 | 24 | 128 | matrix | 900 | 371 | 5.91 MB |
| PrimitiveWetModel | 127 | 24 | 192 | matrix | 600 | 54 | 13.22 MB |

### Primitive Equation, Float32 vs Float64

| Model | NF | T | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | Float32 | 31 | 8 | 2400 | 2252 | 654.00 KB |
| PrimitiveWetModel | Float64 | 31 | 8 | 2400 | 2277 | 654.61 KB |

### Grids

| Model | T | L | Grid | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 63 | 8 | FullGaussianGrid | 96 | 1200 | 662 | 2.40 MB |
| PrimitiveWetModel | 63 | 8 | FullClenshawGrid | 95 | 1200 | 666 | 2.38 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralGaussianGrid | 96 | 1200 | 662 | 2.31 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralClenshawGrid | 95 | 1200 | 617 | 2.29 MB |
| PrimitiveWetModel | 63 | 8 | HEALPixGrid | 95 | 1200 | 666 | 2.24 MB |
| PrimitiveWetModel | 63 | 8 | OctaHEALPixGrid | 95 | 1200 | 664 | 2.27 MB |

### Number of vertical layers

| Model | T | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 4 | 2400 | 2583 | 580.27 KB |
| PrimitiveWetModel | 31 | 8 | 2400 | 2435 | 654.00 KB |
| PrimitiveWetModel | 31 | 12 | 2400 | 2421 | 727.73 KB |
| PrimitiveWetModel | 31 | 16 | 2400 | 2424 | 801.45 KB |

### PrimitiveDryModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveDryModel | 31 | 8 | true | true | 2400 | 309 | 648.86 KB |
| PrimitiveDryModel | 31 | 8 | true | false | 2400 | 319 | 648.86 KB |
| PrimitiveDryModel | 31 | 8 | false | true | 2400 | 355 | 648.86 KB |

### PrimitiveWetModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 8 | true | true | 2400 | 2588 | 654.00 KB |
| PrimitiveWetModel | 31 | 8 | true | false | 2400 | 2770 | 654.00 KB |
| PrimitiveWetModel | 31 | 8 | false | true | 2400 | 7803 | 654.00 KB |

### Individual dynamics functions


#### PrimitiveWetModel | Float32 | T31 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| --- | --- | --- | --- |
| pressure_gradient_flux! | 1.327 ms| 182.47 KiB| 6126 |
| linear_virtual_temperature! | 12.710 μs| 2.77 KiB| 52 |
| geopotential! | 19.130 μs| 3.73 KiB| 101 |
| vertical_integration! | 27.450 μs| 7.17 KiB| 141 |
| surface_pressure_tendency! | N/A| N/A| N/A |
| vertical_velocity! | 24.670 μs| 8.45 KiB| 191 |
| linear_pressure_gradient! | 12.480 μs| 2.34 KiB| 50 |
| vertical_advection! | 24.000 μs| 9.11 KiB| 143 |
| vordiv_tendencies! | N/A| N/A| N/A |
| temperature_tendency! | N/A| N/A| N/A |
| humidity_tendency! | N/A| N/A| N/A |
| bernoulli_potential! | N/A| N/A| N/A |

