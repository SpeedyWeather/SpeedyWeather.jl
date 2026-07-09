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
| 31 | 8 | LT+FFT | 1326 | 908 | 2171 |
| 31 | 8 | MT | 678 | 71 | 4140 |
| 42 | 8 | LT+FFT | 556 | 389 | 1187 |
| 42 | 8 | MT | 157 | 19 | 2118 |
| 63 | 8 | LT+FFT | 146 | 116 | 592 |
| 63 | 8 | MT | 33 | 2.7 | 1016 |
| 85 | 8 | LT+FFT | 63 | 45 | 350 |
| 85 | 8 | MT | 7.4 | 0.7 | 268 |
| 85 | 16 | LT+FFT | 31 | 54 | 1175 |
| 85 | 16 | MT | 5.1 | 0.8 | 572 |
| 85 | 24 | LT+FFT | 15 | 41 | 1209 |
| 85 | 24 | MT | 4.1 | 0.6 | 370 |
| 127 | 8 | LT+FFT | 14 | 12 | 156 |
| 127 | 8 | MT | 0.5 | 0.1 | 58 |
| 127 | 16 | LT+FFT | 7.9 | 16 | 368 |
| 127 | 16 | MT | 0.4 | 0.1 | 82 |
| 127 | 24 | LT+FFT | 4.8 | 12 | 373 |
| 127 | 24 | MT | 0.4 | 0.1 | 54 |
| 170 | 8 | LT+FFT | 5.3 | 4.4 | 72 |
| 170 | 16 | LT+FFT | 2.9 | 6.1 | 141 |
| 170 | 24 | LT+FFT | 1.7 | 4.7 | 167 |
| 255 | 8 | LT+FFT | 1.4 | 1.1 | 24 |
| 255 | 16 | LT+FFT | 0.7 | 1.5 | 43 |
| 255 | 24 | LT+FFT | 0.5 | 1.1 | 43 |

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

Created for SpeedyWeather.jl v0.21.1+DEV on Thu, 09 Jul 2026 17:06:12.

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
| BarotropicModel | 31 | 1 | false | 1800 | 18285 | 764.91 KB |
| ShallowWaterModel | 31 | 1 | false | 2400 | 17906 | 907.93 KB |
| PrimitiveDryModel | 31 | 8 | true | 2400 | 1426 | 5.20 MB |
| PrimitiveWetModel | 31 | 8 | true | 2400 | 905 | 6.15 MB |

### Shallow water model, resolution

| Model | T | L | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 31 | 1 | 48 | 2400 | 18802 | 907.93 KB |
| ShallowWaterModel | 42 | 1 | 64 | 1800 | 7145 | 1.59 MB |
| ShallowWaterModel | 63 | 1 | 96 | 1200 | 2382 | 3.56 MB |
| ShallowWaterModel | 85 | 1 | 128 | 900 | 901 | 6.49 MB |
| ShallowWaterModel | 127 | 1 | 192 | 600 | 232 | 15.33 MB |
| ShallowWaterModel | 170 | 1 | 256 | 450 | 87 | 28.95 MB |
| ShallowWaterModel | 255 | 1 | 384 | 300 | 19 | 72.94 MB |

### Primitive wet model, resolution

| Model | T | L | Rings | Transform | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 8 | 48 | default | 2400 | 908 | 6.15 MB |
| PrimitiveWetModel | 42 | 8 | 64 | default | 1800 | 389 | 10.40 MB |
| PrimitiveWetModel | 63 | 8 | 96 | default | 1200 | 116 | 22.10 MB |
| PrimitiveWetModel | 85 | 8 | 128 | default | 900 | 45 | 38.45 MB |
| PrimitiveWetModel | 127 | 8 | 192 | default | 600 | 12 | 84.59 MB |
| PrimitiveWetModel | 170 | 8 | 256 | default | 450 | 4.4 | 149.99 MB |
| PrimitiveWetModel | 255 | 8 | 384 | default | 300 | 1.1 | 340.17 MB |
| PrimitiveWetModel | 85 | 16 | 128 | default | 900 | 54 | 67.40 MB |
| PrimitiveWetModel | 127 | 16 | 192 | default | 600 | 16 | 147.35 MB |
| PrimitiveWetModel | 170 | 16 | 256 | default | 450 | 6.1 | 259.75 MB |
| PrimitiveWetModel | 255 | 16 | 384 | default | 300 | 1.5 | 582.63 MB |
| PrimitiveWetModel | 85 | 24 | 128 | default | 900 | 41 | 96.39 MB |
| PrimitiveWetModel | 127 | 24 | 192 | default | 600 | 12 | 210.18 MB |
| PrimitiveWetModel | 170 | 24 | 256 | default | 450 | 4.7 | 369.61 MB |
| PrimitiveWetModel | 255 | 24 | 384 | default | 300 | 1.1 | 825.21 MB |
| PrimitiveWetModel | 31 | 8 | 48 | matrix | 2400 | 71 | 48.05 MB |
| PrimitiveWetModel | 42 | 8 | 64 | matrix | 1800 | 19 | 133.77 MB |
| PrimitiveWetModel | 63 | 8 | 96 | matrix | 1200 | 2.7 | 582.56 MB |
| PrimitiveWetModel | 85 | 8 | 128 | matrix | 900 | 0.7 | 1.75 GB |
| PrimitiveWetModel | 127 | 8 | 192 | matrix | 600 | 0.1 | 8.19 GB |
| PrimitiveWetModel | 85 | 16 | 128 | matrix | 900 | 0.8 | 1.78 GB |
| PrimitiveWetModel | 127 | 16 | 192 | matrix | 600 | 0.1 | 8.24 GB |
| PrimitiveWetModel | 85 | 24 | 128 | matrix | 900 | 0.6 | 1.80 GB |
| PrimitiveWetModel | 127 | 24 | 192 | matrix | 600 | 0.1 | 8.30 GB |

### Primitive Equation, Float32 vs Float64

| Model | NF | T | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | Float32 | 31 | 8 | 2400 | 904 | 6.15 MB |
| PrimitiveWetModel | Float64 | 31 | 8 | 2400 | 789 | 11.23 MB |

### Grids

| Model | T | L | Grid | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 63 | 8 | FullGaussianGrid | 96 | 1200 | 76 | 32.04 MB |
| PrimitiveWetModel | 63 | 8 | FullClenshawGrid | 95 | 1200 | 84 | 31.78 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralGaussianGrid | 96 | 1200 | 116 | 22.10 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralClenshawGrid | 95 | 1200 | 108 | 21.82 MB |
| PrimitiveWetModel | 63 | 8 | HEALPixGrid | 95 | 1200 | 170 | 16.24 MB |
| PrimitiveWetModel | 63 | 8 | OctaHEALPixGrid | 95 | 1200 | 126 | 19.68 MB |

### Number of vertical layers

| Model | T | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 4 | 2400 | 1415 | 3.80 MB |
| PrimitiveWetModel | 31 | 8 | 2400 | 901 | 6.15 MB |
| PrimitiveWetModel | 31 | 12 | 2400 | 635 | 8.51 MB |
| PrimitiveWetModel | 31 | 16 | 2400 | 495 | 10.87 MB |

### PrimitiveDryModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveDryModel | 31 | 8 | true | true | 2400 | 1464 | 5.20 MB |
| PrimitiveDryModel | 31 | 8 | true | false | 2400 | 1766 | 5.20 MB |
| PrimitiveDryModel | 31 | 8 | false | true | 2400 | 1603 | 5.20 MB |

### PrimitiveWetModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 8 | true | true | 2400 | 896 | 6.15 MB |
| PrimitiveWetModel | 31 | 8 | true | false | 2400 | 1725 | 6.15 MB |
| PrimitiveWetModel | 31 | 8 | false | true | 2400 | 1026 | 6.15 MB |

### Individual dynamics functions


#### PrimitiveWetModel | Float32 | T31 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| --- | --- | --- | --- |
| pressure_gradient_flux! | 81.283 μs| 52.20 KiB| 593 |
| linear_virtual_temperature! | 3.331 μs| 0 bytes| 0 |
| geopotential! | 11.171 μs| 816 bytes| 13 |
| vertical_integration! | 13.970 μs| 0 bytes| 0 |
| surface_pressure_tendency! | 22.021 μs| 25.41 KiB| 288 |
| vertical_velocity! | 59.412 μs| 0 bytes| 0 |
| linear_pressure_gradient! | 3.076 μs| 0 bytes| 0 |
| vertical_advection! | 167.826 μs| 4.69 KiB| 60 |
| vordiv_tendencies! | 424.156 μs| 249.42 KiB| 711 |
| temperature_tendency! | 552.770 μs| 365.66 KiB| 1012 |
| humidity_tendency! | 529.530 μs| 364.31 KiB| 1002 |
| bernoulli_potential! | 173.087 μs| 120.16 KiB| 327 |

## Architecture: `gpu-nvidia`

Created for SpeedyWeather.jl v0.21.1+DEV on Thu, 09 Jul 2026 16:15:22.

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
- CUDACore: 6.2.0
- GPUArrays: 11.5.8
- GPUCompiler: 1.22.7
- KernelAbstractions: 0.9.42
- CUDA_Driver_jll: 13.3.0+0
- CUDA_Compiler_jll: 0.4.4+1
- CUDA_Runtime_jll: 0.23.0+1
- NVPTX_LLVM_Backend_jll: 22.1.7+1

Toolchain:
- Julia: 1.12.2
- LLVM: 18.1.7

1 device:
  0: NVIDIA H100 80GB HBM3 (sm_90, 65.604 GiB / 79.647 GiB available)
     compiles to sm_90a / PTX 9.3 (LLVM: sm_90a / PTX 9.0)
```


### Models, default setups

| Model | T | L | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| BarotropicModel | 31 | 1 | false | 1800 | 1701 | 470.62 KB |
| ShallowWaterModel | 31 | 1 | false | 2400 | 907 | 473.65 KB |
| PrimitiveDryModel | 31 | 8 | true | 2400 | 300 | 629.74 KB |
| PrimitiveWetModel | 31 | 8 | true | 2400 | 2287 | 634.88 KB |

### Shallow water model, resolution

| Model | T | L | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 31 | 1 | 48 | 2400 | 923 | 473.65 KB |
| ShallowWaterModel | 42 | 1 | 64 | 1800 | 503 | 791.35 KB |
| ShallowWaterModel | 63 | 1 | 96 | 1200 | 216 | 1.68 MB |
| ShallowWaterModel | 85 | 1 | 128 | 900 | 120 | 2.90 MB |
| ShallowWaterModel | 127 | 1 | 192 | 600 | 52 | 6.36 MB |
| ShallowWaterModel | 170 | 1 | 256 | 450 | 29 | 11.15 MB |
| ShallowWaterModel | 255 | 1 | 384 | 300 | 13 | 24.78 MB |

### Primitive wet model, resolution

| Model | T | L | Rings | Transform | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 8 | 48 | default | 2400 | 2171 | 634.88 KB |
| PrimitiveWetModel | 42 | 8 | 64 | default | 1800 | 1187 | 1.06 MB |
| PrimitiveWetModel | 63 | 8 | 96 | default | 1200 | 592 | 2.24 MB |
| PrimitiveWetModel | 85 | 8 | 128 | default | 900 | 350 | 3.87 MB |
| PrimitiveWetModel | 127 | 8 | 192 | default | 600 | 156 | 8.49 MB |
| PrimitiveWetModel | 170 | 8 | 256 | default | 450 | 72 | 14.91 MB |
| PrimitiveWetModel | 255 | 8 | 384 | default | 300 | 24 | 33.15 MB |
| PrimitiveWetModel | 85 | 16 | 128 | default | 900 | 1175 | 4.92 MB |
| PrimitiveWetModel | 127 | 16 | 192 | default | 600 | 368 | 10.85 MB |
| PrimitiveWetModel | 170 | 16 | 256 | default | 450 | 141 | 19.10 MB |
| PrimitiveWetModel | 255 | 16 | 384 | default | 300 | 43 | 42.59 MB |
| PrimitiveWetModel | 85 | 24 | 128 | default | 900 | 1209 | 5.97 MB |
| PrimitiveWetModel | 127 | 24 | 192 | default | 600 | 373 | 13.21 MB |
| PrimitiveWetModel | 170 | 24 | 256 | default | 450 | 167 | 23.30 MB |
| PrimitiveWetModel | 255 | 24 | 384 | default | 300 | 43 | 52.03 MB |
| PrimitiveWetModel | 31 | 8 | 48 | matrix | 2400 | 4140 | 562.25 KB |
| PrimitiveWetModel | 42 | 8 | 64 | matrix | 1800 | 2118 | 960.72 KB |
| PrimitiveWetModel | 63 | 8 | 96 | matrix | 1200 | 1016 | 2.10 MB |
| PrimitiveWetModel | 85 | 8 | 128 | matrix | 900 | 268 | 3.68 MB |
| PrimitiveWetModel | 127 | 8 | 192 | matrix | 600 | 58 | 8.20 MB |
| PrimitiveWetModel | 85 | 16 | 128 | matrix | 900 | 572 | 4.73 MB |
| PrimitiveWetModel | 127 | 16 | 192 | matrix | 600 | 82 | 10.56 MB |
| PrimitiveWetModel | 85 | 24 | 128 | matrix | 900 | 370 | 5.78 MB |
| PrimitiveWetModel | 127 | 24 | 192 | matrix | 600 | 54 | 12.92 MB |

### Primitive Equation, Float32 vs Float64

| Model | NF | T | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | Float32 | 31 | 8 | 2400 | 2006 | 634.88 KB |
| PrimitiveWetModel | Float64 | 31 | 8 | 2400 | 2151 | 635.49 KB |

### Grids

| Model | T | L | Grid | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 63 | 8 | FullGaussianGrid | 96 | 1200 | 615 | 2.33 MB |
| PrimitiveWetModel | 63 | 8 | FullClenshawGrid | 95 | 1200 | 666 | 2.31 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralGaussianGrid | 96 | 1200 | 618 | 2.24 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralClenshawGrid | 95 | 1200 | 625 | 2.22 MB |
| PrimitiveWetModel | 63 | 8 | HEALPixGrid | 95 | 1200 | 668 | 2.17 MB |
| PrimitiveWetModel | 63 | 8 | OctaHEALPixGrid | 95 | 1200 | 657 | 2.20 MB |

### Number of vertical layers

| Model | T | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 4 | 2400 | 2283 | 561.15 KB |
| PrimitiveWetModel | 31 | 8 | 2400 | 2279 | 634.88 KB |
| PrimitiveWetModel | 31 | 12 | 2400 | 2285 | 708.61 KB |
| PrimitiveWetModel | 31 | 16 | 2400 | 2257 | 782.33 KB |

### PrimitiveDryModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveDryModel | 31 | 8 | true | true | 2400 | 285 | 629.74 KB |
| PrimitiveDryModel | 31 | 8 | true | false | 2400 | 292 | 629.74 KB |
| PrimitiveDryModel | 31 | 8 | false | true | 2400 | 314 | 629.74 KB |

### PrimitiveWetModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 8 | true | true | 2400 | 2261 | 634.88 KB |
| PrimitiveWetModel | 31 | 8 | true | false | 2400 | 2597 | 634.88 KB |
| PrimitiveWetModel | 31 | 8 | false | true | 2400 | 7845 | 634.88 KB |

### Individual dynamics functions


#### PrimitiveWetModel | Float32 | T31 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| --- | --- | --- | --- |
| pressure_gradient_flux! | 1.419 ms| 187.58 KiB| 6336 |
| linear_virtual_temperature! | 13.280 μs| 2.77 KiB| 52 |
| geopotential! | 21.820 μs| 4.09 KiB| 108 |
| vertical_integration! | 32.330 μs| 8.83 KiB| 163 |
| surface_pressure_tendency! | N/A| N/A| N/A |
| vertical_velocity! | 26.720 μs| 8.45 KiB| 191 |
| linear_pressure_gradient! | 13.090 μs| 2.34 KiB| 50 |
| vertical_advection! | 32.170 μs| 10.72 KiB| 170 |
| vordiv_tendencies! | N/A| N/A| N/A |
| temperature_tendency! | N/A| N/A| N/A |
| humidity_tendency! | N/A| N/A| N/A |
| bernoulli_potential! | N/A| N/A| N/A |

