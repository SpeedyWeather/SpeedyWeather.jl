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
julia --project=. manual_benchmarking.jl reactant-cpu   # Reactant on CPU
julia --project=. manual_benchmarking.jl reactant-gpu   # Reactant on CUDA GPU
```

Each run updates only its own architecture's section in this `README.md`; results for other architectures are preserved via `benchmark_results.json`.

## Overview: PrimitiveWet resolution across architectures

Simulated years per wallclock day (SYPD) for the `PrimitiveWetModel` resolution sweep, one column per architecture. Each (T, L) configuration is reported for both the standard Legendre transform and fast Fourier transform (LT+FFT) and the single matrix transform (MT). Empty cells mean the architecture has not yet been benchmarked or that suite was skipped. Comparison figures across architectures are available on the documentation's `Benchmarks` page.

| T | L | Transform | cpu-arm | cpu-x86 | gpu-nvidia |
| --- | --- | --- | --- | --- | --- |
| 32 | 8 | LT+FFT | 1400 | 856 | 5276 |
| 32 | 8 | MT | 757 | 107 | 5569 |
| 43 | 8 | LT+FFT | 564 | 370 | 3115 |
| 43 | 8 | MT | 272 | 28 | 3724 |
| 64 | 8 | LT+FFT | 147 | 107 | 1243 |
| 64 | 8 | MT | 48 | 3.7 | 1184 |
| 86 | 8 | LT+FFT | 57 | 39 | 657 |
| 86 | 8 | MT | 13 | 0.9 | 305 |
| 86 | 16 | LT+FFT | 51 | 49 | 546 |
| 86 | 16 | MT | 19 | 1.1 | 138 |
| 86 | 24 | LT+FFT | 48 | 38 | 536 |
| 86 | 24 | MT | 11 | 0.9 | 77 |
| 128 | 8 | LT+FFT | 15 | 11 | 265 |
| 128 | 8 | MT | 1.5 | 0.1 | 33 |
| 128 | 16 | LT+FFT | 21 | 15 | 241 |
| 128 | 16 | MT | 2.1 | 0.2 | 14 |
| 128 | 24 | LT+FFT | 15 | 11 | 218 |
| 128 | 24 | MT | 1.7 | 0.1 | 9.3 |
| 171 | 8 | LT+FFT | 5.5 | 3.9 | 138 |
| 171 | 16 | LT+FFT | 7.7 | 5.6 | 123 |
| 171 | 24 | LT+FFT | 4.3 | 4.5 | 110 |
| 256 | 8 | LT+FFT | 1.4 | 1.0 | 53 |
| 256 | 16 | LT+FFT | 1.9 | 1.3 | 44 |
| 256 | 24 | LT+FFT | 1.5 | 1.0 | 38 |

## Architecture: `cpu-arm`

Created for SpeedyWeather.jl v0.21.1+DEV on Tue, 21 Jul 2026 17:31:54.

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

| Model | T | L | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| BarotropicModel | 31 | 1 | false | 1800 | 46826 | 780.58 KB |
| ShallowWaterModel | 31 | 1 | false | 2400 | 35514 | 962.86 KB |
| PrimitiveDryModel | 31 | 8 | true | 2400 | 2222 | 5.27 MB |
| PrimitiveWetModel | 31 | 8 | true | 2400 | 1435 | 6.22 MB |

### Shallow water model, resolution

| Model | T | L | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 31 | 1 | 48 | 2400 | 18732 | 962.86 KB |
| ShallowWaterModel | 42 | 1 | 64 | 1800 | 15114 | 1.68 MB |
| ShallowWaterModel | 63 | 1 | 96 | 1200 | 3442 | 3.77 MB |
| ShallowWaterModel | 85 | 1 | 128 | 900 | 1658 | 6.84 MB |
| ShallowWaterModel | 127 | 1 | 192 | 600 | 392 | 16.12 MB |
| ShallowWaterModel | 170 | 1 | 256 | 450 | 133 | 30.33 MB |
| ShallowWaterModel | 255 | 1 | 384 | 300 | 28 | 76.02 MB |

### Primitive wet model, resolution

| Model | T | L | Rings | Transform | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 8 | 48 | default | 2400 | 1400 | 6.22 MB |
| PrimitiveWetModel | 42 | 8 | 64 | default | 1800 | 564 | 10.51 MB |
| PrimitiveWetModel | 63 | 8 | 96 | default | 1200 | 147 | 22.34 MB |
| PrimitiveWetModel | 85 | 8 | 128 | default | 900 | 57 | 38.87 MB |
| PrimitiveWetModel | 127 | 8 | 192 | default | 600 | 15 | 85.50 MB |
| PrimitiveWetModel | 170 | 8 | 256 | default | 450 | 5.5 | 151.57 MB |
| PrimitiveWetModel | 255 | 8 | 384 | default | 300 | 1.4 | 343.69 MB |
| PrimitiveWetModel | 85 | 16 | 128 | default | 900 | 51 | 67.81 MB |
| PrimitiveWetModel | 127 | 16 | 192 | default | 600 | 21 | 148.26 MB |
| PrimitiveWetModel | 170 | 16 | 256 | default | 450 | 7.7 | 261.34 MB |
| PrimitiveWetModel | 255 | 16 | 384 | default | 300 | 1.9 | 586.14 MB |
| PrimitiveWetModel | 85 | 24 | 128 | default | 900 | 48 | 96.80 MB |
| PrimitiveWetModel | 127 | 24 | 192 | default | 600 | 15 | 211.09 MB |
| PrimitiveWetModel | 170 | 24 | 256 | default | 450 | 4.3 | 371.20 MB |
| PrimitiveWetModel | 255 | 24 | 384 | default | 300 | 1.5 | 828.73 MB |
| PrimitiveWetModel | 31 | 8 | 48 | matrix | 2400 | 757 | 48.11 MB |
| PrimitiveWetModel | 42 | 8 | 64 | matrix | 1800 | 272 | 133.88 MB |
| PrimitiveWetModel | 63 | 8 | 96 | matrix | 1200 | 48 | 582.80 MB |
| PrimitiveWetModel | 85 | 8 | 128 | matrix | 900 | 13 | 1.75 GB |
| PrimitiveWetModel | 127 | 8 | 192 | matrix | 600 | 1.5 | 8.19 GB |
| PrimitiveWetModel | 85 | 16 | 128 | matrix | 900 | 19 | 1.78 GB |
| PrimitiveWetModel | 127 | 16 | 192 | matrix | 600 | 2.1 | 8.24 GB |
| PrimitiveWetModel | 85 | 24 | 128 | matrix | 900 | 11 | 1.80 GB |
| PrimitiveWetModel | 127 | 24 | 192 | matrix | 600 | 1.7 | 8.30 GB |

### Primitive Equation, Float32 vs Float64

| Model | NF | T | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | Float32 | 31 | 8 | 2400 | 1227 | 6.22 MB |
| PrimitiveWetModel | Float64 | 31 | 8 | 2400 | 1232 | 11.35 MB |

### Grids

| Model | T | L | Grid | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 63 | 8 | FullGaussianGrid | 96 | 1200 | 80 | 32.40 MB |
| PrimitiveWetModel | 63 | 8 | FullClenshawGrid | 95 | 1200 | 109 | 32.13 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralGaussianGrid | 96 | 1200 | 158 | 22.34 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralClenshawGrid | 95 | 1200 | 146 | 22.06 MB |
| PrimitiveWetModel | 63 | 8 | HEALPixGrid | 95 | 1200 | 219 | 16.42 MB |
| PrimitiveWetModel | 63 | 8 | OctaHEALPixGrid | 95 | 1200 | 171 | 19.89 MB |

### Number of vertical layers

| Model | T | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 4 | 2400 | 2324 | 3.87 MB |
| PrimitiveWetModel | 31 | 8 | 2400 | 1404 | 6.22 MB |
| PrimitiveWetModel | 31 | 12 | 2400 | 1003 | 8.58 MB |
| PrimitiveWetModel | 31 | 16 | 2400 | 771 | 10.94 MB |

### PrimitiveDryModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveDryModel | 31 | 8 | true | true | 2400 | 2368 | 5.27 MB |
| PrimitiveDryModel | 31 | 8 | true | false | 2400 | 4046 | 5.27 MB |
| PrimitiveDryModel | 31 | 8 | false | true | 2400 | 2614 | 5.27 MB |

### PrimitiveWetModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 8 | true | true | 2400 | 1423 | 6.22 MB |
| PrimitiveWetModel | 31 | 8 | true | false | 2400 | 3055 | 6.22 MB |
| PrimitiveWetModel | 31 | 8 | false | true | 2400 | 1558 | 6.22 MB |

### Individual dynamics functions


#### PrimitiveWetModel | Float32 | T31 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| --- | --- | --- | --- |
| pressure_gradient_flux! | 39.666 μs| 31.98 KiB| 200 |
| linear_virtual_temperature! | 2.056 μs| 0 bytes| 0 |
| geopotential! | 7.562 μs| 384 bytes| 6 |
| vertical_integration! | 14.500 μs| 0 bytes| 0 |
| surface_pressure_tendency! | 11.875 μs| 15.66 KiB| 96 |
| vertical_velocity! | 23.333 μs| 0 bytes| 0 |
| linear_pressure_gradient! | 2.065 μs| 0 bytes| 0 |
| vertical_advection! | 112.500 μs| 2.44 KiB| 32 |
| vordiv_tendencies! | 224.375 μs| 231.83 KiB| 284 |
| temperature_tendency! | 291.583 μs| 344.77 KiB| 401 |
| humidity_tendency! | 279.375 μs| 344.09 KiB| 396 |
| bernoulli_potential! | 93.667 μs| 114.30 KiB| 129 |

## Architecture: `cpu-x86`

Created for SpeedyWeather.jl v0.21.1+DEV on Tue, 21 Jul 2026 16:59:17.

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
| BarotropicModel | 31 | 1 | false | 1800 | 10654 | 780.69 KB |
| ShallowWaterModel | 31 | 1 | false | 2400 | 19394 | 962.97 KB |
| PrimitiveDryModel | 31 | 8 | true | 2400 | 1467 | 5.27 MB |
| PrimitiveWetModel | 31 | 8 | true | 2400 | 819 | 6.22 MB |

### Shallow water model, resolution

| Model | T | L | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 31 | 1 | 48 | 2400 | 19308 | 962.97 KB |
| ShallowWaterModel | 42 | 1 | 64 | 1800 | 8616 | 1.68 MB |
| ShallowWaterModel | 63 | 1 | 96 | 1200 | 2297 | 3.77 MB |
| ShallowWaterModel | 85 | 1 | 128 | 900 | 960 | 6.84 MB |
| ShallowWaterModel | 127 | 1 | 192 | 600 | 231 | 16.12 MB |
| ShallowWaterModel | 170 | 1 | 256 | 450 | 86 | 30.33 MB |
| ShallowWaterModel | 255 | 1 | 384 | 300 | 18 | 76.02 MB |

### Primitive wet model, resolution

| Model | T | L | Rings | Transform | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 8 | 48 | default | 2400 | 856 | 6.22 MB |
| PrimitiveWetModel | 42 | 8 | 64 | default | 1800 | 370 | 10.51 MB |
| PrimitiveWetModel | 63 | 8 | 96 | default | 1200 | 107 | 22.34 MB |
| PrimitiveWetModel | 85 | 8 | 128 | default | 900 | 39 | 38.87 MB |
| PrimitiveWetModel | 127 | 8 | 192 | default | 600 | 11 | 85.50 MB |
| PrimitiveWetModel | 170 | 8 | 256 | default | 450 | 3.9 | 151.57 MB |
| PrimitiveWetModel | 255 | 8 | 384 | default | 300 | 1.0 | 343.69 MB |
| PrimitiveWetModel | 85 | 16 | 128 | default | 900 | 49 | 67.81 MB |
| PrimitiveWetModel | 127 | 16 | 192 | default | 600 | 15 | 148.26 MB |
| PrimitiveWetModel | 170 | 16 | 256 | default | 450 | 5.6 | 261.34 MB |
| PrimitiveWetModel | 255 | 16 | 384 | default | 300 | 1.3 | 586.14 MB |
| PrimitiveWetModel | 85 | 24 | 128 | default | 900 | 38 | 96.80 MB |
| PrimitiveWetModel | 127 | 24 | 192 | default | 600 | 11 | 211.09 MB |
| PrimitiveWetModel | 170 | 24 | 256 | default | 450 | 4.5 | 371.20 MB |
| PrimitiveWetModel | 255 | 24 | 384 | default | 300 | 1.0 | 828.73 MB |
| PrimitiveWetModel | 31 | 8 | 48 | matrix | 2400 | 107 | 48.11 MB |
| PrimitiveWetModel | 42 | 8 | 64 | matrix | 1800 | 28 | 133.88 MB |
| PrimitiveWetModel | 63 | 8 | 96 | matrix | 1200 | 3.7 | 582.80 MB |
| PrimitiveWetModel | 85 | 8 | 128 | matrix | 900 | 0.9 | 1.75 GB |
| PrimitiveWetModel | 127 | 8 | 192 | matrix | 600 | 0.1 | 8.19 GB |
| PrimitiveWetModel | 85 | 16 | 128 | matrix | 900 | 1.1 | 1.78 GB |
| PrimitiveWetModel | 127 | 16 | 192 | matrix | 600 | 0.2 | 8.24 GB |
| PrimitiveWetModel | 85 | 24 | 128 | matrix | 900 | 0.9 | 1.80 GB |
| PrimitiveWetModel | 127 | 24 | 192 | matrix | 600 | 0.1 | 8.30 GB |

### Primitive Equation, Float32 vs Float64

| Model | NF | T | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | Float32 | 31 | 8 | 2400 | 850 | 6.22 MB |
| PrimitiveWetModel | Float64 | 31 | 8 | 2400 | 767 | 11.35 MB |

### Grids

| Model | T | L | Grid | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 63 | 8 | FullGaussianGrid | 96 | 1200 | 71 | 32.40 MB |
| PrimitiveWetModel | 63 | 8 | FullClenshawGrid | 95 | 1200 | 74 | 32.13 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralGaussianGrid | 96 | 1200 | 100 | 22.34 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralClenshawGrid | 95 | 1200 | 100 | 22.06 MB |
| PrimitiveWetModel | 63 | 8 | HEALPixGrid | 95 | 1200 | 163 | 16.42 MB |
| PrimitiveWetModel | 63 | 8 | OctaHEALPixGrid | 95 | 1200 | 119 | 19.89 MB |

### Number of vertical layers

| Model | T | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 4 | 2400 | 1412 | 3.87 MB |
| PrimitiveWetModel | 31 | 8 | 2400 | 851 | 6.22 MB |
| PrimitiveWetModel | 31 | 12 | 2400 | 603 | 8.58 MB |
| PrimitiveWetModel | 31 | 16 | 2400 | 467 | 10.94 MB |

### PrimitiveDryModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveDryModel | 31 | 8 | true | true | 2400 | 1456 | 5.27 MB |
| PrimitiveDryModel | 31 | 8 | true | false | 2400 | 2293 | 5.27 MB |
| PrimitiveDryModel | 31 | 8 | false | true | 2400 | 1609 | 5.27 MB |

### PrimitiveWetModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 8 | true | true | 2400 | 852 | 6.22 MB |
| PrimitiveWetModel | 31 | 8 | true | false | 2400 | 1773 | 6.22 MB |
| PrimitiveWetModel | 31 | 8 | false | true | 2400 | 951 | 6.22 MB |

### Individual dynamics functions


#### PrimitiveWetModel | Float32 | T31 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| --- | --- | --- | --- |
| pressure_gradient_flux! | 70.633 μs| 31.98 KiB| 200 |
| linear_virtual_temperature! | 3.278 μs| 0 bytes| 0 |
| geopotential! | 10.260 μs| 384 bytes| 6 |
| vertical_integration! | 13.680 μs| 0 bytes| 0 |
| surface_pressure_tendency! | 18.721 μs| 15.66 KiB| 96 |
| vertical_velocity! | 59.542 μs| 0 bytes| 0 |
| linear_pressure_gradient! | 3.179 μs| 0 bytes| 0 |
| vertical_advection! | 160.256 μs| 2.44 KiB| 32 |
| vordiv_tendencies! | 412.104 μs| 218.48 KiB| 284 |
| temperature_tendency! | 532.669 μs| 324.75 KiB| 401 |
| humidity_tendency! | 507.748 μs| 324.08 KiB| 396 |
| bernoulli_potential! | 174.696 μs| 107.62 KiB| 129 |

## Architecture: `gpu-nvidia`

Created for SpeedyWeather.jl v0.22.0+DEV on Tue, 18 Aug 2026 14:54:44.

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
- GPUArrays: 11.5.11
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
  0: NVIDIA H100 80GB HBM3 (sm_90, 65.839 GiB / 79.647 GiB available)
     compiles to sm_90a / PTX 9.3 (LLVM: sm_90a / PTX 9.0)
```


### Models, default setups

| Model | truncation | L | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| BarotropicModel | 32 | 1 | false | 1800 | 2113 | 490.02 KB |
| ShallowWaterModel | 32 | 1 | false | 2400 | 1671 | 493.52 KB |
| PrimitiveDryModel | 32 | 8 | true | 2400 | 355 | 661.63 KB |
| PrimitiveWetModel | 32 | 8 | true | 2400 | 5708 | 666.76 KB |

### Shallow water model, resolution

| Model | truncation | L | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 32 | 1 | 48 | 2400 | 1727 | 493.52 KB |
| ShallowWaterModel | 43 | 1 | 64 | 1800 | 960 | 825.68 KB |
| ShallowWaterModel | 64 | 1 | 96 | 1200 | 425 | 1.75 MB |
| ShallowWaterModel | 86 | 1 | 128 | 900 | 240 | 3.03 MB |
| ShallowWaterModel | 128 | 1 | 192 | 600 | 107 | 6.65 MB |
| ShallowWaterModel | 171 | 1 | 256 | 450 | 59 | 11.68 MB |
| ShallowWaterModel | 256 | 1 | 384 | 300 | 26 | 25.96 MB |

### Primitive wet model, resolution

| Model | truncation | L | Rings | Transform | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | 48 | default | 2400 | 5276 | 666.76 KB |
| PrimitiveWetModel | 43 | 8 | 64 | default | 1800 | 3115 | 1.11 MB |
| PrimitiveWetModel | 64 | 8 | 96 | default | 1200 | 1243 | 2.34 MB |
| PrimitiveWetModel | 86 | 8 | 128 | default | 900 | 657 | 4.04 MB |
| PrimitiveWetModel | 128 | 8 | 192 | default | 600 | 265 | 8.83 MB |
| PrimitiveWetModel | 171 | 8 | 256 | default | 450 | 138 | 15.50 MB |
| PrimitiveWetModel | 256 | 8 | 384 | default | 300 | 53 | 34.43 MB |
| PrimitiveWetModel | 86 | 16 | 128 | default | 900 | 546 | 5.08 MB |
| PrimitiveWetModel | 128 | 16 | 192 | default | 600 | 241 | 11.19 MB |
| PrimitiveWetModel | 171 | 16 | 256 | default | 450 | 123 | 19.69 MB |
| PrimitiveWetModel | 256 | 16 | 384 | default | 300 | 44 | 43.87 MB |
| PrimitiveWetModel | 86 | 24 | 128 | default | 900 | 536 | 6.13 MB |
| PrimitiveWetModel | 128 | 24 | 192 | default | 600 | 218 | 13.55 MB |
| PrimitiveWetModel | 171 | 24 | 256 | default | 450 | 110 | 23.89 MB |
| PrimitiveWetModel | 256 | 24 | 384 | default | 300 | 38 | 53.31 MB |
| PrimitiveWetModel | 32 | 8 | 48 | matrix | 2400 | 5569 | 581.83 KB |
| PrimitiveWetModel | 43 | 8 | 64 | matrix | 1800 | 3724 | 994.76 KB |
| PrimitiveWetModel | 64 | 8 | 96 | matrix | 1200 | 1184 | 2.17 MB |
| PrimitiveWetModel | 86 | 8 | 128 | matrix | 900 | 305 | 3.81 MB |
| PrimitiveWetModel | 128 | 8 | 192 | matrix | 600 | 33 | 8.50 MB |
| PrimitiveWetModel | 86 | 16 | 128 | matrix | 900 | 138 | 4.86 MB |
| PrimitiveWetModel | 128 | 16 | 192 | matrix | 600 | 14 | 10.86 MB |
| PrimitiveWetModel | 86 | 24 | 128 | matrix | 900 | 77 | 5.91 MB |
| PrimitiveWetModel | 128 | 24 | 192 | matrix | 600 | 9.3 | 13.22 MB |

### Primitive Equation, Float32 vs Float64

| Model | NF | truncation | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | Float32 | 32 | 8 | 2400 | 5243 | 666.76 KB |
| PrimitiveWetModel | Float64 | 32 | 8 | 2400 | 4757 | 667.37 KB |

### Grids

| Model | truncation | L | Grid | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 64 | 8 | FullGaussianGrid | 96 | 1200 | 1288 | 2.43 MB |
| PrimitiveWetModel | 64 | 8 | FullClenshawGrid | 127 | 1200 | 1135 | 4.17 MB |
| PrimitiveWetModel | 64 | 8 | OctahedralGaussianGrid | 96 | 1200 | 1279 | 2.34 MB |
| PrimitiveWetModel | 64 | 8 | OctahedralClenshawGrid | 127 | 1200 | 911 | 4.01 MB |
| PrimitiveWetModel | 64 | 8 | HEALPixGrid | 127 | 1200 | 1153 | 3.93 MB |
| PrimitiveWetModel | 64 | 8 | OctaHEALPixGrid | 127 | 1200 | 938 | 3.98 MB |

### Number of vertical layers

| Model | truncation | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 4 | 2400 | 5476 | 593.04 KB |
| PrimitiveWetModel | 32 | 8 | 2400 | 5300 | 666.76 KB |
| PrimitiveWetModel | 32 | 12 | 2400 | 4935 | 740.49 KB |
| PrimitiveWetModel | 32 | 16 | 2400 | 4874 | 814.22 KB |

### PrimitiveDryModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveDryModel | 32 | 8 | true | true | 2400 | 324 | 661.63 KB |
| PrimitiveDryModel | 32 | 8 | true | false | 2400 | 330 | 661.63 KB |
| PrimitiveDryModel | 32 | 8 | false | true | 2400 | 350 | 661.63 KB |

### PrimitiveWetModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | true | true | 2400 | 5431 | 666.76 KB |
| PrimitiveWetModel | 32 | 8 | true | false | 2400 | 5955 | 666.76 KB |
| PrimitiveWetModel | 32 | 8 | false | true | 2400 | 7898 | 666.76 KB |

### Individual dynamics functions


#### PrimitiveWetModel | Float32 | T31 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| --- | --- | --- | --- |
| pressure_gradient_flux! | 226.718 μs| 14.30 KiB| 323 |
| linear_virtual_temperature! | 13.490 μs| 2.77 KiB| 52 |
| geopotential! | 19.740 μs| 3.73 KiB| 101 |
| vertical_integration! | 29.120 μs| 7.17 KiB| 141 |
| surface_pressure_tendency! | N/A| N/A| N/A |
| vertical_velocity! | 25.250 μs| 8.45 KiB| 191 |
| linear_pressure_gradient! | 13.110 μs| 2.34 KiB| 50 |
| vertical_advection! | 25.509 μs| 9.09 KiB| 142 |
| vordiv_tendencies! | N/A| N/A| N/A |
| temperature_tendency! | N/A| N/A| N/A |
| humidity_tendency! | N/A| N/A| N/A |
| bernoulli_potential! | N/A| N/A| N/A |

