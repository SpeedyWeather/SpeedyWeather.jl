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
| 31 | 8 | LT+FFT | 1326 | 884 | 1495 |
| 31 | 8 | MT | 678 | 67 | 1684 |
| 42 | 8 | LT+FFT | 556 | 389 | 990 |
| 42 | 8 | MT | 157 | 18 | 429 |
| 63 | 8 | LT+FFT | 146 | 113 | 405 |
| 63 | 8 | MT | 33 | 2.5 | 175 |
| 85 | 8 | LT+FFT | 63 | 27 | 225 |
| 85 | 8 | MT | 7.4 | 0.6 | 63 |
| 85 | 16 | LT+FFT | 31 | 22 | 236 |
| 85 | 16 | MT | 5.1 | 0.3 | 62 |
| 85 | 24 | LT+FFT | 15 | 13 | 244 |
| 85 | 24 | MT | 4.1 | 0.2 | 62 |
| 127 | 8 | LT+FFT | 14 | 11 | 74 |
| 127 | 8 | MT | 0.5 | 0.1 | 19 |
| 127 | 16 | LT+FFT | 7.9 | 5.6 | 74 |
| 127 | 16 | MT | 0.4 | 0.0 | 18 |
| 127 | 24 | LT+FFT | 4.8 | 3.6 | 74 |
| 127 | 24 | MT | 0.4 | 0.0 | 18 |
| 170 | 8 | LT+FFT | 5.3 | 4.0 | 34 |
| 170 | 16 | LT+FFT | 2.9 | 2.0 | 31 |
| 170 | 24 | LT+FFT | 1.7 | 1.4 | 25 |
| 255 | 8 | LT+FFT | 1.4 | 1.0 | 12 |
| 255 | 16 | LT+FFT | 0.7 | 0.5 | 10 |
| 255 | 24 | LT+FFT | 0.5 | 0.3 | 8.4 |

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

Created for SpeedyWeather.jl v0.21.1+DEV on Tue, 30 Jun 2026 17:08:38.

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
| BarotropicModel | 31 | 1 | false | 1800 | 20342 | 746.49 KB |
| ShallowWaterModel | 31 | 1 | false | 2400 | 15269 | 881.09 KB |
| PrimitiveDryModel | 31 | 8 | true | 2400 | 1493 | 4.70 MB |
| PrimitiveWetModel | 31 | 8 | true | 2400 | 882 | 5.38 MB |

### Shallow water model, resolution

| Model | T | L | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 31 | 1 | 48 | 2400 | 15174 | 881.09 KB |
| ShallowWaterModel | 42 | 1 | 64 | 1800 | 7282 | 1.54 MB |
| ShallowWaterModel | 63 | 1 | 96 | 1200 | 2201 | 3.47 MB |
| ShallowWaterModel | 85 | 1 | 128 | 900 | 907 | 6.32 MB |
| ShallowWaterModel | 127 | 1 | 192 | 600 | 183 | 14.97 MB |
| ShallowWaterModel | 170 | 1 | 256 | 450 | 89 | 28.31 MB |
| ShallowWaterModel | 255 | 1 | 384 | 300 | 18 | 71.53 MB |

### Primitive wet model, resolution

| Model | T | L | Rings | Transform | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 8 | 48 | default | 2400 | 884 | 5.38 MB |
| PrimitiveWetModel | 42 | 8 | 64 | default | 1800 | 389 | 9.10 MB |
| PrimitiveWetModel | 63 | 8 | 96 | default | 1200 | 113 | 19.36 MB |
| PrimitiveWetModel | 85 | 8 | 128 | default | 900 | 27 | 33.68 MB |
| PrimitiveWetModel | 127 | 8 | 192 | default | 600 | 11 | 74.25 MB |
| PrimitiveWetModel | 170 | 8 | 256 | default | 450 | 4.0 | 131.86 MB |
| PrimitiveWetModel | 255 | 8 | 384 | default | 300 | 1.0 | 300.09 MB |
| PrimitiveWetModel | 85 | 16 | 128 | default | 900 | 22 | 57.92 MB |
| PrimitiveWetModel | 127 | 16 | 192 | default | 600 | 5.6 | 126.80 MB |
| PrimitiveWetModel | 170 | 16 | 256 | default | 450 | 2.0 | 223.74 MB |
| PrimitiveWetModel | 255 | 16 | 384 | default | 300 | 0.5 | 502.99 MB |
| PrimitiveWetModel | 85 | 24 | 128 | default | 900 | 13 | 82.21 MB |
| PrimitiveWetModel | 127 | 24 | 192 | default | 600 | 3.6 | 179.42 MB |
| PrimitiveWetModel | 170 | 24 | 256 | default | 450 | 1.4 | 315.70 MB |
| PrimitiveWetModel | 255 | 24 | 384 | default | 300 | 0.3 | 706.03 MB |
| PrimitiveWetModel | 31 | 8 | 48 | matrix | 2400 | 67 | 47.28 MB |
| PrimitiveWetModel | 42 | 8 | 64 | matrix | 1800 | 18 | 132.47 MB |
| PrimitiveWetModel | 63 | 8 | 96 | matrix | 1200 | 2.5 | 579.81 MB |
| PrimitiveWetModel | 85 | 8 | 128 | matrix | 900 | 0.6 | 1.75 GB |
| PrimitiveWetModel | 127 | 8 | 192 | matrix | 600 | 0.1 | 8.18 GB |
| PrimitiveWetModel | 85 | 16 | 128 | matrix | 900 | 0.3 | 1.77 GB |
| PrimitiveWetModel | 127 | 16 | 192 | matrix | 600 | 0.0 | 8.22 GB |
| PrimitiveWetModel | 85 | 24 | 128 | matrix | 900 | 0.2 | 1.79 GB |
| PrimitiveWetModel | 127 | 24 | 192 | matrix | 600 | 0.0 | 8.26 GB |

### Primitive Equation, Float32 vs Float64

| Model | NF | T | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | Float32 | 31 | 8 | 2400 | 874 | 5.38 MB |
| PrimitiveWetModel | Float64 | 31 | 8 | 2400 | 792 | 9.69 MB |

### Grids

| Model | T | L | Grid | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 63 | 8 | FullGaussianGrid | 96 | 1200 | 76 | 28.09 MB |
| PrimitiveWetModel | 63 | 8 | FullClenshawGrid | 95 | 1200 | 53 | 27.86 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralGaussianGrid | 96 | 1200 | 111 | 19.36 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralClenshawGrid | 95 | 1200 | 103 | 19.11 MB |
| PrimitiveWetModel | 63 | 8 | HEALPixGrid | 95 | 1200 | 154 | 14.14 MB |
| PrimitiveWetModel | 63 | 8 | OctaHEALPixGrid | 95 | 1200 | 68 | 17.21 MB |

### Number of vertical layers

| Model | T | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 4 | 2400 | 1290 | 3.41 MB |
| PrimitiveWetModel | 31 | 8 | 2400 | 872 | 5.38 MB |
| PrimitiveWetModel | 31 | 12 | 2400 | 627 | 7.36 MB |
| PrimitiveWetModel | 31 | 16 | 2400 | 520 | 9.35 MB |

### PrimitiveDryModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveDryModel | 31 | 8 | true | true | 2400 | 1477 | 4.70 MB |
| PrimitiveDryModel | 31 | 8 | true | false | 2400 | 2219 | 4.70 MB |
| PrimitiveDryModel | 31 | 8 | false | true | 2400 | 1902 | 4.70 MB |

### PrimitiveWetModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 8 | true | true | 2400 | 663 | 5.38 MB |
| PrimitiveWetModel | 31 | 8 | true | false | 2400 | 1873 | 5.38 MB |
| PrimitiveWetModel | 31 | 8 | false | true | 2400 | 1096 | 5.38 MB |

### Individual dynamics functions


#### PrimitiveWetModel | Float32 | T31 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| --- | --- | --- | --- |
| pressure_gradient_flux! | 62.492 μs| 49.08 KiB| 593 |
| linear_virtual_temperature! | 3.180 μs| 0 bytes| 0 |
| geopotential! | 11.040 μs| 816 bytes| 13 |
| vertical_integration! | 13.970 μs| 0 bytes| 0 |
| surface_pressure_tendency! | 26.411 μs| 25.41 KiB| 288 |
| vertical_velocity! | 62.792 μs| 0 bytes| 0 |
| linear_pressure_gradient! | 2.889 μs| 0 bytes| 0 |
| vertical_advection! | 166.656 μs| 3.94 KiB| 60 |
| vordiv_tendencies! | 403.772 μs| 247.45 KiB| 711 |
| temperature_tendency! | 543.987 μs| 362.72 KiB| 1012 |
| humidity_tendency! | 513.336 μs| 361.47 KiB| 1002 |
| bernoulli_potential! | 161.035 μs| 119.38 KiB| 327 |

## Architecture: `gpu-nvidia`

Created for SpeedyWeather.jl v0.21.1+DEV on Tue, 30 Jun 2026 16:41:06.

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
- cuBLAS: 13.5.1
- cuSPARSE: 12.8.1
- cuSOLVER: 12.2.2
- cuFFT: 12.3.0
- cuRAND: 10.4.3
- CUPTI: 2026.2.0 (API 13.3.0)
- NVML: 13.0.0+580.126.9

Julia packages: 
- CUDACore: 6.2.0
- GPUArrays: 11.5.8
- GPUCompiler: 1.22.7
- KernelAbstractions: 0.9.42
- CUDA_Driver_jll: 13.3.0+0
- CUDA_Compiler_jll: 0.4.4+0
- CUDA_Runtime_jll: 0.23.0+0
- NVPTX_LLVM_Backend_jll: 22.1.7+1

Toolchain:
- Julia: 1.12.2
- LLVM: 18.1.7

1 device:
  0: NVIDIA H100 80GB HBM3 (sm_90, 67.720 GiB / 79.647 GiB available)
     compiles to sm_90a / PTX 9.3 (LLVM: sm_90a / PTX 9.0)
```


### Models, default setups

| Model | T | L | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| BarotropicModel | 31 | 1 | false | 1800 | 1332 | 446.14 KB |
| ShallowWaterModel | 31 | 1 | false | 2400 | 882 | 448.94 KB |
| PrimitiveDryModel | 31 | 8 | true | 2400 | 310 | 604.45 KB |
| PrimitiveWetModel | 31 | 8 | true | 2400 | 1670 | 609.31 KB |

### Shallow water model, resolution

| Model | T | L | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 31 | 1 | 48 | 2400 | 880 | 448.94 KB |
| ShallowWaterModel | 42 | 1 | 64 | 1800 | 500 | 758.70 KB |
| ShallowWaterModel | 63 | 1 | 96 | 1200 | 219 | 1.63 MB |
| ShallowWaterModel | 85 | 1 | 128 | 900 | 123 | 2.84 MB |
| ShallowWaterModel | 127 | 1 | 192 | 600 | 57 | 6.26 MB |
| ShallowWaterModel | 170 | 1 | 256 | 450 | 31 | 11.03 MB |
| ShallowWaterModel | 255 | 1 | 384 | 300 | 15 | 24.59 MB |

### Primitive wet model, resolution

| Model | T | L | Rings | Transform | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 8 | 48 | default | 2400 | 1495 | 609.31 KB |
| PrimitiveWetModel | 42 | 8 | 64 | default | 1800 | 990 | 1.02 MB |
| PrimitiveWetModel | 63 | 8 | 96 | default | 1200 | 405 | 2.19 MB |
| PrimitiveWetModel | 85 | 8 | 128 | default | 900 | 225 | 3.81 MB |
| PrimitiveWetModel | 127 | 8 | 192 | default | 600 | 74 | 8.39 MB |
| PrimitiveWetModel | 170 | 8 | 256 | default | 450 | 34 | 14.78 MB |
| PrimitiveWetModel | 255 | 8 | 384 | default | 300 | 12 | 32.96 MB |
| PrimitiveWetModel | 85 | 16 | 128 | default | 900 | 236 | 4.85 MB |
| PrimitiveWetModel | 127 | 16 | 192 | default | 600 | 74 | 10.75 MB |
| PrimitiveWetModel | 170 | 16 | 256 | default | 450 | 31 | 18.97 MB |
| PrimitiveWetModel | 255 | 16 | 384 | default | 300 | 10 | 42.40 MB |
| PrimitiveWetModel | 85 | 24 | 128 | default | 900 | 244 | 5.90 MB |
| PrimitiveWetModel | 127 | 24 | 192 | default | 600 | 74 | 13.11 MB |
| PrimitiveWetModel | 170 | 24 | 256 | default | 450 | 25 | 23.17 MB |
| PrimitiveWetModel | 255 | 24 | 384 | default | 300 | 8.4 | 51.83 MB |
| PrimitiveWetModel | 31 | 8 | 48 | matrix | 2400 | 1684 | 560.65 KB |
| PrimitiveWetModel | 42 | 8 | 64 | matrix | 1800 | 429 | 959.12 KB |
| PrimitiveWetModel | 63 | 8 | 96 | matrix | 1200 | 175 | 2.09 MB |
| PrimitiveWetModel | 85 | 8 | 128 | matrix | 900 | 63 | 3.68 MB |
| PrimitiveWetModel | 127 | 8 | 192 | matrix | 600 | 19 | 8.20 MB |
| PrimitiveWetModel | 85 | 16 | 128 | matrix | 900 | 62 | 4.73 MB |
| PrimitiveWetModel | 127 | 16 | 192 | matrix | 600 | 18 | 10.56 MB |
| PrimitiveWetModel | 85 | 24 | 128 | matrix | 900 | 62 | 5.78 MB |
| PrimitiveWetModel | 127 | 24 | 192 | matrix | 600 | 18 | 12.92 MB |

### Primitive Equation, Float32 vs Float64

| Model | NF | T | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | Float32 | 31 | 8 | 2400 | 1486 | 609.31 KB |
| PrimitiveWetModel | Float64 | 31 | 8 | 2400 | 1365 | 609.92 KB |

### Grids

| Model | T | L | Grid | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 63 | 8 | FullGaussianGrid | 96 | 1200 | 397 | 2.28 MB |
| PrimitiveWetModel | 63 | 8 | FullClenshawGrid | 95 | 1200 | 388 | 2.26 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralGaussianGrid | 96 | 1200 | 400 | 2.19 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralClenshawGrid | 95 | 1200 | 393 | 2.17 MB |
| PrimitiveWetModel | 63 | 8 | HEALPixGrid | 95 | 1200 | 397 | 2.12 MB |
| PrimitiveWetModel | 63 | 8 | OctaHEALPixGrid | 95 | 1200 | 404 | 2.15 MB |

### Number of vertical layers

| Model | T | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 4 | 2400 | 1371 | 535.58 KB |
| PrimitiveWetModel | 31 | 8 | 2400 | 1498 | 609.31 KB |
| PrimitiveWetModel | 31 | 12 | 2400 | 1501 | 683.04 KB |
| PrimitiveWetModel | 31 | 16 | 2400 | 1657 | 756.76 KB |

### PrimitiveDryModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveDryModel | 31 | 8 | true | true | 2400 | 299 | 604.45 KB |
| PrimitiveDryModel | 31 | 8 | true | false | 2400 | 295 | 604.45 KB |
| PrimitiveDryModel | 31 | 8 | false | true | 2400 | 321 | 604.45 KB |

### PrimitiveWetModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 8 | true | true | 2400 | 1477 | 609.31 KB |
| PrimitiveWetModel | 31 | 8 | true | false | 2400 | 1691 | 609.31 KB |
| PrimitiveWetModel | 31 | 8 | false | true | 2400 | 6830 | 609.31 KB |

### Individual dynamics functions


#### PrimitiveWetModel | Float32 | T31 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| --- | --- | --- | --- |
| pressure_gradient_flux! | 1.427 ms| 192.58 KiB| 6472 |
| linear_virtual_temperature! | 13.340 μs| 2.77 KiB| 52 |
| geopotential! | 21.680 μs| 4.09 KiB| 108 |
| vertical_integration! | 32.090 μs| 8.83 KiB| 163 |
| surface_pressure_tendency! | 552.028 μs| 95.73 KiB| 3105 |
| vertical_velocity! | 27.050 μs| 8.45 KiB| 191 |
| linear_pressure_gradient! | 13.089 μs| 2.34 KiB| 50 |
| vertical_advection! | 31.859 μs| 9.22 KiB| 170 |
| vordiv_tendencies! | 298.509 μs| 28.28 KiB| 482 |
| temperature_tendency! | 435.258 μs| 30.94 KiB| 635 |
| humidity_tendency! | 434.638 μs| 27.83 KiB| 545 |
| bernoulli_potential! | 151.889 μs| 14.17 KiB| 304 |

