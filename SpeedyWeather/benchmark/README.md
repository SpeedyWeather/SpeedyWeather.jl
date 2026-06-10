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
| 31 | 8 | LT+FFT | 1484 | 959 | 1137 |
| 31 | 8 | MT | 484 | 71 | 1480 |
| 42 | 8 | LT+FFT | 662 | 413 | 656 |
| 42 | 8 | MT | 173 | 19 | 278 |
| 63 | 8 | LT+FFT | 188 | 115 | 380 |
| 63 | 8 | MT | 34 | 2.8 | 148 |
| 85 | 8 | LT+FFT | 65 | 54 | 213 |
| 85 | 8 | MT | 8.7 | 0.7 | 76 |
| 85 | 16 | LT+FFT | 42 | 28 | 242 |
| 85 | 16 | MT | 7.1 | 0.4 | 77 |
| 85 | 24 | LT+FFT | 23 | 17 | 231 |
| 85 | 24 | MT | 5.2 | 0.3 | 77 |
| 127 | 8 | LT+FFT | 18 | 13 | 64 |
| 127 | 8 | MT | 0.6 | 0.1 | 18 |
| 127 | 16 | LT+FFT | 10.0 | 7.4 | 71 |
| 127 | 16 | MT | 1.1 | 0.1 | 18 |
| 127 | 24 | LT+FFT | 6.2 | 4.5 | 67 |
| 127 | 24 | MT | 0.9 | 0.0 | 17 |
| 170 | 8 | LT+FFT | 7.0 | 5.1 | 33 |
| 170 | 16 | LT+FFT | 3.7 | 2.8 | 29 |
| 170 | 24 | LT+FFT | 2.2 | 1.7 | 27 |
| 255 | 8 | LT+FFT | 1.8 | 1.3 | 14 |
| 255 | 16 | LT+FFT | 0.9 | 0.7 | 11 |
| 255 | 24 | LT+FFT | 0.6 | 0.4 | 8.7 |

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

Created for SpeedyWeather.jl v0.20.3 on Sat, 06 Jun 2026 21:45:16.

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
| BarotropicModel | 31 | 1 | false | 2400 | 28443 | 688.47 KB |
| ShallowWaterModel | 31 | 1 | false | 2400 | 16131 | 822.99 KB |
| PrimitiveDryModel | 31 | 8 | true | 2400 | 1046 | 4.16 MB |
| PrimitiveWetModel | 31 | 8 | true | 2400 | 888 | 4.84 MB |

### Shallow water model, resolution

| Model | T | L | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 31 | 1 | 48 | 2400 | 15640 | 822.99 KB |
| ShallowWaterModel | 42 | 1 | 64 | 1800 | 7227 | 1.44 MB |
| ShallowWaterModel | 63 | 1 | 96 | 1200 | 2208 | 3.25 MB |
| ShallowWaterModel | 85 | 1 | 128 | 900 | 912 | 5.94 MB |
| ShallowWaterModel | 127 | 1 | 192 | 600 | 186 | 14.14 MB |
| ShallowWaterModel | 170 | 1 | 256 | 450 | 91 | 26.86 MB |
| ShallowWaterModel | 255 | 1 | 384 | 300 | 20 | 68.30 MB |

### Primitive wet model, resolution

| Model | T | L | Rings | Transform | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 8 | 48 | default | 2400 | 959 | 4.84 MB |
| PrimitiveWetModel | 42 | 8 | 64 | default | 1800 | 413 | 8.18 MB |
| PrimitiveWetModel | 63 | 8 | 96 | default | 1200 | 115 | 17.42 MB |
| PrimitiveWetModel | 85 | 8 | 128 | default | 900 | 54 | 30.35 MB |
| PrimitiveWetModel | 127 | 8 | 192 | default | 600 | 13 | 67.03 MB |
| PrimitiveWetModel | 170 | 8 | 256 | default | 450 | 5.1 | 119.25 MB |
| PrimitiveWetModel | 255 | 8 | 384 | default | 300 | 1.3 | 272.27 MB |
| PrimitiveWetModel | 85 | 16 | 128 | default | 900 | 28 | 51.65 MB |
| PrimitiveWetModel | 127 | 16 | 192 | default | 600 | 7.4 | 113.22 MB |
| PrimitiveWetModel | 170 | 16 | 256 | default | 450 | 2.8 | 200.00 MB |
| PrimitiveWetModel | 255 | 16 | 384 | default | 300 | 0.7 | 450.63 MB |
| PrimitiveWetModel | 85 | 24 | 128 | default | 900 | 17 | 73.00 MB |
| PrimitiveWetModel | 127 | 24 | 192 | default | 600 | 4.5 | 159.48 MB |
| PrimitiveWetModel | 170 | 24 | 256 | default | 450 | 1.7 | 280.85 MB |
| PrimitiveWetModel | 255 | 24 | 384 | default | 300 | 0.4 | 629.13 MB |
| PrimitiveWetModel | 31 | 8 | 48 | matrix | 2400 | 71 | 47.07 MB |
| PrimitiveWetModel | 42 | 8 | 64 | matrix | 1800 | 19 | 132.12 MB |
| PrimitiveWetModel | 63 | 8 | 96 | matrix | 1200 | 2.8 | 579.07 MB |
| PrimitiveWetModel | 85 | 8 | 128 | matrix | 900 | 0.7 | 1.74 GB |
| PrimitiveWetModel | 127 | 8 | 192 | matrix | 600 | 0.1 | 8.17 GB |
| PrimitiveWetModel | 85 | 16 | 128 | matrix | 900 | 0.4 | 1.76 GB |
| PrimitiveWetModel | 127 | 16 | 192 | matrix | 600 | 0.1 | 8.22 GB |
| PrimitiveWetModel | 85 | 24 | 128 | matrix | 900 | 0.3 | 1.78 GB |
| PrimitiveWetModel | 127 | 24 | 192 | matrix | 600 | 0.0 | 8.26 GB |

### Primitive Equation, Float32 vs Float64

| Model | NF | T | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | Float32 | 31 | 8 | 2400 | 717 | 4.84 MB |
| PrimitiveWetModel | Float64 | 31 | 8 | 2400 | 787 | 9.00 MB |

### Grids

| Model | T | L | Grid | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 63 | 8 | FullGaussianGrid | 96 | 1200 | 82 | 25.36 MB |
| PrimitiveWetModel | 63 | 8 | FullClenshawGrid | 95 | 1200 | 83 | 25.14 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralGaussianGrid | 96 | 1200 | 115 | 17.42 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralClenshawGrid | 95 | 1200 | 74 | 17.18 MB |
| PrimitiveWetModel | 63 | 8 | HEALPixGrid | 95 | 1200 | 188 | 12.73 MB |
| PrimitiveWetModel | 63 | 8 | OctaHEALPixGrid | 95 | 1200 | 75 | 15.50 MB |

### Number of vertical layers

| Model | T | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 4 | 2400 | 1450 | 3.11 MB |
| PrimitiveWetModel | 31 | 8 | 2400 | 948 | 4.84 MB |
| PrimitiveWetModel | 31 | 12 | 2400 | 639 | 6.58 MB |
| PrimitiveWetModel | 31 | 16 | 2400 | 546 | 8.33 MB |

### PrimitiveDryModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveDryModel | 31 | 8 | true | true | 2400 | 1477 | 4.16 MB |
| PrimitiveDryModel | 31 | 8 | true | false | 2400 | 2247 | 4.16 MB |
| PrimitiveDryModel | 31 | 8 | false | true | 2400 | 1899 | 4.16 MB |

### PrimitiveWetModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 8 | true | true | 2400 | 960 | 4.84 MB |
| PrimitiveWetModel | 31 | 8 | true | false | 2400 | 1682 | 4.84 MB |
| PrimitiveWetModel | 31 | 8 | false | true | 2400 | 1211 | 4.84 MB |

### Individual dynamics functions


#### PrimitiveWetModel | Float32 | T31 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| --- | --- | --- | --- |
| pressure_gradient_flux! | 67.850 μs| 100.30 KiB| 789 |
| linear_virtual_temperature! | 3.405 μs| 0 bytes| 0 |
| geopotential! | 11.580 μs| 5.51 KiB| 22 |
| vertical_integration! | 14.010 μs| 0 bytes| 0 |
| surface_pressure_tendency! | 24.450 μs| 24.66 KiB| 288 |
| vertical_velocity! | 86.940 μs| 346.84 KiB| 12 |
| linear_pressure_gradient! | 3.048 μs| 0 bytes| 0 |
| vertical_advection! | 228.479 μs| 8.69 KiB| 96 |
| vordiv_tendencies! | 395.089 μs| 246.36 KiB| 721 |
| temperature_tendency! | 512.349 μs| 361.98 KiB| 1024 |
| humidity_tendency! | 483.629 μs| 360.62 KiB| 1014 |
| bernoulli_potential! | 172.839 μs| 416.55 KiB| 340 |

## Architecture: `gpu-nvidia`

Created for SpeedyWeather.jl v0.20.3 on Sat, 06 Jun 2026 22:14:51.

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
- runtime 13.2, artifact installation
- driver 580.126.9 for 13.3
- compiler 13.3, artifact installation

CUDA libraries: 
- cuBLAS: 13.4.0
- cuSPARSE: 12.7.10
- cuSOLVER: 12.2.0
- cuFFT: 12.2.0
- cuRAND: 10.4.2
- CUPTI: 2026.1.1 (API 13.2.1)
- NVML: 13.0.0+580.126.9

Julia packages: 
- CUDACore: 6.1.1
- GPUArrays: 11.5.5
- GPUCompiler: 1.17.1
- KernelAbstractions: 0.9.41
- CUDA_Driver_jll: 13.3.0+0
- CUDA_Compiler_jll: 0.4.4+0
- CUDA_Runtime_jll: 0.21.0+1

Toolchain:
- Julia: 1.12.2
- LLVM: 18.1.7

1 device:
  0: NVIDIA H100 80GB HBM3 (sm_90, 47.753 GiB / 79.647 GiB available)
```


### Models, default setups

| Model | T | L | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| BarotropicModel | 31 | 1 | false | 2400 | 1701 | 432.95 KB |
| ShallowWaterModel | 31 | 1 | false | 2400 | 850 | 436.08 KB |
| PrimitiveDryModel | 31 | 8 | true | 2400 | 1505 | 579.99 KB |
| PrimitiveWetModel | 31 | 8 | true | 2400 | 1235 | 585.21 KB |

### Shallow water model, resolution

| Model | T | L | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 31 | 1 | 48 | 2400 | 827 | 436.08 KB |
| ShallowWaterModel | 42 | 1 | 64 | 1800 | 478 | 741.88 KB |
| ShallowWaterModel | 63 | 1 | 96 | 1200 | 198 | 1.61 MB |
| ShallowWaterModel | 85 | 1 | 128 | 900 | 131 | 2.80 MB |
| ShallowWaterModel | 127 | 1 | 192 | 600 | 59 | 6.21 MB |
| ShallowWaterModel | 170 | 1 | 256 | 450 | 33 | 10.96 MB |
| ShallowWaterModel | 255 | 1 | 384 | 300 | 15 | 24.49 MB |

### Primitive wet model, resolution

| Model | T | L | Rings | Transform | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 8 | 48 | default | 2400 | 1137 | 585.21 KB |
| PrimitiveWetModel | 42 | 8 | 64 | default | 1800 | 656 | 991.36 KB |
| PrimitiveWetModel | 63 | 8 | 96 | default | 1200 | 380 | 2.14 MB |
| PrimitiveWetModel | 85 | 8 | 128 | default | 900 | 213 | 3.74 MB |
| PrimitiveWetModel | 127 | 8 | 192 | default | 600 | 64 | 8.30 MB |
| PrimitiveWetModel | 170 | 8 | 256 | default | 450 | 33 | 14.65 MB |
| PrimitiveWetModel | 255 | 8 | 384 | default | 300 | 14 | 32.77 MB |
| PrimitiveWetModel | 85 | 16 | 128 | default | 900 | 242 | 4.79 MB |
| PrimitiveWetModel | 127 | 16 | 192 | default | 600 | 71 | 10.65 MB |
| PrimitiveWetModel | 170 | 16 | 256 | default | 450 | 29 | 18.85 MB |
| PrimitiveWetModel | 255 | 16 | 384 | default | 300 | 11 | 42.21 MB |
| PrimitiveWetModel | 85 | 24 | 128 | default | 900 | 231 | 5.84 MB |
| PrimitiveWetModel | 127 | 24 | 192 | default | 600 | 67 | 13.01 MB |
| PrimitiveWetModel | 170 | 24 | 256 | default | 450 | 27 | 23.04 MB |
| PrimitiveWetModel | 255 | 24 | 384 | default | 300 | 8.7 | 51.64 MB |
| PrimitiveWetModel | 31 | 8 | 48 | matrix | 2400 | 1480 | 561.28 KB |
| PrimitiveWetModel | 42 | 8 | 64 | matrix | 1800 | 278 | 959.74 KB |
| PrimitiveWetModel | 63 | 8 | 96 | matrix | 1200 | 148 | 2.09 MB |
| PrimitiveWetModel | 85 | 8 | 128 | matrix | 900 | 76 | 3.68 MB |
| PrimitiveWetModel | 127 | 8 | 192 | matrix | 600 | 18 | 8.20 MB |
| PrimitiveWetModel | 85 | 16 | 128 | matrix | 900 | 77 | 4.73 MB |
| PrimitiveWetModel | 127 | 16 | 192 | matrix | 600 | 18 | 10.56 MB |
| PrimitiveWetModel | 85 | 24 | 128 | matrix | 900 | 77 | 5.78 MB |
| PrimitiveWetModel | 127 | 24 | 192 | matrix | 600 | 17 | 12.92 MB |

### Primitive Equation, Float32 vs Float64

| Model | NF | T | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | Float32 | 31 | 8 | 2400 | 1097 | 585.21 KB |
| PrimitiveWetModel | Float64 | 31 | 8 | 2400 | 599 | 585.85 KB |

### Grids

| Model | T | L | Grid | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 63 | 8 | FullGaussianGrid | 96 | 1200 | 261 | 2.23 MB |
| PrimitiveWetModel | 63 | 8 | FullClenshawGrid | 95 | 1200 | 386 | 2.21 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralGaussianGrid | 96 | 1200 | 378 | 2.14 MB |
| PrimitiveWetModel | 63 | 8 | OctahedralClenshawGrid | 95 | 1200 | 379 | 2.12 MB |
| PrimitiveWetModel | 63 | 8 | HEALPixGrid | 95 | 1200 | 374 | 2.07 MB |
| PrimitiveWetModel | 63 | 8 | OctaHEALPixGrid | 95 | 1200 | 388 | 2.10 MB |

### Number of vertical layers

| Model | T | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 4 | 2400 | 993 | 511.49 KB |
| PrimitiveWetModel | 31 | 8 | 2400 | 1116 | 585.21 KB |
| PrimitiveWetModel | 31 | 12 | 2400 | 1179 | 658.94 KB |
| PrimitiveWetModel | 31 | 16 | 2400 | 1111 | 732.67 KB |

### PrimitiveDryModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveDryModel | 31 | 8 | true | true | 2400 | 1323 | 579.99 KB |
| PrimitiveDryModel | 31 | 8 | true | false | 2400 | 1353 | 579.99 KB |
| PrimitiveDryModel | 31 | 8 | false | true | 2400 | 3412 | 579.99 KB |

### PrimitiveWetModel: Physics or dynamics only

| Model | T | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 31 | 8 | true | true | 2400 | 1125 | 585.21 KB |
| PrimitiveWetModel | 31 | 8 | true | false | 2400 | 1214 | 585.21 KB |
| PrimitiveWetModel | 31 | 8 | false | true | 2400 | 2735 | 585.21 KB |

### Individual dynamics functions


#### PrimitiveWetModel | Float32 | T31 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| --- | --- | --- | --- |
| pressure_gradient_flux! | 1.499 ms| 191.42 KiB| 6470 |
| linear_virtual_temperature! | 15.550 μs| 2.75 KiB| 72 |
| geopotential! | 36.020 μs| 6.72 KiB| 196 |
| vertical_integration! | 34.330 μs| 7.58 KiB| 170 |
| surface_pressure_tendency! | 617.227 μs| 95.48 KiB| 3097 |
| vertical_velocity! | 59.590 μs| 15.16 KiB| 403 |
| linear_pressure_gradient! | 14.400 μs| 2.28 KiB| 58 |
| vertical_advection! | 32.750 μs| 12.28 KiB| 190 |
| vordiv_tendencies! | 302.219 μs| 26.31 KiB| 450 |
| temperature_tendency! | 443.968 μs| 30.69 KiB| 724 |
| humidity_tendency! | 434.519 μs| 25.92 KiB| 578 |
| bernoulli_potential! | 176.190 μs| 20.30 KiB| 488 |

