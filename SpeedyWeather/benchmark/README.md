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
- Compile: One-off compilation of the time step (Reactant only), not included in the SYPD. Compiled code is cached within a session, so configurations compiled before (e.g. in an earlier suite) compile much faster.

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

| T | L | Transform | cpu-arm | cpu-x86 | gpu-nvidia | gpu-amd | reactant-gpu |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 32 | 8 | LT+FFT | 1308 | 870 | 5961 | 3763 | 28529 |
| 32 | 8 | MT | 772 | 882 | 5853 | 3696 | 19398 |
| 43 | 8 | LT+FFT | 559 | 371 | 4522 | 2149 | 15354 |
| 43 | 8 | MT | 253 | 337 | 4468 | 2830 | 13568 |
| 64 | 8 | LT+FFT | 152 | 108 | 3025 | 828 | 9160 |
| 64 | 8 | MT | 52 | 72 | 3117 | 826 | 6487 |
| 86 | 8 | LT+FFT | 58 | 43 | 1254 | 246 | 5394 |
| 86 | 8 | MT | 14 | 23 | 1260 | 245 | 2012 |
| 86 | 16 | LT+FFT | 29 | 21 | 781 | 202 | 3033 |
| 86 | 16 | MT | 8.8 | 13 | 783 | 150 | 1688 |
| 86 | 24 | LT+FFT | 20 | 28 | 1357 | 196 | 2083 |
| 86 | 24 | MT | 6.1 | 19 | 1358 | 108 | 1248 |
| 128 | 8 | LT+FFT | 14 | 11 | 366 | 86 | 2260 |
| 128 | 8 | MT | 1.2 | — | 238 | 39 | 664 |
| 128 | 16 | LT+FFT | 7.6 | 11 | 612 | 80 | 1359 |
| 128 | 16 | MT | 1.2 | — | 288 | 23 | 550 |
| 128 | 24 | LT+FFT | 5.5 | 7.3 | 553 | 76 | 901 |
| 128 | 24 | MT | 1.0 | — | 255 | 17 | 473 |
| 171 | 8 | LT+FFT | 5.7 | 4.1 | 188 | 44 | 1080 |
| 171 | 8 | MT | — | — | 52 | — | 165 |
| 171 | 16 | LT+FFT | 3.1 | 4.3 | 313 | 41 | 610 |
| 171 | 16 | MT | — | — | 72 | — | 101 |
| 171 | 24 | LT+FFT | 1.9 | 2.9 | 273 | 38 | 396 |
| 171 | 24 | MT | — | — | 63 | — | 81 |
| 256 | 8 | LT+FFT | 1.4 | 1.0 | 69 | 16 | 362 |
| 256 | 16 | LT+FFT | 0.8 | 1.0 | 109 | 13 | 194 |
| 256 | 24 | LT+FFT | 0.5 | 0.7 | 91 | 11 | 134 |

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

Created for SpeedyWeather.jl v0.24.0-DEV on Mon, 19 Oct 2026 19:11:30.

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
| BarotropicModel | 32 | 1 | false | 1800 | 14105 | 725.92 KB |
| ShallowWaterModel | 32 | 1 | false | 2400 | 19255 | 908.28 KB |
| PrimitiveDryModel | 32 | 8 | true | 2400 | 1405 | 4.80 MB |
| PrimitiveWetModel | 32 | 8 | true | 2400 | 850 | 5.85 MB |

### Shallow water model, resolution

| Model | truncation | L | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 32 | 1 | 48 | 2400 | 18400 | 908.28 KB |
| ShallowWaterModel | 43 | 1 | 64 | 1800 | 7955 | 1.59 MB |
| ShallowWaterModel | 64 | 1 | 96 | 1200 | 2204 | 3.57 MB |
| ShallowWaterModel | 86 | 1 | 128 | 900 | 913 | 6.49 MB |
| ShallowWaterModel | 128 | 1 | 192 | 600 | 233 | 15.36 MB |
| ShallowWaterModel | 171 | 1 | 256 | 450 | 85 | 29.00 MB |
| ShallowWaterModel | 256 | 1 | 384 | 300 | 19 | 73.09 MB |

### Primitive wet model, resolution

| Model | truncation | L | Rings | Transform | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | 48 | default | 2400 | 870 | 5.85 MB |
| PrimitiveWetModel | 43 | 8 | 64 | default | 1800 | 371 | 9.86 MB |
| PrimitiveWetModel | 64 | 8 | 96 | default | 1200 | 108 | 20.96 MB |
| PrimitiveWetModel | 86 | 8 | 128 | default | 900 | 43 | 36.46 MB |
| PrimitiveWetModel | 128 | 8 | 192 | default | 600 | 11 | 80.27 MB |
| PrimitiveWetModel | 171 | 8 | 256 | default | 450 | 4.1 | 142.42 MB |
| PrimitiveWetModel | 256 | 8 | 384 | default | 300 | 1.0 | 323.46 MB |
| PrimitiveWetModel | 86 | 16 | 128 | default | 900 | 21 | 62.29 MB |
| PrimitiveWetModel | 128 | 16 | 192 | default | 600 | 11 | 136.27 MB |
| PrimitiveWetModel | 171 | 16 | 256 | default | 450 | 4.3 | 240.35 MB |
| PrimitiveWetModel | 256 | 16 | 384 | default | 300 | 1.0 | 539.78 MB |
| PrimitiveWetModel | 86 | 24 | 128 | default | 900 | 28 | 88.16 MB |
| PrimitiveWetModel | 128 | 24 | 192 | default | 600 | 7.3 | 192.34 MB |
| PrimitiveWetModel | 171 | 24 | 256 | default | 450 | 2.9 | 338.38 MB |
| PrimitiveWetModel | 256 | 24 | 384 | default | 300 | 0.7 | 756.24 MB |
| PrimitiveWetModel | 32 | 8 | 48 | matrix | 2400 | 882 | 34.29 MB |
| PrimitiveWetModel | 43 | 8 | 64 | matrix | 1800 | 337 | 93.00 MB |
| PrimitiveWetModel | 64 | 8 | 96 | matrix | 1200 | 72 | 396.46 MB |
| PrimitiveWetModel | 86 | 8 | 128 | matrix | 900 | 23 | 1.18 GB |
| PrimitiveWetModel | 86 | 16 | 128 | matrix | 900 | 13 | 1.21 GB |
| PrimitiveWetModel | 86 | 24 | 128 | matrix | 900 | 19 | 1.23 GB |

### Primitive Equation, Float32 vs Float64

| Model | NF | truncation | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | Float32 | 32 | 8 | 2400 | 865 | 5.85 MB |
| PrimitiveWetModel | Float64 | 32 | 8 | 2400 | 728 | 11.16 MB |

### Grids

| Model | truncation | L | Grid | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 64 | 8 | FullGaussianGrid | 96 | 1200 | 70 | 30.74 MB |
| PrimitiveWetModel | 64 | 8 | FullClenshawGrid | 127 | 1200 | 42 | 51.35 MB |
| PrimitiveWetModel | 64 | 8 | OctahedralGaussianGrid | 96 | 1200 | 108 | 20.96 MB |
| PrimitiveWetModel | 64 | 8 | OctahedralClenshawGrid | 127 | 1200 | 59 | 32.92 MB |
| PrimitiveWetModel | 64 | 8 | HEALPixGrid | 127 | 1200 | 90 | 24.29 MB |
| PrimitiveWetModel | 64 | 8 | OctaHEALPixGrid | 127 | 1200 | 59 | 30.19 MB |

### Number of vertical layers

| Model | truncation | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 4 | 2400 | 1398 | 3.74 MB |
| PrimitiveWetModel | 32 | 8 | 2400 | 870 | 5.85 MB |
| PrimitiveWetModel | 32 | 12 | 2400 | 609 | 7.95 MB |
| PrimitiveWetModel | 32 | 16 | 2400 | 478 | 10.07 MB |

### PrimitiveDryModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveDryModel | 32 | 8 | true | true | 2400 | 1402 | 4.80 MB |
| PrimitiveDryModel | 32 | 8 | true | false | 2400 | 2313 | 4.80 MB |
| PrimitiveDryModel | 32 | 8 | false | true | 2400 | 1483 | 4.80 MB |

### PrimitiveWetModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | true | true | 2400 | 874 | 5.85 MB |
| PrimitiveWetModel | 32 | 8 | true | false | 2400 | 1782 | 5.85 MB |
| PrimitiveWetModel | 32 | 8 | false | true | 2400 | 965 | 5.85 MB |

### Individual dynamics functions


#### PrimitiveWetModel | Float32 | T31 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| --- | --- | --- | --- |
| pressure_gradient_flux! | 16.109 μs| 0 bytes| 0 |
| linear_virtual_temperature! | 3.374 μs| 0 bytes| 0 |
| geopotential! | 9.780 μs| 384 bytes| 6 |
| vertical_integration! | 13.770 μs| 0 bytes| 0 |
| surface_pressure_tendency! | 18.400 μs| 15.66 KiB| 96 |
| vertical_velocity! | 58.340 μs| 0 bytes| 0 |
| linear_pressure_gradient! | 3.034 μs| 0 bytes| 0 |
| vertical_advection! | 160.239 μs| 2.44 KiB| 32 |
| vordiv_tendencies! | 416.997 μs| 218.62 KiB| 284 |
| temperature_tendency! | 522.466 μs| 324.75 KiB| 401 |
| humidity_tendency! | 497.646 μs| 324.08 KiB| 396 |
| bernoulli_potential! | 163.429 μs| 107.61 KiB| 129 |

## Architecture: `gpu-nvidia`

Created for SpeedyWeather.jl v0.24.0-DEV on Tue, 20 Oct 2026 16:34:44.

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
- GPUArrays: 11.5.16
- GPUCompiler: 1.23.0
- KernelAbstractions: 0.9.45
- CUDA_Driver_jll: 13.3.1+0
- CUDA_Compiler_jll: 0.4.4+1
- CUDA_Runtime_jll: 0.23.0+1
- NVPTX_LLVM_Backend_jll: 22.1.7+1

Toolchain:
- Julia: 1.12.2
- LLVM: 18.1.7

1 device:
  0: NVIDIA H100 80GB HBM3 (sm_90, 77.766 GiB / 79.647 GiB available)
     compiles to sm_90a / PTX 9.3 (LLVM: sm_90a / PTX 9.0)
```


### Models, default setups

| Model | truncation | L | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| BarotropicModel | 32 | 1 | false | 1800 | 16667 | 207.93 KB |
| ShallowWaterModel | 32 | 1 | false | 2400 | 11335 | 211.51 KB |
| PrimitiveDryModel | 32 | 8 | true | 2400 | 7786 | 355.38 KB |
| PrimitiveWetModel | 32 | 8 | true | 2400 | 6187 | 361.68 KB |

### Shallow water model, resolution

| Model | truncation | L | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 32 | 1 | 48 | 2400 | 11601 | 211.51 KB |
| ShallowWaterModel | 43 | 1 | 64 | 1800 | 8517 | 355.38 KB |
| ShallowWaterModel | 64 | 1 | 96 | 1200 | 4469 | 766.00 KB |
| ShallowWaterModel | 86 | 1 | 128 | 900 | 1652 | 1.34 MB |
| ShallowWaterModel | 128 | 1 | 192 | 600 | 296 | 3.31 MB |
| ShallowWaterModel | 171 | 1 | 256 | 450 | 132 | 5.72 MB |
| ShallowWaterModel | 256 | 1 | 384 | 300 | 52 | 12.50 MB |

### Primitive wet model, resolution

| Model | truncation | L | Rings | Transform | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | 48 | default | 2400 | 5961 | 361.68 KB |
| PrimitiveWetModel | 43 | 8 | 64 | default | 1800 | 4522 | 605.90 KB |
| PrimitiveWetModel | 64 | 8 | 96 | default | 1200 | 3025 | 1.30 MB |
| PrimitiveWetModel | 86 | 8 | 128 | default | 900 | 1254 | 2.28 MB |
| PrimitiveWetModel | 128 | 8 | 192 | default | 600 | 366 | 5.45 MB |
| PrimitiveWetModel | 171 | 8 | 256 | default | 450 | 188 | 9.47 MB |
| PrimitiveWetModel | 256 | 8 | 384 | default | 300 | 69 | 20.87 MB |
| PrimitiveWetModel | 86 | 16 | 128 | default | 900 | 781 | 3.33 MB |
| PrimitiveWetModel | 128 | 16 | 192 | default | 600 | 612 | 7.81 MB |
| PrimitiveWetModel | 171 | 16 | 256 | default | 450 | 313 | 13.67 MB |
| PrimitiveWetModel | 256 | 16 | 384 | default | 300 | 109 | 30.31 MB |
| PrimitiveWetModel | 86 | 24 | 128 | default | 900 | 1357 | 4.38 MB |
| PrimitiveWetModel | 128 | 24 | 192 | default | 600 | 553 | 10.17 MB |
| PrimitiveWetModel | 171 | 24 | 256 | default | 450 | 273 | 17.86 MB |
| PrimitiveWetModel | 256 | 24 | 384 | default | 300 | 91 | 39.75 MB |
| PrimitiveWetModel | 32 | 8 | 48 | matrix | 2400 | 5853 | 361.68 KB |
| PrimitiveWetModel | 43 | 8 | 64 | matrix | 1800 | 4468 | 605.90 KB |
| PrimitiveWetModel | 64 | 8 | 96 | matrix | 1200 | 3117 | 1.30 MB |
| PrimitiveWetModel | 86 | 8 | 128 | matrix | 900 | 1260 | 2.28 MB |
| PrimitiveWetModel | 128 | 8 | 192 | matrix | 600 | 238 | 5.07 MB |
| PrimitiveWetModel | 171 | 8 | 256 | matrix | 450 | 52 | 8.97 MB |
| PrimitiveWetModel | 86 | 16 | 128 | matrix | 900 | 783 | 3.33 MB |
| PrimitiveWetModel | 128 | 16 | 192 | matrix | 600 | 288 | 7.43 MB |
| PrimitiveWetModel | 171 | 16 | 256 | matrix | 450 | 72 | 13.16 MB |
| PrimitiveWetModel | 86 | 24 | 128 | matrix | 900 | 1358 | 4.38 MB |
| PrimitiveWetModel | 128 | 24 | 192 | matrix | 600 | 255 | 9.79 MB |
| PrimitiveWetModel | 171 | 24 | 256 | matrix | 450 | 63 | 17.36 MB |

### Primitive Equation, Float32 vs Float64

| Model | NF | truncation | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | Float32 | 32 | 8 | 2400 | 5915 | 361.68 KB |
| PrimitiveWetModel | Float64 | 32 | 8 | 2400 | 4853 | 362.28 KB |

### Grids

| Model | truncation | L | Grid | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 64 | 8 | FullGaussianGrid | 96 | 1200 | 2167 | 1.30 MB |
| PrimitiveWetModel | 64 | 8 | FullClenshawGrid | 127 | 1200 | 1668 | 2.26 MB |
| PrimitiveWetModel | 64 | 8 | OctahedralGaussianGrid | 96 | 1200 | 3117 | 1.30 MB |
| PrimitiveWetModel | 64 | 8 | OctahedralClenshawGrid | 127 | 1200 | 2409 | 2.26 MB |
| PrimitiveWetModel | 64 | 8 | HEALPixGrid | 127 | 1200 | 2852 | 2.26 MB |
| PrimitiveWetModel | 64 | 8 | OctaHEALPixGrid | 127 | 1200 | 2631 | 2.26 MB |

### Number of vertical layers

| Model | truncation | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 4 | 2400 | 6284 | 287.95 KB |
| PrimitiveWetModel | 32 | 8 | 2400 | 6290 | 361.68 KB |
| PrimitiveWetModel | 32 | 12 | 2400 | 6216 | 435.41 KB |
| PrimitiveWetModel | 32 | 16 | 2400 | 6441 | 509.13 KB |

### PrimitiveDryModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveDryModel | 32 | 8 | true | true | 2400 | 8171 | 355.38 KB |
| PrimitiveDryModel | 32 | 8 | true | false | 2400 | 9289 | 355.38 KB |
| PrimitiveDryModel | 32 | 8 | false | true | 2400 | 13432 | 355.38 KB |

### PrimitiveWetModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | true | true | 2400 | 6249 | 361.68 KB |
| PrimitiveWetModel | 32 | 8 | true | false | 2400 | 8004 | 361.68 KB |
| PrimitiveWetModel | 32 | 8 | false | true | 2400 | 9098 | 361.68 KB |

### Individual dynamics functions


#### PrimitiveWetModel | Float32 | T31 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| --- | --- | --- | --- |
| pressure_gradient_flux! | 13.480 μs| 4.25 KiB| 59 |
| linear_virtual_temperature! | 12.770 μs| 2.77 KiB| 52 |
| geopotential! | 18.980 μs| 3.73 KiB| 101 |
| vertical_integration! | 27.290 μs| 7.73 KiB| 143 |
| surface_pressure_tendency! | N/A| N/A| N/A |
| vertical_velocity! | 24.570 μs| 9.03 KiB| 191 |
| linear_pressure_gradient! | 12.550 μs| 2.34 KiB| 50 |
| vertical_advection! | 24.191 μs| 9.09 KiB| 142 |
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

## Architecture: `reactant-gpu`

Created for SpeedyWeather.jl v0.24.0-DEV on Tue, 20 Oct 2026 19:33:02.

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
- GPUArrays: 11.5.16
- GPUCompiler: 1.23.0
- KernelAbstractions: 0.9.45
- CUDA_Driver_jll: 13.3.1+0
- CUDA_Compiler_jll: 0.4.4+1
- CUDA_Runtime_jll: 0.23.0+1
- NVPTX_LLVM_Backend_jll: 22.1.7+1

Toolchain:
- Julia: 1.12.2
- LLVM: 18.1.7

1 device:
  0: NVIDIA H100 80GB HBM3 (sm_90, 19.124 GiB / 79.647 GiB available)
     compiles to sm_90a / PTX 9.3 (LLVM: sm_90a / PTX 9.0)
```

Reactant backend: `gpu` (via `Reactant.set_default_backend("gpu")`), Reactant v0.2.295, Reactant_jll v0.0.424+0.


### Models, default setups

| Model | truncation | L | Physics | Δt | SYPD | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- |
| BarotropicModel | 32 | 1 | false | 1800 | 172596 | 1.81 MB | 30.7 s |
| ShallowWaterModel | 32 | 1 | false | 2400 | 128561 | 2.00 MB | 38.5 s |
| PrimitiveDryModel | 32 | 8 | true | 2400 | 36319 | 5.60 MB | 176.7 s |
| PrimitiveWetModel | 32 | 8 | true | 2400 | 27700 | 6.65 MB | 201.0 s |

### Shallow water model, resolution

| Model | truncation | L | Rings | Δt | SYPD | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 32 | 1 | 48 | 2400 | 115159 | 2.00 MB | 2.1 s |
| ShallowWaterModel | 43 | 1 | 64 | 1786 | 81850 | 4.28 MB | 9.9 s |
| ShallowWaterModel | 64 | 1 | 96 | 1200 | 53711 | 12.85 MB | 13.8 s |
| ShallowWaterModel | 86 | 1 | 128 | 893 | 30469 | 29.04 MB | 16.8 s |
| ShallowWaterModel | 128 | 1 | 192 | 600 | 12661 | 91.85 MB | 10.8 s |
| ShallowWaterModel | 171 | 1 | 256 | 449 | 6477 | 211.92 MB | 11.5 s |
| ShallowWaterModel | 256 | 1 | 384 | 300 | 2042 | 693.70 MB | 13.2 s |

### Primitive wet model, resolution

| Model | truncation | L | Rings | Transform | Δt | SYPD | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | 48 | default | 2400 | 28529 | 6.65 MB | 15.9 s |
| PrimitiveWetModel | 43 | 8 | 64 | default | 1786 | 15354 | 12.05 MB | 42.9 s |
| PrimitiveWetModel | 64 | 8 | 96 | default | 1200 | 9160 | 29.16 MB | 44.5 s |
| PrimitiveWetModel | 86 | 8 | 128 | default | 893 | 5394 | 57.10 MB | 48.7 s |
| PrimitiveWetModel | 128 | 8 | 192 | default | 600 | 2260 | 152.52 MB | 41.0 s |
| PrimitiveWetModel | 171 | 8 | 256 | default | 449 | 1080 | 317.84 MB | 42.4 s |
| PrimitiveWetModel | 256 | 8 | 384 | default | 300 | 362 | 927.31 MB | 43.9 s |
| PrimitiveWetModel | 86 | 16 | 128 | default | 893 | 3033 | 80.75 MB | 182.6 s |
| PrimitiveWetModel | 128 | 16 | 192 | default | 600 | 1359 | 203.69 MB | 50.3 s |
| PrimitiveWetModel | 171 | 16 | 256 | default | 449 | 610 | 407.23 MB | 49.7 s |
| PrimitiveWetModel | 256 | 16 | 384 | default | 300 | 194 | 1.12 GB | 52.8 s |
| PrimitiveWetModel | 86 | 24 | 128 | default | 893 | 2083 | 104.45 MB | 223.7 s |
| PrimitiveWetModel | 128 | 24 | 192 | default | 600 | 901 | 254.93 MB | 87.6 s |
| PrimitiveWetModel | 171 | 24 | 256 | default | 449 | 396 | 496.72 MB | 85.2 s |
| PrimitiveWetModel | 256 | 24 | 384 | default | 300 | 134 | 1.32 GB | 90.1 s |
| PrimitiveWetModel | 32 | 8 | 48 | matrix | 2400 | 19398 | 33.85 MB | 173.6 s |
| PrimitiveWetModel | 43 | 8 | 64 | matrix | 1786 | 13568 | 92.21 MB | 19.7 s |
| PrimitiveWetModel | 64 | 8 | 96 | matrix | 1200 | 6487 | 394.69 MB | 24.2 s |
| PrimitiveWetModel | 86 | 8 | 128 | matrix | 893 | 2012 | 1.18 GB | 26.2 s |
| PrimitiveWetModel | 128 | 8 | 192 | matrix | 600 | 664 | 5.48 GB | 21.1 s |
| PrimitiveWetModel | 171 | 8 | 256 | matrix | 449 | 165 | 16.82 GB | 20.2 s |
| PrimitiveWetModel | 86 | 16 | 128 | matrix | 893 | 1688 | 1.20 GB | 168.2 s |
| PrimitiveWetModel | 128 | 16 | 192 | matrix | 600 | 550 | 5.54 GB | 34.5 s |
| PrimitiveWetModel | 171 | 16 | 256 | matrix | 449 | 101 | 16.92 GB | 34.4 s |
| PrimitiveWetModel | 86 | 24 | 128 | matrix | 893 | 1248 | 1.23 GB | 206.5 s |
| PrimitiveWetModel | 128 | 24 | 192 | matrix | 600 | 473 | 5.59 GB | 71.8 s |
| PrimitiveWetModel | 171 | 24 | 256 | matrix | 449 | 81 | 17.02 GB | 69.8 s |

### Primitive Equation, Float32 vs Float64

| Model | NF | truncation | L | Δt | SYPD | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | Float32 | 32 | 8 | 2400 | 27930 | 6.65 MB | 17.4 s |
| PrimitiveWetModel | Float64 | 32 | 8 | 2400 | 15317 | 13.14 MB | 239.2 s |

### Grids

| Model | truncation | L | Grid | Rings | Δt | SYPD | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 64 | 8 | FullGaussianGrid | 96 | 1200 | 8196 | 38.17 MB | 198.0 s |
| PrimitiveWetModel | 64 | 8 | FullClenshawGrid | 127 | 1200 | 7340 | 64.20 MB | 198.1 s |
| PrimitiveWetModel | 64 | 8 | OctahedralGaussianGrid | 96 | 1200 | 9192 | 29.16 MB | 20.4 s |
| PrimitiveWetModel | 64 | 8 | OctahedralClenshawGrid | 127 | 1200 | 7969 | 46.81 MB | 197.6 s |
| PrimitiveWetModel | 64 | 8 | HEALPixGrid | 127 | 1200 | 9363 | 29.36 MB | 197.2 s |
| PrimitiveWetModel | 64 | 8 | OctaHEALPixGrid | 127 | 1200 | 7968 | 43.11 MB | 198.6 s |

### Number of vertical layers

| Model | truncation | L | Δt | SYPD | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 4 | 2400 | 43490 | 4.70 MB | 153.8 s |
| PrimitiveWetModel | 32 | 8 | 2400 | 28368 | 6.65 MB | 16.9 s |
| PrimitiveWetModel | 32 | 12 | 2400 | 14667 | 8.59 MB | 167.6 s |
| PrimitiveWetModel | 32 | 16 | 2400 | 10686 | 10.55 MB | 52.7 s |

### PrimitiveDryModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveDryModel | 32 | 8 | true | true | 2400 | 35907 | 5.60 MB | 15.8 s |
| PrimitiveDryModel | 32 | 8 | true | false | 2400 | 36877 | 5.60 MB | 10.5 s |
| PrimitiveDryModel | 32 | 8 | false | true | 2400 | 52151 | 5.60 MB | 11.7 s |

### PrimitiveWetModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | true | true | 2400 | 28702 | 6.65 MB | 17.6 s |
| PrimitiveWetModel | 32 | 8 | true | false | 2400 | 44588 | 6.65 MB | 11.7 s |
| PrimitiveWetModel | 32 | 8 | false | true | 2400 | 39791 | 6.65 MB | 15.1 s |

### Individual dynamics functions


#### PrimitiveWetModel | Float32 | T31 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| --- | --- | --- | --- |
| pressure_gradient_flux! | 89.260 μs| 3.70 KiB| 25 |
| linear_virtual_temperature! | 102.341 μs| 4.58 KiB| 37 |
| geopotential! | 88.850 μs| 6.00 KiB| 36 |
| vertical_integration! | 93.960 μs| 16.34 KiB| 99 |
| surface_pressure_tendency! | 96.209 μs| 5.67 KiB| 68 |
| vertical_velocity! | 98.779 μs| 4.09 KiB| 32 |
| linear_pressure_gradient! | 73.200 μs| 3.77 KiB| 25 |
| vertical_advection! | 101.799 μs| 4.20 KiB| 53 |
| vordiv_tendencies! | 120.599 μs| 8.98 KiB| 102 |
| temperature_tendency! | 112.849 μs| 6.30 KiB| 80 |
| humidity_tendency! | 108.949 μs| 6.30 KiB| 80 |
| bernoulli_potential! | 105.619 μs| 7.02 KiB| 79 |

