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

Reactant runs use the `MatrixSpectralTransform` (configurations above T150 are too large for it and skipped), no progress meter and no output.

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
| 32 | 8 | LT+FFT | 1308 | 862 | 5439 | 3763 | — |
| 32 | 8 | MT | 772 | 887 | 5343 | 3696 | 16563 |
| 43 | 8 | LT+FFT | 559 | 368 | 4199 | 2149 | — |
| 43 | 8 | MT | 253 | 308 | 3986 | 2830 | 11888 |
| 64 | 8 | LT+FFT | 152 | 107 | 2822 | 828 | — |
| 64 | 8 | MT | 52 | 63 | 2804 | 826 | 5718 |
| 86 | 8 | LT+FFT | 58 | 41 | 1087 | 246 | — |
| 86 | 8 | MT | 14 | 19 | 1085 | 245 | 1324 |
| 86 | 16 | LT+FFT | 29 | 21 | 749 | 202 | — |
| 86 | 16 | MT | 8.8 | 11 | 746 | 150 | 1105 |
| 86 | 24 | LT+FFT | 20 | 14 | 639 | 196 | — |
| 86 | 24 | MT | 6.1 | 8.6 | 640 | 108 | 955 |
| 128 | 8 | LT+FFT | 14 | 11 | 271 | 86 | — |
| 128 | 8 | MT | 1.2 | 3.3 | 199 | 39 | 466 |
| 128 | 16 | LT+FFT | 7.6 | 5.7 | 243 | 80 | — |
| 128 | 16 | MT | 1.2 | 2.1 | 138 | 23 | 410 |
| 128 | 24 | LT+FFT | 5.5 | 3.6 | 227 | 76 | — |
| 128 | 24 | MT | 1.0 | 1.6 | 118 | 17 | 372 |
| 171 | 8 | LT+FFT | 5.7 | 3.9 | 141 | 44 | — |
| 171 | 16 | LT+FFT | 3.1 | 2.1 | 126 | 41 | — |
| 171 | 24 | LT+FFT | 1.9 | 1.4 | 113 | 38 | — |
| 256 | 8 | LT+FFT | 1.4 | 1.0 | 54 | 16 | — |
| 256 | 16 | LT+FFT | 0.8 | 0.5 | 45 | 13 | — |
| 256 | 24 | LT+FFT | 0.5 | 0.3 | 39 | 11 | — |

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

Created for SpeedyWeather.jl v0.23.0 on Tue, 06 Oct 2026 14:40:04.

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
| BarotropicModel | 32 | 1 | false | 1800 | 24493 | 725.64 KB |
| ShallowWaterModel | 32 | 1 | false | 2400 | 18290 | 908.00 KB |
| PrimitiveDryModel | 32 | 8 | true | 2400 | 1338 | 4.78 MB |
| PrimitiveWetModel | 32 | 8 | true | 2400 | 866 | 5.82 MB |

### Shallow water model, resolution

| Model | truncation | L | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 32 | 1 | 48 | 2400 | 18696 | 908.00 KB |
| ShallowWaterModel | 43 | 1 | 64 | 1800 | 8417 | 1.59 MB |
| ShallowWaterModel | 64 | 1 | 96 | 1200 | 2171 | 3.57 MB |
| ShallowWaterModel | 86 | 1 | 128 | 900 | 904 | 6.49 MB |
| ShallowWaterModel | 128 | 1 | 192 | 600 | 232 | 15.36 MB |
| ShallowWaterModel | 171 | 1 | 256 | 450 | 85 | 29.00 MB |
| ShallowWaterModel | 256 | 1 | 384 | 300 | 18 | 73.09 MB |

### Primitive wet model, resolution

| Model | truncation | L | Rings | Transform | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | 48 | default | 2400 | 862 | 5.82 MB |
| PrimitiveWetModel | 43 | 8 | 64 | default | 1800 | 368 | 9.82 MB |
| PrimitiveWetModel | 64 | 8 | 96 | default | 1200 | 107 | 20.88 MB |
| PrimitiveWetModel | 86 | 8 | 128 | default | 900 | 41 | 36.31 MB |
| PrimitiveWetModel | 128 | 8 | 192 | default | 600 | 11 | 79.95 MB |
| PrimitiveWetModel | 171 | 8 | 256 | default | 450 | 3.9 | 141.86 MB |
| PrimitiveWetModel | 256 | 8 | 384 | default | 300 | 1.0 | 322.22 MB |
| PrimitiveWetModel | 86 | 16 | 128 | default | 900 | 21 | 62.14 MB |
| PrimitiveWetModel | 128 | 16 | 192 | default | 600 | 5.7 | 135.95 MB |
| PrimitiveWetModel | 171 | 16 | 256 | default | 450 | 2.1 | 239.79 MB |
| PrimitiveWetModel | 256 | 16 | 384 | default | 300 | 0.5 | 538.55 MB |
| PrimitiveWetModel | 86 | 24 | 128 | default | 900 | 14 | 88.01 MB |
| PrimitiveWetModel | 128 | 24 | 192 | default | 600 | 3.6 | 192.01 MB |
| PrimitiveWetModel | 171 | 24 | 256 | default | 450 | 1.4 | 337.81 MB |
| PrimitiveWetModel | 256 | 24 | 384 | default | 300 | 0.3 | 755.00 MB |
| PrimitiveWetModel | 32 | 8 | 48 | matrix | 2400 | 887 | 34.26 MB |
| PrimitiveWetModel | 43 | 8 | 64 | matrix | 1800 | 308 | 92.95 MB |
| PrimitiveWetModel | 64 | 8 | 96 | matrix | 1200 | 63 | 396.37 MB |
| PrimitiveWetModel | 86 | 8 | 128 | matrix | 900 | 19 | 1.18 GB |
| PrimitiveWetModel | 128 | 8 | 192 | matrix | 600 | 3.3 | 5.49 GB |
| PrimitiveWetModel | 86 | 16 | 128 | matrix | 900 | 11 | 1.21 GB |
| PrimitiveWetModel | 128 | 16 | 192 | matrix | 600 | 2.1 | 5.55 GB |
| PrimitiveWetModel | 86 | 24 | 128 | matrix | 900 | 8.6 | 1.23 GB |
| PrimitiveWetModel | 128 | 24 | 192 | matrix | 600 | 1.6 | 5.60 GB |

### Primitive Equation, Float32 vs Float64

| Model | NF | truncation | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | Float32 | 32 | 8 | 2400 | 870 | 5.82 MB |
| PrimitiveWetModel | Float64 | 32 | 8 | 2400 | 695 | 11.04 MB |

### Grids

| Model | truncation | L | Grid | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 64 | 8 | FullGaussianGrid | 96 | 1200 | 71 | 30.59 MB |
| PrimitiveWetModel | 64 | 8 | FullClenshawGrid | 127 | 1200 | 42 | 51.09 MB |
| PrimitiveWetModel | 64 | 8 | OctahedralGaussianGrid | 96 | 1200 | 108 | 20.88 MB |
| PrimitiveWetModel | 64 | 8 | OctahedralClenshawGrid | 127 | 1200 | 58 | 32.77 MB |
| PrimitiveWetModel | 64 | 8 | HEALPixGrid | 127 | 1200 | 87 | 24.19 MB |
| PrimitiveWetModel | 64 | 8 | OctaHEALPixGrid | 127 | 1200 | 54 | 30.06 MB |

### Number of vertical layers

| Model | truncation | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 4 | 2400 | 1310 | 3.72 MB |
| PrimitiveWetModel | 32 | 8 | 2400 | 855 | 5.82 MB |
| PrimitiveWetModel | 32 | 12 | 2400 | 610 | 7.93 MB |
| PrimitiveWetModel | 32 | 16 | 2400 | 469 | 10.04 MB |

### PrimitiveDryModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveDryModel | 32 | 8 | true | true | 2400 | 1353 | 4.78 MB |
| PrimitiveDryModel | 32 | 8 | true | false | 2400 | 2276 | 4.78 MB |
| PrimitiveDryModel | 32 | 8 | false | true | 2400 | 1499 | 4.78 MB |

### PrimitiveWetModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | true | true | 2400 | 861 | 5.82 MB |
| PrimitiveWetModel | 32 | 8 | true | false | 2400 | 1733 | 5.82 MB |
| PrimitiveWetModel | 32 | 8 | false | true | 2400 | 973 | 5.82 MB |

### Individual dynamics functions


#### PrimitiveWetModel | Float32 | T31 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| --- | --- | --- | --- |
| pressure_gradient_flux! | 69.340 μs| 31.98 KiB| 200 |
| linear_virtual_temperature! | 3.401 μs| 0 bytes| 0 |
| geopotential! | 9.990 μs| 384 bytes| 6 |
| vertical_integration! | 13.590 μs| 0 bytes| 0 |
| surface_pressure_tendency! | 18.490 μs| 15.66 KiB| 96 |
| vertical_velocity! | 59.720 μs| 0 bytes| 0 |
| linear_pressure_gradient! | 3.099 μs| 0 bytes| 0 |
| vertical_advection! | 160.399 μs| 2.44 KiB| 32 |
| vordiv_tendencies! | 403.267 μs| 218.48 KiB| 284 |
| temperature_tendency! | 515.327 μs| 324.75 KiB| 401 |
| humidity_tendency! | 497.697 μs| 324.08 KiB| 396 |
| bernoulli_potential! | 162.949 μs| 107.61 KiB| 129 |

## Architecture: `gpu-nvidia`

Created for SpeedyWeather.jl v0.24.0-DEV on Wed, 14 Oct 2026 17:44:31.

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
  0: NVIDIA H100 80GB HBM3 (sm_90, 77.874 GiB / 79.647 GiB available)
     compiles to sm_90a / PTX 9.3 (LLVM: sm_90a / PTX 9.0)
```


### Models, default setups

| Model | truncation | L | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| BarotropicModel | 32 | 1 | false | 1800 | 15501 | 207.93 KB |
| ShallowWaterModel | 32 | 1 | false | 2400 | 11584 | 211.51 KB |
| PrimitiveDryModel | 32 | 8 | true | 2400 | 7050 | 355.81 KB |
| PrimitiveWetModel | 32 | 8 | true | 2400 | 5594 | 362.10 KB |

### Shallow water model, resolution

| Model | truncation | L | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 32 | 1 | 48 | 2400 | 11181 | 211.51 KB |
| ShallowWaterModel | 43 | 1 | 64 | 1800 | 8585 | 355.38 KB |
| ShallowWaterModel | 64 | 1 | 96 | 1200 | 4636 | 766.00 KB |
| ShallowWaterModel | 86 | 1 | 128 | 900 | 1649 | 1.34 MB |
| ShallowWaterModel | 128 | 1 | 192 | 600 | 296 | 3.27 MB |
| ShallowWaterModel | 171 | 1 | 256 | 450 | 129 | 5.66 MB |
| ShallowWaterModel | 256 | 1 | 384 | 300 | 52 | 12.40 MB |

### Primitive wet model, resolution

| Model | truncation | L | Rings | Transform | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | 48 | default | 2400 | 5439 | 362.10 KB |
| PrimitiveWetModel | 43 | 8 | 64 | default | 1800 | 4199 | 606.33 KB |
| PrimitiveWetModel | 64 | 8 | 96 | default | 1200 | 2822 | 1.30 MB |
| PrimitiveWetModel | 86 | 8 | 128 | default | 900 | 1087 | 2.28 MB |
| PrimitiveWetModel | 128 | 8 | 192 | default | 600 | 271 | 5.45 MB |
| PrimitiveWetModel | 171 | 8 | 256 | default | 450 | 141 | 9.47 MB |
| PrimitiveWetModel | 256 | 8 | 384 | default | 300 | 54 | 20.87 MB |
| PrimitiveWetModel | 86 | 16 | 128 | default | 900 | 749 | 3.33 MB |
| PrimitiveWetModel | 128 | 16 | 192 | default | 600 | 243 | 7.81 MB |
| PrimitiveWetModel | 171 | 16 | 256 | default | 450 | 126 | 13.67 MB |
| PrimitiveWetModel | 256 | 16 | 384 | default | 300 | 45 | 30.31 MB |
| PrimitiveWetModel | 86 | 24 | 128 | default | 900 | 639 | 4.38 MB |
| PrimitiveWetModel | 128 | 24 | 192 | default | 600 | 227 | 10.17 MB |
| PrimitiveWetModel | 171 | 24 | 256 | default | 450 | 113 | 17.86 MB |
| PrimitiveWetModel | 256 | 24 | 384 | default | 300 | 39 | 39.75 MB |
| PrimitiveWetModel | 32 | 8 | 48 | matrix | 2400 | 5343 | 362.10 KB |
| PrimitiveWetModel | 43 | 8 | 64 | matrix | 1800 | 3986 | 606.33 KB |
| PrimitiveWetModel | 64 | 8 | 96 | matrix | 1200 | 2804 | 1.30 MB |
| PrimitiveWetModel | 86 | 8 | 128 | matrix | 900 | 1085 | 2.28 MB |
| PrimitiveWetModel | 128 | 8 | 192 | matrix | 600 | 199 | 5.07 MB |
| PrimitiveWetModel | 86 | 16 | 128 | matrix | 900 | 746 | 3.33 MB |
| PrimitiveWetModel | 128 | 16 | 192 | matrix | 600 | 138 | 7.43 MB |
| PrimitiveWetModel | 86 | 24 | 128 | matrix | 900 | 640 | 4.38 MB |
| PrimitiveWetModel | 128 | 24 | 192 | matrix | 600 | 118 | 9.79 MB |

### Primitive Equation, Float32 vs Float64

| Model | NF | truncation | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | Float32 | 32 | 8 | 2400 | 5188 | 362.10 KB |
| PrimitiveWetModel | Float64 | 32 | 8 | 2400 | 4565 | 362.70 KB |

### Grids

| Model | truncation | L | Grid | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 64 | 8 | FullGaussianGrid | 96 | 1200 | 2195 | 1.30 MB |
| PrimitiveWetModel | 64 | 8 | FullClenshawGrid | 127 | 1200 | 1514 | 2.26 MB |
| PrimitiveWetModel | 64 | 8 | OctahedralGaussianGrid | 96 | 1200 | 2960 | 1.30 MB |
| PrimitiveWetModel | 64 | 8 | OctahedralClenshawGrid | 127 | 1200 | 2211 | 2.26 MB |
| PrimitiveWetModel | 64 | 8 | HEALPixGrid | 127 | 1200 | 2600 | 2.26 MB |
| PrimitiveWetModel | 64 | 8 | OctaHEALPixGrid | 127 | 1200 | 2331 | 2.26 MB |

### Number of vertical layers

| Model | truncation | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 4 | 2400 | 5626 | 288.37 KB |
| PrimitiveWetModel | 32 | 8 | 2400 | 5778 | 362.10 KB |
| PrimitiveWetModel | 32 | 12 | 2400 | 5507 | 435.83 KB |
| PrimitiveWetModel | 32 | 16 | 2400 | 5749 | 509.56 KB |

### PrimitiveDryModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveDryModel | 32 | 8 | true | true | 2400 | 7597 | 355.81 KB |
| PrimitiveDryModel | 32 | 8 | true | false | 2400 | 8611 | 355.81 KB |
| PrimitiveDryModel | 32 | 8 | false | true | 2400 | 13468 | 355.81 KB |

### PrimitiveWetModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | true | true | 2400 | 5669 | 362.10 KB |
| PrimitiveWetModel | 32 | 8 | true | false | 2400 | 7744 | 362.10 KB |
| PrimitiveWetModel | 32 | 8 | false | true | 2400 | 9399 | 362.10 KB |

### Individual dynamics functions


#### PrimitiveWetModel | Float32 | T31 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| --- | --- | --- | --- |
| pressure_gradient_flux! | 71.928 μs| 13.95 KiB| 307 |
| linear_virtual_temperature! | 13.370 μs| 2.77 KiB| 52 |
| geopotential! | 19.740 μs| 3.73 KiB| 101 |
| vertical_integration! | 28.859 μs| 7.20 KiB| 143 |
| surface_pressure_tendency! | N/A| N/A| N/A |
| vertical_velocity! | 25.849 μs| 8.45 KiB| 191 |
| linear_pressure_gradient! | 13.109 μs| 2.34 KiB| 50 |
| vertical_advection! | 25.590 μs| 9.09 KiB| 142 |
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

Created for SpeedyWeather.jl v0.24.0-DEV on Fri, 16 Oct 2026 10:49:04.

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
  0: NVIDIA H100 80GB HBM3 (sm_90, 19.151 GiB / 79.647 GiB available)
     compiles to sm_90a / PTX 9.3 (LLVM: sm_90a / PTX 9.0)
```

Reactant backend: `gpu` (via `Reactant.set_default_backend("gpu")`), Reactant v0.2.293, Reactant_jll v0.0.421+0.


### Models, default setups

| Model | truncation | L | Physics | Δt | SYPD | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- |
| BarotropicModel | 32 | 1 | false | 1800 | 119071 | 28.73 MB | 29.6 s |
| ShallowWaterModel | 32 | 1 | false | 2400 | 99985 | 28.92 MB | 37.7 s |
| PrimitiveDryModel | 32 | 8 | true | 2400 | 28865 | 32.78 MB | 179.2 s |
| PrimitiveWetModel | 32 | 8 | true | 2400 | 18000 | 33.82 MB | 196.4 s |

### Shallow water model, resolution

| Model | truncation | L | Rings | Δt | SYPD | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 32 | 1 | 48 | 2400 | 89477 | 28.92 MB | 1.9 s |
| ShallowWaterModel | 43 | 1 | 64 | 1786 | 40055 | 83.94 MB | 9.8 s |
| ShallowWaterModel | 64 | 1 | 96 | 1200 | 16759 | 377.31 MB | 13.8 s |
| ShallowWaterModel | 86 | 1 | 128 | 893 | 2804 | 1.15 GB | 17.1 s |
| ShallowWaterModel | 128 | 1 | 192 | 600 | 893 | 5.42 GB | 10.5 s |

### Primitive wet model, resolution

| Model | truncation | L | Rings | Δt | SYPD | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | 48 | 2400 | 16563 | 33.82 MB | 14.5 s |
| PrimitiveWetModel | 43 | 8 | 64 | 1786 | 11888 | 92.17 MB | 43.8 s |
| PrimitiveWetModel | 64 | 8 | 96 | 1200 | 5718 | 394.60 MB | 44.2 s |
| PrimitiveWetModel | 86 | 8 | 128 | 893 | 1324 | 1.18 GB | 48.6 s |
| PrimitiveWetModel | 128 | 8 | 192 | 600 | 466 | 5.48 GB | 40.7 s |
| PrimitiveWetModel | 86 | 16 | 128 | 893 | 1105 | 1.20 GB | 183.5 s |
| PrimitiveWetModel | 128 | 16 | 192 | 600 | 410 | 5.54 GB | 48.9 s |
| PrimitiveWetModel | 86 | 24 | 128 | 893 | 955 | 1.23 GB | 222.5 s |
| PrimitiveWetModel | 128 | 24 | 192 | 600 | 372 | 5.59 GB | 86.5 s |

### Primitive Equation, Float32 vs Float64

| Model | NF | truncation | L | Δt | SYPD | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | Float32 | 32 | 8 | 2400 | 16477 | 33.82 MB | 15.0 s |
| PrimitiveWetModel | Float64 | 32 | 8 | 2400 | 7988 | 67.52 MB | 242.0 s |

### Grids

| Model | truncation | L | Grid | Rings | Δt | SYPD | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 64 | 8 | FullGaussianGrid | 96 | 1200 | 4514 | 661.18 MB | 198.1 s |
| PrimitiveWetModel | 64 | 8 | FullClenshawGrid | 127 | 1200 | 3205 | 1.16 GB | 198.8 s |
| PrimitiveWetModel | 64 | 8 | OctahedralGaussianGrid | 96 | 1200 | 5702 | 394.60 MB | 19.1 s |
| PrimitiveWetModel | 64 | 8 | OctahedralClenshawGrid | 127 | 1200 | 3597 | 660.61 MB | 197.1 s |
| PrimitiveWetModel | 64 | 8 | HEALPixGrid | 127 | 1200 | 5701 | 442.45 MB | 197.2 s |
| PrimitiveWetModel | 64 | 8 | OctaHEALPixGrid | 127 | 1200 | 4878 | 588.27 MB | 201.1 s |

### Number of vertical layers

| Model | truncation | L | Δt | SYPD | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 4 | 2400 | 30540 | 31.72 MB | 156.1 s |
| PrimitiveWetModel | 32 | 8 | 2400 | 16444 | 33.82 MB | 14.7 s |
| PrimitiveWetModel | 32 | 12 | 2400 | 17013 | 35.93 MB | 166.3 s |
| PrimitiveWetModel | 32 | 16 | 2400 | 13105 | 38.05 MB | 53.5 s |

### PrimitiveDryModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveDryModel | 32 | 8 | true | true | 2400 | 24610 | 32.78 MB | 10.0 s |
| PrimitiveDryModel | 32 | 8 | true | false | 2400 | 28786 | 32.78 MB | 7.5 s |
| PrimitiveDryModel | 32 | 8 | false | true | 2400 | 47684 | 32.78 MB | 10.4 s |

### PrimitiveWetModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | true | true | 2400 | 16482 | 33.82 MB | 14.3 s |
| PrimitiveWetModel | 32 | 8 | true | false | 2400 | 24543 | 33.82 MB | 8.0 s |
| PrimitiveWetModel | 32 | 8 | false | true | 2400 | 32011 | 33.82 MB | 13.9 s |

### Individual dynamics functions


#### PrimitiveWetModel | Float32 | T31 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| --- | --- | --- | --- |
| pressure_gradient_flux! | 105.590 μs| 15.19 KiB| 80 |
| linear_virtual_temperature! | 96.940 μs| 4.55 KiB| 38 |
| geopotential! | 99.170 μs| 6.09 KiB| 38 |
| vertical_integration! | 99.380 μs| 16.47 KiB| 104 |
| surface_pressure_tendency! | 129.700 μs| 14.56 KiB| 83 |
| vertical_velocity! | 99.660 μs| 4.19 KiB| 34 |
| linear_pressure_gradient! | 86.250 μs| 3.81 KiB| 26 |
| vertical_advection! | 99.410 μs| 3.92 KiB| 52 |
| vordiv_tendencies! | 136.050 μs| 18.69 KiB| 122 |
| temperature_tendency! | 129.820 μs| 15.83 KiB| 97 |
| humidity_tendency! | 123.390 μs| 15.83 KiB| 97 |
| bernoulli_potential! | 107.970 μs| 15.95 KiB| 95 |

