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
| 32 | 8 | LT+FFT | 1308 | 870 | 5906 | 3763 | 8310 |
| 32 | 8 | MT | 772 | 882 | 5779 | 3696 | 19072 |
| 43 | 8 | LT+FFT | 559 | 371 | 4435 | 2149 | 4379 |
| 43 | 8 | MT | 253 | 337 | 4399 | 2830 | 13266 |
| 64 | 8 | LT+FFT | 152 | 108 | 2978 | 828 | 1972 |
| 64 | 8 | MT | 52 | 72 | 3057 | 826 | 6776 |
| 86 | 8 | LT+FFT | 58 | 43 | 1239 | 246 | 1001 |
| 86 | 8 | MT | 14 | 23 | 1245 | 245 | 2009 |
| 86 | 16 | LT+FFT | 29 | 21 | 777 | 202 | 825 |
| 86 | 16 | MT | 8.8 | 13 | 774 | 150 | 1735 |
| 86 | 24 | LT+FFT | 20 | 28 | 1350 | 196 | 728 |
| 86 | 24 | MT | 6.1 | 19 | 1351 | 108 | 1249 |
| 128 | 8 | LT+FFT | 14 | 11 | 364 | 86 | 355 |
| 128 | 8 | MT | 1.2 | — | 237 | 39 | 658 |
| 128 | 16 | LT+FFT | 7.6 | 11 | 608 | 80 | 313 |
| 128 | 16 | MT | 1.2 | — | 288 | 23 | 554 |
| 128 | 24 | LT+FFT | 5.5 | 7.3 | 561 | 76 | 269 |
| 128 | 24 | MT | 1.0 | — | 254 | 17 | 477 |
| 171 | 8 | LT+FFT | 5.7 | 4.1 | 187 | 44 | 186 |
| 171 | 8 | MT | — | — | 52 | — | 164 |
| 171 | 16 | LT+FFT | 3.1 | 4.3 | 312 | 41 | 164 |
| 171 | 16 | MT | — | — | 72 | — | 101 |
| 171 | 24 | LT+FFT | 1.9 | 2.9 | 272 | 38 | 140 |
| 171 | 24 | MT | — | — | 62 | — | 80 |
| 256 | 8 | LT+FFT | 1.4 | 1.0 | 69 | 16 | 83 |
| 256 | 16 | LT+FFT | 0.8 | 1.0 | 108 | 13 | 68 |
| 256 | 24 | LT+FFT | 0.5 | 0.7 | 91 | 11 | 57 |

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

Created for SpeedyWeather.jl v0.24.0-DEV on Mon, 19 Oct 2026 17:01:22.

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
| BarotropicModel | 32 | 1 | false | 1800 | 16695 | 207.93 KB |
| ShallowWaterModel | 32 | 1 | false | 2400 | 11303 | 211.51 KB |
| PrimitiveDryModel | 32 | 8 | true | 2400 | 7845 | 355.38 KB |
| PrimitiveWetModel | 32 | 8 | true | 2400 | 6137 | 361.68 KB |

### Shallow water model, resolution

| Model | truncation | L | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 32 | 1 | 48 | 2400 | 11560 | 211.51 KB |
| ShallowWaterModel | 43 | 1 | 64 | 1800 | 8575 | 355.38 KB |
| ShallowWaterModel | 64 | 1 | 96 | 1200 | 4425 | 766.00 KB |
| ShallowWaterModel | 86 | 1 | 128 | 900 | 1640 | 1.34 MB |
| ShallowWaterModel | 128 | 1 | 192 | 600 | 292 | 3.31 MB |
| ShallowWaterModel | 171 | 1 | 256 | 450 | 130 | 5.72 MB |
| ShallowWaterModel | 256 | 1 | 384 | 300 | 52 | 12.50 MB |

### Primitive wet model, resolution

| Model | truncation | L | Rings | Transform | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | 48 | default | 2400 | 5906 | 361.68 KB |
| PrimitiveWetModel | 43 | 8 | 64 | default | 1800 | 4435 | 605.90 KB |
| PrimitiveWetModel | 64 | 8 | 96 | default | 1200 | 2978 | 1.30 MB |
| PrimitiveWetModel | 86 | 8 | 128 | default | 900 | 1239 | 2.28 MB |
| PrimitiveWetModel | 128 | 8 | 192 | default | 600 | 364 | 5.45 MB |
| PrimitiveWetModel | 171 | 8 | 256 | default | 450 | 187 | 9.47 MB |
| PrimitiveWetModel | 256 | 8 | 384 | default | 300 | 69 | 20.87 MB |
| PrimitiveWetModel | 86 | 16 | 128 | default | 900 | 777 | 3.33 MB |
| PrimitiveWetModel | 128 | 16 | 192 | default | 600 | 608 | 7.81 MB |
| PrimitiveWetModel | 171 | 16 | 256 | default | 450 | 312 | 13.67 MB |
| PrimitiveWetModel | 256 | 16 | 384 | default | 300 | 108 | 30.31 MB |
| PrimitiveWetModel | 86 | 24 | 128 | default | 900 | 1350 | 4.38 MB |
| PrimitiveWetModel | 128 | 24 | 192 | default | 600 | 561 | 10.17 MB |
| PrimitiveWetModel | 171 | 24 | 256 | default | 450 | 272 | 17.86 MB |
| PrimitiveWetModel | 256 | 24 | 384 | default | 300 | 91 | 39.75 MB |
| PrimitiveWetModel | 32 | 8 | 48 | matrix | 2400 | 5779 | 361.68 KB |
| PrimitiveWetModel | 43 | 8 | 64 | matrix | 1800 | 4399 | 605.90 KB |
| PrimitiveWetModel | 64 | 8 | 96 | matrix | 1200 | 3057 | 1.30 MB |
| PrimitiveWetModel | 86 | 8 | 128 | matrix | 900 | 1245 | 2.28 MB |
| PrimitiveWetModel | 128 | 8 | 192 | matrix | 600 | 237 | 5.07 MB |
| PrimitiveWetModel | 171 | 8 | 256 | matrix | 450 | 52 | 8.97 MB |
| PrimitiveWetModel | 86 | 16 | 128 | matrix | 900 | 774 | 3.33 MB |
| PrimitiveWetModel | 128 | 16 | 192 | matrix | 600 | 288 | 7.43 MB |
| PrimitiveWetModel | 171 | 16 | 256 | matrix | 450 | 72 | 13.16 MB |
| PrimitiveWetModel | 86 | 24 | 128 | matrix | 900 | 1351 | 4.38 MB |
| PrimitiveWetModel | 128 | 24 | 192 | matrix | 600 | 254 | 9.79 MB |
| PrimitiveWetModel | 171 | 24 | 256 | matrix | 450 | 62 | 17.36 MB |

### Primitive Equation, Float32 vs Float64

| Model | NF | truncation | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | Float32 | 32 | 8 | 2400 | 5838 | 361.68 KB |
| PrimitiveWetModel | Float64 | 32 | 8 | 2400 | 4917 | 362.28 KB |

### Grids

| Model | truncation | L | Grid | Rings | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 64 | 8 | FullGaussianGrid | 96 | 1200 | 2134 | 1.30 MB |
| PrimitiveWetModel | 64 | 8 | FullClenshawGrid | 127 | 1200 | 1641 | 2.26 MB |
| PrimitiveWetModel | 64 | 8 | OctahedralGaussianGrid | 96 | 1200 | 3131 | 1.30 MB |
| PrimitiveWetModel | 64 | 8 | OctahedralClenshawGrid | 127 | 1200 | 2333 | 2.26 MB |
| PrimitiveWetModel | 64 | 8 | HEALPixGrid | 127 | 1200 | 2813 | 2.26 MB |
| PrimitiveWetModel | 64 | 8 | OctaHEALPixGrid | 127 | 1200 | 2582 | 2.26 MB |

### Number of vertical layers

| Model | truncation | L | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 4 | 2400 | 6249 | 287.95 KB |
| PrimitiveWetModel | 32 | 8 | 2400 | 6242 | 361.68 KB |
| PrimitiveWetModel | 32 | 12 | 2400 | 6140 | 435.41 KB |
| PrimitiveWetModel | 32 | 16 | 2400 | 6457 | 509.13 KB |

### PrimitiveDryModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveDryModel | 32 | 8 | true | true | 2400 | 8143 | 355.38 KB |
| PrimitiveDryModel | 32 | 8 | true | false | 2400 | 9248 | 355.38 KB |
| PrimitiveDryModel | 32 | 8 | false | true | 2400 | 13333 | 355.38 KB |

### PrimitiveWetModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | true | true | 2400 | 6198 | 361.68 KB |
| PrimitiveWetModel | 32 | 8 | true | false | 2400 | 7898 | 361.68 KB |
| PrimitiveWetModel | 32 | 8 | false | true | 2400 | 9015 | 361.68 KB |

### Individual dynamics functions


#### PrimitiveWetModel | Float32 | T31 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| --- | --- | --- | --- |
| pressure_gradient_flux! | 13.810 μs| 4.27 KiB| 60 |
| linear_virtual_temperature! | 12.918 μs| 2.77 KiB| 52 |
| geopotential! | 19.348 μs| 3.73 KiB| 101 |
| vertical_integration! | 27.687 μs| 7.75 KiB| 144 |
| surface_pressure_tendency! | N/A| N/A| N/A |
| vertical_velocity! | 25.181 μs| 9.03 KiB| 191 |
| linear_pressure_gradient! | 12.820 μs| 2.34 KiB| 50 |
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

## Architecture: `reactant-gpu`

Created for SpeedyWeather.jl v0.24.0-DEV on Mon, 19 Oct 2026 20:39:18.

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
  0: NVIDIA H100 80GB HBM3 (sm_90, 19.059 GiB / 79.647 GiB available)
     compiles to sm_90a / PTX 9.3 (LLVM: sm_90a / PTX 9.0)
```

Reactant backend: `gpu` (via `Reactant.set_default_backend("gpu")`), Reactant v0.2.295, Reactant_jll v0.0.424+0.


### Models, default setups

| Model | truncation | L | Physics | Δt | SYPD | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- |
| BarotropicModel | 32 | 1 | false | 1800 | 5230 | 657.20 KB | 32.9 s |
| ShallowWaterModel | 32 | 1 | false | 2400 | 4460 | 839.45 KB | 41.5 s |
| PrimitiveDryModel | 32 | 8 | true | 2400 | 9333 | 5.82 MB | 178.0 s |
| PrimitiveWetModel | 32 | 8 | true | 2400 | 8398 | 6.87 MB | 194.6 s |

### Shallow water model, resolution

| Model | truncation | L | Rings | Δt | SYPD | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 32 | 1 | 48 | 2400 | 4375 | 839.45 KB | 5.2 s |
| ShallowWaterModel | 43 | 1 | 64 | 1786 | 2395 | 1.48 MB | 14.2 s |
| ShallowWaterModel | 64 | 1 | 96 | 1200 | 962 | 3.44 MB | 21.1 s |
| ShallowWaterModel | 86 | 1 | 128 | 893 | 475 | 6.49 MB | 27.8 s |
| ShallowWaterModel | 128 | 1 | 192 | 600 | 190 | 16.30 MB | 34.3 s |
| ShallowWaterModel | 171 | 1 | 256 | 449 | 101 | 32.44 MB | 55.3 s |
| ShallowWaterModel | 256 | 1 | 384 | 300 | 43 | 89.21 MB | 145.7 s |

### Primitive wet model, resolution

| Model | truncation | L | Rings | Transform | Δt | SYPD | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | 48 | default | 2400 | 8310 | 6.87 MB | 19.4 s |
| PrimitiveWetModel | 43 | 8 | 64 | default | 1786 | 4379 | 11.61 MB | 46.4 s |
| PrimitiveWetModel | 64 | 8 | 96 | default | 1200 | 1972 | 24.82 MB | 48.4 s |
| PrimitiveWetModel | 86 | 8 | 128 | default | 893 | 1001 | 43.38 MB | 52.9 s |
| PrimitiveWetModel | 128 | 8 | 192 | default | 600 | 355 | 96.42 MB | 55.1 s |
| PrimitiveWetModel | 171 | 8 | 256 | default | 449 | 186 | 172.55 MB | 62.7 s |
| PrimitiveWetModel | 256 | 8 | 384 | default | 300 | 83 | 398.88 MB | 87.6 s |
| PrimitiveWetModel | 86 | 16 | 128 | default | 893 | 825 | 77.14 MB | 191.1 s |
| PrimitiveWetModel | 128 | 16 | 192 | default | 600 | 313 | 169.81 MB | 64.5 s |
| PrimitiveWetModel | 171 | 16 | 256 | default | 449 | 164 | 301.02 MB | 71.6 s |
| PrimitiveWetModel | 256 | 16 | 384 | default | 300 | 68 | 683.04 MB | 97.9 s |
| PrimitiveWetModel | 86 | 24 | 128 | default | 893 | 728 | 110.95 MB | 230.3 s |
| PrimitiveWetModel | 128 | 24 | 192 | default | 600 | 269 | 243.28 MB | 103.6 s |
| PrimitiveWetModel | 171 | 24 | 256 | default | 449 | 140 | 429.58 MB | 107.1 s |
| PrimitiveWetModel | 256 | 24 | 384 | default | 300 | 57 | 967.32 MB | 130.8 s |
| PrimitiveWetModel | 32 | 8 | 48 | matrix | 2400 | 19072 | 33.85 MB | 169.7 s |
| PrimitiveWetModel | 43 | 8 | 64 | matrix | 1786 | 13266 | 92.21 MB | 18.9 s |
| PrimitiveWetModel | 64 | 8 | 96 | matrix | 1200 | 6776 | 394.69 MB | 23.6 s |
| PrimitiveWetModel | 86 | 8 | 128 | matrix | 893 | 2009 | 1.18 GB | 27.4 s |
| PrimitiveWetModel | 128 | 8 | 192 | matrix | 600 | 658 | 5.48 GB | 21.8 s |
| PrimitiveWetModel | 171 | 8 | 256 | matrix | 449 | 164 | 16.82 GB | 21.3 s |
| PrimitiveWetModel | 86 | 16 | 128 | matrix | 893 | 1735 | 1.20 GB | 169.5 s |
| PrimitiveWetModel | 128 | 16 | 192 | matrix | 600 | 554 | 5.54 GB | 35.5 s |
| PrimitiveWetModel | 171 | 16 | 256 | matrix | 449 | 101 | 16.92 GB | 36.0 s |
| PrimitiveWetModel | 86 | 24 | 128 | matrix | 893 | 1249 | 1.23 GB | 209.7 s |
| PrimitiveWetModel | 128 | 24 | 192 | matrix | 600 | 477 | 5.59 GB | 72.8 s |
| PrimitiveWetModel | 171 | 24 | 256 | matrix | 449 | 80 | 17.02 GB | 72.2 s |

### Primitive Equation, Float32 vs Float64

| Model | NF | truncation | L | Δt | SYPD | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | Float32 | 32 | 8 | 2400 | 8151 | 6.87 MB | 20.2 s |
| PrimitiveWetModel | Float64 | 32 | 8 | 2400 | 6392 | 13.58 MB | 239.2 s |

### Grids

| Model | truncation | L | Grid | Rings | Δt | SYPD | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 64 | 8 | FullGaussianGrid | 96 | 1200 | 8257 | 34.17 MB | 196.0 s |
| PrimitiveWetModel | 64 | 8 | FullClenshawGrid | 127 | 1200 | 6895 | 57.07 MB | 199.0 s |
| PrimitiveWetModel | 64 | 8 | OctahedralGaussianGrid | 96 | 1200 | 1965 | 24.82 MB | 24.9 s |
| PrimitiveWetModel | 64 | 8 | OctahedralClenshawGrid | 127 | 0 | failed: BoundsError | 0 bytes | — |
| PrimitiveWetModel | 64 | 8 | HEALPixGrid | 127 | 1200 | 2805 | 25.83 MB | 201.2 s |
| PrimitiveWetModel | 64 | 8 | OctaHEALPixGrid | 127 | 0 | failed: BoundsError | 0 bytes | — |

### Number of vertical layers

| Model | truncation | L | Δt | SYPD | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 4 | 2400 | 9534 | 4.14 MB | 156.0 s |
| PrimitiveWetModel | 32 | 8 | 2400 | 8309 | 6.87 MB | 19.7 s |
| PrimitiveWetModel | 32 | 12 | 2400 | 6207 | 9.60 MB | 169.8 s |
| PrimitiveWetModel | 32 | 16 | 2400 | 5384 | 12.35 MB | 56.0 s |

### PrimitiveDryModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveDryModel | 32 | 8 | true | true | 2400 | 9001 | 5.82 MB | 17.4 s |
| PrimitiveDryModel | 32 | 8 | true | false | 2400 | 9071 | 5.82 MB | 11.9 s |
| PrimitiveDryModel | 32 | 8 | false | true | 2400 | 5498 | 5.82 MB | 16.7 s |

### PrimitiveWetModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | true | true | 2400 | 8319 | 6.87 MB | 19.5 s |
| PrimitiveWetModel | 32 | 8 | true | false | 2400 | 9415 | 6.87 MB | 13.7 s |
| PrimitiveWetModel | 32 | 8 | false | true | 2400 | 4520 | 6.87 MB | 23.5 s |

### Individual dynamics functions


#### PrimitiveWetModel | Float32 | T31 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| --- | --- | --- | --- |
| pressure_gradient_flux! | 90.602 μs| 3.70 KiB| 25 |
| linear_virtual_temperature! | 105.791 μs| 4.58 KiB| 37 |
| geopotential! | 104.591 μs| 6.00 KiB| 36 |
| vertical_integration! | 105.261 μs| 16.34 KiB| 99 |
| surface_pressure_tendency! | 305.854 μs| 8.67 KiB| 90 |
| vertical_velocity! | 105.521 μs| 4.09 KiB| 32 |
| linear_pressure_gradient! | 86.228 μs| 3.77 KiB| 25 |
| vertical_advection! | 105.218 μs| 4.20 KiB| 53 |
| vordiv_tendencies! | 408.446 μs| 12.06 KiB| 125 |
| temperature_tendency! | 535.461 μs| 9.30 KiB| 102 |
| humidity_tendency! | 560.861 μs| 9.30 KiB| 102 |
| bernoulli_potential! | 266.021 μs| 10.09 KiB| 102 |

