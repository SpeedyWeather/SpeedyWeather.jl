# Benchmarks

Performance benchmarks for SpeedyWeather.jl, collected across multiple architectures. Each architecture's results live in its own section below; the overview table at the top compares the headline PrimitiveWet resolution sweep across all archs that have been benchmarked so far.

All simulations are benchmarked over several seconds (wallclock time) without output. Benchmarking excludes initialization and is started just before the main time loop and finishes right after. By default each configuration is timed up to 5 times from the same initial state, as long as all its timed runs together take at most 30 s, so short runs are timed 5 times and runs longer than 15 s once; `--repeats N` times every configuration exactly N times. Each architecture's section states its settings. SYPD is the median of these runs, and the Spread column shows how much they differ. Results still vary between processes and machines, e.g. by ±30% on a laptop, so the following presents a rough idea of how fast a SpeedyWeather simulation will run, and how much memory it requires.

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
- SYPD: Speed of simulation, simulated years per wallclock day, median of the timed runs.
- Runs: Number of timed runs, see the settings of each architecture.
- Spread: (slowest − fastest) / median wallclock time of the timed runs, — for a single run.
- Memory: Memory footprint of simulation, variables and constants.
- Compile: One-off compilation of the time step (Reactant only), not included in the SYPD. Compiled code is cached within a session, so configurations compiled before (e.g. in an earlier suite) compile much faster.
- Steps, Forward, Gradient (reverse-mode benchmarks only): wallclock time of that many time steps forward, and of the reverse-mode gradient of a scalar loss after them with respect to the initial state (incl. its forward pass), with the Checkpointing scheme of the time loop. Gradient/Forward is their ratio.

### Running the benchmarks

Reproduce the benchmark suite by running, from `SpeedyWeather/benchmark`:

```
julia --project=. manual_benchmarking.jl                # CPU (auto-labelled cpu-arm or cpu-x86)
julia --project=. manual_benchmarking.jl gpu            # CUDA GPU
julia --project=. manual_benchmarking.jl amdgpu         # AMDGPU (HIP graphs forced on)
julia --project=. manual_benchmarking.jl reactant-cpu   # Reactant on CPU
julia --project=. manual_benchmarking.jl reactant-gpu   # Reactant on CUDA GPU
julia --project=. manual_benchmarking.jl reactant-tpu   # Reactant on TPU
```

Options and flags, in any order (`manual_benchmarking.jl --help` lists them): `--length X`, e.g. `manual_benchmarking.jl gpu --length 10`, multiplies the number of time steps of the timed runs of the PrimitiveWet resolution sweep (the overview below) for more robust timings; the other suites are not affected. `--repeats N` times every configuration exactly N times instead of up to 5 times within 30 s. `--suites benchmark201,benchmark900` restricts the run to specific suites (default: all); a partial run replaces that architecture's stored results with just the suites it ran. `--reverse` (Reactant only, e.g. `manual_benchmarking.jl reactant-gpu --reverse`) additionally runs the reverse-mode benchmarks, which take hours. `--hi-res` (GPU and TPU only) adds T342 and T512 to the PrimitiveWet resolution sweep. `--SPEEDY` additionally runs the PrimitiveWet model over the resolutions of the resolution sweep in a configuration closer to SPEEDY.f90 (seasonal prescribed climatologies of sea surface temperature, land surface temperature and soil moisture, constant surface roughness, wind slowdown and gusts).

Each run replaces its own architecture's section in this `README.md` (including suites added by flags in an earlier run); results for other architectures are preserved via `benchmark_results.json`.

## Overview: PrimitiveWet resolution across architectures

Simulated years per wallclock day (SYPD) for the `PrimitiveWetModel` resolution sweep, one column per architecture. Each (T, L) configuration is reported for both the standard Legendre transform and fast Fourier transform (LT+FFT) and the single matrix transform (MT). Empty cells mean the architecture has not yet been benchmarked or that suite was skipped. Comparison figures across architectures are available on the documentation's `Benchmarks` page.

| T | L | Transform | cpu-arm | cpu-x86 | gpu-nvidia | gpu-amd | reactant-cpu | reactant-gpu |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 32 | 8 | LT+FFT | 1308 | 806 | 5961 | 3763 | 909 | 30033 |
| 32 | 8 | MT | 772 | 834 | 5853 | 3696 | 880 | 18378 |
| 43 | 8 | LT+FFT | 559 | 347 | 4522 | 2149 | 458 | 17233 |
| 43 | 8 | MT | 253 | 317 | 4468 | 2830 | 353 | 15171 |
| 64 | 8 | LT+FFT | 152 | 101 | 3025 | 828 | 159 | 9588 |
| 64 | 8 | MT | 52 | 68 | 3117 | 826 | 97 | 7139 |
| 86 | 8 | LT+FFT | 58 | 41 | 1254 | 246 | 68 | 5962 |
| 86 | 8 | MT | 14 | 22 | 1260 | 245 | 34 | 2829 |
| 86 | 16 | LT+FFT | 29 | 21 | 781 | 202 | 40 | 3324 |
| 86 | 16 | MT | 8.8 | 13 | 783 | 150 | 21 | 1799 |
| 86 | 24 | LT+FFT | 20 | 14 | 1357 | 196 | 22 | 2217 |
| 86 | 24 | MT | 6.1 | 9.1 | 1358 | 108 | 14 | 1282 |
| 128 | 8 | LT+FFT | 14 | 11 | 366 | 86 | 23 | 2077 |
| 128 | 8 | MT | 1.2 | — | 238 | 39 | — | 676 |
| 128 | 16 | LT+FFT | 7.6 | 5.7 | 612 | 80 | 13 | 1444 |
| 128 | 16 | MT | 1.2 | — | 288 | 23 | — | 573 |
| 128 | 24 | LT+FFT | 5.5 | 3.6 | 553 | 76 | 6.6 | 870 |
| 128 | 24 | MT | 1.0 | — | 255 | 17 | — | 485 |
| 171 | 8 | LT+FFT | 5.7 | 4.1 | 188 | 44 | 8.9 | 1121 |
| 171 | 8 | MT | — | — | 52 | — | — | 166 |
| 171 | 16 | LT+FFT | 3.1 | 2.1 | 313 | 41 | 4.9 | 583 |
| 171 | 16 | MT | — | — | 72 | — | — | 98 |
| 171 | 24 | LT+FFT | 1.9 | 1.4 | 273 | 38 | 2.4 | 373 |
| 171 | 24 | MT | — | — | 63 | — | — | 78 |
| 256 | 8 | LT+FFT | 1.4 | 1.0 | 69 | 16 | 2.4 | 366 |
| 256 | 16 | LT+FFT | 0.8 | 0.5 | 109 | 13 | 1.3 | 185 |
| 256 | 24 | LT+FFT | 0.5 | 0.3 | 91 | 11 | 0.7 | 124 |
| 342 | 8 | LT+FFT | — | — | — | — | — | 125 |
| 342 | 16 | LT+FFT | — | — | — | — | — | 72 |
| 342 | 24 | LT+FFT | — | — | — | — | — | 42 |
| 512 | 8 | LT+FFT | — | — | — | — | — | 33 |
| 512 | 16 | LT+FFT | — | — | — | — | — | 19 |
| 512 | 24 | LT+FFT | — | — | — | — | — | 12 |

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

Created for SpeedyWeather.jl v0.24.0-DEV on Thu, 05 Nov 2026 17:11:22.

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


Settings: 10 timed runs per configuration (`--repeats 10`), the PrimitiveWet resolution sweep 2× longer (`--length 2`).

### Models, default setups

| Model | truncation | L | Physics | Δt | SYPD | Runs | Spread | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| BarotropicModel | 32 | 1 | false | 1800 | 26958 | 10 | 2% | 725.97 KB |
| ShallowWaterModel | 32 | 1 | false | 2400 | 18843 | 10 | 3% | 908.33 KB |
| PrimitiveDryModel | 32 | 8 | true | 2400 | 1238 | 10 | 1% | 5.04 MB |
| PrimitiveWetModel | 32 | 8 | true | 2400 | 811 | 10 | 2% | 6.15 MB |

### Shallow water model, resolution

| Model | truncation | L | Rings | Δt | SYPD | Runs | Spread | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 32 | 1 | 48 | 2400 | 18670 | 10 | 8% | 908.33 KB |
| ShallowWaterModel | 43 | 1 | 64 | 1800 | 7808 | 10 | 7% | 1.59 MB |
| ShallowWaterModel | 64 | 1 | 96 | 1200 | 2091 | 10 | 1% | 3.57 MB |
| ShallowWaterModel | 86 | 1 | 128 | 900 | 868 | 10 | 2% | 6.49 MB |
| ShallowWaterModel | 128 | 1 | 192 | 600 | 225 | 10 | 2% | 15.36 MB |
| ShallowWaterModel | 171 | 1 | 256 | 450 | 83 | 10 | 3% | 29.00 MB |
| ShallowWaterModel | 256 | 1 | 384 | 300 | 18 | 10 | 1% | 73.09 MB |

### Primitive wet model, resolution

| Model | truncation | L | Rings | Transform | Δt | SYPD | Runs | Spread | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | 48 | default | 2400 | 806 | 10 | 1% | 6.15 MB |
| PrimitiveWetModel | 43 | 8 | 64 | default | 1800 | 347 | 10 | 1% | 10.37 MB |
| PrimitiveWetModel | 64 | 8 | 96 | default | 1200 | 101 | 10 | 1% | 22.02 MB |
| PrimitiveWetModel | 86 | 8 | 128 | default | 900 | 41 | 10 | 1% | 38.26 MB |
| PrimitiveWetModel | 128 | 8 | 192 | default | 600 | 11 | 10 | 1% | 84.15 MB |
| PrimitiveWetModel | 171 | 8 | 256 | default | 450 | 4.1 | 10 | 2% | 149.16 MB |
| PrimitiveWetModel | 256 | 8 | 384 | default | 300 | 1.0 | 10 | 1% | 338.28 MB |
| PrimitiveWetModel | 86 | 16 | 128 | default | 900 | 21 | 10 | 1% | 64.09 MB |
| PrimitiveWetModel | 128 | 16 | 192 | default | 600 | 5.7 | 10 | 1% | 140.14 MB |
| PrimitiveWetModel | 171 | 16 | 256 | default | 450 | 2.1 | 10 | 3% | 247.09 MB |
| PrimitiveWetModel | 256 | 16 | 384 | default | 300 | 0.5 | 10 | 1% | 554.60 MB |
| PrimitiveWetModel | 86 | 24 | 128 | default | 900 | 14 | 10 | 18% | 89.96 MB |
| PrimitiveWetModel | 128 | 24 | 192 | default | 600 | 3.6 | 10 | 11% | 196.21 MB |
| PrimitiveWetModel | 171 | 24 | 256 | default | 450 | 1.4 | 10 | 2% | 345.11 MB |
| PrimitiveWetModel | 256 | 24 | 384 | default | 300 | 0.3 | 10 | 1% | 771.06 MB |
| PrimitiveWetModel | 32 | 8 | 48 | matrix | 2400 | 834 | 10 | 3% | 34.60 MB |
| PrimitiveWetModel | 43 | 8 | 64 | matrix | 1800 | 317 | 10 | 4% | 93.50 MB |
| PrimitiveWetModel | 64 | 8 | 96 | matrix | 1200 | 68 | 10 | 4% | 397.51 MB |
| PrimitiveWetModel | 86 | 8 | 128 | matrix | 900 | 22 | 10 | 4% | 1.18 GB |
| PrimitiveWetModel | 86 | 16 | 128 | matrix | 900 | 13 | 10 | 2% | 1.21 GB |
| PrimitiveWetModel | 86 | 24 | 128 | matrix | 900 | 9.1 | 10 | 3% | 1.24 GB |

### Primitive Equation, Float32 vs Float64

| Model | NF | truncation | L | Δt | SYPD | Runs | Spread | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | Float32 | 32 | 8 | 2400 | 796 | 10 | 8% | 6.15 MB |
| PrimitiveWetModel | Float64 | 32 | 8 | 2400 | 693 | 10 | 2% | 11.77 MB |

### Grids

| Model | truncation | L | Grid | Rings | Δt | SYPD | Runs | Spread | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 64 | 8 | FullGaussianGrid | 96 | 1200 | 69 | 10 | 2% | 32.51 MB |
| PrimitiveWetModel | 64 | 8 | FullClenshawGrid | 127 | 1200 | 41 | 10 | 1% | 54.47 MB |
| PrimitiveWetModel | 64 | 8 | OctahedralGaussianGrid | 96 | 1200 | 104 | 10 | 1% | 22.02 MB |
| PrimitiveWetModel | 64 | 8 | OctahedralClenshawGrid | 127 | 1200 | 58 | 10 | 1% | 34.69 MB |
| PrimitiveWetModel | 64 | 8 | HEALPixGrid | 127 | 1200 | 468 | 10 | 4% | 25.47 MB |
| PrimitiveWetModel | 64 | 8 | OctaHEALPixGrid | 127 | 1200 | 59 | 10 | 8% | 31.76 MB |

### Number of vertical layers

| Model | truncation | L | Δt | SYPD | Runs | Spread | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 4 | 2400 | 1150 | 10 | 3% | 4.05 MB |
| PrimitiveWetModel | 32 | 8 | 2400 | 818 | 10 | 1% | 6.15 MB |
| PrimitiveWetModel | 32 | 12 | 2400 | 599 | 10 | 4% | 8.26 MB |
| PrimitiveWetModel | 32 | 16 | 2400 | 477 | 10 | 1% | 10.37 MB |

### PrimitiveDryModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Runs | Spread | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveDryModel | 32 | 8 | true | true | 2400 | 1254 | 10 | 2% | 5.04 MB |
| PrimitiveDryModel | 32 | 8 | true | false | 2400 | 2322 | 10 | 2% | 5.04 MB |
| PrimitiveDryModel | 32 | 8 | false | true | 2400 | 1444 | 10 | 1% | 5.04 MB |

### PrimitiveWetModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Runs | Spread | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | true | true | 2400 | 808 | 10 | 1% | 6.15 MB |
| PrimitiveWetModel | 32 | 8 | true | false | 2400 | 1777 | 10 | 2% | 6.15 MB |
| PrimitiveWetModel | 32 | 8 | false | true | 2400 | 981 | 10 | 1% | 6.15 MB |

### Individual dynamics functions


#### PrimitiveWetModel | Float32 | T32 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| --- | --- | --- | --- |
| pressure_gradient_flux! | 16.599 μs| 0 bytes| 0 |
| linear_virtual_temperature! | 3.583 μs| 0 bytes| 0 |
| geopotential! | 9.960 μs| 384 bytes| 6 |
| vertical_integration! | 13.710 μs| 0 bytes| 0 |
| surface_pressure_tendency! | 18.359 μs| 15.66 KiB| 96 |
| vertical_velocity! | 58.509 μs| 0 bytes| 0 |
| linear_pressure_gradient! | 3.255 μs| 0 bytes| 0 |
| vertical_advection! | 159.958 μs| 2.44 KiB| 32 |
| vordiv_tendencies! | 409.305 μs| 218.62 KiB| 284 |
| temperature_tendency! | 518.913 μs| 324.75 KiB| 401 |
| humidity_tendency! | 497.194 μs| 324.08 KiB| 396 |
| bernoulli_potential! | 163.948 μs| 107.61 KiB| 129 |

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

## Architecture: `reactant-cpu`

Created for SpeedyWeather.jl v0.24.0-DEV on Thu, 05 Nov 2026 21:41:23.

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
  JULIA_DEPOT_PATH = /p/tmp/maxgelbr/julia-depot-cpu:/home/maxgelbr/.julia:
  LD_LIBRARY_PATH = /usr/local/lib:/usr/local/lib:
```

Reactant backend: `cpu` (via `Reactant.set_default_backend("cpu")`), 1 device(s): cpu (one used), Reactant v0.2.305, Reactant_jll v0.0.441+0.


Settings: 10 timed runs per configuration (`--repeats 10`), the PrimitiveWet resolution sweep 2× longer (`--length 2`).

### Models, default setups

| Model | truncation | L | Physics | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| BarotropicModel | 32 | 1 | false | 1800 | 9815 | 10 | 3% | 1.81 MB | 29.0 s |
| ShallowWaterModel | 32 | 1 | false | 2400 | 9375 | 10 | 1% | 2.00 MB | 35.1 s |
| PrimitiveDryModel | 32 | 8 | true | 2400 | 1265 | 10 | 5% | 5.84 MB | 178.9 s |
| PrimitiveWetModel | 32 | 8 | true | 2400 | 856 | 10 | 3% | 6.95 MB | 213.3 s |

### Shallow water model, resolution

| Model | truncation | L | Rings | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 32 | 1 | 48 | 2400 | 9154 | 10 | 6% | 2.00 MB | 1.5 s |
| ShallowWaterModel | 43 | 1 | 64 | 1786 | 3510 | 10 | 4% | 4.28 MB | 7.4 s |
| ShallowWaterModel | 64 | 1 | 96 | 1200 | 959 | 10 | 22% | 12.85 MB | 7.6 s |
| ShallowWaterModel | 86 | 1 | 128 | 893 | 394 | 10 | 6% | 29.04 MB | 7.8 s |
| ShallowWaterModel | 128 | 1 | 192 | 600 | 92 | 10 | 3% | 91.85 MB | 8.4 s |
| ShallowWaterModel | 171 | 1 | 256 | 449 | 48 | 10 | 4% | 211.92 MB | 9.0 s |
| ShallowWaterModel | 256 | 1 | 384 | 300 | 9.0 | 10 | 4% | 693.70 MB | 10.1 s |

### Primitive wet model, resolution

| Model | truncation | L | Rings | Transform | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | 48 | default | 2400 | 909 | 10 | 29% | 6.95 MB | 9.0 s |
| PrimitiveWetModel | 43 | 8 | 64 | default | 1786 | 458 | 10 | 16% | 12.56 MB | 30.9 s |
| PrimitiveWetModel | 64 | 8 | 96 | default | 1200 | 159 | 10 | 30% | 30.21 MB | 30.5 s |
| PrimitiveWetModel | 86 | 8 | 128 | default | 893 | 68 | 10 | 10% | 58.89 MB | 31.1 s |
| PrimitiveWetModel | 128 | 8 | 192 | default | 600 | 23 | 10 | 6% | 156.39 MB | 29.7 s |
| PrimitiveWetModel | 171 | 8 | 256 | default | 449 | 8.9 | 10 | 3% | 324.57 MB | 30.6 s |
| PrimitiveWetModel | 256 | 8 | 384 | default | 300 | 2.4 | 10 | 1% | 942.13 MB | 34.6 s |
| PrimitiveWetModel | 86 | 16 | 128 | default | 893 | 40 | 10 | 3% | 82.55 MB | 181.4 s |
| PrimitiveWetModel | 128 | 16 | 192 | default | 600 | 13 | 10 | 10% | 207.56 MB | 34.3 s |
| PrimitiveWetModel | 171 | 16 | 256 | default | 449 | 4.9 | 10 | 5% | 413.97 MB | 36.3 s |
| PrimitiveWetModel | 256 | 16 | 384 | default | 300 | 1.3 | 10 | 6% | 1.14 GB | 45.9 s |
| PrimitiveWetModel | 86 | 24 | 128 | default | 893 | 22 | 10 | 2% | 106.25 MB | 218.2 s |
| PrimitiveWetModel | 128 | 24 | 192 | default | 600 | 6.6 | 10 | 5% | 258.80 MB | 68.4 s |
| PrimitiveWetModel | 171 | 24 | 256 | default | 449 | 2.4 | 10 | 6% | 503.46 MB | 65.5 s |
| PrimitiveWetModel | 256 | 24 | 384 | default | 300 | 0.7 | 10 | 9% | 1.34 GB | 76.2 s |
| PrimitiveWetModel | 32 | 8 | 48 | matrix | 2400 | 880 | 10 | 4% | 34.15 MB | 189.6 s |
| PrimitiveWetModel | 43 | 8 | 64 | matrix | 1786 | 353 | 10 | 4% | 92.72 MB | 13.7 s |
| PrimitiveWetModel | 64 | 8 | 96 | matrix | 1200 | 97 | 10 | 2% | 395.75 MB | 13.6 s |
| PrimitiveWetModel | 86 | 8 | 128 | matrix | 893 | 34 | 10 | 2% | 1.18 GB | 14.2 s |
| PrimitiveWetModel | 86 | 16 | 128 | matrix | 893 | 21 | 10 | 2% | 1.21 GB | 168.7 s |
| PrimitiveWetModel | 86 | 24 | 128 | matrix | 893 | 14 | 10 | 2% | 1.23 GB | 205.7 s |

### Primitive Equation, Float32 vs Float64

| Model | NF | truncation | L | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | Float32 | 32 | 8 | 2400 | 834 | 10 | 1% | 6.95 MB | 8.9 s |
| PrimitiveWetModel | Float64 | 32 | 8 | 2400 | 586 | 10 | 5% | 13.75 MB | 252.7 s |

### Grids

| Model | truncation | L | Grid | Rings | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 64 | 8 | FullGaussianGrid | 96 | 1200 | 90 | 10 | 2% | 39.94 MB | 211.5 s |
| PrimitiveWetModel | 64 | 8 | FullClenshawGrid | 127 | 1200 | 75 | 10 | 2% | 67.33 MB | 212.4 s |
| PrimitiveWetModel | 64 | 8 | OctahedralGaussianGrid | 96 | 1200 | 150 | 10 | 2% | 30.21 MB | 10.1 s |
| PrimitiveWetModel | 64 | 8 | OctahedralClenshawGrid | 127 | 1200 | 98 | 10 | 1% | 48.59 MB | 211.5 s |
| PrimitiveWetModel | 64 | 8 | HEALPixGrid | 127 | 1200 | 107 | 10 | 5% | 30.54 MB | 211.1 s |
| PrimitiveWetModel | 64 | 8 | OctaHEALPixGrid | 127 | 1200 | 71 | 10 | 62% | 44.68 MB | 209.9 s |

### Number of vertical layers

| Model | truncation | L | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 4 | 2400 | 1441 | 10 | 6% | 5.01 MB | 163.5 s |
| PrimitiveWetModel | 32 | 8 | 2400 | 841 | 10 | 5% | 6.95 MB | 10.5 s |
| PrimitiveWetModel | 32 | 12 | 2400 | 712 | 10 | 5% | 8.90 MB | 169.1 s |
| PrimitiveWetModel | 32 | 16 | 2400 | 496 | 10 | 2% | 10.85 MB | 33.0 s |

### PrimitiveDryModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveDryModel | 32 | 8 | true | true | 2400 | 1207 | 10 | 1% | 5.84 MB | 8.4 s |
| PrimitiveDryModel | 32 | 8 | true | false | 2400 | 1375 | 10 | 2% | 5.84 MB | 4.6 s |
| PrimitiveDryModel | 32 | 8 | false | true | 2400 | 1548 | 10 | 4% | 5.84 MB | 8.1 s |

### PrimitiveWetModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | true | true | 2400 | 830 | 10 | 2% | 6.95 MB | 9.1 s |
| PrimitiveWetModel | 32 | 8 | true | false | 2400 | 1176 | 10 | 4% | 6.95 MB | 5.2 s |
| PrimitiveWetModel | 32 | 8 | false | true | 2400 | 1068 | 10 | 1% | 6.95 MB | 11.2 s |

### Individual dynamics functions


#### PrimitiveWetModel | Float32 | T32 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| --- | --- | --- | --- |
| pressure_gradient_flux! | 42.649 μs| 3.84 KiB| 25 |
| linear_virtual_temperature! | 51.261 μs| 4.72 KiB| 37 |
| geopotential! | 46.931 μs| 6.14 KiB| 36 |
| vertical_integration! | 84.501 μs| 16.48 KiB| 99 |
| surface_pressure_tendency! | 411.971 μs| 5.81 KiB| 68 |
| vertical_velocity! | 45.341 μs| 4.23 KiB| 32 |
| linear_pressure_gradient! | 35.770 μs| 3.95 KiB| 25 |
| vertical_advection! | 526.061 μs| 4.34 KiB| 53 |
| vordiv_tendencies! | 653.265 μs| 9.12 KiB| 102 |
| temperature_tendency! | 943.269 μs| 6.44 KiB| 80 |
| humidity_tendency! | 906.840 μs| 6.44 KiB| 80 |
| bernoulli_potential! | 401.237 μs| 7.16 KiB| 79 |

### Primitive wet model, reverse mode (AD) with Reactant

Wallclock time of the time steps forward and of the reverse-mode gradient with respect to the initial state over the same time steps, minimum over up to 3 runs, compilation excluded. The loss is the squared norm of the vorticity after the time steps.

| truncation | L | Rings | Δt | Steps | Checkpointing | Forward | Gradient | Gradient/Forward | Compile forward | Compile gradient | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 32 | 8 | 48 | 2400 | 50 | Binomial(32) | 0.399 s | 1.61 s | 4.0 | 25.7 s | 277.5 s | ok |
| 32 | 8 | 48 | 2400 | 100 | Binomial(32) | 0.816 s | 3.45 s | 4.2 | 10.2 s | 35.9 s | ok |
| 32 | 8 | 48 | 2400 | 200 | Binomial(32) | 1.6 s | 7.23 s | 4.5 | 12.1 s | 37.9 s | ok |
| 32 | 8 | 48 | 2400 | 2000 | Binomial(32) | 15.9 s | 90.7 s | 5.7 | 12.1 s | 36.9 s | non-finite or zero gradient |
| 64 | 8 | 96 | 1200 | 50 | Binomial(32) | 1.19 s | 4.67 s | 3.9 | 16.2 s | 51.4 s | ok |
| 64 | 8 | 96 | 1200 | 100 | Binomial(32) | 2.31 s | 10.2 s | 4.4 | 11.3 s | 42.3 s | ok |
| 64 | 8 | 96 | 1200 | 200 | Binomial(32) | 4.63 s | 21.4 s | 4.6 | 11.8 s | 41.8 s | ok |
| 64 | 8 | 96 | 1200 | 2000 | Binomial(32) | 47.3 s | 266 s | 5.6 | 11.2 s | 38.9 s | non-finite or zero gradient |
| 128 | 8 | 192 | 600 | 50 | Binomial(32) | 3.78 s | 16 s | 4.2 | 11.7 s | 36.3 s | ok |
| 128 | 8 | 192 | 600 | 100 | Binomial(32) | 7.5 s | 35 s | 4.7 | 15.4 s | 48.9 s | ok |
| 128 | 8 | 192 | 600 | 200 | Binomial(32) | 15.2 s | 73.1 s | 4.8 | 12.2 s | 38.2 s | ok |
| 128 | 8 | 192 | 600 | 2000 | Binomial(32) | 156 s | 939 s | 6.0 | 16.4 s | 40.4 s | non-finite or zero gradient |

### SPEEDY.f90-like configuration

The PrimitiveWet model in a configuration closer to SPEEDY.f90, at the resolutions of the PrimitiveWet resolution sweep with the default transform: seasonal prescribed climatologies of sea surface temperature (`SeasonalOceanClimatology`), land surface temperature (`SeasonalLandTemperature`) and soil moisture (`SeasonalSoilMoisture`) instead of the slab ocean and the land bucket model, constant surface roughness lengths over land and ocean (`ConstantSurfaceRoughness`) instead of the learned ones, and SPEEDY.f90's near-surface wind of 0.95 times the lowermost-layer wind (`ConstantWindSlowdown`) with 5 m/s gusts (`ConstantWindGust`) instead of the learned wind slowdown and Beljaars gusts; sea ice stays thermodynamic, driven by the prescribed sea surface temperature.

| Model | truncation | L | Rings | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | 48 | 2400 | 851 | 10 | 24% | 7.47 MB | 196.4 s |
| PrimitiveWetModel | 43 | 8 | 64 | 1786 | 427 | 10 | 1% | 13.42 MB | 19.9 s |
| PrimitiveWetModel | 64 | 8 | 96 | 1200 | 149 | 10 | 3% | 32.00 MB | 16.9 s |
| PrimitiveWetModel | 86 | 8 | 128 | 893 | 64 | 10 | 2% | 61.96 MB | 22.0 s |
| PrimitiveWetModel | 128 | 8 | 192 | 600 | 22 | 10 | 3% | 163.01 MB | 18.6 s |
| PrimitiveWetModel | 171 | 8 | 256 | 449 | 8.5 | 10 | 5% | 336.08 MB | 24.7 s |
| PrimitiveWetModel | 256 | 8 | 384 | 300 | 2.3 | 10 | 3% | 967.45 MB | 27.4 s |
| PrimitiveWetModel | 86 | 16 | 128 | 893 | 40 | 10 | 2% | 85.61 MB | 162.2 s |
| PrimitiveWetModel | 128 | 16 | 192 | 600 | 12 | 10 | 7% | 214.17 MB | 33.6 s |
| PrimitiveWetModel | 171 | 16 | 256 | 449 | 4.6 | 10 | 4% | 425.47 MB | 31.4 s |
| PrimitiveWetModel | 256 | 16 | 384 | 300 | 1.3 | 10 | 8% | 1.16 GB | 38.0 s |
| PrimitiveWetModel | 86 | 24 | 128 | 893 | 21 | 10 | 2% | 109.32 MB | 189.5 s |
| PrimitiveWetModel | 128 | 24 | 192 | 600 | 6.4 | 10 | 3% | 265.41 MB | 51.3 s |
| PrimitiveWetModel | 171 | 24 | 256 | 449 | 2.4 | 10 | 3% | 514.96 MB | 52.7 s |
| PrimitiveWetModel | 256 | 24 | 384 | 300 | 0.7 | 10 | 6% | 1.36 GB | 57.3 s |

## Architecture: `reactant-gpu`

Created for SpeedyWeather.jl v0.24.0-DEV on Thu, 05 Nov 2026 19:02:34.

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
  0: NVIDIA H100 80GB HBM3 (sm_90, 18.929 GiB / 79.647 GiB available)
     compiles to sm_90a / PTX 9.3 (LLVM: sm_90a / PTX 9.0)
```

Reactant backend: `gpu` (via `Reactant.set_default_backend("gpu")`), 1 device(s): NVIDIA H100 80GB HBM3 (one used), Reactant v0.2.305, Reactant_jll v0.0.441+0.


Settings: 10 timed runs per configuration (`--repeats 10`), the PrimitiveWet resolution sweep 2× longer (`--length 2`).

### Models, default setups

| Model | truncation | L | Physics | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| BarotropicModel | 32 | 1 | false | 1800 | 141752 | 10 | 2% | 1.81 MB | 28.0 s |
| ShallowWaterModel | 32 | 1 | false | 2400 | 133520 | 10 | 1% | 2.00 MB | 32.2 s |
| PrimitiveDryModel | 32 | 8 | true | 2400 | 37923 | 10 | 1% | 5.84 MB | 181.6 s |
| PrimitiveWetModel | 32 | 8 | true | 2400 | 29429 | 10 | 11% | 6.95 MB | 214.0 s |

### Shallow water model, resolution

| Model | truncation | L | Rings | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 32 | 1 | 48 | 2400 | 134189 | 10 | 1% | 2.00 MB | 1.5 s |
| ShallowWaterModel | 43 | 1 | 64 | 1786 | 93225 | 10 | 1% | 4.28 MB | 7.4 s |
| ShallowWaterModel | 64 | 1 | 96 | 1200 | 57069 | 10 | 1% | 12.85 MB | 7.7 s |
| ShallowWaterModel | 86 | 1 | 128 | 893 | 34800 | 10 | 5% | 29.04 MB | 7.6 s |
| ShallowWaterModel | 128 | 1 | 192 | 600 | 13188 | 10 | 1% | 91.85 MB | 7.9 s |
| ShallowWaterModel | 171 | 1 | 256 | 449 | 6432 | 10 | 1% | 211.92 MB | 8.4 s |
| ShallowWaterModel | 256 | 1 | 384 | 300 | 1977 | 10 | 1% | 693.70 MB | 10.1 s |

### Primitive wet model, resolution

| Model | truncation | L | Rings | Transform | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | 48 | default | 2400 | 30033 | 10 | 4% | 6.95 MB | 14.2 s |
| PrimitiveWetModel | 43 | 8 | 64 | default | 1786 | 17233 | 10 | 4% | 12.56 MB | 35.3 s |
| PrimitiveWetModel | 64 | 8 | 96 | default | 1200 | 9588 | 10 | 1% | 30.21 MB | 32.8 s |
| PrimitiveWetModel | 86 | 8 | 128 | default | 893 | 5962 | 10 | 2% | 58.89 MB | 33.9 s |
| PrimitiveWetModel | 128 | 8 | 192 | default | 600 | 2077 | 10 | 2% | 156.39 MB | 34.3 s |
| PrimitiveWetModel | 171 | 8 | 256 | default | 449 | 1121 | 10 | 1% | 324.57 MB | 34.1 s |
| PrimitiveWetModel | 256 | 8 | 384 | default | 300 | 366 | 10 | 1% | 942.13 MB | 36.5 s |
| PrimitiveWetModel | 342 | 8 | 512 | default | 225 | 125 | 10 | 1% | 2.06 GB | 40.3 s |
| PrimitiveWetModel | 512 | 8 | 768 | default | 150 | 33 | 10 | 0% | 6.36 GB | 49.0 s |
| PrimitiveWetModel | 86 | 16 | 128 | default | 893 | 3324 | 10 | 5% | 82.55 MB | 186.7 s |
| PrimitiveWetModel | 128 | 16 | 192 | default | 600 | 1444 | 10 | 4% | 207.56 MB | 40.7 s |
| PrimitiveWetModel | 171 | 16 | 256 | default | 449 | 583 | 10 | 2% | 413.97 MB | 41.8 s |
| PrimitiveWetModel | 256 | 16 | 384 | default | 300 | 185 | 10 | 1% | 1.14 GB | 44.1 s |
| PrimitiveWetModel | 342 | 16 | 512 | default | 225 | 72 | 10 | 1% | 2.41 GB | 47.8 s |
| PrimitiveWetModel | 512 | 16 | 768 | default | 150 | 19 | 10 | 0% | 7.14 GB | 57.4 s |
| PrimitiveWetModel | 86 | 24 | 128 | default | 893 | 2217 | 10 | 4% | 106.25 MB | 215.9 s |
| PrimitiveWetModel | 128 | 24 | 192 | default | 600 | 870 | 10 | 2% | 258.80 MB | 71.6 s |
| PrimitiveWetModel | 171 | 24 | 256 | default | 449 | 373 | 10 | 1% | 503.46 MB | 71.1 s |
| PrimitiveWetModel | 256 | 24 | 384 | default | 300 | 124 | 10 | 1% | 1.34 GB | 75.0 s |
| PrimitiveWetModel | 342 | 24 | 512 | default | 225 | 42 | 10 | 1% | 2.76 GB | 81.1 s |
| PrimitiveWetModel | 512 | 24 | 768 | default | 150 | 12 | 10 | 0% | 7.91 GB | 93.1 s |
| PrimitiveWetModel | 32 | 8 | 48 | matrix | 2400 | 18378 | 10 | 1% | 34.16 MB | 192.4 s |
| PrimitiveWetModel | 43 | 8 | 64 | matrix | 1786 | 15171 | 10 | 1% | 92.72 MB | 17.4 s |
| PrimitiveWetModel | 64 | 8 | 96 | matrix | 1200 | 7139 | 10 | 1% | 395.75 MB | 17.5 s |
| PrimitiveWetModel | 86 | 8 | 128 | matrix | 893 | 2829 | 10 | 1% | 1.18 GB | 18.7 s |
| PrimitiveWetModel | 128 | 8 | 192 | matrix | 600 | 676 | 10 | 2% | 5.48 GB | 17.6 s |
| PrimitiveWetModel | 171 | 8 | 256 | matrix | 449 | 166 | 10 | 1% | 16.83 GB | 19.4 s |
| PrimitiveWetModel | 86 | 16 | 128 | matrix | 893 | 1799 | 10 | 5% | 1.21 GB | 181.5 s |
| PrimitiveWetModel | 128 | 16 | 192 | matrix | 600 | 573 | 10 | 2% | 5.54 GB | 29.0 s |
| PrimitiveWetModel | 171 | 16 | 256 | matrix | 449 | 98 | 10 | 0% | 16.93 GB | 31.6 s |
| PrimitiveWetModel | 86 | 24 | 128 | matrix | 893 | 1282 | 10 | 3% | 1.23 GB | 209.6 s |
| PrimitiveWetModel | 128 | 24 | 192 | matrix | 600 | 485 | 10 | 1% | 5.60 GB | 59.6 s |
| PrimitiveWetModel | 171 | 24 | 256 | matrix | 449 | 78 | 10 | 0% | 17.03 GB | 59.8 s |

### Primitive Equation, Float32 vs Float64

| Model | NF | truncation | L | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | Float32 | 32 | 8 | 2400 | 29232 | 10 | 7% | 6.95 MB | 15.2 s |
| PrimitiveWetModel | Float64 | 32 | 8 | 2400 | 19510 | 10 | 7% | 13.75 MB | 254.8 s |

### Grids

| Model | truncation | L | Grid | Rings | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 64 | 8 | FullGaussianGrid | 96 | 1200 | 8563 | 10 | 3% | 39.94 MB | 213.5 s |
| PrimitiveWetModel | 64 | 8 | FullClenshawGrid | 127 | 1200 | 7324 | 10 | 2% | 67.33 MB | 212.6 s |
| PrimitiveWetModel | 64 | 8 | OctahedralGaussianGrid | 96 | 1200 | 9388 | 10 | 2% | 30.21 MB | 13.1 s |
| PrimitiveWetModel | 64 | 8 | OctahedralClenshawGrid | 127 | 1200 | 8223 | 10 | 3% | 48.59 MB | 215.5 s |
| PrimitiveWetModel | 64 | 8 | HEALPixGrid | 127 | 1200 | 9413 | 10 | 2% | 30.54 MB | 216.1 s |
| PrimitiveWetModel | 64 | 8 | OctaHEALPixGrid | 127 | 1200 | 8700 | 10 | 9% | 44.68 MB | 213.0 s |

### Number of vertical layers

| Model | truncation | L | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 4 | 2400 | 41631 | 10 | 1% | 5.01 MB | 166.3 s |
| PrimitiveWetModel | 32 | 8 | 2400 | 29442 | 10 | 4% | 6.95 MB | 14.7 s |
| PrimitiveWetModel | 32 | 12 | 2400 | 21825 | 10 | 4% | 8.90 MB | 183.1 s |
| PrimitiveWetModel | 32 | 16 | 2400 | 10836 | 10 | 1% | 10.85 MB | 42.8 s |

### PrimitiveDryModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveDryModel | 32 | 8 | true | true | 2400 | 37890 | 10 | 6% | 5.84 MB | 12.9 s |
| PrimitiveDryModel | 32 | 8 | true | false | 2400 | 45133 | 10 | 5% | 5.84 MB | 11.5 s |
| PrimitiveDryModel | 32 | 8 | false | true | 2400 | 54417 | 10 | 9% | 5.84 MB | 10.1 s |

### PrimitiveWetModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | true | true | 2400 | 28475 | 10 | 9% | 6.95 MB | 14.9 s |
| PrimitiveWetModel | 32 | 8 | true | false | 2400 | 44795 | 10 | 8% | 6.95 MB | 9.3 s |
| PrimitiveWetModel | 32 | 8 | false | true | 2400 | 42615 | 10 | 1% | 6.95 MB | 13.9 s |

### Individual dynamics functions


#### PrimitiveWetModel | Float32 | T32 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| --- | --- | --- | --- |
| pressure_gradient_flux! | 81.330 μs| 3.84 KiB| 25 |
| linear_virtual_temperature! | 98.097 μs| 4.72 KiB| 37 |
| geopotential! | 89.250 μs| 6.14 KiB| 36 |
| vertical_integration! | 96.429 μs| 16.48 KiB| 99 |
| surface_pressure_tendency! | 85.950 μs| 5.81 KiB| 68 |
| vertical_velocity! | 88.260 μs| 4.23 KiB| 32 |
| linear_pressure_gradient! | 72.979 μs| 3.95 KiB| 25 |
| vertical_advection! | 96.720 μs| 4.34 KiB| 53 |
| vordiv_tendencies! | 115.410 μs| 9.12 KiB| 102 |
| temperature_tendency! | 109.359 μs| 6.44 KiB| 80 |
| humidity_tendency! | 104.050 μs| 6.44 KiB| 80 |
| bernoulli_potential! | 97.539 μs| 7.16 KiB| 79 |

### Primitive wet model, reverse mode (AD) with Reactant

Wallclock time of the time steps forward and of the reverse-mode gradient with respect to the initial state over the same time steps, minimum over up to 3 runs, compilation excluded. The loss is the squared norm of the vorticity after the time steps.

| truncation | L | Rings | Δt | Steps | Checkpointing | Forward | Gradient | Gradient/Forward | Compile forward | Compile gradient | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 32 | 8 | 48 | 2400 | 50 | Binomial(32) | 0.0124 s | 0.0914 s | 7.3 | 29.6 s | 298.4 s | ok |
| 32 | 8 | 48 | 2400 | 100 | Binomial(32) | 0.0218 s | 0.177 s | 8.1 | 16.5 s | 55.4 s | ok |
| 32 | 8 | 48 | 2400 | 200 | Binomial(32) | 0.0442 s | 0.335 s | 7.6 | 16.2 s | 57.3 s | ok |
| 32 | 8 | 48 | 2400 | 2000 | Binomial(32) | 0.412 s | 3.56 s | 8.6 | 15.9 s | 56.6 s | non-finite or zero gradient |
| 64 | 8 | 96 | 1200 | 50 | Binomial(32) | 0.0191 s | 0.121 s | 6.4 | 15.4 s | 58.3 s | ok |
| 64 | 8 | 96 | 1200 | 100 | Binomial(32) | 0.0362 s | 0.237 s | 6.6 | 15.8 s | 57.6 s | ok |
| 64 | 8 | 96 | 1200 | 200 | Binomial(32) | 0.0703 s | 0.47 s | 6.7 | 15.5 s | 57.9 s | ok |
| 64 | 8 | 96 | 1200 | 2000 | Binomial(32) | 0.663 s | 5.23 s | 7.9 | 15.6 s | 57.8 s | non-finite or zero gradient |
| 128 | 8 | 192 | 600 | 50 | Binomial(32) | 0.0387 s | 0.202 s | 5.2 | 17.4 s | 59.3 s | ok |
| 128 | 8 | 192 | 600 | 100 | Binomial(32) | 0.0743 s | 0.408 s | 5.5 | 17.1 s | 58.2 s | ok |
| 128 | 8 | 192 | 600 | 200 | Binomial(32) | 0.146 s | 0.82 s | 5.6 | 17.2 s | 58.2 s | ok |
| 128 | 8 | 192 | 600 | 2000 | Binomial(32) | 1.42 s | 9.45 s | 6.6 | 17.2 s | 58.6 s | non-finite or zero gradient |
| 256 | 8 | 384 | 300 | 50 | Binomial(32) | 0.119 s | 0.58 s | 4.9 | 19.6 s | 61.4 s | ok |
| 256 | 8 | 384 | 300 | 100 | Binomial(32) | 0.234 s | 1.22 s | 5.2 | 19.6 s | 60.1 s | ok |
| 256 | 8 | 384 | 300 | 200 | Binomial(32) | 0.464 s | 2.44 s | 5.2 | 19.6 s | 60.3 s | ok |
| 256 | 8 | 384 | 300 | 2000 | Binomial(32) | 4.62 s | 29.1 s | 6.3 | 19.6 s | 60.0 s | non-finite or zero gradient |

### SPEEDY.f90-like configuration

The PrimitiveWet model in a configuration closer to SPEEDY.f90, at the resolutions of the PrimitiveWet resolution sweep with the default transform: seasonal prescribed climatologies of sea surface temperature (`SeasonalOceanClimatology`), land surface temperature (`SeasonalLandTemperature`) and soil moisture (`SeasonalSoilMoisture`) instead of the slab ocean and the land bucket model, constant surface roughness lengths over land and ocean (`ConstantSurfaceRoughness`) instead of the learned ones, and SPEEDY.f90's near-surface wind of 0.95 times the lowermost-layer wind (`ConstantWindSlowdown`) with 5 m/s gusts (`ConstantWindGust`) instead of the learned wind slowdown and Beljaars gusts; sea ice stays thermodynamic, driven by the prescribed sea surface temperature.

| Model | truncation | L | Rings | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | 48 | 2400 | 28623 | 10 | 2% | 7.47 MB | 199.1 s |
| PrimitiveWetModel | 43 | 8 | 64 | 1786 | 16833 | 10 | 1% | 13.42 MB | 21.0 s |
| PrimitiveWetModel | 64 | 8 | 96 | 1200 | 9181 | 10 | 2% | 32.01 MB | 19.9 s |
| PrimitiveWetModel | 86 | 8 | 128 | 893 | 5679 | 10 | 6% | 61.96 MB | 21.6 s |
| PrimitiveWetModel | 128 | 8 | 192 | 600 | 1949 | 10 | 5% | 163.01 MB | 21.7 s |
| PrimitiveWetModel | 171 | 8 | 256 | 449 | 1064 | 10 | 5% | 336.08 MB | 21.2 s |
| PrimitiveWetModel | 256 | 8 | 384 | 300 | 356 | 10 | 2% | 967.45 MB | 23.8 s |
| PrimitiveWetModel | 342 | 8 | 512 | 225 | 124 | 10 | 1% | 2.11 GB | 27.0 s |
| PrimitiveWetModel | 512 | 8 | 768 | 150 | 33 | 10 | 1% | 6.46 GB | 35.6 s |
| PrimitiveWetModel | 86 | 16 | 128 | 893 | 3088 | 10 | 7% | 85.61 MB | 167.1 s |
| PrimitiveWetModel | 128 | 16 | 192 | 600 | 1375 | 10 | 7% | 214.17 MB | 30.1 s |
| PrimitiveWetModel | 171 | 16 | 256 | 449 | 544 | 10 | 5% | 425.47 MB | 30.3 s |
| PrimitiveWetModel | 256 | 16 | 384 | 300 | 180 | 10 | 1% | 1.16 GB | 32.9 s |
| PrimitiveWetModel | 342 | 16 | 512 | 225 | 69 | 10 | 1% | 2.45 GB | 37.0 s |
| PrimitiveWetModel | 512 | 16 | 768 | 150 | 19 | 10 | 0% | 7.24 GB | 46.9 s |
| PrimitiveWetModel | 86 | 24 | 128 | 893 | 2071 | 10 | 7% | 109.32 MB | 200.7 s |
| PrimitiveWetModel | 128 | 24 | 192 | 600 | 823 | 10 | 7% | 265.41 MB | 62.0 s |
| PrimitiveWetModel | 171 | 24 | 256 | 449 | 359 | 10 | 3% | 514.96 MB | 61.7 s |
| PrimitiveWetModel | 256 | 24 | 384 | 300 | 120 | 10 | 3% | 1.36 GB | 64.5 s |
| PrimitiveWetModel | 342 | 24 | 512 | 225 | 41 | 10 | 1% | 2.80 GB | 69.6 s |
| PrimitiveWetModel | 512 | 24 | 768 | 150 | 12 | 10 | 1% | 8.01 GB | 82.6 s |

