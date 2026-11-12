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
julia --project=. manual_benchmarking.jl                       # CPU (auto-labelled cpu-arm or cpu-x86)
julia --project=. manual_benchmarking.jl gpu                   # CUDA GPU
julia --project=. manual_benchmarking.jl amdgpu                # AMDGPU (HIP graphs forced on)
julia --project=. manual_benchmarking.jl metal                 # Apple GPU (Metal)
julia --project=. manual_benchmarking.jl reactant-cpu          # Reactant on CPU
julia --project=. manual_benchmarking.jl reactant-nvidia-gpu   # Reactant on NVIDIA GPU (CUDA)
julia --project=. manual_benchmarking.jl reactant-amd-gpu      # Reactant on AMD GPU (ROCm)
julia --project=. manual_benchmarking.jl reactant-tpu          # Reactant on TPU
julia --project=. manual_benchmarking.jl --debug               # quick regression check, not stored: PrimitiveWet LT+FFT at L8: T32, T43, T64, T86, T128; L16: T32, T128
julia regression/regression.jl check                           # debug mode on main, latest release and latest benchmarked revision
```

Options and flags, in any order (`manual_benchmarking.jl --help` lists them): `--length X`, e.g. `manual_benchmarking.jl gpu --length 10`, multiplies the number of time steps of the timed runs of the PrimitiveWet resolution sweep (the overview below) for more robust timings; the other suites are not affected. `--repeats N` times every configuration exactly N times instead of up to 5 times within 30 s. `--suites benchmark201,benchmark900` restricts the run to specific suites (default: all); a partial run replaces that architecture's stored results with just the suites it ran. `--reverse` (Reactant only, e.g. `manual_benchmarking.jl reactant-nvidia-gpu --reverse`) additionally runs the reverse-mode benchmarks, which take hours. `--hi-res` (GPU and TPU only) adds T342 and T512 to the PrimitiveWet resolution sweep. `--SPEEDY` additionally runs the PrimitiveWet model over the resolutions of the resolution sweep in a configuration closer to SPEEDY.f90 (seasonal prescribed climatologies of sea surface temperature, land surface temperature and soil moisture, constant surface roughness, wind slowdown and gusts).

Each run replaces its own architecture's section in this `README.md` (including suites added by flags in an earlier run); results for other architectures are preserved via `benchmark_results.json`.

A pull request labelled `benchmark-performance` is benchmarked against main with the same debug mode on every push, and the comparison is posted as a comment.

## Overview: PrimitiveWet resolution across architectures

Simulated years per wallclock day (SYPD) for the `PrimitiveWetModel` resolution sweep, one column per architecture. Each (T, L) configuration is reported for both the standard Legendre transform and fast Fourier transform (LT+FFT) and the single matrix transform (MT). Empty cells mean the architecture has not yet been benchmarked or that suite was skipped. Comparison figures across architectures are available on the documentation's `Benchmarks` page.

| T | L | Transform | cpu-arm | cpu-x86 | gpu-nvidia | gpu-amd | reactant-cpu | reactant-nvidia-gpu |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 32 | 8 | LT+FFT | 1308 | 796 | 8349 | 3763 | 1210 | 33951 |
| 32 | 8 | MT | 772 | 841 | 8249 | 3696 | 1203 | 21090 |
| 43 | 8 | LT+FFT | 559 | 343 | 6088 | 2149 | 584 | 20532 |
| 43 | 8 | MT | 253 | 306 | 6025 | 2830 | 464 | 17556 |
| 64 | 8 | LT+FFT | 152 | 100 | 3354 | 828 | 195 | 11068 |
| 64 | 8 | MT | 52 | 65 | 3355 | 826 | 120 | 8023 |
| 86 | 8 | LT+FFT | 58 | 40 | 1269 | 246 | 87 | 6583 |
| 86 | 8 | MT | 14 | 21 | 1283 | 245 | 39 | 2975 |
| 86 | 16 | LT+FFT | 29 | 21 | 799 | 202 | 52 | 3415 |
| 86 | 16 | MT | 8.8 | 12 | 796 | 150 | 25 | 1882 |
| 86 | 24 | LT+FFT | 20 | 14 | 700 | 196 | 27 | 2440 |
| 86 | 24 | MT | 6.1 | 9.0 | 699 | 108 | 16 | 1326 |
| 128 | 8 | LT+FFT | 14 | 11 | 423 | 86 | 29 | 2383 |
| 128 | 8 | MT | 1.2 | — | 239 | 39 | — | 692 |
| 128 | 16 | LT+FFT | 7.6 | 5.5 | 363 | 80 | 15 | 1339 |
| 128 | 16 | MT | 1.2 | — | 145 | 23 | — | 577 |
| 128 | 24 | LT+FFT | 5.5 | 3.6 | 321 | 76 | 7.7 | 885 |
| 128 | 24 | MT | 1.0 | — | 128 | 17 | — | 493 |
| 171 | 8 | LT+FFT | 5.7 | 4.1 | 215 | 44 | 10 | 1178 |
| 171 | 8 | MT | — | — | 52 | — | — | 168 |
| 171 | 16 | LT+FFT | 3.1 | 2.1 | 179 | 41 | 5.5 | 645 |
| 171 | 16 | MT | — | — | 36 | — | — | 97 |
| 171 | 24 | LT+FFT | 1.9 | 1.4 | 155 | 38 | 2.7 | 372 |
| 171 | 24 | MT | — | — | 31 | — | — | 78 |
| 256 | 8 | LT+FFT | 1.4 | 1.0 | 78 | 16 | 2.6 | 374 |
| 256 | 16 | LT+FFT | 0.8 | 0.5 | 61 | 13 | 1.5 | 185 |
| 256 | 24 | LT+FFT | 0.5 | 0.3 | 50 | 11 | 0.8 | 121 |
| 342 | 8 | LT+FFT | — | — | 36 | — | — | 130 |
| 342 | 16 | LT+FFT | — | — | 26 | — | — | 71 |
| 342 | 24 | LT+FFT | — | — | 21 | — | — | 42 |
| 512 | 8 | LT+FFT | — | — | 11 | — | — | 32 |
| 512 | 16 | LT+FFT | — | — | 7.6 | — | — | 17 |
| 512 | 24 | LT+FFT | — | — | 5.7 | — | — | 12 |

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

Created for SpeedyWeather.jl v0.24.0-DEV on Mon, 09 Nov 2026 12:50:51.

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
| BarotropicModel | 32 | 1 | false | 1800 | 25968 | 10 | 1% | 725.97 KB |
| ShallowWaterModel | 32 | 1 | false | 2400 | 18232 | 10 | 6% | 908.33 KB |
| PrimitiveDryModel | 32 | 8 | true | 2400 | 1245 | 10 | 1% | 4.97 MB |
| PrimitiveWetModel | 32 | 8 | true | 2400 | 803 | 10 | 1% | 6.08 MB |

### Shallow water model, resolution

| Model | truncation | L | Rings | Δt | SYPD | Runs | Spread | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 32 | 1 | 48 | 2400 | 18634 | 10 | 1% | 908.33 KB |
| ShallowWaterModel | 43 | 1 | 64 | 1800 | 7847 | 10 | 4% | 1.59 MB |
| ShallowWaterModel | 64 | 1 | 96 | 1200 | 2095 | 10 | 2% | 3.57 MB |
| ShallowWaterModel | 86 | 1 | 128 | 900 | 866 | 10 | 0% | 6.49 MB |
| ShallowWaterModel | 128 | 1 | 192 | 600 | 225 | 10 | 1% | 15.36 MB |
| ShallowWaterModel | 171 | 1 | 256 | 450 | 83 | 10 | 1% | 29.00 MB |
| ShallowWaterModel | 256 | 1 | 384 | 300 | 17 | 10 | 2% | 73.09 MB |

### Primitive wet model, resolution

| Model | truncation | L | Rings | Transform | Δt | SYPD | Runs | Spread | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | 48 | default | 2400 | 796 | 10 | 0% | 6.08 MB |
| PrimitiveWetModel | 43 | 8 | 64 | default | 1800 | 343 | 10 | 0% | 10.24 MB |
| PrimitiveWetModel | 64 | 8 | 96 | default | 1200 | 100 | 10 | 0% | 21.74 MB |
| PrimitiveWetModel | 86 | 8 | 128 | default | 900 | 40 | 10 | 1% | 37.77 MB |
| PrimitiveWetModel | 128 | 8 | 192 | default | 600 | 11 | 10 | 0% | 83.07 MB |
| PrimitiveWetModel | 171 | 8 | 256 | default | 450 | 4.1 | 10 | 0% | 147.25 MB |
| PrimitiveWetModel | 256 | 8 | 384 | default | 300 | 1.0 | 10 | 1% | 334.04 MB |
| PrimitiveWetModel | 86 | 16 | 128 | default | 900 | 21 | 10 | 1% | 63.11 MB |
| PrimitiveWetModel | 128 | 16 | 192 | default | 600 | 5.5 | 10 | 0% | 138.00 MB |
| PrimitiveWetModel | 171 | 16 | 256 | default | 450 | 2.1 | 10 | 1% | 243.28 MB |
| PrimitiveWetModel | 256 | 16 | 384 | default | 300 | 0.5 | 10 | 3% | 546.12 MB |
| PrimitiveWetModel | 86 | 24 | 128 | default | 900 | 14 | 10 | 1% | 88.49 MB |
| PrimitiveWetModel | 128 | 24 | 192 | default | 600 | 3.6 | 10 | 2% | 192.99 MB |
| PrimitiveWetModel | 171 | 24 | 256 | default | 450 | 1.4 | 10 | 1% | 339.40 MB |
| PrimitiveWetModel | 256 | 24 | 384 | default | 300 | 0.3 | 10 | 1% | 758.33 MB |
| PrimitiveWetModel | 32 | 8 | 48 | matrix | 2400 | 841 | 10 | 2% | 34.52 MB |
| PrimitiveWetModel | 43 | 8 | 64 | matrix | 1800 | 306 | 10 | 2% | 93.38 MB |
| PrimitiveWetModel | 64 | 8 | 96 | matrix | 1200 | 65 | 10 | 4% | 397.23 MB |
| PrimitiveWetModel | 86 | 8 | 128 | matrix | 900 | 21 | 10 | 4% | 1.18 GB |
| PrimitiveWetModel | 86 | 16 | 128 | matrix | 900 | 12 | 10 | 4% | 1.21 GB |
| PrimitiveWetModel | 86 | 24 | 128 | matrix | 900 | 9.0 | 10 | 5% | 1.23 GB |

### Primitive Equation, Float32 vs Float64

| Model | NF | truncation | L | Δt | SYPD | Runs | Spread | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | Float32 | 32 | 8 | 2400 | 802 | 10 | 1% | 6.08 MB |
| PrimitiveWetModel | Float64 | 32 | 8 | 2400 | 668 | 10 | 1% | 11.62 MB |

### Grids

| Model | truncation | L | Grid | Rings | Δt | SYPD | Runs | Spread | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 64 | 8 | FullGaussianGrid | 96 | 1200 | 67 | 10 | 0% | 32.24 MB |
| PrimitiveWetModel | 64 | 8 | FullClenshawGrid | 127 | 1200 | 40 | 10 | 1% | 54.20 MB |
| PrimitiveWetModel | 64 | 8 | OctahedralGaussianGrid | 96 | 1200 | 102 | 10 | 0% | 21.74 MB |
| PrimitiveWetModel | 64 | 8 | OctahedralClenshawGrid | 127 | 1200 | 56 | 10 | 2% | 34.42 MB |
| PrimitiveWetModel | 64 | 8 | HEALPixGrid | 127 | 1200 | 450 | 10 | 1% | 25.20 MB |
| PrimitiveWetModel | 64 | 8 | OctaHEALPixGrid | 127 | 1200 | 65 | 10 | 1% | 31.49 MB |

### Number of vertical layers

| Model | truncation | L | Δt | SYPD | Runs | Spread | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 4 | 2400 | 1162 | 10 | 2% | 4.01 MB |
| PrimitiveWetModel | 32 | 8 | 2400 | 801 | 10 | 0% | 6.08 MB |
| PrimitiveWetModel | 32 | 12 | 2400 | 580 | 10 | 1% | 8.15 MB |
| PrimitiveWetModel | 32 | 16 | 2400 | 457 | 10 | 1% | 10.23 MB |

### PrimitiveDryModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Runs | Spread | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveDryModel | 32 | 8 | true | true | 2400 | 1232 | 10 | 1% | 4.97 MB |
| PrimitiveDryModel | 32 | 8 | true | false | 2400 | 2236 | 10 | 2% | 4.97 MB |
| PrimitiveDryModel | 32 | 8 | false | true | 2400 | 1410 | 10 | 2% | 4.97 MB |

### PrimitiveWetModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Runs | Spread | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | true | true | 2400 | 801 | 10 | 0% | 6.08 MB |
| PrimitiveWetModel | 32 | 8 | true | false | 2400 | 1719 | 10 | 1% | 6.08 MB |
| PrimitiveWetModel | 32 | 8 | false | true | 2400 | 957 | 10 | 1% | 6.08 MB |

### Individual dynamics functions


#### PrimitiveWetModel | Float32 | T32 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| --- | --- | --- | --- |
| pressure_gradient_flux! | 0.000 ns| 0 bytes| 0 |
| linear_virtual_temperature! | 0.000 ns| 0 bytes| 0 |
| geopotential! | 0.000 ns| 0 bytes| 0 |
| vertical_integration! | 0.000 ns| 0 bytes| 0 |
| surface_pressure_tendency! | 0.000 ns| 0 bytes| 0 |
| vertical_velocity! | 0.000 ns| 0 bytes| 0 |
| linear_pressure_gradient! | 0.000 ns| 0 bytes| 0 |
| vertical_advection! | 0.000 ns| 0 bytes| 0 |
| vordiv_tendencies! | 0.000 ns| 0 bytes| 0 |
| temperature_tendency! | 0.000 ns| 0 bytes| 0 |
| humidity_tendency! | 0.000 ns| 0 bytes| 0 |
| bernoulli_potential! | 0.000 ns| 0 bytes| 0 |

## Architecture: `gpu-nvidia`

Created for SpeedyWeather.jl v0.24.0-DEV on Mon, 09 Nov 2026 09:59:05.

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
- runtime 13.4.0, artifact installation
- driver 580.126.9 for CUDA 13.5, forward-compatible
- compiler 13.5.37, artifact installation

CUDA libraries: 
- cuBLAS: 13.8.0
- cuSPARSE: 12.8.6
- cuSOLVER: 12.3.5
- cuFFT: 12.4.0
- cuRAND: 10.4.4
- CUPTI: 2026.3.1 (API 13.4.1)
- NVML: 13.0.0+580.126.9

Julia packages: 
- CUDACore: 6.4.1
- GPUArrays: 11.5.16
- GPUCompiler: 2.10.0
- KernelAbstractions: 0.9.45
- CUDA_Driver_jll: 13.5.0+0
- CUDA_Compiler_jll: 0.6.4+0
- CUDA_Runtime_jll: 0.25.0+3
- NVPTX_LLVM_Backend_jll: 23.1.2+2

Toolchain:
- Julia: 1.12.2
- LLVM: 18.1.7

1 device:
  0: NVIDIA H100 80GB HBM3 (sm_90, 75.628 GiB / 79.647 GiB available)
     compiles to sm_90a / PTX 9.4 (LLVM: sm_90a / PTX 9.2)
```


Settings: 10 timed runs per configuration (`--repeats 10`), the PrimitiveWet resolution sweep 2× longer (`--length 2`).

### Models, default setups

| Model | truncation | L | Physics | Δt | SYPD | Runs | Spread | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| BarotropicModel | 32 | 1 | false | 1800 | 21657 | 10 | 1% | 217.07 KB |
| ShallowWaterModel | 32 | 1 | false | 2400 | 16298 | 10 | 4% | 221.44 KB |
| PrimitiveDryModel | 32 | 8 | true | 2400 | 10705 | 10 | 6% | 375.98 KB |
| PrimitiveWetModel | 32 | 8 | true | 2400 | 7860 | 10 | 3% | 385.97 KB |

### Shallow water model, resolution

| Model | truncation | L | Rings | Δt | SYPD | Runs | Spread | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 32 | 1 | 48 | 2400 | 15995 | 10 | 1% | 221.44 KB |
| ShallowWaterModel | 43 | 1 | 64 | 1800 | 11880 | 10 | 2% | 365.32 KB |
| ShallowWaterModel | 64 | 1 | 96 | 1200 | 4591 | 10 | 0% | 775.94 KB |
| ShallowWaterModel | 86 | 1 | 128 | 900 | 1640 | 10 | 0% | 1.35 MB |
| ShallowWaterModel | 128 | 1 | 192 | 600 | 298 | 10 | 0% | 3.44 MB |
| ShallowWaterModel | 171 | 1 | 256 | 450 | 156 | 10 | 0% | 5.89 MB |
| ShallowWaterModel | 256 | 1 | 384 | 300 | 63 | 10 | 0% | 12.75 MB |

### Primitive wet model, resolution

| Model | truncation | L | Rings | Transform | Δt | SYPD | Runs | Spread | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | 48 | default | 2400 | 8349 | 10 | 2% | 385.97 KB |
| PrimitiveWetModel | 43 | 8 | 64 | default | 1800 | 6088 | 10 | 2% | 630.19 KB |
| PrimitiveWetModel | 64 | 8 | 96 | default | 1200 | 3354 | 10 | 0% | 1.33 MB |
| PrimitiveWetModel | 86 | 8 | 128 | default | 900 | 1269 | 10 | 0% | 2.30 MB |
| PrimitiveWetModel | 128 | 8 | 192 | default | 600 | 423 | 10 | 0% | 5.61 MB |
| PrimitiveWetModel | 171 | 8 | 256 | default | 450 | 215 | 10 | 0% | 9.68 MB |
| PrimitiveWetModel | 256 | 8 | 384 | default | 300 | 78 | 10 | 0% | 21.17 MB |
| PrimitiveWetModel | 342 | 8 | 512 | default | 225 | 36 | 10 | 0% | 37.11 MB |
| PrimitiveWetModel | 512 | 8 | 768 | default | 150 | 11 | 10 | 0% | 82.37 MB |
| PrimitiveWetModel | 86 | 16 | 128 | default | 900 | 799 | 10 | 1% | 3.35 MB |
| PrimitiveWetModel | 128 | 16 | 192 | default | 600 | 363 | 10 | 0% | 7.97 MB |
| PrimitiveWetModel | 171 | 16 | 256 | default | 450 | 179 | 10 | 0% | 13.87 MB |
| PrimitiveWetModel | 256 | 16 | 384 | default | 300 | 61 | 10 | 0% | 30.60 MB |
| PrimitiveWetModel | 342 | 16 | 512 | default | 225 | 26 | 10 | 0% | 53.89 MB |
| PrimitiveWetModel | 512 | 16 | 768 | default | 150 | 7.6 | 10 | 0% | 120.12 MB |
| PrimitiveWetModel | 86 | 24 | 128 | default | 900 | 700 | 10 | 0% | 4.40 MB |
| PrimitiveWetModel | 128 | 24 | 192 | default | 600 | 321 | 10 | 0% | 10.33 MB |
| PrimitiveWetModel | 171 | 24 | 256 | default | 450 | 155 | 10 | 0% | 18.07 MB |
| PrimitiveWetModel | 256 | 24 | 384 | default | 300 | 50 | 10 | 0% | 40.04 MB |
| PrimitiveWetModel | 342 | 24 | 512 | default | 225 | 21 | 10 | 0% | 70.67 MB |
| PrimitiveWetModel | 512 | 24 | 768 | default | 150 | 5.7 | 10 | 0% | 157.87 MB |
| PrimitiveWetModel | 32 | 8 | 48 | matrix | 2400 | 8249 | 10 | 1% | 385.97 KB |
| PrimitiveWetModel | 43 | 8 | 64 | matrix | 1800 | 6025 | 10 | 1% | 630.19 KB |
| PrimitiveWetModel | 64 | 8 | 96 | matrix | 1200 | 3355 | 10 | 0% | 1.33 MB |
| PrimitiveWetModel | 86 | 8 | 128 | matrix | 900 | 1283 | 10 | 1% | 2.30 MB |
| PrimitiveWetModel | 128 | 8 | 192 | matrix | 600 | 239 | 10 | 0% | 5.09 MB |
| PrimitiveWetModel | 171 | 8 | 256 | matrix | 450 | 52 | 10 | 0% | 8.99 MB |
| PrimitiveWetModel | 86 | 16 | 128 | matrix | 900 | 796 | 10 | 0% | 3.35 MB |
| PrimitiveWetModel | 128 | 16 | 192 | matrix | 600 | 145 | 10 | 0% | 7.45 MB |
| PrimitiveWetModel | 171 | 16 | 256 | matrix | 450 | 36 | 10 | 0% | 13.19 MB |
| PrimitiveWetModel | 86 | 24 | 128 | matrix | 900 | 699 | 10 | 0% | 4.40 MB |
| PrimitiveWetModel | 128 | 24 | 192 | matrix | 600 | 128 | 10 | 0% | 9.81 MB |
| PrimitiveWetModel | 171 | 24 | 256 | matrix | 450 | 31 | 10 | 0% | 17.38 MB |

### Primitive Equation, Float32 vs Float64

| Model | NF | truncation | L | Δt | SYPD | Runs | Spread | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | Float32 | 32 | 8 | 2400 | 8588 | 10 | 5% | 385.97 KB |
| PrimitiveWetModel | Float64 | 32 | 8 | 2400 | 8207 | 10 | 1% | 387.15 KB |

### Grids

| Model | truncation | L | Grid | Rings | Δt | SYPD | Runs | Spread | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 64 | 8 | FullGaussianGrid | 96 | 1200 | 2508 | 10 | 1% | 1.33 MB |
| PrimitiveWetModel | 64 | 8 | FullClenshawGrid | 127 | 1200 | 1703 | 10 | 1% | 2.29 MB |
| PrimitiveWetModel | 64 | 8 | OctahedralGaussianGrid | 96 | 1200 | 3315 | 10 | 1% | 1.33 MB |
| PrimitiveWetModel | 64 | 8 | OctahedralClenshawGrid | 127 | 1200 | 2528 | 10 | 1% | 2.29 MB |
| PrimitiveWetModel | 64 | 8 | HEALPixGrid | 127 | 1200 | 13283 | 10 | 16% | 2.29 MB |
| PrimitiveWetModel | 64 | 8 | OctaHEALPixGrid | 127 | 1200 | 2752 | 10 | 0% | 2.29 MB |

### Number of vertical layers

| Model | truncation | L | Δt | SYPD | Runs | Spread | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 4 | 2400 | 8527 | 10 | 2% | 312.24 KB |
| PrimitiveWetModel | 32 | 8 | 2400 | 8601 | 10 | 1% | 385.97 KB |
| PrimitiveWetModel | 32 | 12 | 2400 | 8371 | 10 | 6% | 459.70 KB |
| PrimitiveWetModel | 32 | 16 | 2400 | 8556 | 10 | 1% | 533.42 KB |

### PrimitiveDryModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Runs | Spread | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveDryModel | 32 | 8 | true | true | 2400 | 10838 | 10 | 1% | 375.98 KB |
| PrimitiveDryModel | 32 | 8 | true | false | 2400 | 13811 | 10 | 0% | 375.98 KB |
| PrimitiveDryModel | 32 | 8 | false | true | 2400 | 15410 | 10 | 0% | 375.98 KB |

### PrimitiveWetModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Runs | Spread | Memory|
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | true | true | 2400 | 8611 | 10 | 1% | 385.97 KB |
| PrimitiveWetModel | 32 | 8 | true | false | 2400 | 12417 | 10 | 1% | 385.97 KB |
| PrimitiveWetModel | 32 | 8 | false | true | 2400 | 11534 | 10 | 1% | 385.97 KB |

### Individual dynamics functions


#### PrimitiveWetModel | Float32 | T32 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| --- | --- | --- | --- |
| pressure_gradient_flux! | 0.000 ns| 0 bytes| 0 |
| linear_virtual_temperature! | 0.000 ns| 0 bytes| 0 |
| geopotential! | 0.000 ns| 0 bytes| 0 |
| vertical_integration! | 0.000 ns| 0 bytes| 0 |
| surface_pressure_tendency! | 0.000 ns| 0 bytes| 0 |
| vertical_velocity! | 0.000 ns| 0 bytes| 0 |
| linear_pressure_gradient! | 0.000 ns| 0 bytes| 0 |
| vertical_advection! | 0.000 ns| 0 bytes| 0 |
| vordiv_tendencies! | 0.000 ns| 0 bytes| 0 |
| temperature_tendency! | 0.000 ns| 0 bytes| 0 |
| humidity_tendency! | 0.000 ns| 0 bytes| 0 |
| bernoulli_potential! | 0.000 ns| 0 bytes| 0 |

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

Created for SpeedyWeather.jl v0.24.0-DEV on Mon, 09 Nov 2026 16:38:37.

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

Reactant backend: `cpu` (via `Reactant.set_default_backend("cpu")`), 1 device(s): cpu (one used), Reactant v0.2.307, Reactant_jll v0.0.442+0.


Settings: 10 timed runs per configuration (`--repeats 10`), the PrimitiveWet resolution sweep 2× longer (`--length 2`).

### Models, default setups

| Model | truncation | L | Physics | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| BarotropicModel | 32 | 1 | false | 1800 | 12873 | 10 | 1% | 1.81 MB | 28.0 s |
| ShallowWaterModel | 32 | 1 | false | 2400 | 11636 | 10 | 5% | 2.00 MB | 38.1 s |
| PrimitiveDryModel | 32 | 8 | true | 2400 | 1776 | 10 | 5% | 5.77 MB | 176.1 s |
| PrimitiveWetModel | 32 | 8 | true | 2400 | 1201 | 10 | 4% | 6.88 MB | 213.1 s |

### Shallow water model, resolution

| Model | truncation | L | Rings | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 32 | 1 | 48 | 2400 | 11771 | 10 | 6% | 2.00 MB | 1.8 s |
| ShallowWaterModel | 43 | 1 | 64 | 1786 | 4413 | 10 | 3% | 4.28 MB | 8.3 s |
| ShallowWaterModel | 64 | 1 | 96 | 1200 | 1173 | 10 | 6% | 12.85 MB | 6.9 s |
| ShallowWaterModel | 86 | 1 | 128 | 893 | 529 | 10 | 5% | 29.04 MB | 7.3 s |
| ShallowWaterModel | 128 | 1 | 192 | 600 | 105 | 10 | 1% | 91.85 MB | 9.2 s |
| ShallowWaterModel | 171 | 1 | 256 | 449 | 60 | 10 | 2% | 211.92 MB | 9.9 s |
| ShallowWaterModel | 256 | 1 | 384 | 300 | 9.7 | 10 | 1% | 693.70 MB | 11.0 s |

### Primitive wet model, resolution

| Model | truncation | L | Rings | Transform | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | 48 | default | 2400 | 1210 | 10 | 3% | 6.88 MB | 10.5 s |
| PrimitiveWetModel | 43 | 8 | 64 | default | 1786 | 584 | 10 | 2% | 12.43 MB | 33.5 s |
| PrimitiveWetModel | 64 | 8 | 96 | default | 1200 | 195 | 10 | 7% | 29.94 MB | 27.5 s |
| PrimitiveWetModel | 86 | 8 | 128 | default | 893 | 87 | 10 | 2% | 58.40 MB | 28.2 s |
| PrimitiveWetModel | 128 | 8 | 192 | default | 600 | 29 | 10 | 1% | 155.32 MB | 27.6 s |
| PrimitiveWetModel | 171 | 8 | 256 | default | 449 | 10 | 10 | 3% | 322.67 MB | 28.7 s |
| PrimitiveWetModel | 256 | 8 | 384 | default | 300 | 2.6 | 10 | 6% | 937.89 MB | 30.6 s |
| PrimitiveWetModel | 86 | 16 | 128 | default | 893 | 52 | 10 | 4% | 81.57 MB | 187.2 s |
| PrimitiveWetModel | 128 | 16 | 192 | default | 600 | 15 | 10 | 6% | 205.42 MB | 32.2 s |
| PrimitiveWetModel | 171 | 16 | 256 | default | 449 | 5.5 | 10 | 2% | 410.16 MB | 33.1 s |
| PrimitiveWetModel | 256 | 16 | 384 | default | 300 | 1.5 | 10 | 11% | 1.13 GB | 35.7 s |
| PrimitiveWetModel | 86 | 24 | 128 | default | 893 | 27 | 10 | 4% | 104.78 MB | 255.8 s |
| PrimitiveWetModel | 128 | 24 | 192 | default | 600 | 7.7 | 10 | 5% | 255.58 MB | 61.3 s |
| PrimitiveWetModel | 171 | 24 | 256 | default | 449 | 2.7 | 10 | 3% | 497.74 MB | 62.8 s |
| PrimitiveWetModel | 256 | 24 | 384 | default | 300 | 0.8 | 10 | 3% | 1.32 GB | 66.0 s |
| PrimitiveWetModel | 32 | 8 | 48 | matrix | 2400 | 1203 | 10 | 2% | 34.08 MB | 195.1 s |
| PrimitiveWetModel | 43 | 8 | 64 | matrix | 1786 | 464 | 10 | 1% | 92.59 MB | 12.4 s |
| PrimitiveWetModel | 64 | 8 | 96 | matrix | 1200 | 120 | 10 | 1% | 395.47 MB | 12.9 s |
| PrimitiveWetModel | 86 | 8 | 128 | matrix | 893 | 39 | 10 | 2% | 1.18 GB | 13.9 s |
| PrimitiveWetModel | 86 | 16 | 128 | matrix | 893 | 25 | 10 | 10% | 1.20 GB | 176.0 s |
| PrimitiveWetModel | 86 | 24 | 128 | matrix | 893 | 16 | 10 | 2% | 1.23 GB | 204.1 s |

### Primitive Equation, Float32 vs Float64

| Model | NF | truncation | L | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | Float32 | 32 | 8 | 2400 | 1188 | 10 | 3% | 6.88 MB | 8.6 s |
| PrimitiveWetModel | Float64 | 32 | 8 | 2400 | 823 | 10 | 2% | 13.61 MB | 279.1 s |

### Grids

| Model | truncation | L | Grid | Rings | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 64 | 8 | FullGaussianGrid | 96 | 1200 | 110 | 10 | 2% | 39.66 MB | 218.3 s |
| PrimitiveWetModel | 64 | 8 | FullClenshawGrid | 127 | 1200 | 92 | 10 | 3% | 67.05 MB | 238.2 s |
| PrimitiveWetModel | 64 | 8 | OctahedralGaussianGrid | 96 | 1200 | 193 | 10 | 2% | 29.94 MB | 9.4 s |
| PrimitiveWetModel | 64 | 8 | OctahedralClenshawGrid | 127 | 1200 | 127 | 10 | 1% | 48.31 MB | 217.5 s |
| PrimitiveWetModel | 64 | 8 | HEALPixGrid | 127 | 1200 | 126 | 10 | 10% | 30.27 MB | 215.4 s |
| PrimitiveWetModel | 64 | 8 | OctaHEALPixGrid | 127 | 1200 | 68 | 10 | 3% | 44.41 MB | 236.2 s |

### Number of vertical layers

| Model | truncation | L | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 4 | 2400 | 2440 | 10 | 6% | 4.97 MB | 171.8 s |
| PrimitiveWetModel | 32 | 8 | 2400 | 1312 | 10 | 10% | 6.88 MB | 8.4 s |
| PrimitiveWetModel | 32 | 12 | 2400 | 1015 | 10 | 7% | 8.79 MB | 177.6 s |
| PrimitiveWetModel | 32 | 16 | 2400 | 700 | 10 | 4% | 10.71 MB | 36.7 s |

### PrimitiveDryModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveDryModel | 32 | 8 | true | true | 2400 | 1851 | 10 | 6% | 5.77 MB | 6.1 s |
| PrimitiveDryModel | 32 | 8 | true | false | 2400 | 2142 | 10 | 3% | 5.77 MB | 5.5 s |
| PrimitiveDryModel | 32 | 8 | false | true | 2400 | 2480 | 10 | 7% | 5.77 MB | 7.9 s |

### PrimitiveWetModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | true | true | 2400 | 1218 | 10 | 2% | 6.88 MB | 11.0 s |
| PrimitiveWetModel | 32 | 8 | true | false | 2400 | 1734 | 10 | 4% | 6.88 MB | 6.1 s |
| PrimitiveWetModel | 32 | 8 | false | true | 2400 | 1553 | 10 | 3% | 6.88 MB | 10.6 s |

### Individual dynamics functions


#### PrimitiveWetModel | Float32 | T32 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| --- | --- | --- | --- |
| pressure_gradient_flux! | 0.000 ns| 0 bytes| 0 |
| linear_virtual_temperature! | 0.000 ns| 0 bytes| 0 |
| geopotential! | 0.000 ns| 0 bytes| 0 |
| vertical_integration! | 0.000 ns| 0 bytes| 0 |
| surface_pressure_tendency! | 0.000 ns| 0 bytes| 0 |
| vertical_velocity! | 0.000 ns| 0 bytes| 0 |
| linear_pressure_gradient! | 0.000 ns| 0 bytes| 0 |
| vertical_advection! | 0.000 ns| 0 bytes| 0 |
| vordiv_tendencies! | 0.000 ns| 0 bytes| 0 |
| temperature_tendency! | 0.000 ns| 0 bytes| 0 |
| humidity_tendency! | 0.000 ns| 0 bytes| 0 |
| bernoulli_potential! | 0.000 ns| 0 bytes| 0 |

### Primitive wet model, reverse mode (AD) with Reactant

Wallclock time of the time steps forward and of the reverse-mode gradient with respect to the initial state over the same time steps, minimum over up to 3 runs, compilation excluded. The loss is the squared norm of the vorticity after the time steps.

| truncation | L | Rings | Δt | Steps | Checkpointing | Forward | Gradient | Gradient/Forward | Compile forward | Compile gradient | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 32 | 8 | 48 | 2400 | 50 | Binomial(32) | 0.278 s | 1.24 s | 4.5 | 34.1 s | 297.1 s | ok |
| 32 | 8 | 48 | 2400 | 100 | Binomial(32) | 0.557 s | 2.69 s | 4.8 | 11.8 s | 46.1 s | ok |
| 32 | 8 | 48 | 2400 | 200 | Binomial(32) | 1.04 s | 5.62 s | 5.4 | 9.5 s | 33.8 s | ok |
| 32 | 8 | 48 | 2400 | 2000 | Binomial(32) | 11.1 s | 68.2 s | 6.1 | 9.6 s | 33.7 s | non-finite or zero gradient |
| 64 | 8 | 96 | 1200 | 50 | Binomial(32) | 0.878 s | 3.81 s | 4.3 | 13.6 s | 46.8 s | ok |
| 64 | 8 | 96 | 1200 | 100 | Binomial(32) | 1.74 s | 8.28 s | 4.8 | 14.0 s | 48.2 s | ok |
| 64 | 8 | 96 | 1200 | 200 | Binomial(32) | 3.48 s | 17 s | 4.9 | 13.8 s | 48.8 s | ok |
| 64 | 8 | 96 | 1200 | 2000 | Binomial(32) | 34.8 s | 206 s | 5.9 | 10.7 s | 35.8 s | non-finite or zero gradient |
| 128 | 8 | 192 | 600 | 50 | Binomial(32) | 2.94 s | 13.8 s | 4.7 | 11.6 s | 35.3 s | ok |
| 128 | 8 | 192 | 600 | 100 | Binomial(32) | 5.83 s | 29.5 s | 5.1 | 15.1 s | 47.5 s | ok |
| 128 | 8 | 192 | 600 | 200 | Binomial(32) | 11.6 s | 61.8 s | 5.3 | 15.2 s | 47.9 s | ok |
| 128 | 8 | 192 | 600 | 2000 | Binomial(32) | 117 s | 744 s | 6.4 | 15.2 s | 49.0 s | non-finite or zero gradient |

### SPEEDY.f90-like configuration

The PrimitiveWet model in a configuration closer to SPEEDY.f90, at the resolutions of the PrimitiveWet resolution sweep with the default transform: seasonal prescribed climatologies of sea surface temperature (`SeasonalOceanClimatology`), land surface temperature (`SeasonalLandTemperature`) and soil moisture (`SeasonalSoilMoisture`) instead of the slab ocean and the land bucket model, constant surface roughness lengths over land and ocean (`ConstantSurfaceRoughness`) instead of the learned ones, and SPEEDY.f90's near-surface wind of 0.95 times the lowermost-layer wind (`ConstantWindSlowdown`) with 5 m/s gusts (`ConstantWindGust`) instead of the learned wind slowdown and Beljaars gusts; sea ice stays thermodynamic, driven by the prescribed sea surface temperature.

| Model | truncation | L | Rings | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | 48 | 2400 | 1228 | 10 | 4% | 7.40 MB | 217.0 s |
| PrimitiveWetModel | 43 | 8 | 64 | 1786 | 591 | 10 | 3% | 13.29 MB | 17.7 s |
| PrimitiveWetModel | 64 | 8 | 96 | 1200 | 205 | 10 | 4% | 31.73 MB | 14.8 s |
| PrimitiveWetModel | 86 | 8 | 128 | 893 | 86 | 10 | 5% | 61.47 MB | 15.9 s |
| PrimitiveWetModel | 128 | 8 | 192 | 600 | 28 | 10 | 5% | 161.93 MB | 16.4 s |
| PrimitiveWetModel | 171 | 8 | 256 | 449 | 9.7 | 10 | 4% | 334.17 MB | 20.6 s |
| PrimitiveWetModel | 256 | 8 | 384 | 300 | 2.5 | 10 | 4% | 963.21 MB | 19.4 s |
| PrimitiveWetModel | 86 | 16 | 128 | 893 | 49 | 10 | 6% | 84.63 MB | 165.4 s |
| PrimitiveWetModel | 128 | 16 | 192 | 600 | 15 | 10 | 7% | 212.03 MB | 22.7 s |
| PrimitiveWetModel | 171 | 16 | 256 | 449 | 5.4 | 10 | 3% | 421.66 MB | 23.6 s |
| PrimitiveWetModel | 256 | 16 | 384 | 300 | 1.4 | 10 | 4% | 1.16 GB | 26.6 s |
| PrimitiveWetModel | 86 | 24 | 128 | 893 | 26 | 10 | 4% | 107.85 MB | 231.9 s |
| PrimitiveWetModel | 128 | 24 | 192 | 600 | 7.6 | 10 | 3% | 262.19 MB | 70.1 s |
| PrimitiveWetModel | 171 | 24 | 256 | 449 | 2.6 | 10 | 2% | 509.25 MB | 67.9 s |
| PrimitiveWetModel | 256 | 24 | 384 | 300 | 0.8 | 10 | 11% | 1.35 GB | 69.4 s |

## Architecture: `reactant-nvidia-gpu`

Created for SpeedyWeather.jl v0.24.0-DEV on Mon, 09 Nov 2026 14:45:29.

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
- runtime 13.4.0, artifact installation
- driver 580.126.9 for CUDA 13.5, forward-compatible
- compiler 13.5.37, artifact installation

CUDA libraries: 
- cuBLAS: 13.8.0
- cuSPARSE: 12.8.6
- cuSOLVER: 12.3.5
- cuFFT: 12.4.0
- cuRAND: 10.4.4
- CUPTI: 2026.3.1 (API 13.4.1)
- NVML: 13.0.0+580.126.9

Julia packages: 
- CUDACore: 6.4.1
- GPUArrays: 11.5.16
- GPUCompiler: 2.10.0
- KernelAbstractions: 0.9.45
- CUDA_Driver_jll: 13.5.0+0
- CUDA_Compiler_jll: 0.6.4+0
- CUDA_Runtime_jll: 0.25.0+3
- NVPTX_LLVM_Backend_jll: 23.1.2+2

Toolchain:
- Julia: 1.12.2
- LLVM: 18.1.7

1 device:
  0: NVIDIA H100 80GB HBM3 (sm_90, 18.934 GiB / 79.647 GiB available)
     compiles to sm_90a / PTX 9.4 (LLVM: sm_90a / PTX 9.2)
```

Reactant backend: `gpu` (via `Reactant.set_default_backend("gpu")`), 1 device(s): NVIDIA H100 80GB HBM3 (one used), Reactant v0.2.307, Reactant_jll v0.0.442+0.


Settings: 10 timed runs per configuration (`--repeats 10`), the PrimitiveWet resolution sweep 2× longer (`--length 2`).

### Models, default setups

| Model | truncation | L | Physics | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| BarotropicModel | 32 | 1 | false | 1800 | 155202 | 10 | 13% | 1.81 MB | 28.5 s |
| ShallowWaterModel | 32 | 1 | false | 2400 | 148129 | 10 | 1% | 2.00 MB | 33.4 s |
| PrimitiveDryModel | 32 | 8 | true | 2400 | 35862 | 10 | 2% | 5.77 MB | 184.8 s |
| PrimitiveWetModel | 32 | 8 | true | 2400 | 33880 | 10 | 2% | 6.88 MB | 222.2 s |

### Shallow water model, resolution

| Model | truncation | L | Rings | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 32 | 1 | 48 | 2400 | 147160 | 10 | 14% | 2.00 MB | 1.6 s |
| ShallowWaterModel | 43 | 1 | 64 | 1786 | 105532 | 10 | 1% | 4.28 MB | 7.6 s |
| ShallowWaterModel | 64 | 1 | 96 | 1200 | 62888 | 10 | 2% | 12.85 MB | 7.9 s |
| ShallowWaterModel | 86 | 1 | 128 | 893 | 38207 | 10 | 4% | 29.04 MB | 7.8 s |
| ShallowWaterModel | 128 | 1 | 192 | 600 | 14612 | 10 | 1% | 91.85 MB | 8.1 s |
| ShallowWaterModel | 171 | 1 | 256 | 449 | 6884 | 10 | 0% | 211.92 MB | 8.9 s |
| ShallowWaterModel | 256 | 1 | 384 | 300 | 2050 | 10 | 1% | 693.70 MB | 10.6 s |

### Primitive wet model, resolution

| Model | truncation | L | Rings | Transform | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | 48 | default | 2400 | 33951 | 10 | 3% | 6.88 MB | 15.0 s |
| PrimitiveWetModel | 43 | 8 | 64 | default | 1786 | 20532 | 10 | 2% | 12.43 MB | 33.8 s |
| PrimitiveWetModel | 64 | 8 | 96 | default | 1200 | 11068 | 10 | 3% | 29.94 MB | 32.0 s |
| PrimitiveWetModel | 86 | 8 | 128 | default | 893 | 6583 | 10 | 3% | 58.40 MB | 32.7 s |
| PrimitiveWetModel | 128 | 8 | 192 | default | 600 | 2383 | 10 | 3% | 155.32 MB | 33.9 s |
| PrimitiveWetModel | 171 | 8 | 256 | default | 449 | 1178 | 10 | 2% | 322.67 MB | 33.6 s |
| PrimitiveWetModel | 256 | 8 | 384 | default | 300 | 374 | 10 | 1% | 937.89 MB | 34.6 s |
| PrimitiveWetModel | 342 | 8 | 512 | default | 225 | 130 | 10 | 1% | 2.05 GB | 39.2 s |
| PrimitiveWetModel | 512 | 8 | 768 | default | 150 | 32 | 10 | 0% | 6.35 GB | 50.0 s |
| PrimitiveWetModel | 86 | 16 | 128 | default | 893 | 3415 | 10 | 3% | 81.57 MB | 198.8 s |
| PrimitiveWetModel | 128 | 16 | 192 | default | 600 | 1339 | 10 | 3% | 205.42 MB | 41.0 s |
| PrimitiveWetModel | 171 | 16 | 256 | default | 449 | 645 | 10 | 2% | 410.16 MB | 41.7 s |
| PrimitiveWetModel | 256 | 16 | 384 | default | 300 | 185 | 10 | 1% | 1.13 GB | 45.8 s |
| PrimitiveWetModel | 342 | 16 | 512 | default | 225 | 71 | 10 | 0% | 2.39 GB | 49.0 s |
| PrimitiveWetModel | 512 | 16 | 768 | default | 150 | 17 | 10 | 0% | 7.10 GB | 59.3 s |
| PrimitiveWetModel | 86 | 24 | 128 | default | 893 | 2440 | 10 | 9% | 104.78 MB | 229.2 s |
| PrimitiveWetModel | 128 | 24 | 192 | default | 600 | 885 | 10 | 3% | 255.58 MB | 72.5 s |
| PrimitiveWetModel | 171 | 24 | 256 | default | 449 | 372 | 10 | 2% | 497.74 MB | 71.4 s |
| PrimitiveWetModel | 256 | 24 | 384 | default | 300 | 121 | 10 | 1% | 1.32 GB | 75.1 s |
| PrimitiveWetModel | 342 | 24 | 512 | default | 225 | 42 | 10 | 1% | 2.73 GB | 82.3 s |
| PrimitiveWetModel | 512 | 24 | 768 | default | 150 | 12 | 10 | 1% | 7.86 GB | 94.9 s |
| PrimitiveWetModel | 32 | 8 | 48 | matrix | 2400 | 21090 | 10 | 1% | 34.08 MB | 199.5 s |
| PrimitiveWetModel | 43 | 8 | 64 | matrix | 1786 | 17556 | 10 | 1% | 92.59 MB | 17.7 s |
| PrimitiveWetModel | 64 | 8 | 96 | matrix | 1200 | 8023 | 10 | 4% | 395.47 MB | 17.7 s |
| PrimitiveWetModel | 86 | 8 | 128 | matrix | 893 | 2975 | 10 | 1% | 1.18 GB | 18.2 s |
| PrimitiveWetModel | 128 | 8 | 192 | matrix | 600 | 692 | 10 | 5% | 5.48 GB | 17.7 s |
| PrimitiveWetModel | 171 | 8 | 256 | matrix | 449 | 168 | 10 | 3% | 16.83 GB | 19.6 s |
| PrimitiveWetModel | 86 | 16 | 128 | matrix | 893 | 1882 | 10 | 4% | 1.20 GB | 190.3 s |
| PrimitiveWetModel | 128 | 16 | 192 | matrix | 600 | 577 | 10 | 2% | 5.54 GB | 29.8 s |
| PrimitiveWetModel | 171 | 16 | 256 | matrix | 449 | 97 | 10 | 1% | 16.92 GB | 31.7 s |
| PrimitiveWetModel | 86 | 24 | 128 | matrix | 893 | 1326 | 10 | 3% | 1.23 GB | 218.9 s |
| PrimitiveWetModel | 128 | 24 | 192 | matrix | 600 | 493 | 10 | 3% | 5.59 GB | 60.0 s |
| PrimitiveWetModel | 171 | 24 | 256 | matrix | 449 | 78 | 10 | 1% | 17.02 GB | 60.0 s |

### Primitive Equation, Float32 vs Float64

| Model | NF | truncation | L | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | Float32 | 32 | 8 | 2400 | 31808 | 10 | 15% | 6.88 MB | 16.3 s |
| PrimitiveWetModel | Float64 | 32 | 8 | 2400 | 15820 | 10 | 2% | 13.61 MB | 259.9 s |

### Grids

| Model | truncation | L | Grid | Rings | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 64 | 8 | FullGaussianGrid | 96 | 1200 | 9596 | 10 | 3% | 39.66 MB | 221.0 s |
| PrimitiveWetModel | 64 | 8 | FullClenshawGrid | 127 | 1200 | 8176 | 10 | 8% | 67.05 MB | 219.6 s |
| PrimitiveWetModel | 64 | 8 | OctahedralGaussianGrid | 96 | 1200 | 10708 | 10 | 3% | 29.94 MB | 13.4 s |
| PrimitiveWetModel | 64 | 8 | OctahedralClenshawGrid | 127 | 1200 | 9377 | 10 | 4% | 48.31 MB | 224.2 s |
| PrimitiveWetModel | 64 | 8 | HEALPixGrid | 127 | 1200 | 10960 | 10 | 9% | 30.27 MB | 222.8 s |
| PrimitiveWetModel | 64 | 8 | OctaHEALPixGrid | 127 | 1200 | 9890 | 10 | 6% | 44.41 MB | 221.1 s |

### Number of vertical layers

| Model | truncation | L | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 4 | 2400 | 41984 | 10 | 5% | 4.97 MB | 174.7 s |
| PrimitiveWetModel | 32 | 8 | 2400 | 32811 | 10 | 14% | 6.88 MB | 15.5 s |
| PrimitiveWetModel | 32 | 12 | 2400 | 26430 | 10 | 2% | 8.79 MB | 197.3 s |
| PrimitiveWetModel | 32 | 16 | 2400 | 12063 | 10 | 19% | 10.71 MB | 42.3 s |

### PrimitiveDryModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveDryModel | 32 | 8 | true | true | 2400 | 35864 | 10 | 12% | 5.77 MB | 15.8 s |
| PrimitiveDryModel | 32 | 8 | true | false | 2400 | 43463 | 10 | 15% | 5.77 MB | 13.6 s |
| PrimitiveDryModel | 32 | 8 | false | true | 2400 | 61428 | 10 | 19% | 5.77 MB | 10.3 s |

### PrimitiveWetModel: Physics or dynamics only

| Model | truncation | L | Dynamics | Physics | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | true | true | 2400 | 32525 | 10 | 16% | 6.88 MB | 15.9 s |
| PrimitiveWetModel | 32 | 8 | true | false | 2400 | 44148 | 10 | 22% | 6.88 MB | 10.3 s |
| PrimitiveWetModel | 32 | 8 | false | true | 2400 | 47853 | 10 | 23% | 6.88 MB | 13.9 s |

### Individual dynamics functions


#### PrimitiveWetModel | Float32 | T32 L8 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| --- | --- | --- | --- |
| pressure_gradient_flux! | 0.000 ns| 0 bytes| 0 |
| linear_virtual_temperature! | 0.000 ns| 0 bytes| 0 |
| geopotential! | 0.000 ns| 0 bytes| 0 |
| vertical_integration! | 0.000 ns| 0 bytes| 0 |
| surface_pressure_tendency! | 0.000 ns| 0 bytes| 0 |
| vertical_velocity! | 0.000 ns| 0 bytes| 0 |
| linear_pressure_gradient! | 0.000 ns| 0 bytes| 0 |
| vertical_advection! | 0.000 ns| 0 bytes| 0 |
| vordiv_tendencies! | 0.000 ns| 0 bytes| 0 |
| temperature_tendency! | 0.000 ns| 0 bytes| 0 |
| humidity_tendency! | 0.000 ns| 0 bytes| 0 |
| bernoulli_potential! | 0.000 ns| 0 bytes| 0 |

### Primitive wet model, reverse mode (AD) with Reactant

Wallclock time of the time steps forward and of the reverse-mode gradient with respect to the initial state over the same time steps, minimum over up to 3 runs, compilation excluded. The loss is the squared norm of the vorticity after the time steps.

| truncation | L | Rings | Δt | Steps | Checkpointing | Forward | Gradient | Gradient/Forward | Compile forward | Compile gradient | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 32 | 8 | 48 | 2400 | 50 | Binomial(32) | 0.0111 s | 0.0765 s | 6.9 | 40.4 s | 318.6 s | ok |
| 32 | 8 | 48 | 2400 | 100 | Binomial(32) | 0.0206 s | 0.14 s | 6.8 | 15.7 s | 55.0 s | ok |
| 32 | 8 | 48 | 2400 | 200 | Binomial(32) | 0.0378 s | 0.269 s | 7.1 | 15.8 s | 55.5 s | ok |
| 32 | 8 | 48 | 2400 | 2000 | Binomial(32) | 0.362 s | 2.93 s | 8.1 | 15.8 s | 55.4 s | non-finite or zero gradient |
| 64 | 8 | 96 | 1200 | 50 | Binomial(32) | 0.0168 s | 0.0948 s | 5.6 | 15.5 s | 53.8 s | ok |
| 64 | 8 | 96 | 1200 | 100 | Binomial(32) | 0.0308 s | 0.182 s | 5.9 | 15.5 s | 52.6 s | ok |
| 64 | 8 | 96 | 1200 | 200 | Binomial(32) | 0.0591 s | 0.353 s | 6.0 | 15.5 s | 53.1 s | ok |
| 64 | 8 | 96 | 1200 | 2000 | Binomial(32) | 0.574 s | 3.94 s | 6.9 | 15.5 s | 53.3 s | non-finite or zero gradient |
| 128 | 8 | 192 | 600 | 50 | Binomial(32) | 0.0363 s | 0.172 s | 4.7 | 17.2 s | 59.4 s | ok |
| 128 | 8 | 192 | 600 | 100 | Binomial(32) | 0.0697 s | 0.348 s | 5.0 | 17.1 s | 57.2 s | ok |
| 128 | 8 | 192 | 600 | 200 | Binomial(32) | 0.137 s | 0.696 s | 5.1 | 17.0 s | 57.0 s | ok |
| 128 | 8 | 192 | 600 | 2000 | Binomial(32) | 1.35 s | 8.03 s | 6.0 | 17.0 s | 57.1 s | non-finite or zero gradient |
| 256 | 8 | 384 | 300 | 50 | Binomial(32) | 0.112 s | 0.509 s | 4.5 | 20.0 s | 61.0 s | ok |
| 256 | 8 | 384 | 300 | 100 | Binomial(32) | 0.221 s | 1.07 s | 4.8 | 19.8 s | 59.2 s | ok |
| 256 | 8 | 384 | 300 | 200 | Binomial(32) | 0.439 s | 2.2 s | 5.0 | 19.9 s | 59.1 s | ok |
| 256 | 8 | 384 | 300 | 2000 | Binomial(32) | 4.37 s | 26.2 s | 6.0 | 19.9 s | 59.2 s | non-finite or zero gradient |

### SPEEDY.f90-like configuration

The PrimitiveWet model in a configuration closer to SPEEDY.f90, at the resolutions of the PrimitiveWet resolution sweep with the default transform: seasonal prescribed climatologies of sea surface temperature (`SeasonalOceanClimatology`), land surface temperature (`SeasonalLandTemperature`) and soil moisture (`SeasonalSoilMoisture`) instead of the slab ocean and the land bucket model, constant surface roughness lengths over land and ocean (`ConstantSurfaceRoughness`) instead of the learned ones, and SPEEDY.f90's near-surface wind of 0.95 times the lowermost-layer wind (`ConstantWindSlowdown`) with 5 m/s gusts (`ConstantWindGust`) instead of the learned wind slowdown and Beljaars gusts; sea ice stays thermodynamic, driven by the prescribed sea surface temperature.

| Model | truncation | L | Rings | Δt | SYPD | Runs | Spread | Memory| Compile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PrimitiveWetModel | 32 | 8 | 48 | 2400 | 32344 | 10 | 6% | 7.40 MB | 206.8 s |
| PrimitiveWetModel | 43 | 8 | 64 | 1786 | 20314 | 10 | 27% | 13.29 MB | 21.5 s |
| PrimitiveWetModel | 64 | 8 | 96 | 1200 | 10895 | 10 | 67% | 31.73 MB | 20.1 s |
| PrimitiveWetModel | 86 | 8 | 128 | 893 | 6306 | 10 | 5% | 61.47 MB | 20.2 s |
| PrimitiveWetModel | 128 | 8 | 192 | 600 | 2172 | 10 | 9% | 161.93 MB | 22.0 s |
| PrimitiveWetModel | 171 | 8 | 256 | 449 | 1112 | 10 | 3% | 334.17 MB | 21.7 s |
| PrimitiveWetModel | 256 | 8 | 384 | 300 | 362 | 10 | 4% | 963.21 MB | 23.6 s |
| PrimitiveWetModel | 342 | 8 | 512 | 225 | 126 | 10 | 2% | 2.10 GB | 26.4 s |
| PrimitiveWetModel | 512 | 8 | 768 | 150 | 31 | 10 | 1% | 6.45 GB | 35.6 s |
| PrimitiveWetModel | 86 | 16 | 128 | 893 | 3187 | 10 | 22% | 84.63 MB | 175.4 s |
| PrimitiveWetModel | 128 | 16 | 192 | 600 | 1255 | 10 | 5% | 212.03 MB | 30.6 s |
| PrimitiveWetModel | 171 | 16 | 256 | 449 | 616 | 10 | 8% | 421.66 MB | 30.0 s |
| PrimitiveWetModel | 256 | 16 | 384 | 300 | 177 | 10 | 1% | 1.16 GB | 33.8 s |
| PrimitiveWetModel | 342 | 16 | 512 | 225 | 70 | 10 | 1% | 2.44 GB | 36.9 s |
| PrimitiveWetModel | 512 | 16 | 768 | 150 | 17 | 10 | 1% | 7.20 GB | 46.8 s |
| PrimitiveWetModel | 86 | 24 | 128 | 893 | 2224 | 10 | 18% | 107.85 MB | 208.3 s |
| PrimitiveWetModel | 128 | 24 | 192 | 600 | 835 | 10 | 11% | 262.19 MB | 62.7 s |
| PrimitiveWetModel | 171 | 24 | 256 | 449 | 358 | 10 | 3% | 509.25 MB | 61.4 s |
| PrimitiveWetModel | 256 | 24 | 384 | 300 | 118 | 10 | 2% | 1.35 GB | 64.1 s |
| PrimitiveWetModel | 342 | 24 | 512 | 225 | 40 | 10 | 1% | 2.78 GB | 70.6 s |
| PrimitiveWetModel | 512 | 24 | 768 | 150 | 12 | 10 | 0% | 7.96 GB | 82.1 s |

