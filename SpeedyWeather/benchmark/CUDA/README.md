# Benchmarks

created for SpeedyWeather.jl v0.13.0 on Thu, 23 Jan 2025 14:57:16. 

### Machine details

All benchmark simulatons were run on a server with access to an NVIDIA A100 GPU.
```julia
julia> CUDA.versioninfo()
CUDA runtime 12.6, artifact installation
CUDA driver 12.6
NVIDIA driver 535.216.3

CUDA libraries: 
- CUBLAS: 12.6.4
- CURAND: 10.3.7
- CUFFT: 11.3.0
- CUSOLVER: 11.7.1
- CUSPARSE: 12.5.4
- CUPTI: 2024.3.2 (API 24.0.0)
- NVML: 12.0.0+535.216.3

Julia packages: 
- CUDA: 5.5.2
- CUDA_Driver_jll: 0.10.4+0
- CUDA_Runtime_jll: 0.15.5+0

Toolchain:
- Julia: 1.10.4
- LLVM: 15.0.7

1 device:
  0: NVIDIA A40 (sm_86, 10.562 MiB / 44.988 GiB available)
```

### Explanation

Abbreviations in the tables below are as follows, omitted columns use defaults.
- NF: Number format, default: Float32
- T: Spectral resolution, maximum degree of spherical harmonics, default: T31
- L: Number of vertical layers, default: 8 (for 3D models)
- Grid: Horizontal grid, default: OctahedralGaussianGrid
- Rings: Grid-point resolution, number of latitude rings pole to pole
- Model: whether run on a CPU or GPU, default: true
### Running the benchmarks

The benchmark suite here can be reproduced by executing:

```> julia manual_benchmarking.jl```

inside `the SpeedyWeather.jl/benchmark` folder. It will create this `README.md` which can be pushed to the repository for updates or comparison.
## Transform benchmarks, CPU


### CPU() | Float32 | T31 L64 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| - | - | - | - |
| forward_legendre | 421.191 μs| 0 bytes| 0 |
| inverse_legendre | 412.735 μs| 0 bytes| 0 |
| forward_fourier | 521.278 μs| 6.00 KiB| 96 |
| inverse_fourier | 509.406 μs| 6.75 KiB| 96 |

### CPU() | Float32 | T127 L64 | OctahedralGaussianGrid | 192 Rings

| Function | Time | Memory | Allocations |
| - | - | - | - |
| forward_legendre | 19.356 ms| 0 bytes| 0 |
| inverse_legendre | 19.522 ms| 0 bytes| 0 |
| forward_fourier | 10.908 ms| 24.00 KiB| 384 |
| inverse_fourier | 10.829 ms| 27.00 KiB| 384 |

### CPU() | Float32 | T511 L64 | OctahedralGaussianGrid | 768 Rings

| Function | Time | Memory | Allocations |
| - | - | - | - |
| forward_legendre | 1.909 s| 0 bytes| 0 |
| inverse_legendre | 1.841 s| 0 bytes| 0 |
| forward_fourier | 237.968 ms| 96.00 KiB| 1536 |
| inverse_fourier | 237.545 ms| 108.00 KiB| 1536 |

## Transform benchmarks, GPU


### GPU() | Float32 | T31 L64 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| - | - | - | - |
| forward_legendre | 64.542 μs| 3.31 KiB| 115 |
| inverse_legendre | 34.054 μs| 3.95 KiB| 127 |
| forward_fourier | 862.168 μs| 144.47 KiB| 4994 |
| inverse_fourier | 1.630 ms| 291.47 KiB| 10754 |

### GPU() | Float32 | T127 L64 | OctahedralGaussianGrid | 192 Rings

| Function | Time | Memory | Allocations |
| - | - | - | - |
| forward_legendre | 2.096 ms| 3.33 KiB| 116 |
| inverse_legendre | 449.895 μs| 3.97 KiB| 128 |
| forward_fourier | 3.641 ms| 576.47 KiB| 19970 |
| inverse_fourier | 6.768 ms| 1.14 MiB| 43010 |

### GPU() | Float32 | T511 L64 | OctahedralGaussianGrid | 768 Rings

| Function | Time | Memory | Allocations |
| - | - | - | - |
| forward_legendre | 137.463 ms| 3.34 KiB| 117 |
| inverse_legendre | 145.007 ms| 4.00 KiB| 130 |
| forward_fourier | 15.543 ms| 2.25 MiB| 79874 |
| inverse_fourier | 28.117 ms| 4.55 MiB| 172034 |


## Benchmark graphs

Times and speedup graphs for the benchmarks are shown below.

![benchmark times](https://github.com/user-attachments/assets/f6ec4605-fd89-4e87-b4fb-dae0f50981e4 "Times") 
![benchmark times](https://github.com/user-attachments/assets/7422a253-cc3e-4fc9-a259-e276223e7a93 "Speedup") 
