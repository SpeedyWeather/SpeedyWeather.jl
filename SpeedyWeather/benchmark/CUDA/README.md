# Benchmarks

created for SpeedyWeather.jl v0.17.4 on Wed, 21 Jan 2026 09:42:23. 

### Machine details

All benchmark simulatons were run on a server with access to an NVIDIA A100 GPU.
```julia
julia> CUDA.versioninfo()
CUDA toolchain: 
- runtime 13.0, artifact installation
- driver 565.57.1 for 13.1
- compiler 13.1

CUDA libraries: 
- CUBLAS: 13.1.0
- CURAND: 10.4.0
- CUFFT: 12.0.0
- CUSOLVER: 12.0.4
- CUSPARSE: 12.6.3
- CUPTI: 2025.3.1 (API 13.0.1)
- NVML: 12.0.0+565.57.1

Julia packages: 
- CUDA: 5.9.6
- GPUArrays: 11.3.4
- GPUCompiler: 1.8.0
- KernelAbstractions: 0.9.39
- CUDA_Driver_jll: 13.1.0+2
- CUDA_Compiler_jll: 0.4.1+1
- CUDA_Runtime_jll: 0.19.2+0

Toolchain:
- Julia: 1.11.7
- LLVM: 16.0.6

1 device:
  0: NVIDIA A40 (sm_86, 21.010 GiB / 44.988 GiB available)
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


### CPU | Float32 | T31 L64 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| - | - | - | - |
| forward_legendre | 528.031 μs| 0 bytes| 0 |
| inverse_legendre | 404.255 μs| 0 bytes| 0 |
| forward_fourier | 425.887 μs| 831.38 KiB| 336 |
| inverse_fourier | 467.349 μs| 806.62 KiB| 336 |

### CPU | Float32 | T127 L64 | OctahedralGaussianGrid | 192 Rings

| Function | Time | Memory | Allocations |
| - | - | - | - |
| forward_legendre | 24.857 ms| 0 bytes| 0 |
| inverse_legendre | 18.349 ms| 0 bytes| 0 |
| forward_fourier | 10.304 ms| 10.00 MiB| 1344 |
| inverse_fourier | 10.671 ms| 9.90 MiB| 1344 |

### CPU | Float32 | T511 L64 | OctahedralGaussianGrid | 768 Rings

| Function | Time | Memory | Allocations |
| - | - | - | - |
| forward_legendre | 2.187 s| 0 bytes| 0 |
| inverse_legendre | 1.838 s| 0 bytes| 0 |
| forward_fourier | 217.900 ms| 148.00 MiB| 5376 |
| inverse_fourier | 224.381 ms| 147.61 MiB| 5376 |

## PrimitiveWet benchmarks, GPU


### GPU(CUDABackend(false, true)) | Float32 | T31 L16 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| - | - | - | - |
| parameterization_tendencies | 277.716 μs| 123.27 KiB| 2105 |
| dynamics_tendencies | 12.857 ms| 2.73 MiB| 85686 |
| implicit_correction | 88.323 μs| 15.70 KiB| 438 |

### GPU(CUDABackend(false, true)) | Float32 | T63 L16 | OctahedralGaussianGrid | 96 Rings

| Function | Time | Memory | Allocations |
| - | - | - | - |
| parameterization_tendencies | 312.098 μs| 123.27 KiB| 2105 |
| dynamics_tendencies | 25.767 ms| 5.30 MiB| 167142 |
| implicit_correction | 93.580 μs| 15.70 KiB| 438 |

### GPU(CUDABackend(false, true)) | Float32 | T127 L16 | OctahedralGaussianGrid | 192 Rings

| Function | Time | Memory | Allocations |
| - | - | - | - |
| parameterization_tendencies | 433.499 μs| 123.27 KiB| 2105 |
| dynamics_tendencies | 52.539 ms| 10.44 MiB| 329924 |
| implicit_correction | 102.023 μs| 15.70 KiB| 438 |

## Transform benchmarks, GPU


### GPU(CUDABackend(false, true)) | Float32 | T31 L64 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| - | - | - | - |
| forward_legendre | 60.701 μs| 3.44 KiB| 139 |
| inverse_legendre | 32.268 μs| 5.64 KiB| 175 |
| forward_fourier | 1.049 ms| 251.94 KiB| 7388 |
| inverse_fourier | 1.300 ms| 263.94 KiB| 8060 |

### GPU(CUDABackend(false, true)) | Float32 | T127 L64 | OctahedralGaussianGrid | 192 Rings

| Function | Time | Memory | Allocations |
| - | - | - | - |
| forward_legendre | 1.899 ms| 3.47 KiB| 141 |
| inverse_legendre | 408.431 μs| 5.67 KiB| 177 |
| forward_fourier | 4.426 ms| 1007.69 KiB| 29548 |
| inverse_fourier | 5.466 ms| 1.03 MiB| 32236 |

### GPU(CUDABackend(false, true)) | Float32 | T511 L64 | OctahedralGaussianGrid | 768 Rings

| Function | Time | Memory | Allocations |
| - | - | - | - |
| forward_legendre | 123.369 ms| 3.53 KiB| 145 |
| inverse_legendre | 177.554 ms| 5.84 KiB| 188 |
| forward_fourier | 24.554 ms| 3.97 MiB| 120253 |
| inverse_fourier | 24.137 ms| 4.16 MiB| 131005 |

## PrimitiveWet benchmarks, CPU


### CPU | Float32 | T31 L16 | OctahedralGaussianGrid | 48 Rings

| Function | Time | Memory | Allocations |
| - | - | - | - |
| parameterization_tendencies | 5.381 ms| 688 bytes| 22 |
| dynamics_tendencies | 3.580 ms| 2.82 MiB| 4463 |
| implicit_correction | 380.150 μs| 5.48 KiB| 79 |

### CPU | Float32 | T63 L16 | OctahedralGaussianGrid | 96 Rings

| Function | Time | Memory | Allocations |
| - | - | - | - |
| parameterization_tendencies | 16.521 ms| 688 bytes| 22 |
| dynamics_tendencies | 17.454 ms| 9.32 MiB| 8543 |
| implicit_correction | 1.401 ms| 5.48 KiB| 79 |

### CPU | Float32 | T127 L16 | OctahedralGaussianGrid | 192 Rings

| Function | Time | Memory | Allocations |
| - | - | - | - |
| parameterization_tendencies | 57.824 ms| 688 bytes| 22 |
| dynamics_tendencies | 103.156 ms| 33.46 MiB| 16703 |
| implicit_correction | 5.607 ms| 5.48 KiB| 79 |


## Benchmark graphs

Times and speedup graphs for the benchmarks are shown below.

![benchmark times](benchmark_times.png "Times") 
![benchmark times](benchmark_speedup.png "Speedup") 
