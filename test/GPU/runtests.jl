using CUDA
using SpeedyWeather
using Adapt
using Test
using KernelAbstractions

# transforms
include("spectral_transform.jl")

# kernels
include("kernels_GPU.jl")