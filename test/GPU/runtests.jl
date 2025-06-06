using CUDA
using SpeedyWeather
using Adapt
using Test
using KernelAbstractions

# transforms
include("spectral_transform.jl")

# kernels
include("kernels_GPU.jl")

# test if the models run on GPU 
include("barotropic.jl")
