using CUDA
using SpeedyWeather
using Adapt
using Test
using KernelAbstractions

# transforms
include("spectral_transform.jl")

# kernels
include("kernels_GPU.jl")

# interpolation
include("interpolate.jl")

# test if the models run on GPU 
include("barotropic.jl")
include("primitive_dry.jl")

