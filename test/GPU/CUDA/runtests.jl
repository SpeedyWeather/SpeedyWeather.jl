using CUDA
using SpeedyWeather
using Adapt
using Test
using KernelAbstractions

# ARCHITECTURE / DEVICE HANDLING 
include("architecture.jl")

# SPECTRAL TRANSFORMS
include("spectral_transform.jl")
