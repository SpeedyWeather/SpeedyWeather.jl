using CUDA
using SpeedyWeather
using Adapt
using Test
using KernelAbstractions

# ARCHITECTURE / DEVICE HANDLING 
include("architecture.jl")

# KERNEL LAUNCHING AND UTILS
include("kernels_GPU.jl")

# SPECTRAL TRANSFORMS
include("spectral_transform.jl")

# INTERPOLATION OF RINGGRIDS
include("interpolate.jl")

# FULL MODELS
include("barotropic.jl")
include("primitive_dry.jl")

