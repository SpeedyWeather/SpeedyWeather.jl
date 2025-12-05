using SpeedyWeather
using Adapt
using Test
using KernelAbstractions

function load_gpu_package()
    gpu_backend = nothing
    try
        @eval using AMDGPU
        gpu_backend = :AMDGPU
    catch
    end
    if gpu_backend === nothing
        try
            @eval using CUDA
            gpu_backend = :CUDA
        catch
        end
    end
    if gpu_backend === nothing
        try
            @eval using Metal
            gpu_backend = :Metal
        catch
        end
    end
    if gpu_backend === nothing
	throw(ErrorException("No compatible GPU backend found. Neither CUDA, AMDGPU, nor Metal is available. Please ensure that a supported GPU and the corresponding Julia package are installed."))
    end
    return gpu_backend
end

gpu_backend = load_gpu_package()

# KERNEL LAUNCHING AND UTILS
include("kernels_GPU.jl")

# SPECTRAL TRANSFORMS
include("spectral_transform.jl")

# INTERPOLATION OF RINGGRIDS
include("interpolate.jl")

# SET FUNCTIONS, GPU SPECIFIC 
include("set.jl")

# VERTICAL, GPU SPECIFIC 
include("vertical_integration.jl")

# FULL MODELS
include("barotropic.jl")
#include("primitive_dry.jl")

if gpu_backend === :CUDA

    include("CUDA/architecture.jl")

elseif gpu_backend === :AMDGPU

    include("AMDGPU/architecture.jl")

elseif gpu_backend === :Metal
    include("MetalGPU/metal.jl")
end
