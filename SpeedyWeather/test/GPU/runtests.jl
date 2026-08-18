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

# TEMPORARY: while pinning down the AMDGPU InvalidIRError, only the MRE crash
# reproducer below is run so buildkite jobs stay fast and focused. All other
# GPU test includes are disabled here on purpose -- restore them once the
# crash is root-caused.
include("mre_amdgpu_crash.jl")

# Companion bisection tiers (tier4-tier8), added on top of the tier1-3 crash
# reproducer above without modifying it -- see that file's header comment.
include("mre_amdgpu_bisect.jl")

# # KERNEL LAUNCHING AND UTILS
# include("kernels_GPU.jl")
#
# # SPECTRAL TRANSFORMS
# include("spectral_transform.jl")
#
# # INTERPOLATION OF RINGGRIDS
# include("interpolate.jl")
#
# # SET FUNCTIONS, GPU SPECIFIC
# include("set.jl")
#
# # VERTICAL, GPU SPECIFIC
# include("vertical_integration.jl")
#
# # FULL MODELS
# include("barotropic.jl")
# include("shallowwater.jl")
# include("primitive_wet.jl")
#
# if gpu_backend === :CUDA
#
#     include("CUDA/architecture.jl")
#
#     # CUDA-GRAPHS ACCELERATED FOURIER TRANSFORM (CUDA-only feature)
#     include("cuda_graphs.jl")
#
#     # REACTANT ON GPU (currently only works tested with CUDA)
#     #include("reactant.jl")
#
# elseif gpu_backend === :AMDGPU
#
#     include("AMDGPU/architecture.jl")
#
# elseif gpu_backend === :Metal
#     include("MetalGPU/metal.jl")
# end
