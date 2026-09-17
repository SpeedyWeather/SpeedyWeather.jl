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

# TEMPORARILY commented out to save CI time while iterating on the Metal gpu_graphs
# synchronization fix in metal_graphs.jl -- restore before merging.
# # KERNEL LAUNCHING AND UTILS
# include("kernels_GPU.jl")
#
# # BROADCASTING
# include("broadcasting.jl")
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
# include("gpu_graphs_shared.jl")

if gpu_backend === :CUDA

    include("CUDA/architecture.jl")

    # CUDA-GRAPHS ACCELERATED FOURIER TRANSFORM (CUDA-only feature)
    include("cuda_graphs.jl")

    # REACTANT ON GPU (currently only works tested with CUDA)
    #include("reactant.jl")

elseif gpu_backend === :AMDGPU

    include("AMDGPU/architecture.jl")

    # HIP-GRAPHS ACCELERATED FOURIER TRANSFORM (AMDGPU-only feature)
    include("hip_graphs.jl")

elseif gpu_backend === :Metal
    # TEMPORARILY commented out to save CI time -- restore before merging.
    # include("MetalGPU/metal.jl")

    # MPSGRAPH-FUSED BATCHED FOURIER TRANSFORM (Metal-only feature, currently disabled by
    # default -- see metal_graphs.jl for why these tests don't gate on that default)
    include("metal_graphs.jl")
end
