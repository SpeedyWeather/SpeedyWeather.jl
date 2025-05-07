module Architectures 
    
    import KernelAbstractions 

    export AbstractArchitecture
    export CPU, GPU, CUDAGPU
    export array_type, on_architecture, architecture, device, convert_to_device

    """
    AbstractArchitecture

    Abstract supertype for architectures supported by SpeedyWeather.
    """
    abstract type AbstractArchitecture end

    """
    CPU <: AbstractArchitecture

    Run SpeedyWeather on one CPU node.
    """
    struct CPU <: AbstractArchitecture end

    """
        GPU(device)

    Return a GPU architecture using `device`.
    `device` defauls to CUDA.CUDABackend(always_inline=true)
    """
    struct GPU{D} <: AbstractArchitecture
        device :: D
    end

    # defined here so that it can be extended in SpeedyWeatherCUDAExt
    function CUDAGPU end
    function GPU end 

    #####
    ##### These methods are extended in SpeedyWeatherCUDAExt
    #####

    device(a::CPU) = KernelAbstractions.CPU()
    device(a::GPU) = a.device

    architecture() = nothing
    architecture(::Number) = nothing
    architecture(::Array) = CPU()
    architecture(::Type{<:Array}) = CPU()
    architecture(a::SubArray) = architecture(parent(a))

    array_type(::CPU) = Array
    array_type(::Type{CPU}) = Array

    # Fallback
    on_architecture(arch, a) = a

    # Tupled implementation
    on_architecture(arch::AbstractArchitecture, t::Tuple) = Tuple(on_architecture(arch, elem) for elem in t)
    on_architecture(arch::AbstractArchitecture, nt::NamedTuple) = NamedTuple{keys(nt)}(on_architecture(arch, Tuple(nt)))

    # On architecture for array types
    on_architecture(::CPU, a::Array) = a
    on_architecture(::CPU, a::BitArray) = a
    on_architecture(::CPU, a::SubArray{<:Any, <:Any, <:Array}) = a
    on_architecture(::CPU, a::StepRangeLen) = a

    # Convert arguments to GPU-compatible types
    @inline convert_to_device(arch, args)  = args
    @inline convert_to_device(::CPU, args) = args
  
end 
