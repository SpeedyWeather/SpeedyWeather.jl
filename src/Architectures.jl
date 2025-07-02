module Architectures 
    
    import KernelAbstractions 

    export AbstractArchitecture
    export CPU, CPUStatic, GPU, CUDAGPU
    export array_type, on_architecture, architecture, device
    export convert_to_device, ismatching, compatible_array_types
    export synchronize 
    
    """
        AbstractArchitecture

    Abstract supertype for architectures supported by SpeedyWeather.
    """
    abstract type AbstractArchitecture end

    """
        AbstractCPU <: AbstractArchitecture

    Abstract supertype for CPU architectures supported by SpeedyWeather.
    """
    abstract type AbstractCPU <: AbstractArchitecture end

    """
        CPU <: AbstractCPU

    Run SpeedyWeather on one CPU node.
    """
    struct CPU{D} <: AbstractCPU
        device::D
    end

    CPU() = CPU(KernelAbstractions.CPU())
    CPUStatic() = CPU(KernelAbstractions.CPU(; static=true))

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

    device(a::CPU) = a.device
    device(a::GPU) = a.device

    """
        architecture(x)

    Return the default architecture that's associated with the input `x`.
    """
    architecture() = nothing
    architecture(::Number) = nothing
    architecture(::Array) = CPU()
    architecture(::Type{<:Array}) = CPU()
    architecture(a::SubArray) = architecture(parent(a))

    """
        architecture(a::AbstractArchitecture)

    Return the architecture of `a`.
    """
    architecture(a::AbstractArchitecture) = a

    """
        array_type(arch::AbstractArchitecture)

    Return the array type that's used with the architecture `arch`.
    """
    array_type(::CPU) = Array
    array_type(::Type{<:CPU}) = Array

    """
        compatible_array_types(arch::AbstractArchitecture)

    Return the array types that are compatible with the architecture `arch`.
    This includes the regular `array_type(arch)` as well as any other array types
    that are compatible with the architecture, e.g. device arrays on GPU.
    """
    compatible_array_types(arch::Type{<:AbstractArchitecture}) = (array_type(arch),) # fallback to array_type
    compatible_array_types(arch::AbstractArchitecture) = (array_type(arch),) # fallback to array_type

    ismatching(arch::Type{<:AbstractArchitecture}, array_T::Type{<:AbstractArray}) = any(arch_array -> array_T <: arch_array, compatible_array_types(arch))
    ismatching(arch::AbstractArchitecture, array_T::Type{<:AbstractArray}) = ismatching(typeof(arch), array_T)
    ismatching(arch::AbstractArchitecture, array::AbstractArray) = ismatching(arch, typeof(array))

    # Fallback
    """
        on_architecture(arch::AbstractArchitecture, a)

    Return `a`, but on the architecture `arch`. 
    """
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

    KernelAbstractions.synchronize(arch::AbstractArchitecture) = KernelAbstractions.synchronize(arch.device)
end 
