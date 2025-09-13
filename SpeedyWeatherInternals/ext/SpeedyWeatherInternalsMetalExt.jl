module SpeedyWeatherInternalsMetalExt

    import  Metal: Metal, MtlArray, MtlDeviceArray
    import  SpeedyWeatherInternals.Architectures: Architectures, 
            GPU, CPU, MetalGPU, array_type, architecture, on_architecture, architecture,
            compatible_array_types, nonparametric_type


    const MtlGPU = GPU{MetalBackend}

    # DEVICE SETUP FOR METAL
    # extend functions from Architectures
    Architectures.array_type(::MtlGPU) = MtlArray
    Architectures.array_type(::Type{<:MtlGPU}) = MtlArray
    Architectures.array_type(::GPU, NF::Type, N::Int) = MtlArray{NF, N, Metal.PrivateStorage}

    Architectures.compatible_array_types(::MtlGPU) = (MtlArray, MtlDeviceArray)
    Architectures.compatible_array_types(::Type{<:MtlGPU}) = (MtlArray, MtlDeviceArray)

    Architectures.nonparametric_type(::Type{<:MtlArray}) = MtlArray

    Architectures.MetalGPU() = GPU(Metal.MetalBackend())
    Architectures.GPU() = MetalGPU() # default to Metal

    Architectures.architecture(::MtlArray) = MetalGPU()
    Architectures.architecture(::Type{<:MtlArray}) = MetalGPU()
    Architectures.architecture(::Type{<:MtlDeviceArray}) = MetalGPU()

    Architectures.on_architecture(::CPU, a::MtlArray) = Array(a)
    Architectures.on_architecture(::MtlGPU, a::MtlArray) = a
    Architectures.on_architecture(::CPU, a::SubArray{<:Any, <:Any, <:MtlArray}) = Array(a)

    Architectures.on_architecture(::MtlGPU, a::Array) = MtlArray(a)
    Architectures.on_architecture(::MtlGPU, a::BitArray) = MtlArray(a)
    Architectures.on_architecture(::MtlGPU, a::SubArray{<:Any, <:Any, <:MtlArray}) = a
    Architectures.on_architecture(::MtlGPU, a::SubArray{<:Any, <:Any, <:Array}) = MtlArray(a)
    Architectures.on_architecture(::MtlGPU, a::StepRangeLen) = a

    @inline Architectures.convert_to_device(::MtlGPU, args) = Metal.mtlconvert(args)
    @inline Architectures.convert_to_device(::MtlGPU, args::Tuple) = map(Metal.mtlconvert, args)
end
