module SpeedyWeatherInternalsMetalExt

    import  Metal: Metal, MtlArray, MtlDeviceArray
    import  SpeedyWeatherInternals.Architectures: Architectures, 
            GPU, CPU, MetalGPU, array_type, architecture, on_architecture, architecture,
            compatible_array_types, nonparametric_type


    # DEVICE SETUP FOR METAL
    # extend functions from Architectures
    Architectures.array_type(::GPU) = MtlArray
    Architectures.array_type(::Type{<:GPU}) = MtlArray
    Architectures.array_type(::GPU, NF::Type, N::Int) = MtlArray{NF, N, Metal.PrivateStorage}

    Architectures.compatible_array_types(::GPU) = (MtlArray, MtlDeviceArray)
    Architectures.compatible_array_types(::Type{<:GPU}) = (MtlArray, MtlDeviceArray)

    Architectures.nonparametric_type(::Type{<:MtlArray}) = MtlArray

    Architectures.MetalGPU() = GPU(Metal.MetalBackend())
    Architectures.GPU() = MetalGPU() # default to Metal

    Architectures.architecture(::MtlArray) = MetalGPU()
    Architectures.architecture(::Type{<:MtlArray}) = MetalGPU()
    Architectures.architecture(::Type{<:MtlDeviceArray}) = MetalGPU()

    Architectures.on_architecture(::CPU, a::MtlArray) = Array(a)
    Architectures.on_architecture(::GPU, a::MtlArray) = a
    Architectures.on_architecture(::CPU, a::SubArray{<:Any, <:Any, <:MtlArray}) = Array(a)

    Architectures.on_architecture(::GPU, a::Array) = MtlArray(a)
    Architectures.on_architecture(::GPU, a::BitArray) = MtlArray(a)
    Architectures.on_architecture(::GPU, a::SubArray{<:Any, <:Any, <:MtlArray}) = a
    Architectures.on_architecture(::GPU, a::SubArray{<:Any, <:Any, <:Array}) = MtlArray(a)
    Architectures.on_architecture(::GPU, a::StepRangeLen) = a

    @inline Architectures.convert_to_device(::GPU, args) = Metal.mtlconvert(args)
    @inline Architectures.convert_to_device(::GPU, args::Tuple) = map(Metal.mtlconvert, args)
end
