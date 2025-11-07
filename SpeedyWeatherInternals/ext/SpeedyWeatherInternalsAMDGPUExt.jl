module SpeedyWeatherInternalsAMDGPUExt

    import AMDGPU: ROCArray, ROCDeviceArray, ROCBackend, rocconvert, Runtime
    import SpeedyWeatherInternals.Architectures: Architectures, GPU, CPU, ROCGPU, array_type, architecture, on_architecture, architecture, compatible_array_types, nonparametric_type

    # DEVICE SETUP FOR AMDGPU
    # extend functions from Architectures
    Architectures.array_type(::GPU) = ROCArray
    Architectures.array_type(::Type{<:GPU}) = ROCArray
    Architectures.array_type(::GPU, NF::Type, N::Int) = ROCArray{NF, N, Runtime.Mem.HIPBuffer}

    Architectures.compatible_array_types(::GPU) = (ROCArray, ROCDeviceArray)
    Architectures.compatible_array_types(::Type{<:GPU}) = (ROCArray, ROCDeviceArray)

    Architectures.nonparametric_type(::Type{<:ROCArray}) = ROCArray

    Architectures.ROCGPU() = GPU(ROCBackend())
    Architectures.GPU() = ROCGPU()

    Architectures.architecture(::ROCArray) = ROCGPU()
    Architectures.architecture(::Type{<:ROCArray}) = ROCGPU()
    Architectures.architecture(::Type{<:ROCDeviceArray}) = ROCGPU()

    Architectures.on_architecture(::CPU, a::ROCArray) = Array(a)
    Architectures.on_architecture(::GPU, a::ROCArray) = a
    Architectures.on_architecture(::CPU, a::SubArray{<:Any, <:Any, <:ROCArray}) = Array(a)

    Architectures.on_architecture(::GPU, a::Array) = ROCArray(a)
    Architectures.on_architecture(::GPU, a::BitArray) = ROCArray(a)
    Architectures.on_architecture(::GPU, a::SubArray{<:Any, <:Any, <:ROCArray}) = a
    Architectures.on_architecture(::GPU, a::SubArray{<:Any, <:Any, <:Array}) = ROCArray(a)
    Architectures.on_architecture(::GPU, a::StepRangeLen) = a

    @inline Architectures.convert_to_device(::GPU, args) = rocconvert(args)
    @inline Architectures.convert_to_device(::GPU, args::Tuple) = map(rocconvert, args)

end
