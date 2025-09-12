module SpeedyWeatherInternalsJLArraysExt

    using JLArrays
    import SpeedyWeatherInternals.Architectures: Architectures, ismatching, CPU, GPU, architecture, array_type, compatible_array_types, nonparametric_type

    # make JLArrays compatible with standard GPU Architecture
    Architectures.ismatching(arch::GPU, array_type::Type{<:JLArray}) = true
    Architectures.ismatching(arch::Type{<:GPU}, array_type::Type{<:JLArray}) = true
    Architectures.ismatching(arch::GPU, array_type::Type{<:JLArrays.JLDeviceArray}) = true
    Architectures.ismatching(arch::Type{<:GPU}, array_type::Type{<:JLArrays.JLDeviceArray}) = true

    Architectures.architecture(::JLArray) = GPU(JLArrays.JLBackend())
    Architectures.architecture(::Type{<:JLArray}) = GPU(JLArrays.JLBackend())
    Architectures.architecture(::Type{<:JLArrays.JLDeviceArray}) = GPU(JLArrays.JLBackend())

    Architectures.on_architecture(::CPU, a::JLArray) = Array(a)
    Architectures.on_architecture(::GPU{JLBackend}, a::JLArray) = a
    Architectures.on_architecture(::CPU, a::SubArray{<:Any, <:Any, <:JLArray}) = Array(a)

    Architectures.on_architecture(::GPU{JLBackend}, a::Array) = JLArray(a)
    Architectures.on_architecture(::GPU{JLBackend}, a::BitArray) = JLArray(a)
    Architectures.on_architecture(::GPU{JLBackend}, a::SubArray{<:Any, <:Any, <:Array}) = JLArray(a)
    Architectures.on_architecture(::GPU{JLBackend}, a::SubArray{<:Any, <:Any, <:JLArray}) = a
    Architectures.on_architecture(::GPU{JLBackend}, a::StepRangeLen) = a
    
    Architectures.array_type(::Type{GPU{JLBackend}}) = JLArray
    Architectures.array_type(::GPU{JLBackend}) = JLArray
    Architectures.array_type(::GPU{JLBackend}, NF::Type, N::Int) = JLArray{NF, N}
    Architectures.compatible_array_types(::GPU) = (JLArray, JLArrays.JLDeviceArray)
    Architectures.nonparametric_type(::Type{<:JLArray}) = JLArray

end
