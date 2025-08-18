module SpeedyWeatherJLArraysExt

using SpeedyWeather, JLArrays
import SpeedyWeather.Architectures: Architectures, ismatching, GPU, architecture, array_type, compatible_array_types, nonparametric_type

# make JLArrays compatible with standard GPU Architecture
Architectures.ismatching(arch::GPU, array_type::Type{<:JLArray}) = true
Architectures.ismatching(arch::Type{<:GPU}, array_type::Type{<:JLArray}) = true
Architectures.ismatching(arch::GPU, array_type::Type{<:JLArrays.JLDeviceArray}) = true
Architectures.ismatching(arch::Type{<:GPU}, array_type::Type{<:JLArrays.JLDeviceArray}) = true
Architectures.architecture(::Type{<:JLArray}) = GPU(JLArrays.JLBackend())
Architectures.array_type(::Type{GPU{JLBackend}}) = JLArray
Architectures.array_type(::GPU{JLBackend}) = JLArray
Architectures.array_type(::GPU{JLBackend}, NF::Type, N::Int) = JLArray{NF, N}
Architectures.compatible_array_types(::GPU) = (JLArray, JLArrays.JLDeviceArray)
Architectures.nonparametric_type(::Type{<:JLArray}) = JLArray

end # module