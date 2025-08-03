module SpeedyWeatherJLArraysExt

using SpeedyWeather, JLArrays
import SpeedyWeather.Architectures: ismatching, GPU, architecture, array_type, compatible_array_types, nonparametric_type

# make JLArrays compatible with standard GPU Architecture
ismatching(arch::GPU, array_type::Type{<:JLArray}) = true
ismatching(arch::Type{<:GPU}, array_type::Type{<:JLArray}) = true
ismatching(arch::GPU, array_type::Type{<:JLArrays.JLDeviceArray}) = true
ismatching(arch::Type{<:GPU}, array_type::Type{<:JLArrays.JLDeviceArray}) = true
architecture(::Type{<:JLArray}) = GPU(JLArrays.JLBackend())
array_type(::Type{GPU{JLBackend}}) = JLArray
array_type(::GPU{JLBackend}) = JLArray
array_type(::GPU{JLBackend}, NF::Type, N::Int) = JLArray{NF, N}
compatible_array_types(::GPU) = (JLArray, JLArrays.JLDeviceArray)
nonparametric_type(::Type{<:JLArray}) = JLArray

end # module