module SpeedyWeatherJLArraysExt

using SpeedyWeather, JLArrays
import SpeedyWeather.Architectures: ismatching, GPU, architecture, array_type, compatible_array_types

# for RingGrids and LowerTriangularArrays:
# every Array needs this method to strip away the parameters
RingGrids.nonparametric_type(::Type{<:JLArray}) = JLArray
LowerTriangularArrays.nonparametric_type(::Type{<:JLArray}) = JLArray

# make JLArrays compatible with standard GPU Architecture
ismatching(arch::GPU, array_type::Type{<:JLArray}) = true
ismatching(arch::Type{<:GPU}, array_type::Type{<:JLArray}) = true
ismatching(arch::GPU, array_type::Type{<:JLArrays.JLDeviceArray}) = true
ismatching(arch::Type{<:GPU}, array_type::Type{<:JLArrays.JLDeviceArray}) = true
architecture(::Type{<:JLArray}) = GPU(JLArrays.JLBackend())
array_type(::Type{GPU{JLBackend}}) = JLArray
array_type(::GPU{JLBackend}) = JLArray
compatible_array_types(::GPU) = (JLArray, JLArrays.JLDeviceArray)

end # module