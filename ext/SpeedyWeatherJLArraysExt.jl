module SpeedyWeatherJLArraysExt

using SpeedyWeather, JLArrays

# for RingGrids and LowerTriangularArrays:
# every Array needs this method to strip away the parameters
RingGrids.nonparametric_type(::Type{<:JLArray}) = JLArray
LowerTriangularArrays.nonparametric_type(::Type{<:JLArray}) = JLArray

end # module