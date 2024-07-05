module SpeedyWeatherJLArraysExt

using SpeedyWeather, JLArrays

# for RingGrids and LowerTriangularMatrices:
# every Array needs this method to strip away the parameters
RingGrids.nonparametric_type(::Type{<:JLArray}) = JLArray
LowerTriangularMatrices.nonparametric_type(::Type{<:JLArray}) = JLArray

end # module