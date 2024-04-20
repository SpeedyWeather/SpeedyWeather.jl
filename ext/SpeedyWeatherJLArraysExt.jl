module SpeedyWeatherJLArraysExt

using SpeedyWeather, JLArrays

# for RingGrids, every Array needs this method to strip away the parameters
SpeedyWeather.RingGrids.nonparametric_type(::Type{JLArray}) = JLArray

end # module