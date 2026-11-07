module SpeedyWeatherTerrariumExt

using SpeedyWeather
using Terrarium

using Adapt
using Dates
using DocStringExtensions

include("coupling.jl")
include("output.jl")
include("adapt.jl")

end # module
