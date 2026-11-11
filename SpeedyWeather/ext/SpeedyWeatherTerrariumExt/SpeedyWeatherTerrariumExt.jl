module SpeedyWeatherTerrariumExt

using SpeedyWeather
using Terrarium

using Dates
using DocStringExtensions
import SpeedyWeather: TYPEDSIGNATURES

include("coupling.jl")
include("output.jl")

end # module
