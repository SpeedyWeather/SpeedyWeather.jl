module SpeedyWeatherTerrariumExt

using SpeedyWeather
using Terrarium

using Dates
using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using SpeedyWeather: TYPEDSIGNATURES

include("coupling.jl")
include("output.jl")

end # module
