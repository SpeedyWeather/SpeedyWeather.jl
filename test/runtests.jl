using SpeedyWeather
using Test

include("utility_functions.jl")
include("spectral_transform.jl")
include("spectral_gradients.jl")
include("diffusion.jl")
include("time_stepping.jl")
include("initialize_from_rest.jl")
include("run_speedy.jl")

# PHYSICS
include("humidity.jl")
include("large_scale_condensation.jl")
