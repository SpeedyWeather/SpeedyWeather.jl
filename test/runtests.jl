using SpeedyWeather
using Test

# GENERAL
include("utility_functions.jl")
include("lower_triangular_matrix.jl")

# SPECTRAL TRANSFORM
include("spectral_transform.jl")
include("spectral_gradients.jl")

# DYNAMICS
include("diffusion.jl")
include("time_stepping.jl")

# PHYSICS
# include("humidity.jl")
# include("large_scale_condensation.jl")

# INITIALIZATION AND INTEGRATION
include("initialize_prognostics.jl")
include("initialize_from_rest.jl")
include("run_speedy.jl")
include("model_hierarchy.jl")