using SpeedyWeather
using Test

# GENERAL
include("utility_functions.jl")
include("lower_triangular_matrix.jl")
include("grids.jl")
include("set_vars.jl")

# GPU/KERNELABSTRACTIONS
include("kernelabstractions.jl")

# SPECTRAL TRANSFORM
include("spectral_transform.jl")
include("spectral_gradients.jl")

# DYNAMICS
include("diffusion.jl")
include("time_stepping.jl")

# VERTICAL LEVELS
include("vertical_levels.jl")

# PHYSICS
include("column_variables.jl")
include("thermodynamics.jl")
include("large_scale_condensation.jl")
include("convection.jl")
include("longwave_radiation.jl")

# INITIALIZATION AND INTEGRATION
include("initialize.jl")
include("run_speedy.jl")
include("run_speedy_with_output.jl")
