using SpeedyWeather
using Test

# GENERAL
include("utility_functions.jl")
include("lower_triangular_matrix.jl")
include("grids.jl")

# GPU/KERNELABSTRACTIONS
include("kernelabstractions.jl")

# SPECTRAL TRANSFORM
include("spectral_transform.jl")
include("spectral_gradients.jl")

# DYNAMICS
include("diffusion.jl")
include("time_stepping.jl")

# PHYSICS
include("column_variables.jl")
include("thermodynamics.jl")
include("large_scale_condensation.jl")
include("convection.jl")

# INITIALIZATION AND INTEGRATION
include("initialize.jl")
include("run_speedy.jl")
