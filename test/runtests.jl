using SpeedyWeather
using Test

# GENERAL
include("utility_functions.jl")
include("lower_triangular_matrix.jl")

# GPU/KERNELABSTRACTIONS 
include("kernelabstractions.jl")

# AD/ENZYME 
include("enzyme.jl")

# SPECTRAL TRANSFORM
include("spectral_transform.jl")
include("spectral_gradients.jl")

# DYNAMICS
include("diffusion.jl")
include("time_stepping.jl")

# PHYSICS (tests need to be updated)
# include("humidity.jl")
# include("large_scale_condensation.jl")

# INITIALIZATION AND INTEGRATION
include("initialize.jl")
include("run_speedy.jl")