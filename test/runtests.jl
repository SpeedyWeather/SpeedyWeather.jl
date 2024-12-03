using SpeedyWeather
using Test

FLAG_EXTENDED_TESTS = "extended_tests" in ARGS ? true : false

if FLAG_EXTENDED_TESTS 
    @info "Running extended test suite"
end 

# GENERAL
include("utility_functions.jl")
include("dates.jl")
include("lower_triangular_matrix.jl")
include("grids.jl")
include("geodesics.jl")
include("interpolation.jl")
include("set.jl")

# GPU/KERNELABSTRACTIONS
include("kernelabstractions.jl")

# SPECTRAL TRANSFORM
include("spectral_transform.jl")
include("spectral_gradients.jl")
include("spectrum.jl")
include("spectral_transform_ad_rules.jl") 

# DYNAMICS
include("diffusion.jl")
include("time_stepping.jl")
include("vertical_advection.jl")
include("particles.jl")
include("particle_advection.jl")

# VERTICAL LEVELS
include("vertical_coordinates.jl")
include("geopotential.jl")

# PHYSICS
include("column_variables.jl")
include("orography.jl")
include("land_sea_mask.jl")
include("ocean.jl")
# include("thermodynamics.jl")
include("large_scale_condensation.jl")
# include("convection.jl")
include("optical_depth.jl")
include("longwave_radiation.jl")
# include("shortwave_radiation.jl")
include("random_process.jl")

# INITIALIZATION AND INTEGRATION
include("run_speedy.jl")

# EXTENSION
include("extending.jl")
include("callbacks.jl")
include("schedule.jl")

# OUTPUT 
include("netcdf_output.jl")