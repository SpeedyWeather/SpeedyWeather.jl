using SpeedyWeather
using Test

FLAG_EXTENDED_TESTS = "extended_tests" in ARGS ? true : false

if FLAG_EXTENDED_TESTS 
    @info "Running extended test suite"
end 

# GENERAL
include("utility_functions.jl")

# GRIDS
include("grids/grids.jl")
include("grids/geodesics.jl")
include("grids/interpolation.jl")
include("grids/reverse.jl")
include("grids/rotate.jl")

# GPU/KERNELABSTRACTIONS
include("gpu/kernelabstractions.jl")

# SPECTRAL TRANSFORM
include("transforms/lower_triangular_matrix.jl")
include("transforms/spectral_transform.jl")
include("transforms/spectral_gradients.jl")
include("transforms/spectrum.jl")
include("transforms/spectral_transform_ad_rules.jl") 
include("transforms/resolutions.jl")

# DYNAMICS
include("dynamics/diffusion.jl")
include("dynamics/time_stepping.jl")
include("dynamics/vertical_advection.jl")
include("dynamics/particles.jl")
include("dynamics/particle_advection.jl")
include("dynamics/forcing_drag.jl")
include("dynamics/dates.jl")
include("dynamics/set.jl")
include("dynamics/orography.jl")
include("dynamics/run_speedy.jl")

# VERTICAL LEVELS
include("dynamics/vertical_coordinates.jl")
include("dynamics/geopotential.jl")

# PHYSICS
include("physics/column_variables.jl")
include("physics/land_sea_mask.jl")
include("physics/ocean.jl")
# include("thermodynamics.jl")
include("physics/large_scale_condensation.jl")
include("physics/convection.jl")

include("physics/albedo.jl")
include("physics/land.jl")
include("physics/optical_depth.jl")
include("physics/longwave_radiation.jl")
# include("physics/shortwave_radiation.jl")
include("physics/surface_fluxes.jl")
include("physics/random_process.jl")
include("physics/stochastic_physics.jl")

# DIFFERENTIABILITY 
if FLAG_EXTENDED_TESTS
    include("differentiability/runtests.jl")
end

# OUTPUT/EXTENSION
include("output/callbacks.jl")
include("output/schedule.jl")
include("output/netcdf_output.jl")
include("output/jld2_output.jl")
