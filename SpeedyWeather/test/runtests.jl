using SpeedyWeather
using Test

# GENERAL
include("spectral_grid.jl")
include("parameters.jl")

# GPU/KERNELABSTRACTIONS
include("GPU/kernelabstractions.jl")

# DYNAMICS
include("dynamics/run_defaults.jl")
include("dynamics/horizontal_diffusion.jl")
include("dynamics/time_stepping.jl")
include("dynamics/vertical_advection.jl")
include("dynamics/particles.jl")
include("dynamics/particle_advection.jl")
include("dynamics/forcing_drag.jl")
include("dynamics/dates.jl")
include("dynamics/set.jl")
include("dynamics/copy_variables.jl")
include("dynamics/orography.jl")
include("dynamics/initial_conditions.jl")

# VERTICAL LEVELS
include("dynamics/vertical_coordinates.jl")
include("dynamics/geopotential.jl")

# PARAMETERIZATIONS
include("parameterizations/variables.jl")
include("parameterizations/custom_parameterization.jl")
include("parameterizations/zenith.jl")
include("parameterizations/land_sea_mask.jl")
include("parameterizations/ocean_sea_ice.jl")
include("parameterizations/large_scale_condensation.jl")
include("parameterizations/convection.jl")
include("parameterizations/albedo.jl")
include("parameterizations/land.jl")
include("parameterizations/longwave_radiation.jl")
include("parameterizations/shortwave_radiation.jl")
include("parameterizations/surface_fluxes.jl")
include("parameterizations/random_process.jl")
include("parameterizations/stochastic_physics.jl")
include("parameterizations/all_parametrizations.jl")

# OUTPUT/EXTENSION
include("output/callbacks.jl")
include("output/schedule.jl")
include("output/netcdf_output.jl")
include("output/jld2_output.jl")
include("output/feedback.jl")
