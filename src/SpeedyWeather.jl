module SpeedyWeather

# STRUCTURE
using DocStringExtensions

import ConstructionBase: constructorof, getproperties, setproperties

# NUMERICS
import Primes
import Random
import LinearAlgebra: LinearAlgebra, Diagonal
export rotate, rotate!

# GPU, PARALLEL
import Base.Threads: Threads, @threads
import KernelAbstractions
import KernelAbstractions: synchronize
import Adapt: Adapt, adapt, adapt_structure

# INPUT OUTPUT
import TOML
import Dates: Dates, DateTime, Period, Millisecond, Second, Minute, Hour, Day, Week, Month, Year
import Printf: Printf, @sprintf
import Random: randstring
import NCDatasets: NCDatasets, NCDataset, defDim, defVar
import JLD2: jldopen, jldsave, JLDFile 
import CodecZlib
import BitInformation: round, round!
import ProgressMeter

# INTERVALS
using DomainSets.IntervalSets

# to avoid a `using Dates` to pass on DateTime arguments
export DateTime, Millisecond, Second, Minute, Hour, Day, Week, Month, Year, Century, Millenium

# export functions that have many cross-component methods
export initialize!, finalize!

# import device architectures
include("Architectures.jl")
using .Architectures

# export device functions 
export on_architecture, architecture

# import utilities
include("Utils/Utils.jl")
using .Utils

import .Utils: parameters

# export user-facing parameter handling types and methods
export  SpeedyParam, SpeedyParams, parameters, stripparams

# LowerTriangularArrays for spherical harmonics
export  LowerTriangularArrays, 
        LowerTriangularArray,
        LowerTriangularMatrix

export  Spectrum

# indexing styles for LowerTriangularArray/Matrix
export  OneBased, ZeroBased
export  eachmatrix, eachharmonic, eachorder
        
include("LowerTriangularArrays/LowerTriangularArrays.jl")
using .LowerTriangularArrays

# RingGrids
export  RingGrids
export  AbstractGrid, AbstractFullGrid, AbstractReducedGrid
export  AbstractField, AbstractField2D, AbstractField3D
export  Field, Field2D, Field3D,
        FullClenshawField, FullGaussianField,
        FullHEALPixField, FullOctaHEALPixField,
        OctahedralGaussianField, OctahedralClenshawField,
        HEALPixField, OctaHEALPixField,
        OctaminimalGaussianField

export  FullClenshawGrid, FullGaussianGrid,
        FullHEALPixGrid, FullOctaHEALPixGrid,
        OctahedralGaussianGrid, OctahedralClenshawGrid,
        HEALPixGrid, OctaHEALPixGrid,
        OctaminimalGaussianGrid
        
export  eachring, eachlayer, eachgridpoint
export  AnvilInterpolator
export  spherical_distance
export  zonal_mean

include("RingGrids/RingGrids.jl")
using .RingGrids

# SpeedyTransforms
export SpeedyTransforms, SpectralTransform
export transform, transform!
export spectral_truncation, spectral_truncation!
export curl, divergence, curl!, divergence!
export ∇, ∇², ∇⁻², ∇!, ∇²!, ∇⁻²!
export power_spectrum

include("SpeedyTransforms/SpeedyTransforms.jl")
using .SpeedyTransforms
import .SpeedyTransforms: prettymemory

# to be defined in GeoMakie extension
export globe, animate
function globe end
function animate end

# abstract types
include("models/abstract_models.jl")
include("dynamics/abstract_types.jl")
include("physics/abstract_types.jl")

# GEOMETRY CONSTANTS ETC
include("dynamics/vertical_coordinates.jl")
include("dynamics/spectral_grid.jl")
include("dynamics/geometry.jl")
include("dynamics/coriolis.jl")
include("dynamics/planet.jl")
include("dynamics/atmosphere.jl")
include("dynamics/adiabatic_conversion.jl")
include("dynamics/orography.jl")
include("physics/land_sea_mask.jl")

# VARIABLES
include("dynamics/tracers.jl")
include("dynamics/particles.jl")
include("dynamics/clock.jl")
include("dynamics/prognostic_variables.jl")
include("dynamics/set.jl")
include("physics/define_column.jl")
include("dynamics/diagnostic_variables.jl")

# MODEL COMPONENTS
include("dynamics/time_integration.jl")
include("dynamics/forcing.jl")
include("dynamics/drag.jl")
include("dynamics/geopotential.jl")
include("dynamics/virtual_temperature.jl")
include("dynamics/initial_conditions.jl")
include("dynamics/horizontal_diffusion.jl")
include("dynamics/vertical_advection.jl")
include("dynamics/implicit.jl")
include("dynamics/scaling.jl")
include("dynamics/tendencies.jl")
include("dynamics/hole_filling.jl")
include("dynamics/particle_advection.jl")
include("dynamics/random_process.jl")

# PARAMETERIZATIONS
include("physics/albedo.jl")
include("physics/tendencies.jl")
include("physics/column_variables.jl")
include("physics/thermodynamics.jl")
include("physics/boundary_layer.jl")
include("physics/temperature_relaxation.jl")
include("physics/vertical_diffusion.jl")
include("physics/large_scale_condensation.jl")
include("physics/surface_fluxes/surface_fluxes.jl")
include("physics/surface_fluxes/momentum.jl")
include("physics/surface_fluxes/heat.jl")
include("physics/surface_fluxes/moisture.jl")
include("physics/convection.jl")
include("physics/zenith.jl")
include("physics/optical_depth.jl")
include("physics/longwave_radiation.jl")
include("physics/shortwave_radiation.jl")
include("physics/stochastic_physics.jl")

# OCEAN AND LAND
include("physics/ocean.jl")
include("physics/land/land.jl")

# OUTPUT
include("output/schedule.jl")
include("output/feedback.jl")
include("output/netcdf_output.jl")
include("output/restart_file.jl")
include("output/callbacks.jl")
include("output/particle_tracker.jl")
include("output/jld2_output.jl")

# MODELS
include("models/simulation.jl")
include("models/barotropic.jl")
include("models/shallow_water.jl")
include("models/primitive_dry.jl")
include("models/primitive_wet.jl")
include("models/tree.jl")
include("models/set.jl")
end