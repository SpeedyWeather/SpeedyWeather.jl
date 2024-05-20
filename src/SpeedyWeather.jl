module SpeedyWeather

# STRUCTURE
using DocStringExtensions

# NUMERICS
import Primes
import Random
import FastGaussQuadrature
import LinearAlgebra: LinearAlgebra, Diagonal

# GPU, PARALLEL
import Base.Threads: Threads, @threads
import FLoops: FLoops, @floop
import KernelAbstractions
import CUDA: CUDA, CUDAKernels
import Adapt: Adapt, adapt, adapt_structure

# INPUT OUTPUT
import TOML
import Dates: Dates, DateTime, Period, Millisecond, Second, Minute, Hour, Day
import Printf: Printf, @sprintf
import Random: randstring
import NCDatasets: NCDatasets, NCDataset, defDim, defVar
import JLD2: jldopen
import CodecZlib
import BitInformation: round, round!
import UnicodePlots
import ProgressMeter

# to avoid a `using Dates` to pass on DateTime arguments
export DateTime, Second, Minute, Hour, Day

# export functions that have many cross-component methods
export initialize!, finish!

include("utility_functions.jl")

# LowerTriangularMatrices for spherical harmonics
export  LowerTriangularMatrices, 
        LowerTriangularMatrix,
        LowerTriangularArray
        
include("LowerTriangularMatrices/LowerTriangularMatrices.jl")
using .LowerTriangularMatrices

# RingGrids
export RingGrids
export  AbstractGrid, AbstractGridArray,
        AbstractFullGridarray, AbstractReducedGridArray
export  FullClenshawGrid, FullClenshawArray,
        FullGaussianGrid, FullGaussianArray,
        FullHEALPixGrid, FullHEALPixArray,
        FullOctaHEALPixGrid, FullOctaHEALPixArray,
        OctahedralGaussianGrid, OctahedralGaussianArray,
        OctahedralClenshawGrid, OctahedralClenshawArray,
        HEALPixGrid, HEALPixArray,
        OctaHEALPixGrid, OctaHEALPixArray,
        eachring, eachgrid, plot

include("RingGrids/RingGrids.jl")
using .RingGrids

# SpeedyTransforms
export SpeedyTransforms, SpectralTransform
export spectral, gridded, spectral!, gridded!
export spectral_truncation, spectral_truncation!
export curl, divergence, curl!, divergence!
export ∇, ∇², ∇⁻², ∇!, ∇²!, ∇⁻²!

include("SpeedyTransforms/SpeedyTransforms.jl")
using .SpeedyTransforms
import .SpeedyTransforms: prettymemory

# Utility for GPU / KernelAbstractions
include("gpu.jl")                               

# abstract types
include("models/abstract_models.jl")
include("dynamics/abstract_types.jl")
include("output/abstract_types.jl")
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
include("dynamics/particles.jl")
include("dynamics/clock.jl")
include("dynamics/prognostic_variables.jl")
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

# PARAMETERIZATIONS
include("physics/albedo.jl")
include("physics/tendencies.jl")
include("physics/column_variables.jl")
include("physics/thermodynamics.jl")
include("physics/boundary_layer.jl")
include("physics/temperature_relaxation.jl")
include("physics/vertical_diffusion.jl")
include("physics/large_scale_condensation.jl")
include("physics/surface_fluxes.jl")
include("physics/convection.jl")
include("physics/zenith.jl")
include("physics/shortwave_radiation.jl")
include("physics/longwave_radiation.jl")

# OCEAN AND LAND
include("physics/ocean.jl")
include("physics/land.jl")

# OUTPUT
include("output/schedule.jl")
include("output/output.jl")
include("output/feedback.jl")
include("output/plot.jl")
include("output/callbacks.jl")
include("output/particle_tracker.jl")

# MODELS
include("models/simulation.jl")
include("models/barotropic.jl")
include("models/shallow_water.jl")
include("models/primitive_dry.jl")
include("models/primitive_wet.jl")
include("models/tree.jl")
end