module SpeedyWeather

# STRUCTURE
using DocStringExtensions

# NUMERICS
import Random
import FastGaussQuadrature
import LinearAlgebra: LinearAlgebra, Diagonal

# GPU, PARALLEL
import Base.Threads: Threads, @threads
import FLoops: FLoops, @floop
import KernelAbstractions
import CUDA
import CUDAKernels
import Adapt: Adapt, adapt, adapt_structure

# INPUT OUTPUT
import TOML
import Dates: Dates, DateTime
import Printf: @sprintf
import NetCDF: NetCDF, NcFile, NcDim, NcVar
import JLD2: jldopen
import CodecZlib
import BitInformation: round, round!
import UnicodePlots
import ProgressMeter

# EXPORT MONOLITHIC INTERFACE TO SPEEDY
export  run_speedy,
        run_speedy!,
        initialize_speedy,
        initialize!,
        run!

export  NoVerticalCoordinates,
        SigmaCoordinates,
        SigmaPressureCoordinates

# EXPORT MODELS
export  Barotropic,             # abstract
        ShallowWater,
        PrimitiveEquation,
        PrimitiveDry,
        PrimitiveWet

export  Model,
        BarotropicModel,        # concrete
        ShallowWaterModel,
        PrimitiveDryModel,
        PrimitiveWetModel

export  Earth,
        EarthAtmosphere

# EXPORT GRIDS
export  SpectralGrid,
        Geometry

export  LowerTriangularMatrix,
        FullClenshawGrid,
        FullGaussianGrid,
        FullHEALPixGrid,
        FullOctaHEALPixGrid,
        OctahedralGaussianGrid,
        OctahedralClenshawGrid,
        HEALPixGrid,
        OctaHEALPixGrid

export  Leapfrog

# EXPORT OROGRAPHIES
export  NoOrography,
        EarthOrography,
        ZonalRidge

export  HyperDiffusion

# EXPORT INITIAL CONDITIONS
export  StartFromFile,
        StartFromRest,
        ZonalJet,
        ZonalWind,
        StartWithRandomVorticity

# EXPORT TEMPERATURE RELAXATION SCHEMES
export  NoTemperatureRelaxation,
        HeldSuarez,
        JablonowskiRelaxation

# EXPORT BOUNDARY LAYER SCHEMES
export  NoBoundaryLayer,
        LinearDrag

# EXPORT VERTICAL DIFFUSION
export  NoVerticalDiffusion,
        VerticalLaplacian

# Large scale condensation
export  SpeedyCondensation

# EXPORT STRUCTS
export  DynamicsConstants,
        SpectralTransform,
        Boundaries,
        PrognosticVariables,
        DiagnosticVariables,
        ColumnVariables

# EXPORT SPECTRAL FUNCTIONS
export  SpectralTransform,
        spectral,
        gridded,
        spectral_truncation

export  OutputWriter, Feedback
        
include("utility_functions.jl")

# LowerTriangularMatrices for spherical harmonics
export LowerTriangularMatrices
include("LowerTriangularMatrices/LowerTriangularMatrices.jl")
using .LowerTriangularMatrices

# RingGrids
export RingGrids
include("RingGrids/RingGrids.jl")
using .RingGrids

# SpeedyTransforms
export SpeedyTransforms
include("SpeedyTransforms/SpeedyTransforms.jl")
using .SpeedyTransforms

# Utility for GPU / KernelAbstractions
include("gpu.jl")                               

# GEOMETRY CONSTANTS ETC
include("abstract_types.jl")
include("dynamics/vertical_coordinates.jl")
include("dynamics/spectral_grid.jl")
include("dynamics/planets.jl")
include("dynamics/atmospheres.jl")
include("dynamics/constants.jl")
include("dynamics/orography.jl")

# VARIABLES
include("dynamics/prognostic_variables.jl")
include("physics/define_column.jl")
include("dynamics/diagnostic_variables.jl")

# MODEL COMPONENTS
include("dynamics/time_integration.jl")
include("dynamics/forcing.jl")
include("dynamics/geopotential.jl")
include("dynamics/initial_conditions.jl")
include("dynamics/horizontal_diffusion.jl")
include("dynamics/implicit.jl")
include("dynamics/scaling.jl")
include("dynamics/tendencies.jl")
include("dynamics/tendencies_dynamics.jl")

# PARAMETERIZATIONS
include("physics/tendencies.jl")
include("physics/column_variables.jl")
include("physics/thermodynamics.jl")
include("physics/boundary_layer.jl")
include("physics/temperature_relaxation.jl")
include("physics/vertical_diffusion.jl")
include("physics/large_scale_condensation.jl")
include("physics/pretty_printing.jl")

# MODELS
include("dynamics/models.jl")

# # PHYSICS
# include("physics/convection.jl")
# include("physics/longwave_radiation.jl")
# include("physics/shortwave_radiation.jl")

# OUTPUT
include("output/output.jl")                     # defines Output
include("output/feedback.jl")                   # defines Feedback

# INTERFACE
include("run_speedy.jl")
end