module SpeedyWeather

    # STRUCTURE
    import Parameters: @unpack

    # NUMERICS
    import Random
    import FastGaussQuadrature

    # GPU
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

    # EXPORT MAIN INTERFACE TO SPEEDY
    export  run_speedy,
            run_speedy!,
            initialize_speedy

    # EXPORT MODELS
    export  Barotropic,
            BarotropicModel,
            ShallowWater,
            ShallowWaterModel,
            PrimitiveEquation,
            PrimitiveDryCore,
            PrimitiveWetCore,
            PrimitiveDryCoreModel,
            PrimitiveWetCoreModel

    # EXPORT GRIDS
    export  LowerTriangularMatrix,
            FullClenshawGrid,
            FullGaussianGrid,
            FullHEALPixGrid,
            FullOctaHEALPixGrid,
            OctahedralGaussianGrid,
            OctahedralClenshawGrid,
            HEALPixGrid,
            OctaHEALPixGrid

    # EXPORT OROGRAPHIES
    export  NoOrography,
            EarthOrography,
            ZonalRidge

    # EXPORT INITIAL CONDITIONS
    export  StartFromFile,
            StartFromRest,
            ZonalJet,
            ZonalWind,
            StartWithVorticity

    # EXPORT STRUCTS
    export  Parameters,
            DynamicsConstants,
            ParameterizationConstants,
            Geometry,
            SpectralTransform,
            Boundaries,
            PrognosticVariables,
            DiagnosticVariables,
            ColumnVariables

    # EXPORT SPECTRAL FUNCTIONS
    export  spectral,
            gridded,
            spectral_truncation
            
    include("utility_functions.jl")
    include("abstract_types.jl")
    
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
        
    include("gpu.jl")                       # defines utility for GPU / KernelAbstractions
    include("default_parameters.jl")        # defines Parameters
    
    # DYNAMICS
    include("dynamics/constants.jl")                # defines DynamicsConstants
    include("dynamics/geometry.jl")                 # defines Geometry
    include("dynamics/boundaries.jl")               # defines Boundaries
    include("dynamics/define_diffusion.jl")         # defines HorizontalDiffusion
    include("dynamics/define_implicit.jl")          # defines ImplicitShallowWater, ImplicitPrimitiveEq
    include("dynamics/parameter_structs.jl")        # defines GenLogisticCoefs
    include("dynamics/models.jl")                   # defines ModelSetups
    include("dynamics/prognostic_variables.jl")     # defines PrognosticVariables
    include("dynamics/diagnostic_variables.jl")     # defines DiagnosticVariables
    include("dynamics/initial_conditions.jl")
    include("dynamics/scaling.jl")
    include("dynamics/geopotential.jl")
    include("dynamics/tendencies_dynamics.jl")
    include("dynamics/tendencies.jl")
    include("dynamics/implicit.jl")
    include("dynamics/diffusion.jl")
    include("dynamics/time_integration.jl")
    
    # PHYSICS
    include("physics/constants.jl")                 # defines Parameterization Constants
    include("physics/parameter_structs.jl")         # defines MagnusCoefs, RadiationCoefs
    include("physics/column_variables.jl")          # defines ColumnVariables
    include("physics/thermodynamics.jl")
    include("physics/tendencies_parametrizations.jl")
    include("physics/convection.jl")
    include("physics/large_scale_condensation.jl")
    include("physics/longwave_radiation.jl")
    include("physics/shortwave_radiation.jl")
    
    # OUTPUT
    include("output/output.jl")                     # defines Output
    include("output/feedback.jl")                   # defines Feedback
    include("output/pretty_printing.jl")
    
    # INTERFACE
    include("run_speedy.jl")
end