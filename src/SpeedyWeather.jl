module SpeedyWeather

    # STRUCTURE
    import Parameters: @with_kw, @unpack
    import InteractiveUtils: subtypes

    # NUMERICS
    import Random
    import FastGaussQuadrature
    import AssociatedLegendrePolynomials as Legendre
    import AbstractFFTs
    import FFTW
    import GenericFFT
    import Primes
    import LinearAlgebra
    import Statistics: mean

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
            AbstractGrid,
            FullClenshawGrid,
            FullGaussianGrid,
            FullHEALPixGrid,
            OctahedralGaussianGrid,
            OctahedralClenshawGrid,
            HEALPixGrid,
            HEALPix4Grid,
            FullHEALPix4Grid

    # EXPORT INTERPOLATION FOR GRIDS
    export  Interpolator,
            GridGeometry,
            InterpolationLocator

    # EXPORT STRUCTS
    export  Parameters,
            Constants,
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
    include("lower_triangular_matrix.jl")   # defines LowerTriangularMatrix
    include("Grids/Grids.jl")               # defines AbstractGrid and concrete Grid types
    include("gpu.jl")                       # defines utility for GPU / KernelAbstractions

    include("parameter_structs.jl")         # defines 
    include("spectral_truncation.jl")
    include("abstract_types.jl")            # defines ModelSetup, Barotropic, ShallowWater,
                                            # PrimitiveEquation and others

    include("default_parameters.jl")        # defines Parameters
    include("constants.jl")                 # defines Constants
    include("parameterization_constants.jl")# defines Parameterization Constants
    include("geometry.jl")                  # defines Geometry

    include("spectral_transform.jl")        # defines SpectralTransform
    include("spectral_gradients.jl")

    include("boundaries.jl")                # defines Boundaries
    include("define_diffusion.jl")          # defines HorizontalDiffusion
    include("define_implicit.jl")           # defines Implicit
    include("models.jl")                    # defines ModelSetups

    include("prognostic_variables.jl")      # defines PrognosticVariables
    include("diagnostic_variables.jl")      # defines DiagnosticVariables
    include("initial_conditions.jl")
    include("scaling.jl")
    include("geopotential.jl")

    include("run_speedy.jl")
    include("tendencies_dynamics.jl")
    include("tendencies.jl")
    include("implicit_correction.jl")
    include("diffusion.jl")
    include("output.jl")                    # defines Output
    include("feedback.jl")                  # defines Feedback
    
    # PHYSICS
    include("column_variables.jl")          # defines ColumnVariables
    include("thermodynamics.jl")
    include("tendencies_parametrizations.jl")
    include("convection.jl")
    include("large_scale_condensation.jl")
    include("longwave_radiation.jl")
    include("shortwave_radiation.jl")

    include("time_integration.jl")
    include("pretty_printing.jl")
end
