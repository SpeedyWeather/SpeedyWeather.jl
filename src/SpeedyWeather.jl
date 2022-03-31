module SpeedyWeather

    # STRUCTURE
    import Parameters: @with_kw, @unpack
    
    # NUMERICS
    import FastGaussQuadrature
    import AssociatedLegendrePolynomials
    import FFTW
    import LinearAlgebra

    # INPUT OUTPUT 
    import Dates: DateTime
    import Printf: @sprintf
    import NetCDF: NetCDF, NcFile, NcDim, NcVar
    import BitInformation: round, round!

    # EXPORT MAIN INTERFACE TO SPEEDY
    export run_speedy

    # EXPORT STRUCTS
    export Parameters, GenLogisticCoefs,
        GeoSpectral, Boundaries, Constants, Geometry, SpectralTransform,
        PrognosticVariables, DiagnosticVariables
    
    # EXPORT SPECTRAL FUNCTIONS
    export  spectral, gridded,
        spectral!, gridded!,
        spectral_truncation, spectral_truncation!,
        spectral_interpolation, triangular_truncation

    include("utility_functions.jl")
    include("parameter_structs.jl")
    include("spectral_truncation.jl")

    include("parameters.jl")
    include("constants.jl")
    include("geometry.jl")
    include("spectral_transform.jl")
    include("spectral_gradients.jl")

    include("boundaries.jl")
    include("diagnostics.jl")
    include("prognostic_variables.jl")
    include("diagnostic_variables.jl")
    include("geopotential.jl")
    include("horizontal_diffusion.jl")
    include("implicit.jl")
    include("parametrization_tendencies.jl")
    include("dynamics_tendencies.jl")
    include("tendencies.jl")
    include("feedback.jl")
    include("output.jl")
    include("time_integration.jl")
    include("run_speedy.jl")
end