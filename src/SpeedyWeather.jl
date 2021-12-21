module SpeedyWeather

    # STRUCTURE
    import Parameters: @with_kw, @unpack
    
    # NUMERICS
    import FastGaussQuadrature: gausslegendre
    import FFTW
    import LinearAlgebra

    # INPUT OUTPUT 
    import Dates: DateTime
    import Printf: @sprintf
    import NetCDF: NcFile, NcDim, NcVar
    import BitInformation: round, round!

    export run_speedy, 
        Params, GenLogisticCoefs,
        GeoSpectral, Boundaries, Constants, Geometry,
        PrognosticVariables, DiagnosticVariables,
        fourier, fourier_inverse,
        legendre, legendre_inverse,
        spectral, gridded

    include("parameter_structs.jl")
    include("parameters.jl")
    include("constants.jl")
    include("geometry.jl")
    include("spectral_transform.jl")
    include("legendre.jl")
    include("fourier.jl")
    include("boundaries.jl")
    include("diagnostics.jl")
    include("prognostic_variables.jl")
    include("diagnostic_variables.jl")
    include("geopotential.jl")
    include("horizontal_diffusion.jl")
    include("implicit.jl")
    include("tendencies.jl")
    include("feedback.jl")
    include("output.jl")
    include("time_integration.jl")
    include("utils.jl")
    include("run_speedy.jl")
end