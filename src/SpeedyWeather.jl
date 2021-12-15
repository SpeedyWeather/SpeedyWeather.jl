module SpeedyWeather

    using NetCDF, FFTW, LinearAlgebra, Parameters, Dates
    
    import FastGaussQuadrature: gausslegendre

    export run_speedy, 
        Params, GenLogisticCoefs,
        GeoSpectral, Boundaries, Constants, Geometry,
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
    include("geopotential.jl")
    include("horizontal_diffusion.jl")
    include("implicit.jl")
    include("tendencies.jl")
    include("time_stepping.jl")
    include("utils.jl")
    include("run_speedy.jl")
end