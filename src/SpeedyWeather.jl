module SpeedyWeather

    using NetCDF, FFTW, LinearAlgebra, Parameters, Dates, FastGaussQuadrature

    export run_speedy, Params, GeoSpectral, Boundaries,
        fourier, fourier_inverse,
        legendre, legendre_inverse,
        spectral, gridded

    include("parameters.jl")
    include("constants.jl")
    include("geometry.jl")
    include("spectral_transform.jl")
    include("legendre.jl")
    include("fourier.jl")
    include("boundaries.jl")
    include("diagnostics.jl")
    include("prognostics.jl")
    include("geopotential.jl")
    include("horizontal_diffusion.jl")
    include("implicit.jl")
    include("tendencies.jl")

    include("run_speedy.jl")

end
