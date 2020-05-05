module SpeedyWeather

using NetCDF, FFTW, LinearAlgebra, Parameters, Dates,
        FastGaussQuadrature

export run_speedy

include("constants.jl")
include("geometry.jl")
include("fourier.jl")
include("legendre.jl")
include("spectral_transform.jl")

# include("params.jl")
# include("spectral_trans.jl")
# include("prognostics.jl")
# include("input_output.jl")
# include("boundaries.jl")
# include("diagnostics.jl")
# include("geopotential.jl")
# include("horizontal_diffusion.jl")
# include("implicit.jl")
# include("models.jl")

end
