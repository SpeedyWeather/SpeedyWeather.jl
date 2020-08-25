using Parameters, FFTW, LinearAlgebra, NetCDF, FastGaussQuadrature

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
include("run_speedy.jl")

G,B = run_speedy()
