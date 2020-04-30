using Parameters, FFTW, LinearAlgebra, NetCDF

include("src/constants.jl")
include("src/geometry.jl")
include("src/legendre.jl")
include("src/fourier.jl")
include("src/spectral_trans.jl")

function run_speedy(::Type{T}=Float32;      # number format
                    kwargs...               # all additional parameters
                    ) where {T<:AbstractFloat}

    P = Parameter(T=T;kwargs...)
    return RunModel(T,P)
end
