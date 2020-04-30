using Parameters, FFTW, LinearAlgebra, NetCDF

include("constants.jl")
include("geometry.jl")
include("legendre.jl")
include("fourier.jl")
include("spectral_trans.jl")

function run_speedy(::Type{T}=Float32;      # number format
                    kwargs...               # all additional parameters
                    ) where {T<:AbstractFloat}

    constants = Constants()

    return RunModel(T,P)
end

constants = Constants()
geometry = Geometry{Float64}(96,48,8,30,constants)
spectral_trans = SpectralTrans{Float64}(constants,geometry)
boundaries = Boundaries{Float64}(constants,geometry,spectral_trans)
