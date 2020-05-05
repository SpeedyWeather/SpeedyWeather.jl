using Parameters, FFTW, LinearAlgebra, NetCDF, FastGaussQuadrature

include("parameters.jl")
include("constants.jl")
include("geometry.jl")
include("spectral_transform.jl")
include("legendre.jl")
include("fourier.jl")
include("boundaries.jl")
include("diagnostics.jl")

function run_speedy(::Type{T}=Float32;      # number format
                    kwargs...               # all additional parameters
                    ) where {T<:AbstractFloat}

    P = Params(T=T,kwargs...)
    C = Constants{T}(P)
    G = GeoSpectral{T}(P,C)
    B = Boundaries{T}(C,G,S)

    # # TODO
    # Prog = PrognosticVars{T}()
    # Diag = DiagnosticVars{T}()

    #time_stepping!()

    return G,S,B
end
