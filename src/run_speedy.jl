# using Parameters, FFTW, LinearAlgebra, NetCDF, FastGaussQuadrature
#
# include("parameters.jl")
# include("constants.jl")
# include("geometry.jl")
# include("spectral_transform.jl")
# include("legendre.jl")
# include("fourier.jl")
# include("boundaries.jl")
# include("diagnostics.jl")
# include("prognostics.jl")
# include("geopotential.jl")

function run_speedy(::Type{T}=Float64;      # number format
                    kwargs...               # all additional parameters
                    ) where {T<:AbstractFloat}

    P = Params(T=T,kwargs...)
    C = Constants{T}(P)
    G = GeoSpectral{T}(P)
    B = Boundaries{T}(P,G)

    Prog = initial_conditions(P,B,G)

    # # TODO
    # Prog = PrognosticVars{T}()
    # Diag = DiagnosticVars{T}()

    #time_stepping!()

    return G,B
end
