"""
`Prog = run_speedy(NF,kwargs...)`
`Prog = run_speedy(kwargs...)`

Runs SpeedyWeather.jl with number format `NF` and any additional parameters in the keyword arguments
`kwargs...`. Any unspeficied parameters will use the default values as defined in `src/parameters.jl`."""
function run_speedy(::Type{NF}=Float64;         # number format, use Float64 as default
                    kwargs...                   # all additional non-default parameters
                    ) where {NF<:AbstractFloat}

    P = Params(NF=NF,kwargs...)
    C = Constants{NF}(P)
    G = GeoSpectral{NF}(P)
    B = Boundaries{NF}(P,G)

    Prog = initial_conditions(P,B,G)

    # TODO
    # Prog = PrognosticVars{T}()
    # Diag = DiagnosticVars{T}()

    M = ModelSetup(P,C,G,B,...)

    #time_stepping!(Prog,Diag,M)

    return Prog
end