"""
`Prog = run_speedy(NF,kwargs...)`
`Prog = run_speedy(kwargs...)`

Runs SpeedyWeather.jl with number format `NF` and any additional parameters in the keyword arguments
`kwargs...`. Any unspeficied parameters will use the default values as defined in `src/parameters.jl`."""
function run_speedy(::Type{NF}=Float64;         # number format, use Float64 as default
                    kwargs...                   # all additional non-default parameters
                    ) where {NF<:AbstractFloat}

    P = Parameters(NF=NF,kwargs...)
    C = Constants(P)
    G = GeoSpectral(P)
    B = Boundaries(P,G)

    Prog = initial_conditions(P,B,G)
    

    M = Model(P,C,G,B)

    #time_stepping!(Prog,Diag,M)

    return Prog
end



"temporary function, testing out different structure"
function initialize_model(::Type{NF}=Float64;         # number format, use Float64 as default
                          kwargs...                   # all additional non-default parameters
                          ) where {NF<:AbstractFloat}

    P          = Parameters(NF=NF,kwargs...)
    C          = Constants{NF}(P)
    G          = GeoSpectral{NF}(P)
    Geo        = Geometry{NF}(P)
    B          = Boundaries{NF}(P,G)
    
    DiagnosticVars = DiagnosticVariables{NF}(G)
    PrognosticVars = initial_conditions(P,B,G)


    M = Model(P,C,G,Geo,B)
    return PrognosticVars,DiagnosticVars,M
end


"""Struct holding all the model structs """
struct Model
    Parameters
    Constants 
    GeoSpectral
    Geometry
    Boundaries
end
