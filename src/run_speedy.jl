"""
    Prog = run_speedy(NF,kwargs...)

Runs SpeedyWeather.jl with number format `NF` and any additional parameters in the keyword arguments
`kwargs...`. Any unspeficied parameters will use the default values as defined in `src/parameters.jl`."""
function run_speedy(::Type{NF}=Float64;         # number format, use Float64 as default
                    kwargs...                   # all additional non-default parameters
                    ) where {NF<:AbstractFloat}

    # initialize the model setup, precalculate constants
    P = Parameters(NF=NF,kwargs...)
    C = Constants(P)
    G = GeoSpectral(P)
    B = Boundaries(P)
    HD = HorizontalDiffusion(P,C,G,B)
    M = Model(P,C,G,B,HD)               # pack all of the above into a Model struct

    # initialize variables
    Prog = initial_conditions(P,B,G)
    # Diag = DiagnosticVariables(P)

    # start speedy model integration
    # time_stepping!(Prog,Diag,M)

    return Prog
end

"temporary function, testing out different structure"
function initialize_model(::Type{NF}=Float64;         # number format, use Float64 as default
                          kwargs...                   # all additional non-default parameters
                          ) where {NF<:AbstractFloat}

    P   = Parameters(NF=NF,kwargs...)
    C   = Constants(P)
    G   = GeoSpectral(P)
    B   = Boundaries(P)
    HD  = HorizontalDiffusion(P,C,G,B)
    M   = Model(P,C,G,B,HD)

    DiagnosticVars = DiagnosticVariables{NF}(G)
    PrognosticVars = initial_conditions(P,B,G)

    return PrognosticVars,DiagnosticVars,M
end

"""Struct holding all the model structs """
struct Model{NF<:AbstractFloat}
    parameters::Parameters
    constants::Constants{NF}
    geospectral::GeoSpectral{NF}
    boundaries::Boundaries{NF}
    horizontal_diffusion::HorizontalDiffusion{NF}
end
