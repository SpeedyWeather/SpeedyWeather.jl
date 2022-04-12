"""
    model_setup = ModelSetup(   ::Parameters,
                                ::Constants,
                                ::GeoSpectral,
                                ::Boundaries,
                                ::HorizontalDiffusion)

The ModelSetup struct holds all other structs that contain precalculated constants, whether scalars or
arrays that do not change throughout model integration."""
struct ModelSetup{NF<:AbstractFloat}
    parameters::Parameters
    constants::Constants{NF}
    geospectral::GeoSpectral{NF}
    boundaries::Boundaries{NF}
    horizontal_diffusion::HorizontalDiffusion{NF}
end

"""
    prog_vars = run_speedy(NF,kwargs...)

Runs SpeedyWeather.jl with number format `NF` and any additional parameters in the keyword arguments
`kwargs...`. Any unspeficied parameters will use the default values as defined in `src/parameters.jl`."""
function run_speedy(::Type{NF}=Float64;             # number format, use Float64 as default
                    kwargs...                       # all additional non-default parameters
                    ) where {NF<:AbstractFloat}

    # INITALIZE MODEL
    prog_vars,diag_vars,model_setup = initialize_model(NF;kwargs...)

    # START MODEL INTEGRATION
    time_stepping!(prog_vars,diag_vars,model_setup) 
    return prog_vars                                # return prognostic variables when finished
end

"""
    prog_vars, diag_vars, model_setup = initialize_model(NF,kwargs...)
    
Initialize the model by returning
- `prog_vars`, the initial conditions of the prognostic variables
- `diag_vars`, the preallocated the diagnotic variables (initialised to zero)
- `model_setup`, the collected pre-calculated structs that don't change throughout integration:
parametes, constants, geometry, spectral transform, boundaries, diffusion."""
function initialize_model(::Type{NF}=Float64;       # number format, use Float64 as default
                          kwargs...                 # all additional non-default parameters
                          ) where {NF<:AbstractFloat}

    P = Parameters(NF=NF;kwargs...)                 # all model parameters chosen through kwargs
    C = Constants(P)                                # constants used in model integration
    G = GeoSpectral(P)                              # geometry and spectral transform structs
    B = Boundaries(P)                               # arrays for boundary conditions
    H = HorizontalDiffusion(P,C,G,B)                # precomputed arrays for horizontal diffusion
    M = ModelSetup(P,C,G,B,H)                       # pack all of the above into a ModelSetup struct

    prognostic_vars = initial_conditions(P,B,G)     # initialize prognostic variables
    diagnostic_vars = DiagnosticVariables(G)        # preallocate all diagnostic variables with zeros
    
    return prognostic_vars,diagnostic_vars,M
end