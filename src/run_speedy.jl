"""
    progn_vars = run_speedy(NF,kwargs...)

Runs SpeedyWeather.jl with number format `NF` and any additional parameters in the keyword arguments
`kwargs...`. Any unspecified parameters will use the default values as defined in `src/parameters.jl`."""
function run_speedy(::Type{NF}=Float32;             # number format, use Float32 as default
                    kwargs...                       # all additional non-default parameters
                    ) where {NF<:AbstractFloat}

    # INITALIZE MODEL
    progn_vars,diagn_vars,model_setup = initialize_speedy(NF;kwargs...)

    # START MODEL INTEGRATION
    time_stepping!(progn_vars,diagn_vars,model_setup)
    return progn_vars                               # return prognostic variables when finished
end

"""
    progn_vars, diagn_vars, model_setup = initialize_speedy(NF,kwargs...)

Initialize the model by returning
- `progn_vars`, the initial conditions of the prognostic variables
- `diagn_vars`, the preallocated the diagnotic variables (initialised to zero)
- `model_setup`, the collected pre-calculated structs that don't change throughout integration.

The keyword arguments `kwargs` are the same as for `run_speedy`. The `model_setup` contains
fields that hold the parameters, constants, geometry, spectral transform, boundaries and diffusion."""
function initialize_speedy(::Type{NF}=Float32;      # number format, use Float32 as default
                          kwargs...                 # all additional non-default parameters
                          ) where {NF<:AbstractFloat}

    P = Parameters(NF=NF;kwargs...)                 # all model parameters chosen through kwargs
    C = Constants(P)                                # constants used in model integration
    G = GeoSpectral(P)                              # geometry and spectral transform structs
    B = Boundaries(P)                               # arrays for boundary conditions
    H = HorizontalDiffusion(P,C,G,B)                # precomputed arrays for horizontal diffusion
    D = DeviceSetup(CPUDevice())                    # device the model is running on, so far only CPU

    if P.model == :barotropic                       # pack all of the above into a *Model struct
        M = BarotropicModel(P,C,G,H,D)                # typeof(M) is used to dispatch dynamically
    elseif P.model == :shallowwater                 # to the supported model types
        I = Implicit(P,C,G)
        M = ShallowWaterModel(P,C,G,B,H,I,D)
    elseif P.model == :primitive
        M = PrimitiveEquationModel(P,C,G,B,H,D)
    end

    prognostic_vars = initial_conditions(M)         # initialize prognostic variables
    diagnostic_vars = DiagnosticVariables(G)        # preallocate all diagnostic variables with zeros

    return prognostic_vars, diagnostic_vars, M
end
