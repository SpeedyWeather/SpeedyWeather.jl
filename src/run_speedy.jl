"""
    progn_vars = run_speedy(NF,kwargs...)

Runs SpeedyWeather.jl with number format `NF` and any additional parameters in the keyword arguments
`kwargs...`. Any unspecified parameters will use the default values as defined in `src/parameters.jl`."""
function run_speedy(::Type{NF}=Float64;             # default number format
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
function initialize_speedy(::Type{NF}=Float64;      # default number format
                          kwargs...                 # all additional non-default parameters
                          ) where {NF<:AbstractFloat}

    P = Parameters(NF=NF;kwargs...)                 # all model parameters chosen through kwargs
    Random.seed!(P.seed)                            # seed Julia's default RNG for reproducibility
    
    C = Constants(P)                                # constants used in model integration
    G = Geometry(P)                                 # everything grid
    S = SpectralTransform(P)                        # everything spectral transform
    B = Boundaries(P,S)                             # arrays for boundary conditions
    H = HorizontalDiffusion(P,C,G,S,B)              # precomputed arrays for horizontal diffusion
    D = DeviceSetup(CPUDevice())                    # device the model is running on, so far only CPU

    if P.model <: Barotropic                        # pack all of the above into a *Model struct
        M = BarotropicModel(P,C,G,S,H,D)            # typeof(M) is used to dispatch dynamically
    elseif P.model <: ShallowWater                  # to the supported model types
        I = Implicit(P,C,S)
        M = ShallowWaterModel(P,C,G,S,B,H,I,D)
    elseif P.model <: PrimitiveEquation
        I = Implicit(P,C,S)
        M = PrimitiveEquationModel(P,C,G,S,B,H,I,D)
    end

    prognostic_vars = initial_conditions(M)         # initialize prognostic variables
    diagnostic_vars = DiagnosticVariables(G,S)      # preallocate all diagnostic variables with zeros

    return prognostic_vars, diagnostic_vars, M
end

"""
    progn = run_speedy!(progn::PrognosticVariables,
                        diagn::DiagnosticVariables,
                        M::ModelSetup)

Convenience function that can be used in combination with `initialize_speedy(kwargs...)` as

    P,D,M = initialize_speedy(kwargs...)
    # possibly change P, D, M manually
    run_speedy!(P,D,M)

to allow for access to the prognostic/diagnostic variables before the time integration is started."""
function run_speedy!(   progn::PrognosticVariables, # all prognostic variables
                        diagn::DiagnosticVariables, # all pre-allocated diagnostic variables
                        M::ModelSetup,              # all precalculated structs
                        )
    time_stepping!(progn,diagn,M)
    return progn
end