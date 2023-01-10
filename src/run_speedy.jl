"""
    progn_vars = run_speedy(NF,Model;kwargs...)     or
    progn_vars = run_speedy(NF;kwargs...)           or
    progn_vars = run_speedy(Model;kwargs...)

Runs SpeedyWeather.jl with number format `NF` and the model `Model` and any additional parameters
in the keyword arguments `kwargs...`. Any unspecified parameters will use the default values as
defined in `src/default_parameters.jl`."""
function run_speedy(::Type{NF}=DEFAULT_NF,          # default number format
                    ::Type{Model}=DEFAULT_MODEL;    # default model
                    kwargs...                       # all additional non-default parameters
                    ) where {NF<:AbstractFloat,Model<:ModelSetup}

    # INITALIZE MODEL
    progn_vars,diagn_vars,model_setup = initialize_speedy(NF,Model;kwargs...)

    # START MODEL INTEGRATION
    time_stepping!(progn_vars,diagn_vars,model_setup)
    return progn_vars                               # return prognostic variables when finished
end

# if only Model M provided, use default number format NF
run_speedy(::Type{Model};kwargs...) where {Model<:ModelSetup} = run_speedy(DEFAULT_NF,Model;kwargs...)

"""
    progn_vars, diagn_vars, model_setup = initialize_speedy(NF,Model;kwargs...) or
    progn_vars, diagn_vars, model_setup = initialize_speedy(NF,kwargs...)       or
    progn_vars, diagn_vars, model_setup = initialize_speedy(Model,kwargs...)

Initialize the model by returning
- `progn_vars`, the initial conditions of the prognostic variables
- `diagn_vars`, the preallocated the diagnotic variables (initialised to zero)
- `model_setup`, the collected pre-calculated structs that don't change throughout integration.

The keyword arguments `kwargs` are the same as for `run_speedy`. The `model_setup` contains
fields that hold the parameters, constants, geometry, spectral transform, boundaries and diffusion."""
function initialize_speedy( ::Type{NF}=DEFAULT_NF,          # default number format
                            ::Type{Model}=DEFAULT_MODEL;    # default model
                            kwargs...                       # all additional non-default parameters
                            ) where {NF<:AbstractFloat,Model<:ModelSetup}

    ConcreteModel = default_concrete_model(Model)   # pick default concrete type if Model abstract
    P = Parameters{ConcreteModel}(NF=NF;kwargs...)  # all model parameters chosen through kwargs
    Random.seed!(P.seed)                            # seed Julia's default RNG for reproducibility
    
    C = Constants(P)                                # constants used in model integration
    G = Geometry(P)                                 # everything grid
    S = SpectralTransform(P)                        # everything spectral transform
    B = Boundaries(P,S,G)                           # arrays for boundary conditions
    H = HorizontalDiffusion(P,C,G,S,B)              # precomputed arrays for horizontal diffusion
    D = DeviceSetup(CPUDevice())                    # device the model is running on, so far only CPU
    
    if ConcreteModel <: Barotropic                  # pack all of the above into a *Model struct
        M = BarotropicModel(P,C,G,S,H,D)            # typeof(M) is used to dispatch dynamically
    elseif ConcreteModel <: ShallowWater            # to the supported model types
        I = Implicit(P,C,S)                         # precompute arrays for semi-implicit corrections
        M = ShallowWaterModel(P,C,G,S,B,H,I,D)
    elseif ConcreteModel <: PrimitiveDryCore        # no humidity 
        I = Implicit(P,C,S)
        K = ParameterizationConstants(P)
        M = PrimitiveDryCoreModel(P,C,K,G,S,B,H,I,D)
    elseif ConcreteModel <: PrimitiveWetCore        # with humidity
        I = Implicit(P,C,S)
        K = ParameterizationConstants(P)
        M = PrimitiveWetCoreModel(P,C,K,G,S,B,H,I,D)
    end

    prognostic_vars = initial_conditions(M)         # initialize prognostic variables
    diagnostic_vars = DiagnosticVariables(G,S)      # preallocate all diagnostic variables with zeros

    return prognostic_vars, diagnostic_vars, M
end

# if only Model M provided, use default number format NF
initialize_speedy(::Type{Model};kwargs...) where {Model<:ModelSetup} = initialize_speedy(DEFAULT_NF,Model;kwargs...)

"""
    progn = run_speedy!(progn::PrognosticVariables,
                        diagn::DiagnosticVariables,
                        M::ModelSetup)

Convenience function that can be used in combination with `initialize_speedy(args...;kwargs...)` as

    P,D,M = initialize_speedy(kwargs...)
    # possibly change P, D, M manually
    run_speedy!(P,D,M)
    # or investigate D, M afterwards

to allow for access to the prognostic/diagnostic variables before the time integration is started."""
function run_speedy!(   progn::PrognosticVariables, # all prognostic variables
                        diagn::DiagnosticVariables, # all pre-allocated diagnostic variables
                        M::ModelSetup,              # all precalculated structs
                        )
    time_stepping!(progn,diagn,M)
    return progn
end