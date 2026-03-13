export BarotropicModel, initialize!

"""
The BarotropicModel contains all model components needed for the simulation of the
barotropic vorticity equations. To be constructed like

    model = BarotropicModel(spectral_grid; kwargs...)

with `spectral_grid::SpectralGrid` used to initalize all non-default components
passed on as keyword arguments, e.g. `planet=Earth(spectral_grid)`. Fields, representing
model components, are
$(TYPEDFIELDS)"""
@parameterized @kwdef mutable struct BarotropicModel{
        SG,     # <:SpectralGrid
        AR,     # <:AbstractArchitecture,
        GE,     # <:AbstractGeometry,
        PL,     # <:AbstractPlanet,
        AT,     # <:AbstractAtmosphere,
        CO,     # <:AbstractCoriolis,
        FR,     # <:AbstractForcing,
        DR,     # <:AbstractDrag,
        PA,     # <:AbstractParticleAdvection,
        IC,     # <:AbstractInitialConditions,
        RP,     # <:AbstractRandomProcess,
        TS,     # <:AbstractTimeStepper,
        ST,     # <:SpectralTransform{NF},
        IM,     # <:AbstractImplicit,
        HD,     # <:AbstractHorizontalDiffusion,
        OU,     # <:AbstractOutput,
        FB,     # <:AbstractFeedback,
    } <: Barotropic

    spectral_grid::SG
    architecture::AR = spectral_grid.architecture

    # DYNAMICS
    @component geometry::GE = Geometry(spectral_grid)
    @component planet::PL = Earth(spectral_grid)
    @component atmosphere::AT = EarthDryAtmosphere(spectral_grid)
    @component coriolis::CO = Coriolis(spectral_grid)
    @component forcing::FR = KolmogorovFlow(spectral_grid)
    @component drag::DR = LinearVorticityDrag(spectral_grid)
    @component particle_advection::PA = nothing
    @component initial_conditions::IC = InitialConditions(spectral_grid, Barotropic)

    # VARIABLES
    random_process::RP = nothing
    tracers::TRACER_DICT = TRACER_DICT()

    # NUMERICS
    time_stepping::TS = Leapfrog(spectral_grid)
    spectral_transform::ST = SpectralTransform(spectral_grid)
    implicit::IM = nothing
    horizontal_diffusion::HD = HyperDiffusion(spectral_grid)

    # OUTPUT
    output::OU = NetCDFOutput(spectral_grid, Barotropic)
    callbacks::Dict{Symbol, AbstractCallback} = Dict{Symbol, AbstractCallback}()
    feedback::FB = Feedback()
end

prognostic_variables(::Type{<:Barotropic}) = (:vor,)
default_concrete_model(::Type{Barotropic}) = BarotropicModel

parameters(model::Barotropic; kwargs...) = SpeedyParams(
    planet = parameters(model.planet; component = :planet, kwargs...),
    atmosphere = parameters(model.atmosphere; component = :atmosphere, kwargs...),
    forcing = parameters(model.forcing; component = :forcing, kwargs...),
    drag = parameters(model.drag; component = :drag, kwargs...),
)

"""
$(TYPEDSIGNATURES)
Calls all `initialize!` functions for most fields, representing components, of `model`,
except for `model.output` and `model.feedback` which are always called
at in `time_stepping!`."""
function initialize!(model::Barotropic; time::DateTime = DEFAULT_DATE)
    (; spectral_grid) = model

    spectral_grid.nlayers > 1 && @error "Only nlayers=1 supported for BarotropicModel, \
        SpectralGrid with nlayers=$(spectral_grid.nlayers) provided."

    # initialize components
    arch = model.architecture
    @maybe_jit arch initialize!(model.geometry, model)
    @maybe_jit arch initialize!(model.time_stepping, model)
    @maybe_jit arch initialize!(model.coriolis, model)
    @maybe_jit arch initialize!(model.forcing, model)
    @maybe_jit arch initialize!(model.drag, model)
    @maybe_jit arch initialize!(model.horizontal_diffusion, model)
    @maybe_jit arch initialize!(model.random_process, model)
    @maybe_jit arch initialize!(model.particle_advection, model)

    # allocate prognostic and diagnostic variables
    prognostic_variables = PrognosticVariables(model)
    diagnostic_variables = DiagnosticVariables(model)
    # initialize particles (or other non-atmosphere prognostic variables)
    initialize!(prognostic_variables.particles, prognostic_variables, diagnostic_variables, model)

    # set the initial conditions
    @maybe_jit arch initialize!(prognostic_variables, model.initial_conditions, model)
    (; clock) = prognostic_variables
    clock.time = time       # set the current time
    clock.start = time      # and store the start time

    return Simulation(prognostic_variables, diagnostic_variables, model)
end
