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
    
    spectral_grid::SpectralGrid
    architecture::AR = spectral_grid.architecture
    
    # DYNAMICS
    geometry::GE = Geometry(spectral_grid)
    @component planet::PL = Earth(spectral_grid)
    @component atmosphere::AT = EarthAtmosphere(spectral_grid)
    @component coriolis::CO = Coriolis(spectral_grid)
    @component forcing::FR = KolmogorovFlow(spectral_grid)
    @component drag::DR = LinearVorticityDrag(spectral_grid)
    @component particle_advection::PA = NoParticleAdvection()
    @component initial_conditions::IC = InitialConditions(Barotropic)
    
    # VARIABLES
    @component random_process::RP = NoRandomProcess(spectral_grid)
    tracers::TRACER_DICT = TRACER_DICT()

    # NUMERICS
    time_stepping::TS = Leapfrog(spectral_grid)
    spectral_transform::ST = SpectralTransform(spectral_grid)
    implicit::IM = NoImplicit(spectral_grid)
    horizontal_diffusion::HD = HyperDiffusion(spectral_grid)

    # OUTPUT
    output::OU = NetCDFOutput(spectral_grid, Barotropic)
    callbacks::Dict{Symbol, AbstractCallback} = Dict{Symbol, AbstractCallback}()
    feedback::FB = Feedback()
end

prognostic_variables(::Type{<:Barotropic}) = (:vor,)
default_concrete_model(::Type{Barotropic}) = BarotropicModel

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
    initialize!(model.time_stepping, model)
    initialize!(model.coriolis, model)
    initialize!(model.forcing, model)
    initialize!(model.drag, model)
    initialize!(model.horizontal_diffusion, model)
    initialize!(model.random_process, model)

    # initial conditions
    prognostic_variables = PrognosticVariables(spectral_grid, model)
    initialize!(prognostic_variables, model.initial_conditions, model)
    prognostic_variables.clock.time = time       # set the current time
    prognostic_variables.clock.start = time      # and store the start time

    # particle advection
    initialize!(model.particle_advection, model)
    initialize!(prognostic_variables.particles, model)

    diagnostic_variables = DiagnosticVariables(spectral_grid, model)
    return Simulation(prognostic_variables, diagnostic_variables, model)
end