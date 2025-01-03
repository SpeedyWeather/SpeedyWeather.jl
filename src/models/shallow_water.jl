export ShallowWaterModel

"""
The ShallowWaterModel contains all model components needed for the simulation of the
shallow water equations. To be constructed like

    model = ShallowWaterModel(spectral_grid; kwargs...)

with `spectral_grid::SpectralGrid` used to initalize all non-default components
passed on as keyword arguments, e.g. `planet=Earth(spectral_grid)`. Fields, representing
model components, are
$(TYPEDFIELDS)"""
@kwdef mutable struct ShallowWaterModel{
    DS,     # <:DeviceSetup,
    GE,     # <:AbstractGeometry,
    PL,     # <:AbstractPlanet,
    AT,     # <:AbstractAtmosphere,
    CO,     # <:AbstractCoriolis,
    OR,     # <:AbstractOrography,
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
} <: ShallowWater
    
    spectral_grid::SpectralGrid
    device_setup::DS = DeviceSetup(spectral_grid.device)

    # DYNAMICS
    geometry::GE = Geometry(spectral_grid)
    planet::PL = Earth(spectral_grid)
    atmosphere::AT = EarthAtmosphere(spectral_grid)
    coriolis::CO = Coriolis(spectral_grid)
    orography::OR = EarthOrography(spectral_grid)
    forcing::FR = NoForcing()
    drag::DR = NoDrag()
    particle_advection::PA = NoParticleAdvection()
    initial_conditions::IC = InitialConditions(ShallowWater)
    
    # VARIABLES
    random_process::RP = NoRandomProcess()
    tracers::TRACER_DICT = TRACER_DICT()

    # NUMERICS
    time_stepping::TS = Leapfrog(spectral_grid)
    spectral_transform::ST = SpectralTransform(spectral_grid)
    implicit::IM = ImplicitShallowWater(spectral_grid)
    horizontal_diffusion::HD = HyperDiffusion(spectral_grid)

    # OUTPUT
    output::OU = NetCDFOutput(spectral_grid, ShallowWater)
    callbacks::Dict{Symbol, AbstractCallback} = Dict{Symbol, AbstractCallback}()
    feedback::FB = Feedback()
end

prognostic_variables(::Type{<:ShallowWater}) = (:vor, :div, :pres)
default_concrete_model(::Type{ShallowWater}) = ShallowWaterModel

"""
$(TYPEDSIGNATURES)
Calls all `initialize!` functions for most components (=fields) of `model`,
except for `model.output` and `model.feedback` which are always initialized
in `time_stepping!` and `model.implicit` which is done in `first_timesteps!`."""
function initialize!(model::ShallowWater; time::DateTime = DEFAULT_DATE)
    (; spectral_grid) = model

    spectral_grid.nlayers > 1 && @error "Only nlayers=1 supported for ShallowWaterModel, \
                                SpectralGrid with nlayers=$(spectral_grid.nlayers) provided."

    # initialize components
    initialize!(model.time_stepping, model)
    initialize!(model.coriolis, model)
    initialize!(model.orography, model)
    initialize!(model.forcing, model)
    initialize!(model.drag, model)
    initialize!(model.horizontal_diffusion, model)
    # model.implicit is initialized in first_timesteps!
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