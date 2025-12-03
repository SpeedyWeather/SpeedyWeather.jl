export PrimitiveDryModel

"""
The PrimitiveDryModel contains all model components (themselves structs) needed for the
simulation of the primitive equations without humidity. To be constructed like

    model = PrimitiveDryModel(spectral_grid; kwargs...)

with `spectral_grid::SpectralGrid` used to initalize all non-default components
passed on as keyword arguments, e.g. `planet=Earth(spectral_grid)`. Fields, representing
model components, are
$(TYPEDFIELDS)"""
@kwdef mutable struct PrimitiveDryModel{
    AR,     # <:AbstractArchitecture,
    GE,     # <:AbstractGeometry,
    PL,     # <:AbstractPlanet,
    AT,     # <:AbstractAtmosphere,
    CO,     # <:AbstractCoriolis,
    GO,     # <:AbstractGeopotential,
    AC,     # <:AbstractAdiabaticConversion,
    PA,     # <:AbstractParticleAdvection,
    IC,     # <:AbstractInitialConditions,
    FR,     # <:AbstractForcing,
    DR,     # <:AbstractDrag,
    RP,     # <:AbstractRandomProcess,
    OR,     # <:AbstractOrography,
    LS,     # <:AbstractLandSeaMask,
    OC,     # <:AbstractOcean,
    SI,     # <:AbstractSeaIce,
    LA,     # <:AbstractLand,
    ZE,     # <:AbstractZenith,
    AL,     # <:AbstractAlbedo,
    BL,     # <:AbstractBoundaryLayer,
    VD,     # <:AbstractVerticalDiffusion,
    SUT,    # <:AbstractSurfaceThermodynamics,
    SUW,    # <:AbstractSurfaceWind,
    SH,     # <:AbstractSurfaceHeatFlux,
    CV,     # <:AbstractConvection,
    SW,     # <:AbstractShortwave,
    LW,     # <:AbstractLongwave,
    SP,     # <:AbstractStochasticPhysics,
    TS,     # <:AbstractTimeStepper,
    ST,     # <:SpectralTransform{NF},
    IM,     # <:AbstractImplicit,
    HD,     # <:AbstractHorizontalDiffusion,
    VA,     # <:AbstractVerticalAdvection,
    OU,     # <:AbstractOutput,
    FB,     # <:AbstractFeedback,
} <: PrimitiveDry

    spectral_grid::SpectralGrid
    architecture::AR = spectral_grid.architecture
    
    # DYNAMICS
    dynamics::Bool = true
    geometry::GE = Geometry(spectral_grid)
    planet::PL = Earth(spectral_grid)
    atmosphere::AT = EarthAtmosphere(spectral_grid)
    coriolis::CO = Coriolis(spectral_grid)
    geopotential::GO = Geopotential(spectral_grid)
    adiabatic_conversion::AC = AdiabaticConversion(spectral_grid)
    particle_advection::PA = nothing
    initial_conditions::IC = InitialConditions(PrimitiveDry)
    forcing::FR = nothing
    drag::DR = nothing

    # VARIABLES
    random_process::RP = nothing
    tracers::TRACER_DICT = TRACER_DICT()

    # BOUNDARY CONDITIONS
    orography::OR = EarthOrography(spectral_grid)
    land_sea_mask::LS = EarthLandSeaMask(spectral_grid)
    ocean::OC = SlabOcean(spectral_grid)
    sea_ice::SI = ThermodynamicSeaIce(spectral_grid)
    land::LA = DryLandModel(spectral_grid)
    solar_zenith::ZE = WhichZenith(spectral_grid, planet)
    albedo::AL = DefaultAlbedo(spectral_grid)
    
    # PHYSICS/PARAMETERIZATIONS
    physics::Bool = true
    boundary_layer_drag::BL = BulkRichardsonDrag(spectral_grid)
    vertical_diffusion::VD = BulkRichardsonDiffusion(spectral_grid)
    surface_thermodynamics::SUT = SurfaceThermodynamicsConstant(spectral_grid)
    surface_wind::SUW = SurfaceWind(spectral_grid)
    surface_heat_flux::SH = SurfaceHeatFlux(spectral_grid)
    convection::CV = DryBettsMiller(spectral_grid)
    shortwave_radiation::SW = TransparentShortwave(spectral_grid)
    longwave_radiation::LW = JeevanjeeRadiation(spectral_grid)
    stochastic_physics::SP = nothing
    
    # NUMERICS
    time_stepping::TS = Leapfrog(spectral_grid)
    spectral_transform::ST = SpectralTransform(spectral_grid)
    implicit::IM = ImplicitPrimitiveEquation(spectral_grid)
    horizontal_diffusion::HD = HyperDiffusion(spectral_grid)
    vertical_advection::VA = CenteredVerticalAdvection(spectral_grid)
    
    # OUTPUT
    output::OU = NetCDFOutput(spectral_grid, PrimitiveDry)
    callbacks::Dict{Symbol, AbstractCallback} = Dict{Symbol, AbstractCallback}()
    feedback::FB = Feedback()
end

prognostic_variables(::Type{<:PrimitiveDry}) = (:vor, :div, :temp, :pres)
default_concrete_model(::Type{PrimitiveDry}) = PrimitiveDryModel

parameters(model::PrimitiveDry; kwargs...) = SpeedyParams(
    planet = parameters(model.planet; component=:planet, kwargs...),
    atmosphere = parameters(model.atmosphere; component=:atmosphere, kwargs...),
    forcing = parameters(model.forcing; component=:forcing, kwargs...),
    drag = parameters(model.drag; component=:drag, kwargs...),
    boundary_layer_drag = parameters(model.boundary_layer_drag; component=:boundary_layer_drag, kwargs...),
    vertical_diffusion = parameters(model.vertical_diffusion; component=:vertical_diffusion, kwargs...),
    surface_thermodynamics = parameters(model.surface_thermodynamics; component=:surface_thermodynamics, kwargs...),
    surface_wind = parameters(model.surface_wind; component=:surface_wind, kwargs...),
    surface_heat_flux = parameters(model.surface_heat_flux; component=:surface_heat_flux, kwargs...),
    convection = parameters(model.convection; component=:convection, kwargs...),
    shortwave_radiation = parameters(model.shortwave_radiation; component=:shortwave_radiation, kwargs...),
    longwave_radiation = parameters(model.longwave_radiation; component=:longwave_radiation, kwargs...),
    albedo = parameters(model.albedo; component=:albedo, kwargs...),
    ocean = parameters(model.ocean; component=:ocean, kwargs...),
    sea_ice = parameters(model.sea_ice; component=:sea_ice, kwargs...),
    land = parameters(model.land; component=:land, kwargs...),
    solar_zenith = parameters(model.solar_zenith; component=:solar_zenith, kwargs...),
)

"""
$(TYPEDSIGNATURES)
Calls all `initialize!` functions for components of `model`,
except for `model.output` and `model.feedback` which are always called
at in `time_stepping!` and `model.implicit` which is done in `first_timesteps!`."""
function initialize!(model::PrimitiveDry; time::DateTime = DEFAULT_DATE)
    (; spectral_grid) = model

    # NUMERICS (implicit is initialized later)
    initialize!(model.geometry, model)
    initialize!(model.time_stepping, model)
    initialize!(model.horizontal_diffusion, model)

    # DYNAMICS
    initialize!(model.coriolis, model)
    initialize!(model.geopotential, model)
    initialize!(model.adiabatic_conversion, model)
    initialize!(model.random_process, model)
    initialize!(model.forcing, model)
    initialize!(model.drag, model)

    # boundary conditionss
    initialize!(model.orography, model)
    initialize!(model.land_sea_mask, model)
    initialize!(model.ocean, model)
    initialize!(model.sea_ice, model)
    initialize!(model.land, model)
    initialize!(model.solar_zenith, time, model)
    initialize!(model.albedo, model)

    # parameterizations
    initialize!(model.boundary_layer_drag, model)
    initialize!(model.vertical_diffusion, model)
    initialize!(model.convection, model)
    initialize!(model.shortwave_radiation, model)
    initialize!(model.longwave_radiation, model)
    initialize!(model.surface_thermodynamics, model)
    initialize!(model.surface_wind, model)
    initialize!(model.surface_heat_flux, model)
    initialize!(model.stochastic_physics, model)
    initialize!(model.particle_advection, model)

    # allocate prognostic and diagnostic variables
    prognostic_variables = PrognosticVariables(spectral_grid, model)
    diagnostic_variables = DiagnosticVariables(spectral_grid, model)
    
    # initialize non-atmosphere prognostic variables
    (; particles, ocean, land) = prognostic_variables
    initialize!(particles, prognostic_variables, diagnostic_variables, model)
    initialize!(ocean,     prognostic_variables, diagnostic_variables, model)
    initialize!(land,      prognostic_variables, diagnostic_variables, model)

    # set the initial conditions (may overwrite variables set in intialize! ocean/land)
    initialize!(prognostic_variables, model.initial_conditions, model)
    (; clock) = prognostic_variables
    clock.time = time       # set the current time
    clock.start = time      # and store the start time

    # pack prognostic, diagnostic variables and model into a simulation
    return Simulation(prognostic_variables, diagnostic_variables, model)
end