export PrimitiveDryModel

"""
The PrimitiveDryModel contains all model components (themselves structs) needed for the
simulation of the primitive equations without humidity. To be constructed like

    model = PrimitiveDryModel(; spectral_grid, kwargs...)

with `spectral_grid::SpectralGrid` used to initalize all non-default components
passed on as keyword arguments, e.g. `planet=Earth(spectral_grid)`. Fields, representing
model components, are
$(TYPEDFIELDS)"""
Base.@kwdef mutable struct PrimitiveDryModel{
    # TODO add constraints again when we stop supporting julia v1.9
    DS,     # <:DeviceSetup,
    GE,     # <:AbstractGeometry,
    PL,     # <:AbstractPlanet,
    AT,     # <:AbstractAtmosphere,
    CO,     # <:AbstractCoriolis,
    GO,     # <:AbstractGeopotential,
    OR,     # <:AbstractOrography,
    AC,     # <:AbstractAdiabaticConversion,
    PA,     # <:AbstractParticleAdvection,
    IC,     # <:AbstractInitialConditions,
    LS,     # <:AbstractLandSeaMask,
    OC,     # <:AbstractOcean,
    LA,     # <:AbstractLand,
    ZE,     # <:AbstractZenith,
    AL,     # <:AbstractAlbedo,
    BL,     # <:AbstractBoundaryLayer,
    TR,     # <:AbstractTemperatureRelaxation,
    VD,     # <:AbstractVerticalDiffusion,
    SUT,    # <:AbstractSurfaceThermodynamics,
    SUW,    # <:AbstractSurfaceWind,
    SH,     # <:AbstractSurfaceHeatFlux,
    CV,     # <:AbstractConvection,
    SW,     # <:AbstractShortwave,
    LW,     # <:AbstractLongwave,
    TS,     # <:AbstractTimeStepper,
    ST,     # <:SpectralTransform{NF},
    IM,     # <:AbstractImplicit,
    HD,     # <:AbstractHorizontalDiffusion,
    VA,     # <:AbstractVerticalAdvection,
    OW,     # <:AbstractOutputWriter,
    FB,     # <:AbstractFeedback,
} <: PrimitiveDry

    spectral_grid::SpectralGrid
    device_setup::DS = DeviceSetup(CPUDevice())
    
    # DYNAMICS
    dynamics::Bool = true
    geometry::GE = Geometry(spectral_grid)
    planet::PL = Earth(spectral_grid)
    atmosphere::AT = EarthAtmosphere(spectral_grid)
    coriolis::CO = Coriolis(spectral_grid)
    geopotential::GO = Geopotential(spectral_grid)
    adiabatic_conversion::AC = AdiabaticConversion(spectral_grid)
    particle_advection::PA = NoParticleAdvection()
    initial_conditions::IC = InitialConditions(PrimitiveDry)
    
    # BOUNDARY CONDITIONS
    orography::OR = EarthOrography(spectral_grid)
    land_sea_mask::LS = LandSeaMask(spectral_grid)
    ocean::OC = SeasonalOceanClimatology(spectral_grid)
    land::LA = SeasonalLandTemperature(spectral_grid)
    solar_zenith::ZE = WhichZenith(spectral_grid, planet)
    albedo::AL = AlbedoClimatology(spectral_grid)
    
    # PHYSICS/PARAMETERIZATIONS
    physics::Bool = true
    boundary_layer_drag::BL = BulkRichardsonDrag(spectral_grid)
    temperature_relaxation::TR = NoTemperatureRelaxation(spectral_grid)
    vertical_diffusion::VD = BulkRichardsonDiffusion(spectral_grid)
    surface_thermodynamics::SUT = SurfaceThermodynamicsConstant(spectral_grid)
    surface_wind::SUW = SurfaceWind(spectral_grid)
    surface_heat_flux::SH = SurfaceHeatFlux(spectral_grid)
    convection::CV = DryBettsMiller(spectral_grid)
    shortwave_radiation::SW = TransparentShortwave(spectral_grid)
    longwave_radiation::LW = JeevanjeeRadiation(spectral_grid)
    
    # NUMERICS
    time_stepping::TS = Leapfrog(spectral_grid)
    spectral_transform::ST = SpectralTransform(spectral_grid)
    implicit::IM = ImplicitPrimitiveEquation(spectral_grid)
    horizontal_diffusion::HD = HyperDiffusion(spectral_grid)
    vertical_advection::VA = CenteredVerticalAdvection(spectral_grid)
    
    # OUTPUT
    output::OW = OutputWriter(spectral_grid, PrimitiveDry)
    callbacks::Dict{Symbol, AbstractCallback} = Dict{Symbol, AbstractCallback}()
    feedback::FB = Feedback()
end

has(::Type{<:PrimitiveDry}, var_name::Symbol) = var_name in (:vor, :div, :temp, :pres)
default_concrete_model(::Type{PrimitiveDry}) = PrimitiveDryModel

"""
$(TYPEDSIGNATURES)
Calls all `initialize!` functions for components of `model`,
except for `model.output` and `model.feedback` which are always called
at in `time_stepping!` and `model.implicit` which is done in `first_timesteps!`."""
function initialize!(model::PrimitiveDry; time::DateTime = DEFAULT_DATE)
    (; spectral_grid) = model

    # NUMERICS (implicit is initialized later)
    initialize!(model.time_stepping, model)
    initialize!(model.horizontal_diffusion, model)

    # DYNAMICS
    initialize!(model.coriolis, model)
    initialize!(model.geopotential, model)
    initialize!(model.adiabatic_conversion, model)

    # boundary conditionss
    initialize!(model.orography, model)
    initialize!(model.land_sea_mask, model)
    initialize!(model.ocean, model)
    initialize!(model.land, model)
    initialize!(model.solar_zenith, time, model)
    initialize!(model.albedo, model)

    # parameterizations
    initialize!(model.boundary_layer_drag, model)
    initialize!(model.temperature_relaxation, model)
    initialize!(model.vertical_diffusion, model)
    initialize!(model.shortwave_radiation, model)
    initialize!(model.longwave_radiation, model)
    initialize!(model.surface_thermodynamics, model)
    initialize!(model.surface_wind, model)
    initialize!(model.surface_heat_flux, model)

    # initial conditions
    prognostic_variables = PrognosticVariables(spectral_grid, model)
    initialize!(prognostic_variables, model.initial_conditions, model)
    (; clock) = prognostic_variables
    clock.time = time       # set the current time
    clock.start = time      # and store the start time

    # particle advection
    initialize!(model.particle_advection, model)
    initialize!(prognostic_variables.particles, model)

    # initialize ocean and land and synchronize clocks
    initialize!(prognostic_variables.ocean, clock.time, model)
    initialize!(prognostic_variables.land, clock.time, model)

    diagnostic_variables = DiagnosticVariables(spectral_grid, model)
    return Simulation(prognostic_variables, diagnostic_variables, model)
end