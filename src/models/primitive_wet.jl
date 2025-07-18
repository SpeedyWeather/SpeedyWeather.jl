export PrimitiveWetModel

"""
The PrimitiveWetModel contains all model components (themselves structs) needed for the
simulation of the primitive equations with humidity. To be constructed like

    model = PrimitiveWetModel(spectral_grid; kwargs...)

with `spectral_grid::SpectralGrid` used to initalize all non-default components
passed on as keyword arguments, e.g. `planet=Earth(spectral_grid)`. Fields, representing
model components, are
$(TYPEDFIELDS)"""
@kwdef mutable struct PrimitiveWetModel{
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
    LA,     # <:AbstractLand,
    ZE,     # <:AbstractZenith,
    AL,     # <:AbstractAlbedo,
    CC,     # <:AbstractClausiusClapeyron,
    BL,     # <:AbstractBoundaryLayer,
    TR,     # <:AbstractTemperatureRelaxation,
    VD,     # <:AbstractVerticalDiffusion,
    SUT,    # <:AbstractSurfaceThermodynamics,
    SUW,    # <:AbstractSurfaceWind,
    SH,     # <:AbstractSurfaceHeatFlux,
    EV,     # <:AbstractSurfaceEvaporation,
    LSC,    # <:AbstractCondensation,
    CV,     # <:AbstractConvection,
    OD,     # <:AbstractOpticalDepth,
    SW,     # <:AbstractShortwave,
    LW,     # <:AbstractLongwave,
    SP,     # <:AbstractStochasticPhysics,
    TS,     # <:AbstractTimeStepper,
    ST,     # <:SpectralTransform{NF},
    IM,     # <:AbstractImplicit,
    HD,     # <:AbstractHorizontalDiffusion,
    VA,     # <:AbstractVerticalAdvection,
    HF,     # <:AbstractHoleFilling,
    OU,     # <:AbstractOutput,
    FB,     # <:AbstractFeedback,
} <: PrimitiveWet

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
    particle_advection::PA = NoParticleAdvection()
    initial_conditions::IC = InitialConditions(PrimitiveWet)
    forcing::FR = NoForcing()
    drag::DR = NoDrag()  
  
    # VARIABLES
    random_process::RP = NoRandomProcess()
    tracers::TRACER_DICT = TRACER_DICT()
    
    # BOUNDARY CONDITIONS
    orography::OR = EarthOrography(spectral_grid)
    land_sea_mask::LS = EarthLandSeaMask(spectral_grid)
    ocean::OC = SeasonalOceanClimatology(spectral_grid)
    land::LA = LandModel(spectral_grid)
    solar_zenith::ZE = WhichZenith(spectral_grid, planet)
    albedo::AL = DefaultAlbedo(spectral_grid)
    
    # PHYSICS/PARAMETERIZATIONS
    physics::Bool = true
    clausius_clapeyron::CC = ClausiusClapeyron(spectral_grid, atmosphere)
    boundary_layer_drag::BL = BulkRichardsonDrag(spectral_grid)
    temperature_relaxation::TR = NoTemperatureRelaxation(spectral_grid)
    vertical_diffusion::VD = BulkRichardsonDiffusion(spectral_grid)
    surface_thermodynamics::SUT = SurfaceThermodynamicsConstant(spectral_grid)
    surface_wind::SUW = SurfaceWind(spectral_grid)
    surface_heat_flux::SH = SurfaceHeatFlux(spectral_grid)
    surface_evaporation::EV = SurfaceEvaporation(spectral_grid)
    large_scale_condensation::LSC = ImplicitCondensation(spectral_grid)
    convection::CV = SimplifiedBettsMiller(spectral_grid)
    optical_depth::OD = ZeroOpticalDepth(spectral_grid)
    shortwave_radiation::SW = TransparentShortwave(spectral_grid)
    longwave_radiation::LW = JeevanjeeRadiation(spectral_grid)
    stochastic_physics::SP = NoStochasticPhysics()
    
    # NUMERICS
    time_stepping::TS = Leapfrog(spectral_grid)
    spectral_transform::ST = SpectralTransform(spectral_grid)
    implicit::IM = ImplicitPrimitiveEquation(spectral_grid)
    horizontal_diffusion::HD = HyperDiffusion(spectral_grid)
    vertical_advection::VA = CenteredVerticalAdvection(spectral_grid)
    hole_filling::HF = ClipNegatives(spectral_grid)
    
    # OUTPUT
    output::OU = NetCDFOutput(spectral_grid, PrimitiveWet)
    callbacks::Dict{Symbol, AbstractCallback} = Dict{Symbol, AbstractCallback}()
    feedback::FB = Feedback()
end

prognostic_variables(::Type{<:PrimitiveWet}) = (:vor, :div, :temp, :humid, :pres)
default_concrete_model(::Type{PrimitiveWet}) = PrimitiveWetModel
 
"""
$(TYPEDSIGNATURES)
Calls all `initialize!` functions for components of `model`,
except for `model.output` and `model.feedback` which are always called
at in `time_stepping!` and `model.implicit` which is done in `first_timesteps!`."""
function initialize!(model::PrimitiveWet; time::DateTime = DEFAULT_DATE)
    (; spectral_grid) = model

    # NUMERICS (implicit is initialized later)
    initialize!(model.time_stepping, model)
    initialize!(model.horizontal_diffusion, model)

    # DYNAMICS
    initialize!(model.coriolis, model)
    initialize!(model.geopotential, model)
    initialize!(model.adiabatic_conversion, model)
    initialize!(model.random_process, model)
    initialize!(model.forcing, model)
    initialize!(model.drag, model)

    # boundary conditions
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
    initialize!(model.large_scale_condensation, model)
    initialize!(model.convection, model)
    initialize!(model.optical_depth, model)
    initialize!(model.shortwave_radiation, model)
    initialize!(model.longwave_radiation, model)
    initialize!(model.surface_thermodynamics, model)
    initialize!(model.surface_wind, model)
    initialize!(model.surface_heat_flux, model)
    initialize!(model.surface_evaporation, model)
    initialize!(model.stochastic_physics, model)

    # allocate prognostic and diagnostic variables
    prognostic_variables = PrognosticVariables(spectral_grid, model)
    diagnostic_variables = DiagnosticVariables(spectral_grid, model)
    
    # particle advection
    initialize!(model.particle_advection, model)
    initialize!(prognostic_variables.particles, model)
    
    # initialize ocean and land
    initialize!(prognostic_variables.ocean, prognostic_variables, diagnostic_variables, model)
    initialize!(prognostic_variables.land,  prognostic_variables, diagnostic_variables, model)

    # set the initial conditions (may overwrite variables set in intialize! ocean/land)
    initialize!(prognostic_variables, model.initial_conditions, model)
    (; clock) = prognostic_variables
    clock.time = time       # set the current time
    clock.start = time      # and store the start time

    # pack prognostic, diagnostic variables and model into a simulation
    return Simulation(prognostic_variables, diagnostic_variables, model)
end