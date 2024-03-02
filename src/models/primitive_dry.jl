export PrimitiveDryModel

"""
$(SIGNATURES)
The PrimitiveDryModel struct holds all other structs that contain precalculated constants,
whether scalars or arrays that do not change throughout model integration.
$(TYPEDFIELDS)"""
Base.@kwdef mutable struct PrimitiveDryModel{
    NF<:AbstractFloat,
    DS<:DeviceSetup,
    PL<:AbstractPlanet,
    AT<:AbstractAtmosphere,
    CO<:AbstractCoriolis,
    GO<:AbstractGeopotential,
    OR<:AbstractOrography,
    AC<:AbstractAdiabaticConversion,
    IC<:InitialConditions,
    LS<:AbstractLandSeaMask,
    OC<:AbstractOcean,
    LA<:AbstractLand,
    ZE<:AbstractZenith,
    BL<:AbstractBoundaryLayer,
    TR<:AbstractTemperatureRelaxation,
    SE<:AbstractVerticalDiffusion,
    SUT<:AbstractSurfaceThermodynamics,
    SUW<:AbstractSurfaceWind,
    SH<:AbstractSurfaceHeat,
    SW<:AbstractShortwave,
    LW<:AbstractLongwave,
    TS<:AbstractTimeStepper,
    ST<:SpectralTransform{NF},
    IM<:AbstractImplicit,
    HD<:AbstractHorizontalDiffusion,
    VA<:AbstractVerticalAdvection,
    GE<:AbstractGeometry,
    OW<:AbstractOutputWriter,
    FB<:AbstractFeedback,
} <: PrimitiveDry

    spectral_grid::SpectralGrid = SpectralGrid()
    geometry::GE = Geometry(spectral_grid)
    
    # DYNAMICS
    dynamics::Bool = true
    planet::PL = Earth(spectral_grid)
    atmosphere::AT = EarthAtmosphere(spectral_grid)
    coriolis::CO = Coriolis(spectral_grid)
    geopotential::GO = Geopotential(spectral_grid)
    adiabatic_conversion::AC = AdiabaticConversion(spectral_grid)
    initial_conditions::IC = ZonalWind()
    
    # BOUNDARY CONDITIONS
    orography::OR = EarthOrography(spectral_grid)
    land_sea_mask::LS = LandSeaMask(spectral_grid)
    ocean::OC = SeasonalOceanClimatology(spectral_grid)
    land::LA = SeasonalLandTemperature(spectral_grid)
    solar_zenith::ZE = WhichZenith(spectral_grid, planet)
    
    # PHYSICS/PARAMETERIZATIONS
    physics::Bool = true
    boundary_layer_drag::BL = BulkRichardsonDrag(spectral_grid)
    temperature_relaxation::TR = NoTemperatureRelaxation(spectral_grid)
    static_energy_diffusion::SE = NoVerticalDiffusion(spectral_grid)
    surface_thermodynamics::SUT = SurfaceThermodynamicsConstant(spectral_grid)
    surface_wind::SUW = SurfaceWind(spectral_grid)
    surface_heat_flux::SH = SurfaceSensibleHeat(spectral_grid)
    shortwave_radiation::SW = NoShortwave(spectral_grid)
    longwave_radiation::LW = UniformCooling(spectral_grid)
    
    # NUMERICS
    device_setup::DS = DeviceSetup(CPUDevice())
    time_stepping::TS = Leapfrog(spectral_grid)
    spectral_transform::ST = SpectralTransform(spectral_grid)
    implicit::IM = ImplicitPrimitiveEquation(spectral_grid)
    horizontal_diffusion::HD = HyperDiffusion(spectral_grid)
    vertical_advection::VA = CenteredVerticalAdvection(spectral_grid)
    
    # OUTPUT
    output::OW = OutputWriter(spectral_grid, PrimitiveDry)
    callbacks::Dict{Symbol,AbstractCallback} = Dict{Symbol,AbstractCallback}()
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
    (;spectral_grid) = model

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

    # parameterizations
    initialize!(model.boundary_layer_drag, model)
    initialize!(model.temperature_relaxation, model)
    initialize!(model.static_energy_diffusion, model)
    initialize!(model.shortwave_radiation, model)
    initialize!(model.longwave_radiation, model)

    # initial conditions
    prognostic_variables = PrognosticVariables(spectral_grid, model)
    initialize!(prognostic_variables, model.initial_conditions, model)
    (;clock) = prognostic_variables
    clock.time = time       # set the time

    # initialize ocean and land and synchronize clocks
    initialize!(prognostic_variables.ocean, clock.time, model)
    initialize!(prognostic_variables.land, clock.time, model)

    diagnostic_variables = DiagnosticVariables(spectral_grid, model)
    return Simulation(prognostic_variables, diagnostic_variables, model)
end