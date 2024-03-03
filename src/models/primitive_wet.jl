export PrimitiveWetModel

"""
$(SIGNATURES)
The PrimitiveDryModel struct holds all other structs that contain precalculated constants,
whether scalars or arrays that do not change throughout model integration.
$(TYPEDFIELDS)"""
Base.@kwdef mutable struct PrimitiveWetModel{
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
    SO<:AbstractSoil,
    VG<:AbstractVegetation,
    CC<:AbstractClausiusClapeyron,
    BL<:AbstractBoundaryLayer,
    TR<:AbstractTemperatureRelaxation,
    SE<:AbstractVerticalDiffusion,
    HU<:AbstractVerticalDiffusion,
    SUT<:AbstractSurfaceThermodynamics,
    SUW<:AbstractSurfaceWind,
    SH<:AbstractSurfaceHeat,
    EV<:AbstractEvaporation,
    LSC<:AbstractCondensation,
    CV<:AbstractConvection,
    SW<:AbstractShortwave,
    LW<:AbstractLongwave,
    TS<:AbstractTimeStepper,
    ST<:SpectralTransform{NF},
    IM<:AbstractImplicit,
    HD<:AbstractHorizontalDiffusion,
    VA<:AbstractVerticalAdvection,
    HF<:AbstractHoleFilling,
    GE<:AbstractGeometry,
    OW<:AbstractOutputWriter,
    FB<:AbstractFeedback,
} <: PrimitiveWet

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
    soil::SO = SeasonalSoilMoisture(spectral_grid)
    vegetation::VG = VegetationClimatology(spectral_grid)
    
    # PHYSICS/PARAMETERIZATIONS
    physics::Bool = true
    clausius_clapeyron::CC = ClausiusClapeyron(spectral_grid, atmosphere)
    boundary_layer_drag::BL = BulkRichardsonDrag(spectral_grid)
    temperature_relaxation::TR = NoTemperatureRelaxation(spectral_grid)
    static_energy_diffusion::SE = NoVerticalDiffusion(spectral_grid)
    humidity_diffusion::HU = NoVerticalDiffusion(spectral_grid)
    surface_thermodynamics::SUT = SurfaceThermodynamicsConstant(spectral_grid)
    surface_wind::SUW = SurfaceWind(spectral_grid)
    surface_heat_flux::SH = SurfaceSensibleHeat(spectral_grid)
    evaporation::EV = SurfaceEvaporation(spectral_grid)
    large_scale_condensation::LSC = ImplicitCondensation(spectral_grid)
    convection::CV = SimplifiedBettsMiller(spectral_grid)
    shortwave_radiation::SW = NoShortwave(spectral_grid)
    longwave_radiation::LW = UniformCooling(spectral_grid)
    
    # NUMERICS
    device_setup::DS = DeviceSetup(CPUDevice())
    time_stepping::TS = Leapfrog(spectral_grid)
    spectral_transform::ST = SpectralTransform(spectral_grid)
    implicit::IM = ImplicitPrimitiveEquation(spectral_grid)
    horizontal_diffusion::HD = HyperDiffusion(spectral_grid)
    vertical_advection::VA = CenteredVerticalAdvection(spectral_grid)
    hole_filling::HF = ClipNegatives()
    
    # OUTPUT
    output::OW = OutputWriter(spectral_grid, PrimitiveDry)
    callbacks::Dict{Symbol,AbstractCallback} = Dict{Symbol,AbstractCallback}()
    feedback::FB = Feedback()
end
 
has(::Type{<:PrimitiveWet}, var_name::Symbol) = var_name in (:vor, :div, :temp, :pres, :humid)
default_concrete_model(::Type{PrimitiveWet}) = PrimitiveWetModel
 
"""
$(TYPEDSIGNATURES)
Calls all `initialize!` functions for components of `model`,
except for `model.output` and `model.feedback` which are always called
at in `time_stepping!` and `model.implicit` which is done in `first_timesteps!`."""
function initialize!(model::PrimitiveWet;time::DateTime = DEFAULT_DATE)
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
    initialize!(model.soil, model)
    initialize!(model.vegetation, model)
    initialize!(model.solar_zenith, time, model)

    # parameterizations
    initialize!(model.boundary_layer_drag, model)
    initialize!(model.temperature_relaxation, model)
    initialize!(model.static_energy_diffusion, model)
    initialize!(model.humidity_diffusion, model)
    initialize!(model.large_scale_condensation, model)
    initialize!(model.convection, model)
    initialize!(model.shortwave_radiation, model)
    initialize!(model.longwave_radiation, model)

    # initial conditions
    prognostic_variables = PrognosticVariables(spectral_grid,model)
    initialize!(prognostic_variables,model.initial_conditions,model)
    (;clock) = prognostic_variables
    clock.time = time       # set the time

    # initialize ocean and land and synchronize clocks
    initialize!(prognostic_variables.ocean, clock.time, model)
    initialize!(prognostic_variables.land, clock.time, model)

    diagnostic_variables = DiagnosticVariables(spectral_grid,model)
    return Simulation(prognostic_variables,diagnostic_variables,model)
end