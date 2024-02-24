export PrimitiveDryModel

"""
$(SIGNATURES)
The PrimitiveDryModel struct holds all other structs that contain precalculated constants,
whether scalars or arrays that do not change throughout model integration.
$(TYPEDFIELDS)"""
Base.@kwdef mutable struct PrimitiveDryModel{NF<:AbstractFloat, D<:AbstractDevice} <: PrimitiveDry
    spectral_grid::SpectralGrid = SpectralGrid()

    # DYNAMICS
    dynamics::Bool = true
    planet::AbstractPlanet = Earth()
    atmosphere::AbstractAtmosphere = EarthAtmosphere()
    coriolis::AbstractCoriolis = Coriolis(spectral_grid)
    initial_conditions::InitialConditions = ZonalWind()
    orography::AbstractOrography{NF} = EarthOrography(spectral_grid)
    adiabatic_conversion::AbstractAdiabaticConversion = AdiabaticConversion(spectral_grid)

    # BOUNDARY CONDITIONS
    land_sea_mask::AbstractLandSeaMask{NF} = LandSeaMask(spectral_grid)
    ocean::AbstractOcean{NF} = SeasonalOceanClimatology(spectral_grid)
    land::AbstractLand{NF} = SeasonalLandTemperature(spectral_grid)
    solar_zenith::AbstractZenith{NF} = WhichZenith(spectral_grid,planet)

    # PHYSICS/PARAMETERIZATIONS
    physics::Bool = true
    boundary_layer_drag::BoundaryLayerDrag{NF} = BulkRichardsonDrag(spectral_grid)
    temperature_relaxation::TemperatureRelaxation{NF} = NoTemperatureRelaxation(spectral_grid)
    static_energy_diffusion::VerticalDiffusion{NF} = NoVerticalDiffusion(spectral_grid)
    surface_thermodynamics::AbstractSurfaceThermodynamics{NF} = SurfaceThermodynamicsConstant(spectral_grid)
    surface_wind::AbstractSurfaceWind{NF} = SurfaceWind(spectral_grid)
    surface_heat_flux::AbstractSurfaceHeat{NF} = SurfaceSensibleHeat(spectral_grid)
    shortwave_radiation::AbstractShortwave{NF} = NoShortwave(spectral_grid)
    longwave_radiation::AbstractLongwave{NF} = UniformCooling(spectral_grid)

    # NUMERICS
    time_stepping::TimeStepper = Leapfrog(spectral_grid)
    spectral_transform::SpectralTransform{NF} = SpectralTransform(spectral_grid)
    horizontal_diffusion::HorizontalDiffusion = HyperDiffusion(spectral_grid)
    implicit::AbstractImplicit = ImplicitPrimitiveEquation(spectral_grid)
    vertical_advection::VerticalAdvection{NF} = CenteredVerticalAdvection(spectral_grid)
    
    # INTERNALS
    geometry::AbstractGeometry = Geometry(spectral_grid)
    constants::DynamicsConstants{NF} = DynamicsConstants(spectral_grid,planet,atmosphere,geometry)
    device_setup::DeviceSetup{D} = DeviceSetup(CPUDevice())

    # OUTPUT
    output::AbstractOutputWriter = OutputWriter(spectral_grid,PrimitiveDry)
    feedback::AbstractFeedback = Feedback()
end

has(::Type{<:PrimitiveDry}, var_name::Symbol) = var_name in (:vor, :div, :temp, :pres)
default_concrete_model(::Type{PrimitiveDry}) = PrimitiveDryModel

"""
$(TYPEDSIGNATURES)
Calls all `initialize!` functions for components of `model`,
except for `model.output` and `model.feedback` which are always called
at in `time_stepping!` and `model.implicit` which is done in `first_timesteps!`."""
function initialize!(model::PrimitiveDry;time::DateTime = DEFAULT_DATE)
    (;spectral_grid) = model

    # NUMERICS (implicit is initialized later)
    # slightly adjust model time step to be a convenient divisor of output timestep
    initialize!(model.time_stepping, model)
    initialize!(model.horizontal_diffusion, model)

    # DYNAMICS
    initialize!(model.coriolis, model)
    initialize!(model.gepotential, model)
    initialize!(model.adiabatic_conversion, model)

    # boundary conditionss
    initialize!(model.orography, model)
    initialize!(model.land_sea_mask, model)
    initialize!(model.ocean, model)
    initialize!(model.land, model)
    initialize!(model.solar_zenith, time, model)

    # parameterizations
    initialize!(model.boundary_layer_drag,model)
    initialize!(model.temperature_relaxation,model)
    initialize!(model.static_energy_diffusion,model)
    initialize!(model.shortwave_radiation,model)
    initialize!(model.longwave_radiation,model)

    # initial conditions
    prognostic_variables = PrognosticVariables(spectral_grid, model)
    initialize!(prognostic_variables, model.initial_conditions, model)
    (;clock) = prognostic_variables
    clock.time = time       # set the time

    # initialize ocean and land and synchronize clocks
    initialize!(prognostic_variables.ocean, clock.time, model)
    initialize!(prognostic_variables.land, clock.time, model)

    diagnostic_variables = DiagnosticVariables(spectral_grid,model)
    return Simulation(prognostic_variables,diagnostic_variables,model)
end