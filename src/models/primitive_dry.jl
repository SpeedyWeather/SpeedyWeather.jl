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
        SG,     # <:SpectralGrid
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
        SC,     # <:AbstractSurfaceCondition,
        SM,     # <:AbstractSurfaceMomentumFlux,
        SH,     # <:AbstractSurfaceHeatFlux,
        CV,     # <:AbstractConvection,
        SW,     # <:AbstractShortwave,
        LW,     # <:AbstractLongwave,
        SP,     # <:AbstractStochasticPhysics,
        CP,     # <:AbstractParameterization,
        TS,     # <:AbstractTimeStepper,
        ST,     # <:SpectralTransform{NF},
        IM,     # <:AbstractImplicit,
        HD,     # <:AbstractHorizontalDiffusion,
        VA,     # <:AbstractVerticalAdvection,
        OU,     # <:AbstractOutput,
        FB,     # <:AbstractFeedback,
        TS1,    # <:Tuple{Symbol}
        TS2,    # <:Tuple{Symbol}
        TS3,    # <:Tuple{Symbol}
        PV,     # <:Val
    } <: PrimitiveDry

    spectral_grid::SG
    architecture::AR = spectral_grid.architecture

    # DYNAMICS
    dynamics::Bool = true
    geometry::GE = Geometry(spectral_grid)
    planet::PL = Earth(spectral_grid)
    atmosphere::AT = EarthDryAtmosphere(spectral_grid)
    coriolis::CO = Coriolis(spectral_grid)
    geopotential::GO = Geopotential(spectral_grid)
    adiabatic_conversion::AC = AdiabaticConversion(spectral_grid)
    particle_advection::PA = nothing
    initial_conditions::IC = InitialConditions(spectral_grid, PrimitiveDry)
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
    surface_condition::SC = SurfaceCondition(spectral_grid)
    surface_momentum_flux::SM = SurfaceMomentumFlux(spectral_grid)
    surface_heat_flux::SH = SurfaceHeatFlux(spectral_grid)
    convection::CV = BettsMillerDryConvection(spectral_grid)
    shortwave_radiation::SW = TransparentShortwave(spectral_grid)
    longwave_radiation::LW = OneBandGreyLongwave(spectral_grid)
    stochastic_physics::SP = nothing
    custom_parameterization::CP = nothing

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

    # Tuples with symbols or instances of all parameterizations and parameter functions
    # Used to initiliaze variables and for the column-based parameterizations
    model_parameters::TS1 = (
        :architecture, :time_stepping, :orography, :geopotential, :atmosphere,
        :planet, :geometry, :land_sea_mask,
    )
    parameterizations::TS2 = (  # mixing
        :vertical_diffusion, :convection,

        # radiation
        :albedo, :shortwave_radiation, :longwave_radiation,

        # surface fluxes
        :boundary_layer_drag, :surface_condition, :surface_momentum_flux, :surface_heat_flux,

        # perturbations
        :stochastic_physics,
    )
    extra_parameterizations::TS3 = (:solar_zenith, :land, :ocean, :sea_ice)

    # DERIVED
    # used to infer parameterizations at compile-time
    params::PV = Val(parameterizations)
end

prognostic_variables(::Type{<:PrimitiveDry}) = (:vor, :div, :temp, :pres)
default_concrete_model(::Type{PrimitiveDry}) = PrimitiveDryModel

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
    initialize!(model.surface_condition, model)
    initialize!(model.surface_momentum_flux, model)
    initialize!(model.surface_heat_flux, model)
    initialize!(model.stochastic_physics, model)
    initialize!(model.particle_advection, model)

    # allocate prognostic and diagnostic variables
    prognostic_variables = PrognosticVariables(model)
    diagnostic_variables = DiagnosticVariables(model)

    # initialize non-atmosphere prognostic variables
    (; particles, ocean, land) = prognostic_variables
    initialize!(particles, prognostic_variables, diagnostic_variables, model.particle_advection, model)
    initialize!(ocean, prognostic_variables, diagnostic_variables, model.ocean, model)
    initialize!(land, prognostic_variables, diagnostic_variables, model.land, model)

    # set the initial conditions (may overwrite variables set in initialize! ocean/land)
    initialize!(prognostic_variables, model.initial_conditions, model)
    (; clock) = prognostic_variables
    clock.time = time       # set the current time
    clock.start = time      # and store the start time

    # pack prognostic, diagnostic variables and model into a simulation
    return Simulation(prognostic_variables, diagnostic_variables, model)
end

function Adapt.adapt_structure(to, model::PrimitiveDryModel)
    adapt_fields = model.model_parameters
    return NamedTuple{adapt_fields}(
        adapt_structure(to, getfield(model, field)) for field in adapt_fields
    )
end 
