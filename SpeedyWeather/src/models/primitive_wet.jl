export PrimitiveWetModel

"""
The PrimitiveWetModel contains all model components (themselves structs) needed for the
simulation of the primitive equations with humidity. To be constructed like

    model = PrimitiveWetModel(spectral_grid; kwargs...)

with `spectral_grid::SpectralGrid` used to initalize all non-default components
passed on as keyword arguments, e.g. `planet=Earth(spectral_grid)`. Fields, representing
model components, are
$(TYPEDFIELDS)"""
@parameterized @kwdef mutable struct PrimitiveWetModel{
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
        HF,     # <:AbstractSurfaceHumidityFlux,
        LSC,    # <:AbstractCondensation,
        CV,     # <:AbstractConvection,
        SW,     # <:AbstractShortwave,
        LW,     # <:AbstractLongwave,
        SP,     # <:AbstractStochasticPhysics,
        CP,     # <:AbstractParameterization
        TS,     # <:AbstractTimeStepper,
        ST,     # <:SpectralTransform{NF},
        IM,     # <:AbstractImplicit,
        HD,     # <:AbstractHorizontalDiffusion,
        VA,     # <:AbstractVerticalAdvection,
        HO,     # <:AbstractHoleFilling,
        OU,     # <:AbstractOutput,
        FB,     # <:AbstractFeedback,
        TS1,    # <:Tuple{Symbol}
        TS2,    # <:Tuple{Symbol}
        PV,     # <:Val
    } <: PrimitiveWet

    spectral_grid::SG
    architecture::AR = spectral_grid.architecture

    # DYNAMICS
    dynamics::Bool = true
    @component geometry::GE = Geometry(spectral_grid)
    @component planet::PL = Earth(spectral_grid)
    @component atmosphere::AT = EarthAtmosphere(spectral_grid)
    @component coriolis::CO = Coriolis(spectral_grid)
    @component geopotential::GO = Geopotential(spectral_grid)
    @component adiabatic_conversion::AC = AdiabaticConversion(spectral_grid)
    @component particle_advection::PA = nothing
    @component initial_conditions::IC = InitialConditions(spectral_grid, PrimitiveWet)
    @component forcing::FR = nothing
    @component drag::DR = SpeedLimitDrag(spectral_grid)

    # VARIABLES
    random_process::RP = nothing
    tracers::TRACER_DICT = TRACER_DICT()

    # BOUNDARY CONDITIONS
    dynamics_only::Bool = false
    @component orography::OR = EarthOrography(spectral_grid)
    @component land_sea_mask::LS = EarthLandSeaMask(spectral_grid)
    @component ocean::OC = SlabOcean(spectral_grid)
    @component sea_ice::SI = ThermodynamicSeaIce(spectral_grid)
    @component land::LA = LandModel(spectral_grid)

    # PHYSICS/PARAMETERIZATIONS
    @component solar_zenith::ZE = WhichZenith(spectral_grid, planet)
    @component albedo::AL = OceanLandAlbedo(spectral_grid)
    @component boundary_layer_drag::BL = BulkRichardsonDrag(spectral_grid)
    @component vertical_diffusion::VD = BulkRichardsonDiffusion(spectral_grid)
    @component surface_condition::SC = SurfaceCondition(spectral_grid)
    @component surface_momentum_flux::SM = SurfaceMomentumFlux(spectral_grid)
    @component surface_heat_flux::SH = SurfaceHeatFlux(spectral_grid)
    @component surface_humidity_flux::HF = SurfaceHumidityFlux(spectral_grid)
    @component large_scale_condensation::LSC = ImplicitCondensation(spectral_grid)
    @component convection::CV = BettsMillerConvection(spectral_grid)
    @component shortwave_radiation::SW = OneBandShortwave(spectral_grid)
    @component longwave_radiation::LW = OneBandLongwave(spectral_grid)
    @component stochastic_physics::SP = nothing
    @component custom_parameterization::CP = nothing

    # NUMERICS
    time_stepping::TS = Leapfrog(spectral_grid)
    spectral_transform::ST = SpectralTransform(spectral_grid)
    implicit::IM = ImplicitPrimitiveEquation(spectral_grid)
    horizontal_diffusion::HD = HyperDiffusion(spectral_grid)
    vertical_advection::VA = CenteredVerticalAdvection(spectral_grid)
    hole_filling::HO = ClipNegatives(spectral_grid)

    # OUTPUT
    output::OU = NetCDFOutput(spectral_grid, PrimitiveWet)
    callbacks::Dict{Symbol, AbstractCallback} = Dict{Symbol, AbstractCallback}()
    feedback::FB = Feedback()

    # COMPONENTS
    # Tuples with symbols or instances of all parameterizations and parameter functions
    # Used to initiliaze variables and for the column-based parameterizations
    # also determine order in which parameterizations are called
    core_components::TS1 = (
        :architecture, :time_stepping, :orography, :geopotential, :atmosphere,
        :planet, :geometry, :land_sea_mask,
    )
    parameterizations::TS2 = (
        # external forcing
        :solar_zenith,

        # mixing and precipitation
        :vertical_diffusion, :large_scale_condensation, :convection,

        # radiation
        :albedo, :shortwave_radiation, :longwave_radiation,

        # surface fluxes
        :boundary_layer_drag, :surface_condition,
        :surface_momentum_flux, :surface_heat_flux, :surface_humidity_flux,

        # perturbations
        :stochastic_physics,
    )

    # DERIVED
    # used to infer parameterizations at compile-time
    params::PV = Val(parameterizations)
end

function variables(model::PrimitiveWet)
    nsteps = get_prognostic_steps(model.time_stepping)
    return variables(typeof(model), nsteps)
end

"""($TYPEDSIGNATURES) All variables needed for the primitive wet model itself (components excluded)."""
function variables(::Type{<:PrimitiveWet}, nsteps)
    return (
        variables(PrimitiveDry, nsteps)...,

        # Add humidity
        PrognosticVariable(:humid, Spectral4D(nsteps), desc = "Specific humidity", units = "kg/kg"),
        GridVariable(:humid, Grid3D(), desc = "Humidity", units = "kg/kg"),
        GridVariable(:humid_prev, Grid3D(), desc = "Specific humidity at previous time step", units = "kg/kg"),
        TendencyVariable(:humid, Spectral3D(), desc = "Tendency of specific humidity", units = "kg/kg/s"),
        TendencyVariable(:humid, Grid3D(), namespace = :grid, desc = "Tendency of specific humidity on the grid", units = "kg/kg/s"),
    )
end

"""
$(TYPEDSIGNATURES)
Calls all `initialize!` functions for components of `model`,
except for `model.output` and `model.feedback` which are always called
at in `time_stepping!` and `model.implicit` which is done in `first_timesteps!`."""
function initialize!(model::PrimitiveWet; time::DateTime = DEFAULT_DATE)
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

    # boundary conditions
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
    initialize!(model.large_scale_condensation, model)
    initialize!(model.convection, model)
    initialize!(model.shortwave_radiation, model)
    initialize!(model.longwave_radiation, model)
    initialize!(model.surface_condition, model)
    initialize!(model.surface_momentum_flux, model)
    initialize!(model.surface_heat_flux, model)
    initialize!(model.surface_humidity_flux, model)
    initialize!(model.stochastic_physics, model)
    initialize!(model.particle_advection, model)

    # allocate all variables
    variables = Variables(model)

    # set the time first
    (; clock) = variables.prognostic
    set!(clock, time = time, start = time)

    # set all initial conditions for the ocean, sea ice, land then atmosphere
    initialize!(variables, model)

    return Simulation(variables, model)
end

"""$(TYPEDSIGNATURES)
A `model` is adapted to the GPU or CPU by wrapping some (but not all!)
of its fields (determined by `model.core_components`) into a NamedTuple.
Importantly, while accessing fields `model.field` still works as usual,
one cannot use multiple dispatch on the model as a whole, e.g. `::PrimitiveDry`
will not work on GPU-adapted models."""
function Adapt.adapt_structure(to, model::PrimitiveWetModel)
    adapt_fields = model.core_components
    return NamedTuple{adapt_fields}(
        adapt_structure(to, getfield(model, field)) for field in adapt_fields
    )
end
