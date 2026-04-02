export PrimitiveDryModel

"""
The PrimitiveDryModel contains all model components (themselves structs) needed for the
simulation of the primitive equations without humidity. To be constructed like

    model = PrimitiveDryModel(spectral_grid; kwargs...)

with `spectral_grid::SpectralGrid` used to initalize all non-default components
passed on as keyword arguments, e.g. `planet=Earth(spectral_grid)`. Fields, representing
model components, are
$(TYPEDFIELDS)"""
@parameterized @kwdef mutable struct PrimitiveDryModel{
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
        PV,     # <:Val
    } <: PrimitiveDry

    spectral_grid::SG
    architecture::AR = spectral_grid.architecture

    # DYNAMICS
    dynamics::Bool = true
    @component geometry::GE = Geometry(spectral_grid)
    @component planet::PL = Earth(spectral_grid)
    @component atmosphere::AT = EarthDryAtmosphere(spectral_grid)
    @component coriolis::CO = Coriolis(spectral_grid)
    @component geopotential::GO = Geopotential(spectral_grid)
    @component adiabatic_conversion::AC = AdiabaticConversion(spectral_grid)
    @component particle_advection::PA = nothing
    @component initial_conditions::IC = InitialConditions(spectral_grid, PrimitiveDry)
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
    @component land::LA = DryLandModel(spectral_grid)

    # PHYSICS/PARAMETERIZATIONS
    @component solar_zenith::ZE = WhichZenith(spectral_grid, planet)
    @component albedo::AL = OceanLandAlbedo(spectral_grid)
    @component boundary_layer_drag::BL = BulkRichardsonDrag(spectral_grid)
    @component vertical_diffusion::VD = BulkRichardsonDiffusion(spectral_grid)
    @component surface_condition::SC = SurfaceCondition(spectral_grid)
    @component surface_momentum_flux::SM = SurfaceMomentumFlux(spectral_grid)
    @component surface_heat_flux::SH = SurfaceHeatFlux(spectral_grid)
    @component convection::CV = BettsMillerDryConvection(spectral_grid)
    @component shortwave_radiation::SW = OneBandGreyShortwave(spectral_grid)
    @component longwave_radiation::LW = OneBandGreyLongwave(spectral_grid)
    @component stochastic_physics::SP = nothing
    @component custom_parameterization::CP = nothing

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

    # Tuples with symbols pointing at the core components needed to be adapted
    # for availability within the parameterizations
    core_components::TS1 = (
        :architecture, :time_stepping, :orography, :geopotential, :atmosphere,
        :planet, :geometry, :land_sea_mask,
    )

    parameterizations::TS2 = (
        # orbit or external forcing
        :solar_zenith,

        # mixing
        :vertical_diffusion, :convection,

        # radiation
        :albedo, :shortwave_radiation, :longwave_radiation,

        # surface fluxes
        :boundary_layer_drag, :surface_condition, :surface_momentum_flux, :surface_heat_flux,

        # perturbations
        :stochastic_physics,
    )

    # DERIVED
    # used to infer parameterizations at compile-time
    params::PV = Val(parameterizations)
end

function variables(model::PrimitiveDry)
    nsteps = get_prognostic_steps(model.time_stepping)
    return variables(typeof(model), nsteps)
end

"""($TYPEDSIGNATURES) All variables needed for the primitive dry model itself (components excluded)."""
function variables(::Type{<:PrimitiveDry}, nsteps)
    return (
        variables(BarotropicModel, nsteps)...,
        PrognosticVariable(:divergence, Spectral4D(nsteps), desc = "Divergence", units = "1/s"),
        PrognosticVariable(:temperature, Spectral4D(nsteps), desc = "Temperature", units = "K"),
        PrognosticVariable(:pressure, Spectral3D(nsteps), desc = "Logarithm of surface pressure", units = "log(Pa)"),

        GridVariable(:divergence, Grid3D(), desc = "Divergence", units = "1/s"),
        GridVariable(:temperature, Grid3D(), desc = "Temperature", units = "K"),
        GridVariable(:pressure, Grid2D(), desc = "Logarithm of surface pressure", units = ""),
        GridVariable(:divergence_prev, Grid3D(), desc = "Divergence at previous time step", units = "1/s"),
        GridVariable(:temperature_prev, Grid3D(), desc = "Temperature at previous time step", units = "K"),
        GridVariable(:pressure_prev, Grid2D(), desc = "Logarithm of surface pressure at previous time step", units = ""),
        GridVariable(:u_prev, Grid3D(), desc = "Zonal wind at previous time step", units = "m/s"),
        GridVariable(:v_prev, Grid3D(), desc = "Meridional wind at previous time step", units = "m/s"),

        TendencyVariable(:divergence, Spectral3D(), desc = "Tendency of divergence", units = "1/s²"),
        TendencyVariable(:temperature, Spectral3D(), desc = "Tendency of temperature", units = "K/s"),
        TendencyVariable(:pressure, Spectral2D(), desc = "Tendency of surface pressure", units = "log(Pa)/s"),
        TendencyVariable(:divergence, Grid3D(), namespace = :grid, desc = "Tendency of divergence on the grid", units = "1/s²"),
        TendencyVariable(:temperature, Grid3D(), namespace = :grid, desc = "Tendency of temperature on the grid", units = "K/s"),
        TendencyVariable(:pressure, Grid2D(), namespace = :grid, desc = "Tendency of surface pressure on the grid", units = "log(Pa)/s"),

        DynamicsVariable(:dpres_dx, Grid2D(), desc = "Zonal gradient of the logarithm of surface pressure"),
        DynamicsVariable(:dpres_dy, Grid2D(), desc = "Meridional gradient of the logarithm of surface pressure"),
        DynamicsVariable(:pres_flux, Grid3D(), desc = "Pressure gradient flux, (u, v) ⋅ ∇lnp_s"),
        DynamicsVariable(:virtual_temperature, Spectral3D(), desc = "Virtual temperature", units = "K"),
        DynamicsVariable(:u_mean_grid, Grid2D(), desc = "Vertically integrated zonal velocity", units = "m/s"),
        DynamicsVariable(:v_mean_grid, Grid2D(), desc = "Vertically integrated meridional velocity", units = "m/s"),
        DynamicsVariable(:div_mean_grid, Grid2D(), desc = "Vertically integrated divergence", units = "1/s"),
        DynamicsVariable(:div_mean, Spectral2D(), desc = "Vertically integrated divergence", units = "1/s"),
        DynamicsVariable(:div_sum_above, Grid3D(), desc = "Partially vertically integrated divergence, top to layer above", units = "1/s"),
        DynamicsVariable(:pres_flux_sum_above, Grid3D(), desc = "Partially vertically integrated pressure gradient flux, top to layer above"),
        DynamicsVariable(:w, Grid3D(), desc = "Vertical velocity, dσ/dt.", units = "1/s"),

        ScratchVariable(:a, Grid3D(), desc = "Scratch array", namespace = :grid),
        ScratchVariable(:b, Grid3D(), desc = "Scratch array", namespace = :grid),
        ScratchVariable(:a_2D, Spectral2D(), desc = "Scratch array"),
        ScratchVariable(:b_2D, Spectral2D(), desc = "Scratch array"),
        ScratchVariable(:a_2D, Grid2D(), desc = "Scratch array", namespace = :grid),
        ScratchVariable(:b_2D, Grid2D(), desc = "Scratch array", namespace = :grid),
    )
end

"""
$(TYPEDSIGNATURES)
Calls all `initialize!` functions for components of `model`,
except for `model.output` and `model.feedback` which are always called
at in `time_stepping!` and `model.implicit` which is done in `first_timesteps!`."""
function initialize!(model::PrimitiveDry; time::DateTime = DEFAULT_DATE)
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

    # allocate all variables
    variables = Variables(model)

    # set the time first
    (; clock) = variables.prognostic
    set!(clock, time = time, start = time)

    # set all initial conditions for the ocean, seaice, land then atmosphere
    initialize!(variables, model)

    return Simulation(variables, model)
end

"""$(TYPEDSIGNATURES)
A `model` is adapted to the GPU or CPU by wrapping some (but not all!)
of its fields (determined by `model.core_components`) into a NamedTuple.
Importantly, while accessing fields `model.field` still works as usual,
one cannot use multiple dispatch on the model as a whole, e.g. `::PrimitiveDry`
will not work on GPU-adapted models."""
function Adapt.adapt_structure(to, model::PrimitiveDryModel)
    adapt_fields = model.core_components
    return NamedTuple{adapt_fields}(
        adapt_structure(to, getfield(model, field)) for field in adapt_fields
    )
end
