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
        GHG,    # NamedTuple of <:AbstractGreenhouseGas,
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
    @component boundary_layer::BL = BoundaryLayer(spectral_grid)
    @component vertical_diffusion::VD = BulkRichardsonDiffusion(spectral_grid)
    @component surface_momentum_flux::SM = SurfaceMomentumFlux(spectral_grid)
    @component surface_heat_flux::SH = SurfaceHeatFlux(spectral_grid)
    @component convection::CV = BettsMillerDryConvection(spectral_grid)
    @component shortwave_radiation::SW = OneBandGreyShortwave(spectral_grid)
    @component longwave_radiation::LW = OneBandGreyLongwave(spectral_grid)
    @component greenhouse_gases::GHG = (;)
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
        :solar_zenith,              # orbit or external forcing
        :vertical_diffusion,        # mixing
        :convection,
        :albedo,                    # radiation
        :shortwave_radiation,
        :longwave_radiation,
        :boundary_layer,            # surface fluxes
        :surface_momentum_flux,
        :surface_heat_flux,
        :stochastic_physics,        # perturbations
    )

    # DERIVED
    # used to infer parameterizations at compile-time
    params::PV = Val(parameterizations)
end

variables(model::PrimitiveDry) = variables(typeof(model), get_nsteps(model.time_stepping, model))

"""($TYPEDSIGNATURES) All variables needed for the primitive dry model itself (components excluded)."""
function variables(::Type{<:PrimitiveDry}, nsteps = DEFAULT_NSTEPS)
    pg = nsteps.prognostic_grid
    ps = nsteps.prognostic_spectral
    tg = nsteps.tendency_grid
    ts = nsteps.tendency_spectral
    return (
        variables(BarotropicModel, nsteps)...,
        PrognosticVariable(:divergence, SpectralXYZT(ps), desc = "Divergence", units = "1/s", fuse = :prognostic),
        PrognosticVariable(:temperature, SpectralXYZT(ps), desc = "Temperature", units = "K", fuse = :prognostic),
        PrognosticVariable(:pressure, SpectralXYT(ps), desc = "Logarithm of surface pressure", units = "log(Pa)", fuse = :prognostic),

        TendencyVariable(:divergence, SpectralXYZT(ts), desc = "Tendency of divergence", units = "1/s²"), # not fused because computed directly by divergence op
        TendencyVariable(:temperature, SpectralXYZT(ts), desc = "Tendency of temperature", units = "K/s", fuse = :spectral_tendencies),
        TendencyVariable(:pressure, SpectralXYT(ts), desc = "Tendency of surface pressure", units = "log(Pa)/s", fuse = :spectral_tendencies),
        TendencyVariable(:divergence, GridXYZT(tg), namespace = :grid, desc = "Tendency of divergence on the grid", units = "1/s²"),
        TendencyVariable(:temperature, GridXYZT(tg), namespace = :grid, desc = "Tendency of temperature on the grid", units = "K/s", fuse = :grid_tendencies),
        TendencyVariable(:pressure, GridXYT(tg), namespace = :grid, desc = "Tendency of surface pressure on the grid", units = "log(Pa)/s", fuse = :grid_tendencies),
        
        GridVariable(:divergence, GridXYZT(pg), desc = "Divergence", units = "1/s", fuse=:grid),
        GridVariable(:temperature, GridXYZT(pg), desc = "Temperature", units = "K", fuse=:grid),
        GridVariable(:pressure, GridXYT(pg), desc = "Logarithm of surface pressure", units = "log(Pa)", fuse=:grid),
        ParameterizationVariable(:surface_pressure, Grid2D(), desc = "Surface pressure", units = "Pa"),

        DynamicsVariable(:uT_anomaly, GridXYZT(tg), desc = "u*T anomaly intermediate on grid", namespace = :grid, fuse = :grid_tendencies),
        DynamicsVariable(:vT_anomaly, GridXYZT(tg), desc = "v*T anomaly intermediate on grid", namespace = :grid, fuse = :grid_tendencies),
        DynamicsVariable(:uT_anomaly, SpectralXYZT(ts), desc = "u*T anomaly intermediate in spectral space", fuse = :spectral_tendencies),
        DynamicsVariable(:vT_anomaly, SpectralXYZT(ts), desc = "v*T anomaly intermediate in spectral space", fuse = :spectral_tendencies),

        DynamicsVariable(:kinetic_energy, GridXYZT(tg), desc = "Kinetic energy intermediate, ½(u²+v²)", namespace = :grid, fuse = :grid_tendencies),
        DynamicsVariable(:kinetic_energy, SpectralXYZT(ts), desc = "Kinetic energy intermediate in spectral space", fuse = :spectral_tendencies),

        DynamicsVariable(:dpres_dx, Grid2D(), desc = "Zonal gradient of the logarithm of surface pressure", fuse = :dpres_grad),
        DynamicsVariable(:dpres_dy, Grid2D(), desc = "Meridional gradient of the logarithm of surface pressure", fuse = :dpres_grad),
        DynamicsVariable(:dpres_dx_spec, Spectral2D(), desc = "Zonal gradient of lnpₛ in spectral space", fuse = :dpres_grad_spec),
        DynamicsVariable(:dpres_dy_spec, Spectral2D(), desc = "Meridional gradient of lnpₛ in spectral space", fuse = :dpres_grad_spec),
        DynamicsVariable(:pres_flux, GridXYZ(), desc = "Pressure gradient flux, (u, v) ⋅ ∇lnp_s"),
        DynamicsVariable(:virtual_temperature, SpectralXYZ(), desc = "Virtual temperature", units = "K"),
        DynamicsVariable(:u_mean_grid, Grid2D(), desc = "Vertically integrated zonal velocity", units = "m/s"),
        DynamicsVariable(:v_mean_grid, Grid2D(), desc = "Vertically integrated meridional velocity", units = "m/s"),
        DynamicsVariable(:div_mean_grid, Grid2D(), desc = "Vertically integrated divergence", units = "1/s"),
        DynamicsVariable(:div_mean, Spectral2D(), desc = "Vertically integrated divergence", units = "1/s"),
        DynamicsVariable(:div_sum_above, GridXYZ(), desc = "Partially vertically integrated divergence, top to layer above", units = "1/s"),
        DynamicsVariable(:pres_flux_sum_above, GridXYZ(), desc = "Partially vertically integrated pressure gradient flux, top to layer above"),
        DynamicsVariable(:w, GridXYZ(), desc = "Vertical velocity, dσ/dt.", units = "1/s"),

        ScratchVariable(:a, GridXYZ(), desc = "Scratch array", namespace = :grid),
        ScratchVariable(:b, GridXYZ(), desc = "Scratch array", namespace = :grid),
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
    initialize!(model.boundary_layer, model)
    initialize!(model.vertical_diffusion, model)
    initialize!(model.convection, model)
    initialize!(model.shortwave_radiation, model)
    initialize!(model.longwave_radiation, model)
    initialize!(model.greenhouse_gases, model)
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

function reinitialize!(model::PrimitiveDryModel, vars::AbstractVariables)
    reinitialize!(model.implicit, model, vars)
    return nothing
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
