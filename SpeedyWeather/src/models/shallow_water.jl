export ShallowWaterModel

"""The ShallowWaterModel contains all model components needed for the simulation
of the shallow water equations. To be constructed like

    model = ShallowWaterModel(spectral_grid; kwargs...)

with `spectral_grid::SpectralGrid` used to initalize all non-default components
passed on as keyword arguments, e.g. `planet = Earth(spectral_grid)`.
Fields, representing model components, are $(TYPEDFIELDS)"""
@parameterized @kwdef mutable struct ShallowWaterModel{
        SG,     # <:SpectralGrid
        AR,     # <:AbstractArchitecture,
        GE,     # <:AbstractGeometry,
        PL,     # <:AbstractPlanet,
        AT,     # <:AbstractAtmosphere,
        CO,     # <:AbstractCoriolis,
        OR,     # <:AbstractOrography,
        FR,     # <:AbstractForcing,
        DR,     # <:AbstractDrag,
        PA,     # <:AbstractParticleAdvection,
        IC,     # <:AbstractInitialConditions,
        RP,     # <:AbstractRandomProcess,
        TS,     # <:AbstractTimeStepper,
        ST,     # <:SpectralTransform{NF},
        IM,     # <:AbstractImplicit,
        HD,     # <:AbstractHorizontalDiffusion,
        OU,     # <:AbstractOutput,
        FB,     # <:AbstractFeedback,
    } <: ShallowWater

    spectral_grid::SG
    architecture::AR = spectral_grid.architecture

    # DYNAMICS
    @component geometry::GE = Geometry(spectral_grid)
    @component planet::PL = Earth(spectral_grid)
    @component atmosphere::AT = EarthDryAtmosphere(spectral_grid)
    @component coriolis::CO = Coriolis(spectral_grid)
    @component orography::OR = EarthOrography(spectral_grid)
    @component forcing::FR = nothing
    @component drag::DR = nothing
    @component particle_advection::PA = nothing
    @component initial_conditions::IC = InitialConditions(spectral_grid, ShallowWater)

    # VARIABLES
    random_process::RP = nothing
    tracers::TRACER_DICT = TRACER_DICT()

    # NUMERICS
    time_stepping::TS = Leapfrog(spectral_grid)
    spectral_transform::ST = SpectralTransform(spectral_grid)
    implicit::IM = ImplicitShallowWater(spectral_grid)
    horizontal_diffusion::HD = HyperDiffusion(spectral_grid)

    # OUTPUT
    output::OU = NetCDFOutput(spectral_grid, ShallowWater)
    callbacks::Dict{Symbol, AbstractCallback} = Dict{Symbol, AbstractCallback}()
    feedback::FB = Feedback()
end

"""($TYPEDSIGNATURES) All variables needed for the shallow water model itself (components excluded)."""
function variables(model::ShallowWater)
    nsteps = get_nsteps(model.time_stepping, model)
    pg = nsteps.prognostic_grid
    ps = nsteps.prognostic_spectral
    tg = nsteps.tendency_grid
    ts = nsteps.tendency_spectral
    return (
        variables(BarotropicModel, nsteps)...,
        PrognosticVariable(:divergence, SpectralXYZT(ps), desc = "Divergence", units = "1/s", fuse = :prognostic),
        PrognosticVariable(:η, SpectralXYT(ps), desc = "Interface displacement", units = "m", fuse = :prognostic),

        TendencyVariable(:divergence, SpectralXYZT(ts), desc = "Tendency of divergence", units = "1/s²"),
        TendencyVariable(:divergence, GridXYZT(tg), namespace = :grid, desc = "Tendency of divergence on the grid", units = "1/s²"),

        TendencyVariable(:η, SpectralXYT(ps), desc = "Tendency of interface displacement", units = "m/s"),
        TendencyVariable(:η, GridXYT(tg), namespace = :grid, desc = "Tendency of interface displacement on the grid", units = "m/s"),

        GridVariable(:divergence, GridXYZT(pg), desc = "Divergence", units = "1/s", fuse = :grid),
        GridVariable(:η, GridXYT(pg), desc = "Interface displacement", units = "m", fuse = :grid),
        DynamicsVariable(:geopotential, GridXYZ(), desc = "Geopotential", units = "m²/s²"),

        ScratchVariable(:a, GridXYZ(), desc = "Scratch array", namespace = :grid),
        ScratchVariable(:b, GridXYZ(), desc = "Scratch array", namespace = :grid),
    )
end

""" 
$(TYPEDSIGNATURES)
Calls all `initialize!` functions for most components (=fields) of `model`,
except for `model.output` and `model.feedback` which are always initialized
in `time_stepping!` and `model.implicit` which is done in `first_timesteps!`."""
function initialize!(model::ShallowWater; time::DateTime = DEFAULT_DATE)
    (; spectral_grid) = model

    spectral_grid.nlayers > 1 && @error "Only nlayers=1 supported for ShallowWaterModel, \
                                SpectralGrid with nlayers=$(spectral_grid.nlayers) provided."

    # initialize components
    initialize!(model.geometry, model)
    initialize!(model.time_stepping, model)
    initialize!(model.coriolis, model)
    initialize!(model.orography, model)
    initialize!(model.forcing, model)
    initialize!(model.drag, model)
    initialize!(model.horizontal_diffusion, model)
    # model.implicit is initialized in first_timesteps!
    initialize!(model.random_process, model)
    initialize!(model.particle_advection, model)

    # allocate all variables and set initial conditions
    variables = Variables(model)

    # set the time first
    (; clock) = variables.prognostic
    set!(clock, time = time, start = time)

    # now set initial conditions
    initialize!(variables, model)

    return Simulation(variables, model)
end
