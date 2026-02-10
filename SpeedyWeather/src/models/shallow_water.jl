export ShallowWaterModel

"""
The ShallowWaterModel contains all model components needed for the simulation of the
shallow water equations. To be constructed like

    model = ShallowWaterModel(spectral_grid; kwargs...)

with `spectral_grid::SpectralGrid` used to initalize all non-default components
passed on as keyword arguments, e.g. `planet=Earth(spectral_grid)`. Fields, representing
model components, are
$(TYPEDFIELDS)"""
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

function variables(::Type{<:ShallowWater})
    return (variables(BarotropicModel)...,
        PrognosticVariable(:div, Spectral3D(), desc = "Divergence", units = "1/s"),
        PrognosticVariable(:pres, Spectral2D(), desc = "Interface displacement", units = "m"),
        GridVariable(:div, Grid3D(), desc = "Divergence", units = "1/s"),
        GridVariable(:pres, Grid2D(), desc = "Interface displacement", units = "m"),
        TendencyVariable(:div, Spectral3D(), desc = "Tendency of divergence", units = "1/s²"),
        TendencyVariable(:pres, Spectral2D(), desc = "Tendency of interface displacement", units = "m/s"),
        TendencyVariable(:div, Grid3D(), namespace = :grid, desc = "Tendency of divergence on the grid", units = "1/s²"),
        TendencyVariable(:pres, Grid2D(), namespace = :grid, desc = "Tendency of interface displacement on the grid", units = "m/s"),
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
    # initialize!(variables.prognostic, model.initial_conditions, model)

    # set the time
    (; clock) = variables.prognostic
    clock.time = time       # set the current time
    clock.start = time      # and store the start time

    return Simulation(variables, model)
end
