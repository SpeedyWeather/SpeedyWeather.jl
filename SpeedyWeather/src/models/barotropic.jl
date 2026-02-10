export BarotropicModel, initialize!

"""
The BarotropicModel contains all model components needed for the simulation of the
barotropic vorticity equations. To be constructed like

    model = BarotropicModel(spectral_grid; kwargs...)

with `spectral_grid::SpectralGrid` used to initalize all non-default components
passed on as keyword arguments, e.g. `planet=Earth(spectral_grid)`. Fields, representing
model components, are
$(TYPEDFIELDS)"""
@parameterized @kwdef mutable struct BarotropicModel{
        SG,     # <:SpectralGrid
        AR,     # <:AbstractArchitecture,
        GE,     # <:AbstractGeometry,
        PL,     # <:AbstractPlanet,
        AT,     # <:AbstractAtmosphere,
        CO,     # <:AbstractCoriolis,
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
    } <: Barotropic

    spectral_grid::SG
    architecture::AR = spectral_grid.architecture

    # DYNAMICS
    @component geometry::GE = Geometry(spectral_grid)
    @component planet::PL = Earth(spectral_grid)
    @component atmosphere::AT = EarthDryAtmosphere(spectral_grid)
    @component coriolis::CO = Coriolis(spectral_grid)
    @component forcing::FR = KolmogorovFlow(spectral_grid)
    @component drag::DR = LinearVorticityDrag(spectral_grid)
    @component particle_advection::PA = nothing
    @component initial_conditions::IC = InitialConditions(spectral_grid, Barotropic)

    # VARIABLES
    random_process::RP = nothing
    tracers::TRACER_DICT = TRACER_DICT()

    # NUMERICS
    time_stepping::TS = Leapfrog(spectral_grid)
    spectral_transform::ST = SpectralTransform(spectral_grid)
    implicit::IM = nothing
    horizontal_diffusion::HD = HyperDiffusion(spectral_grid)

    # OUTPUT
    output::OU = NetCDFOutput(spectral_grid, Barotropic)
    callbacks::Dict{Symbol, AbstractCallback} = Dict{Symbol, AbstractCallback}()
    feedback::FB = Feedback()
end

function variables(::Type{<:Barotropic})
    return (
        PrognosticVariable(:clock, ClockDim(), desc = "Clock", units = "s"),
        PrognosticVariable(:scale, ScalarDim(1), desc = "Scaling of vor and div in the dynamical core", units = "m"),
        PrognosticVariable(:vor, Spectral3D(), desc = "Relative vorticity", units = "1/s"),
        GridVariable(:vor, Grid3D(), desc = "Relative vorticity", units = "1/s"),
        GridVariable(:u, Grid3D(), desc = "Zonal wind", units = "m/s"),
        GridVariable(:v, Grid3D(), desc = "Meridional wind", units = "m/s"),
        TendencyVariable(:vor, Spectral3D(), desc = "Tendency of relative vorticity", units = "1/s²"),
        TendencyVariable(:vor, Grid3D(), namespace = :grid, desc = "Tendency of relative vorticity on the grid", units = "1/s²"),
        TendencyVariable(:u, Grid3D(), namespace = :grid, desc = "Tendency of zonal wind on the grid", units = "m/s²"),
        TendencyVariable(:v, Grid3D(), namespace = :grid, desc = "Tendency of meridional wind on the grid", units = "m/s²"),
        ScratchVariable(:a, Spectral3D(), desc = "Spectral3D scratch memory", units = "?"),
        ScratchVariable(:b, Spectral3D(), desc = "Spectral3D scratch memory", units = "?"),
    )
end

"""
$(TYPEDSIGNATURES)
Calls all `initialize!` functions for most fields, representing components, of `model`,
except for `model.output` and `model.feedback` which are always called
at in `time_stepping!`."""
function initialize!(model::Barotropic; time::DateTime = DEFAULT_DATE)
    (; spectral_grid) = model

    spectral_grid.nlayers > 1 && @error "Only nlayers=1 supported for BarotropicModel, \
        SpectralGrid with nlayers=$(spectral_grid.nlayers) provided."

    # initialize components
    initialize!(model.geometry, model)
    initialize!(model.time_stepping, model)
    initialize!(model.coriolis, model)
    initialize!(model.forcing, model)
    initialize!(model.drag, model)
    initialize!(model.horizontal_diffusion, model)
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
