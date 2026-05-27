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
end

"""$(TYPEDSIGNATURES)
Outer constructor that intercepts the deprecated `output`, `callbacks` and
`feedback` kwargs. They now live on `Simulation` (pass to `initialize!(model, ...)`)
and emit a deprecation warning here for backward-compatibility."""
function BarotropicModel(
        spectral_grid::SpectralGrid;
        output = nothing,
        callbacks = nothing,
        feedback = nothing,
        kwargs...,
    )
    _warn_deprecated_model_kwargs(:BarotropicModel; output, callbacks, feedback)
    return BarotropicModel(; spectral_grid, kwargs...)
end

function variables(model::Barotropic)
    nsteps = get_prognostic_steps(model.time_stepping)
    return variables(typeof(model), nsteps)
end

"""($TYPEDSIGNATURES) All variables needed for the barotropic model itself (components excluded)."""
function variables(::Type{<:Barotropic}, nsteps)
    return (
        PrognosticVariable(:clock, ClockDim(), desc = "Clock", units = "s"),
        PrognosticVariable(:scale, ScalarDim(1), desc = "Scaling of vor and div in the dynamical core", units = "m"),
        PrognosticVariable(:vorticity, Spectral4D(nsteps), desc = "Relative vorticity", units = "1/s"),

        TendencyVariable(:vorticity, Spectral3D(), desc = "Tendency of relative vorticity", units = "1/s²"),
        TendencyVariable(:vorticity, Grid3D(), namespace = :grid, desc = "Tendency of relative vorticity on the grid", units = "1/s²"),
        TendencyVariable(:u, Grid3D(), namespace = :grid, desc = "Tendency of zonal wind on the grid", units = "m/s²"),
        TendencyVariable(:v, Grid3D(), namespace = :grid, desc = "Tendency of meridional wind on the grid", units = "m/s²"),

        GridVariable(:vorticity, Grid3D(), desc = "Relative vorticity", units = "1/s"),
        GridVariable(:u, Grid3D(), desc = "Zonal wind", units = "m/s"),
        GridVariable(:v, Grid3D(), desc = "Meridional wind", units = "m/s"),

        ScratchVariable(:a, Spectral3D(), desc = "Scratch array", units = "?"),
        ScratchVariable(:b, Spectral3D(), desc = "Scratch array", units = "?"),
    )
end

"""
$(TYPEDSIGNATURES)
Calls all `initialize!` functions for most fields, representing components, of `model`.
The output writers, callbacks and feedback live on the returned `Simulation` and are
initialized just before the first time step inside `time_stepping!`."""
function initialize!(
        model::Barotropic;
        time::DateTime = DEFAULT_DATE,
        output = (),
        callbacks::CALLBACK_DICT = CallbackDict(),
        feedback::AbstractFeedback = Feedback(),
    )
    (; spectral_grid) = model

    spectral_grid.nlayers > 1 && @error "Only nlayers=1 supported for BarotropicModel, \
        SpectralGrid with nlayers=$(spectral_grid.nlayers) provided."

    outputs = output_collection(output)

    # initialize components
    initialize!(model.geometry, model)
    initialize!(model.time_stepping, model; interval = combined_output_interval(outputs))
    initialize!(model.coriolis, model)
    initialize!(model.forcing, model)
    initialize!(model.drag, model)
    initialize!(model.horizontal_diffusion, model)
    initialize!(model.random_process, model)
    initialize!(model.particle_advection, model)

    # allocate all variables
    variables = Variables(model)

    # set the time first
    (; clock) = variables.prognostic
    set!(clock, time = time, start = time)

    # now set initial conditions
    initialize!(variables, model)

    return Simulation(variables, model, outputs, callbacks, feedback)
end
