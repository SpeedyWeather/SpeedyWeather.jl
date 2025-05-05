export Simulation

"""
$(TYPEDSIGNATURES)
Simulation is a container struct to be used with `run!(::Simulation)`.
It contains
$(TYPEDFIELDS)"""
struct Simulation{Model<:AbstractModel, Progn, Diagn} <: AbstractSimulation{Model}
    "define the current state of the model"
    prognostic_variables::Progn

    "contain the tendencies and auxiliary arrays to compute them"
    diagnostic_variables::Diagn

    "all parameters, constant at runtime"
    model::Model
end

function Base.show(io::IO, S::AbstractSimulation)
    println(io, "Simulation{$(model_type(S.model))}")
    println(io, "├ prognostic_variables::PrognosticVariables{...}")
    println(io, "├ diagnostic_variables::DiagnosticVariables{...}")
    print(io,   "└ model::$(model_type(S.model)){...}")
end

unpack(sim::AbstractSimulation) = (sim.prognostic_variables, sim.diagnostic_variables, sim.model)

const DEFAULT_PERIOD = Day(10)
const DEFAULT_TIMESTEPS = -1    # -1 means unspecified = use the period kwarg

export run!

"""
$(TYPEDSIGNATURES)
Run a SpeedyWeather.jl `simulation`. The `simulation.model` is assumed to be initialized."""
function run!(
    simulation::AbstractSimulation;
    period::Period = DEFAULT_PERIOD,
    steps::Int = DEFAULT_TIMESTEPS,
    output::Bool = false,
)
    initialize!(simulation; period, steps, output)      # scaling, initialize output, store initial conditions
    time_stepping!(simulation)                          # run it, yeah!
    finalize!(simulation)                               # unscale, finalize output, write restart file, finalize callbacks             

    # return a UnicodePlot of surface vorticity
    surface_vorticity = simulation.diagnostic_variables.grid.vor_grid[:, end]
    return plot(surface_vorticity, title="Surface relative vorticity [1/s]")
end

"""$(TYPEDSIGNATURES)
Initializes a `simulation`. Scales the variables, initializes
the output, stores initial conditions, initializes the progress meter feedback,
callbacks and performs the first two initial time steps to spin up the
leapfrogging scheme."""
function initialize!(
    simulation::AbstractSimulation;
    period::Period = DEFAULT_PERIOD,
    steps::Int = DEFAULT_TIMESTEPS,
    output::Bool = false,
)
    progn, diagn, model = unpack(simulation)

    # SET THE CLOCK
    (; clock) = progn
    (; time_stepping) = model
    if steps != DEFAULT_TIMESTEPS
        # sets the steps, calculate period from it, store the start date, reset counter
        @assert period == DEFAULT_PERIOD "Period and steps cannot be set simultaneously"
        initialize!(clock, time_stepping, steps)    
    else
        # set period = how long to integrate for, tore the start date, reset counter
        initialize!(clock, time_stepping, period)
    end

    # OUTPUT
    simulation.model.output.active = output                     # enable/disable output

    # SCALING: we use vorticity*radius, divergence*radius in the dynamical core
    scale!(progn, diagn, model.spectral_grid.radius)

    # OUTPUT INITIALISATION AND STORING INITIAL CONDITIONS + FEEDBACK
    # propagate spectral state to grid variables for initial condition output
    lf = 1                                  # use first leapfrog index
    transform!(diagn, progn, lf, model, initialize=true)
    initialize!(progn.particles, progn, diagn, model.particle_advection)
    initialize!(model.output, model.feedback, progn, diagn, model)
    initialize!(model.callbacks, progn, diagn, model)
end

"""$(TYPEDSIGNATURES)
Finalize a `simulation`. Finishes the progress meter, unscales variables,
finalizes the output, writes a restart file and finalizes callbacks."""
function finalize!(simulation::AbstractSimulation)
    progn, diagn, model = unpack(simulation)

    finalize!(model.feedback)                       # finish the progress meter, do first for benchmark accuracy
    unscale!(progn)                                 # undo radius-scaling for vor, div from the dynamical core
    unscale!(diagn)                                 # undo radius-scaling for vor, div from the dynamical core
    finalize!(model.output, simulation)             # possibly post-process output, then close netCDF file
    write_restart_file!(model.output, progn)         # as JLD2 
    finalize!(model.callbacks, progn, diagn, model) # any callbacks to finalize?
end