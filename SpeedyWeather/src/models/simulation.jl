export Simulation

"""$(TYPEDSIGNATURES)
Simulation is a container struct simply wrapping the variables (prognostic and diagnostic)
and the model. It contains $(TYPEDFIELDS)"""
struct Simulation{V, M <: AbstractModel} <: AbstractSimulation{M}
    "All variables"
    variables::V

    "All model components, containing parameters, constant at runtime"
    model::M
end

function Base.show(io::IO, S::AbstractSimulation)
    vsize = prettymemory(Base.summarysize(S.variables))
    Msize = prettymemory(Base.summarysize(S.model))
    Ssize = prettymemory(Base.summarysize(S))
    println(io, styled"{warning:Simulation}", "{Variables{...}, $(model_type(S.model)){...}} ", styled"({note:$Ssize})")
    println(io, "├ ", styled"{info:variables}" * ": Variables{...} " * styled"{note:($vsize)}")
    print(io, "└ ", styled"{info:model}" * ": $(model_type(S.model)){...} " * styled"{note:($Msize)}")
    return nothing
end

unpack(sim::AbstractSimulation) = (sim.variables, sim.model)

const DEFAULT_PERIOD = Day(10)
const DEFAULT_TIMESTEPS = -1    # -1 means unspecified = use the period kwarg

export run!

"""$(TYPEDSIGNATURES)
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
    return unicodeplot(simulation)                      # maybe unicodeplot?
end

# fallback to be extended when UnicodePlots extension is loaded
unicodeplot(x) = nothing

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
    (; variables, model) = simulation
    progn = variables.prognostic

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

    # OUTPUT, enable/disable output
    set!(simulation.model.output, active = output, reset_path = true)

    # SCALING: we use vorticity*radius, divergence*radius in the dynamical core
    scale!(variables, model.planet.radius)

    # OUTPUT INITIALISATION AND STORING INITIAL CONDITIONS + FEEDBACK
    # propagate spectral state to grid variables for initial condition output
    lf = model.time_stepping.first_step_euler ? 1 : 2       # use 2nd leapfrog index when restarting

    # raise a warning if starting with leapfrog but there's zero vorticity
    vor = get_step(progn.vor, lf)
    lf == 2 && all(vor .== 0) && @warn "Vorticity is zero on 2nd leapfrog index though you use it to calculate tendencies." *
        " You may wanted to continue with a leapfrog step without data for it in the 2nd step."

    transform!(variables, lf, model, initialize = true)
    haskey(progn, :particles) && initialize!(variables, progn.particles, model)
    initialize!(model.output, model.feedback, variables, model)
    initialize!(model.callbacks, variables, model)
    return simulation
end

"""$(TYPEDSIGNATURES)
Finalize a `simulation`. Finishes the progress meter, unscales variables,
finalizes the output, writes a restart file and finalizes callbacks."""
function finalize!(simulation::AbstractSimulation)
    (; variables, model) = simulation
    finalize!(model.feedback)                       # finish the progress meter, do first for benchmark accuracy
    unscale!(variables)                             # undo radius-scaling for vor, div from the dynamical core
    finalize!(model.output, simulation)             # possibly post-process output, then close netCDF file
    finalize!(model.callbacks, variables, model)    # any callbacks to finalize?
    return simulation
end
