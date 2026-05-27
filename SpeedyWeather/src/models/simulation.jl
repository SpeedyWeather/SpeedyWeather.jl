export Simulation

"""$(TYPEDSIGNATURES)
`Simulation` is a container struct wrapping the variables (prognostic and diagnostic),
the model, all output writers, callbacks and the runtime feedback. It contains
$(TYPEDFIELDS)"""
struct Simulation{V, M <: AbstractModel, O <: Tuple, F} <: AbstractSimulation{M}
    "All variables"
    variables::V

    "All model components, containing parameters, constant at runtime"
    model::M

    "Tuple of output writers (`AbstractOutput`) attached to this simulation"
    output::O

    "Dictionary of callbacks invoked every time step"
    callbacks::CALLBACK_DICT

    "Runtime feedback (progress meter, NaN detection)"
    feedback::F
end

"""$(TYPEDSIGNATURES)
Minimal constructor for an `AbstractSimulation` with no output writers, an empty
callback dictionary and a default `Feedback`. Used internally and for backward
compatibility, e.g. inside `output!` when only `variables` and `model` are known."""
Simulation(variables, model::AbstractModel) =
    Simulation(variables, model, (), CallbackDict(), Feedback())

function Base.show(io::IO, S::AbstractSimulation)
    vsize = prettymemory(Base.summarysize(S.variables))
    Msize = prettymemory(Base.summarysize(S.model))
    Ssize = prettymemory(Base.summarysize(S))
    println(io, styled"{warning:Simulation}", "{...} ", styled"{note:($Ssize)}")
    println(io, "├ ", styled"{info:variables}" * "::Variables{...} " * styled"{note:($vsize)}")
    println(io, "├ ", styled"{info:model}" * "::$(model_type(S.model)){...} " * styled"{note:($Msize)}")

    # output writers: show one line per attached writer, using only the short
    # (non-parameterized) type name — mirrors how AbstractModel show truncates
    # parameter chains. Falls back to "(none)" when no writer is attached.
    nout = length(S.output)
    if nout == 0
        println(io, "├ ", styled"{info:output}" * "::Tuple{} " * styled"{note:(none)}")
    else
        println(io, "├ ", styled"{info:output}" * "::Tuple of $nout writer(s)")
        for (i, writer) in enumerate(S.output)
            branch = i == nout ? "│ └" : "│ ├"
            short = split(string(typeof(writer)), "{", limit = 2)[1]
            status = writer.active ? "active" : "inactive"
            println(io, "$branch ", styled"{magenta:$short} " * styled"{note:($status)}")
        end
    end

    println(io, "├ ", styled"{info:callbacks}" * "::Dict of $(length(S.callbacks)) callback(s)")
    print(io, "└ ", styled"{info:feedback}" * "::$(nameof(typeof(S.feedback)))")
    return nothing
end

unpack(sim::AbstractSimulation) = (sim.variables, sim.model)

"""$(TYPEDSIGNATURES)
Normalize the `output` argument passed to `initialize!`/`Simulation` to a `Tuple` of
`AbstractOutput` writers. Accepts a single writer, a tuple of writers, or `nothing`.
When more than one writer is given, all writers must share the same `interval` —
the model time step is rounded to a divisor of this single interval, so allowing
mismatched intervals would silently desynchronise the writers' time axes."""
output_collection(::Nothing) = ()
output_collection(o::AbstractOutput) = (o,)
function output_collection(outputs::Tuple)
    length(outputs) <= 1 && return outputs
    first_interval = first(outputs).interval
    for (i, o) in enumerate(outputs)
        o.interval == first_interval || throw(ArgumentError(
            "All output writers must share the same `interval`. " *
                "Writer 1 has interval=$(first_interval) but writer $i has interval=$(o.interval). " *
                "Construct your writers with the same `interval` keyword."
        ))
    end
    return outputs
end

"""$(TYPEDSIGNATURES)
Return the output `interval` used to round the model time step to a divisor. All
attached writers share the same `interval` (enforced by `output_collection`).
Falls back to `DEFAULT_OUTPUT_INTERVAL` when no writers are attached."""
combined_output_interval(::Tuple{}) = Second(DEFAULT_OUTPUT_INTERVAL)
combined_output_interval(outputs::Tuple) = first(outputs).interval

"""$(TYPEDSIGNATURES)
Emit a deprecation warning if any of the deprecated kwargs `output`, `callbacks`,
or `feedback` is passed to a model constructor. These now live on `Simulation`
and are passed to `initialize!(model, output=..., callbacks=..., feedback=...)`."""
function _warn_deprecated_model_kwargs(model_name::Symbol; output, callbacks, feedback)
    isnothing(output) || Base.depwarn(
        "Passing `output` to $(model_name)(...) is deprecated; it is ignored. " *
            "Construct the output writer (e.g. `NetCDFOutput(model)`) and pass it to " *
            "`initialize!(model, output=...)` instead.",
        model_name,
    )
    isnothing(callbacks) || Base.depwarn(
        "Passing `callbacks` to $(model_name)(...) is deprecated; it is ignored. " *
            "Pass to `initialize!(model, callbacks=...)` instead.",
        model_name,
    )
    isnothing(feedback) || Base.depwarn(
        "Passing `feedback` to $(model_name)(...) is deprecated; it is ignored. " *
            "Pass to `initialize!(model, feedback=...)` instead.",
        model_name,
    )
    return nothing
end

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
Initializes a `simulation` (does not include the construction step from `model`).
Scales the variables, initializes the output, stores initial conditions, initializes
the progress meter feedback, callbacks and performs the first two initial time steps
to spin up the leapfrogging scheme."""
function initialize!(
        simulation::AbstractSimulation;
        period::Period = DEFAULT_PERIOD,
        steps::Int = DEFAULT_TIMESTEPS,
        output::Bool = false,
    )
    (; variables, model, callbacks) = simulation
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

    # OUTPUT, enable/disable every attached writer
    for writer in simulation.output
        set!(writer, active = output, reset_path = true)
    end

    # SCALING: we use vorticity*radius, divergence*radius in the dynamical core
    scale_prognostic!(variables, model.planet.radius)

    # OUTPUT INITIALISATION AND STORING INITIAL CONDITIONS + FEEDBACK
    # propagate spectral state to grid variables for initial condition output
    lf = model.time_stepping.first_step_euler ? 1 : 2       # use 2nd leapfrog index when restarting

    # raise a warning if starting with leapfrog but there's zero vorticity
    vor = get_step(progn.vorticity, lf)

    @trace if lf == 2 && all(vor .== 0)
        @warn "Vorticity is zero on 2nd leapfrog index though you use it to calculate tendencies." *
            " You may wanted to continue with a leapfrog step without data for it in the 2nd step."
    end

    # transform variables from spectral to grid (= set the diagnostic variables in the correct initial state)
    transform!(variables, lf, model, initialize = true)
    haskey(progn, :particles) && initialize!(variables, progn.particles, model)     # initialize particle work arrays

    # only initialize output and callbacks just before the simulation starts
    for writer in simulation.output
        initialize!(writer, variables, model, callbacks)
    end
    initialize!(callbacks, simulation)
    return simulation
end

"""$(TYPEDSIGNATURES)
Finalize a `simulation`. Finishes the progress meter, unscales variables,
finalizes the output, writes a restart file and finalizes callbacks."""
function finalize!(simulation::AbstractSimulation)
    (; variables, callbacks, feedback) = simulation
    finalize!(feedback)                             # finish the progress meter, do first for benchmark accuracy
    unscale!(variables)                             # undo radius-scaling for vor, div in the dynamical core
    for writer in simulation.output
        finalize!(writer, simulation)               # possibly post-process output, then close netCDF/Zarr file
    end
    finalize!(callbacks, simulation)                # any callbacks to finalize?
    return simulation
end
