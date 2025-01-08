export Simulation

"""
$(TYPEDSIGNATURES)
Simulation is a container struct to be used with `run!(::Simulation)`.
It contains
$(TYPEDFIELDS)"""
struct Simulation{Model<:AbstractModel} <: AbstractSimulation{Model}
    "define the current state of the model"
    prognostic_variables::PrognosticVariables

    "contain the tendencies and auxiliary arrays to compute them"
    diagnostic_variables::DiagnosticVariables

    "all parameters, constant at runtime"
    model::Model
end

function Base.show(io::IO, S::AbstractSimulation)
    println(io, "Simulation{$(model_type(S.model))}")
    println(io, "├ prognostic_variables::PrognosticVariables{...}")
    println(io, "├ diagnostic_variables::DiagnosticVariables{...}")
    print(io,   "└ model::$(model_type(S.model)){...}")
end

export run!

"""
$(TYPEDSIGNATURES)
Run a SpeedyWeather.jl `simulation`. The `simulation.model` is assumed to be initialized."""
function run!(
    simulation::AbstractSimulation;
    period::Period = Day(10),
    output::Bool = false,
)
    (; clock) = simulation.prognostic_variables

    # CLOCK
    set_period!(clock, period)                          # set how long to integrate for
    initialize!(clock, simulation.model.time_stepping)  # store the start date, reset counter

    # OUTPUT
    simulation.model.output.active = output             # enable/disable output

    # run it, yeah!
    time_stepping!(simulation)
end
