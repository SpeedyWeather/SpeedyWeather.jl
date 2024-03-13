export Simulation

"""
$(TYPEDSIGNATURES)
Simulation is a container struct to be used with `run!(::Simulation)`.
It contains
$(TYPEDFIELDS)"""
struct Simulation{Model<:ModelSetup} <: AbstractSimulation{Model}
    "define the current state of the model"
    prognostic_variables::PrognosticVariables

    "contain the tendencies and auxiliary arrays to compute them"
    diagnostic_variables::DiagnosticVariables

    "all parameters, constant at runtime"
    model::Model
end

function Base.show(io::IO, S::AbstractSimulation)
    println(io, "Simulation{$(model_type(S.model))}")
    println(io, "├ $(typeof(S.prognostic_variables))")
    println(io, "├ $(typeof(S.diagnostic_variables))")
    print(io, "└ model::$(model_type(S.model))")
end

export run!

"""
$(TYPEDSIGNATURES)
Run a SpeedyWeather.jl `simulation`. The `simulation.model` is assumed to be initialized."""
function run!(  simulation::AbstractSimulation;
                period = Day(10),
                output::Bool = false,
                n_days::Union{Nothing, Real} = nothing)
    
    if !isnothing(n_days)
        @warn "run!: n_days keyword is deprecated, use period = Day(n_days) instead."
        period = Day(n_days) 
    end 

    (; prognostic_variables, diagnostic_variables, model) = simulation
    (; clock) = prognostic_variables

    # CLOCK
    clock.period = period                   # set the clock's enddate      
    initialize!(clock, model.time_stepping) # store the start date, reset counter

    # OUTPUT
    model.output.output = output            # enable/disable output

    # run it, yeah!
    time_stepping!(prognostic_variables, diagnostic_variables, model)
end
