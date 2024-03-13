const DEFAULT_DATE = DateTime(2000, 1, 1)

export Clock
"""
Clock struct keeps track of the model time, how many days to integrate for
and how many time steps this takes
$(TYPEDFIELDS)."""
Base.@kwdef mutable struct Clock
    "current model time"
    time::DateTime = DEFAULT_DATE

    "start time of simulation"
    start::DateTime = DEFAULT_DATE
    
    "period to integrate for, set in set_period!(::Clock, ::Dates.Period)"
    period::Second = Second(0)

    "Counting all time steps during simulation"
    timestep_counter::Int = 0

    "number of time steps to integrate for, set in initialize!(::Clock, ::AbstractTimeStepper)"
    n_timesteps::Int = 0

    "Time step"
    Δt::Millisecond = Millisecond(0)
end

function timestep!(clock::Clock, Δt; increase_counter::Bool=true)
    clock.time += Δt
    # the first timestep is a half-step and doesn't count
    clock.timestep_counter += increase_counter  
end

# pretty printing
function Base.show(io::IO, C::Clock)
    println(io, "$(typeof(C))")
    keys = propertynames(C)
    print_fields(io, C, keys)
end

"""
$(TYPEDSIGNATURES)
Initialize the clock with the time step `Δt` in the `time_stepping`."""
function initialize!(clock::Clock, time_stepping::AbstractTimeStepper)
    clock.n_timesteps = ceil(Int, clock.period.value/time_stepping.Δt_sec)
    clock.start = clock.time    # store the start time
    clock.timestep_counter = 0  # reset counter
    clock.Δt = time_stepping.Δt_millisec
    return clock
end

"""
$(TYPEDSIGNATURES)
Create and initialize a clock from `time_stepping`"""
function Clock(time_stepping::AbstractTimeStepper; kwargs...)
    clock = Clock(; kwargs...)
    initialize!(clock, time_stepping)
end