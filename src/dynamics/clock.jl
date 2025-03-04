const DEFAULT_DATE = DateTime(2000, 1, 1)

abstract type AbstractClock end

export Clock
"""
Clock struct keeps track of the model time, how many days to integrate for
and how many time steps this takes
$(TYPEDFIELDS)."""
Base.@kwdef mutable struct Clock <: AbstractClock
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

# copy! 
function Base.copy!(clock::Clock, clock_old::Clock)
    clock.time = clock_old.time 
    clock.start = clock_old.start 
    clock.period = clock_old.period 
    clock.timestep_counter = clock_old.timestep_counter
    clock.n_timesteps = clock_old.n_timesteps
    clock.Δt = clock_old.Δt 
    
    return nothing 
end 

"""
$(TYPEDSIGNATURES)
Initialize the clock with the time step `Δt` in the `time_stepping`."""
function initialize!(
    clock::Clock,
    time_stepping::AbstractTimeStepper;
    n_timesteps::Int = -1,
)
    # if >=0 use to set n_timesteps, otherwise calculate from period
    if n_timesteps >= 0
        clock.n_timesteps = n_timesteps
        clock.period = Second(time_stepping.Δt_sec * n_timesteps)
    else
        clock.n_timesteps = ceil(Int, clock.period.value/time_stepping.Δt_sec)
    end
    
    clock.start = clock.time    # store the start time
    clock.timestep_counter = 0  # reset counter
    clock.Δt = time_stepping.Δt_millisec
    return clock
end

function set_period!(clock::Clock, steps::Integer)
    clock.period = Second(period)
end


"""
$(TYPEDSIGNATURES)
Set the `period` of the clock to a new value. Converts any `Dates.Period` input to `Second`."""
function set_period!(clock::Clock, period::Period)
    clock.period = Second(period)
end

"""
$(TYPEDSIGNATURES)
Set the `period` of the clock to a new value. Converts any `::Real` input to `Day`."""
function set_period!(clock::Clock, period::Real)
    @info "Input $period assumed to have units of days. Use Week($period), Hour($period), Minute($period) otherwise."
    clock.period = Day(period)
end

"""
$(TYPEDSIGNATURES)
Create and initialize a clock from `time_stepping`"""
function Clock(time_stepping::AbstractTimeStepper; kwargs...)
    clock = Clock(; kwargs...)
    initialize!(clock, time_stepping)
end