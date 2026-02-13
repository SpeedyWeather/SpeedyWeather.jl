const DEFAULT_DATE = DateTime(2000, 1, 1)

abstract type AbstractClock end

export Clock
"""
Clock struct keeps track of the model time, how many days to integrate for
and how many time steps this takes.
$(TYPEDFIELDS)"""
@kwdef mutable struct Clock{DT, S, I, MS} <: AbstractClock
    "current model time"
    time::DT = DEFAULT_DATE

    "start time of simulation"
    start::DT = DEFAULT_DATE

    "period to integrate for"
    period::S = Second(0)

    "Counting all time steps during simulation"
    timestep_counter::I = 0

    "number of time steps to integrate for, set in `initialize!(::Clock, ::AbstractTimeStepper)`"
    n_timesteps::I = 0

    "Time step"
    Δt::MS = Millisecond(0)
end

# we don't want to adapt the clock, it has to stay mutable,
# we just explicitly transfer nothing to the kernels
Adapt.adapt_structure(to, ::Clock) = nothing

function timestep!(clock::Clock, Δt; increase_counter::Bool = true)
    clock.time += Δt
    # the first timestep is a half-step and doesn't count
    return clock.timestep_counter += increase_counter
end

# pretty printing
function Base.show(io::IO, C::Clock)
    println(io, "$(typeof(C))")
    keys = propertynames(C)
    return print_fields(io, C, keys)
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

"""$(TYPEDSIGNATURES)
Initialize the clock with the time step `Δt` from `time_stepping`."""
initialize!(clock::Clock, time_stepping::AbstractTimeStepper, args...) =
    initialize!(clock, time_stepping.Δt_millisec, args...)

"""$(TYPEDSIGNATURES)
Initialize the clock with the time step `Δt` and `period` to integrate for."""
function initialize!(clock::Clock, Δt::Period, period::Period)
    clock.Δt = Δt
    clock.period = period
    clock.n_timesteps = ceil(Int, Millisecond(period).value / Millisecond(Δt).value)
    return initialize!(clock)      # set start time, reset counter
end

"""$(TYPEDSIGNATURES)
Initialize the clock with the time step `Δt` and the number of time steps `n_timesteps`."""
function initialize!(clock::Clock, Δt::Period, n_timesteps::Integer)
    clock.Δt = Δt
    clock.n_timesteps = n_timesteps
    clock.period = Δt * n_timesteps
    return initialize!(clock)
end

"""$(TYPEDSIGNATURES)
Initialize the clock with setting the start time and resetting the timestep counter."""
function initialize!(clock::Clock)
    clock.start = clock.time    # store the start time
    clock.timestep_counter = 0  # reset counter
    return clock
end

"""
$(TYPEDSIGNATURES)
Create and initialize a clock from `time_stepping`."""
function Clock(time_stepping::AbstractTimeStepper; kwargs...)
    clock = Clock(; kwargs...)
    return initialize!(clock, time_stepping, clock.period)
end

# conversion rules for floating point -> time types
Dates.Second(x::AbstractFloat) = convert(Second, x)
Dates.Minute(x::AbstractFloat) = Second(60x)
Dates.Hour(x::AbstractFloat) = Minute(60x)
Dates.Day(x::AbstractFloat) = Hour(24x)
Dates.Week(x::AbstractFloat) = Day(7x)
Dates.Month(x::AbstractFloat) = Day(30x)  # approximate
Dates.Year(x::AbstractFloat) = Day(365x) # approximate
Dates.Century(x::AbstractFloat) = Year(100x)
Dates.Millenium(x::AbstractFloat) = Century(10x)

# defined to convert from floats to Second (which require ints by default) via rounding
function Base.convert(::Type{<:Second}, x::AbstractFloat)
    xr = round(Int64, x)
    x == xr || @warn "Rounding and converting $x to $xr for integer seconds."
    return Second(xr)
end

# conversion rule that allows integers to be autmoatically converted into Seconds
function Base.convert(::Type{<:Second}, x::Integer)
    @warn "Input '$x' assumed to have units of seconds. Use Minute($x), Hour($x), or Day($x) otherwise."
    return Second(round(Int64, x))
end

# month conversions (approximate: 1 month = 30 days)
Base.convert(::Type{<:Day}, m::Month) = Day(Hour(m))
Base.convert(::Type{<:Hour}, m::Month) = Hour(Second(m))
Base.convert(::Type{<:Minute}, m::Month) = Minute(Second(m))
function Base.convert(::Type{<:Second}, m::Month)
    return Second(m.value * 30 * 24 * 60 * 60) # approximate
end
Base.convert(::Type{<:Millisecond}, m::Month) = Millisecond(Second(m))

# year conversions (approximate: 1 year = 365 days)
Base.convert(::Type{<:Day}, y::Year) = Day(Hour(y))
Base.convert(::Type{<:Hour}, y::Year) = Hour(Second(y))
Base.convert(::Type{<:Minute}, y::Year) = Minute(Second(y))

function Base.convert(::Type{<:Second}, y::Year)
    return Second(y.value * 365 * 24 * 60 * 60) # approximate
end

Base.convert(::Type{<:Millisecond}, y::Year) = Millisecond(Second(y))
