const DEFAULT_DATE = DateTime(2000, 1, 1)

abstract type AbstractClock end

export Clock
"""
Clock struct keeps track of the model time, how many days to integrate for
and how many time steps this takes.
$(TYPEDFIELDS)"""
Base.@kwdef mutable struct Clock <: AbstractClock
    "current model time"
    time::DateTime = DEFAULT_DATE

    "start time of simulation"
    start::DateTime = DEFAULT_DATE

    "period to integrate for"
    period::Second = Second(0)

    "Counting all time steps during simulation"
    timestep_counter::Int = 0

    "number of time steps to integrate for, set in `initialize!(::Clock, ::AbstractTimeStepper)`"
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

"""$(TYPEDSIGNATURES)
Initialize the clock with the time step `Δt` from `time_stepping`."""
initialize!(clock::Clock, time_stepping::AbstractTimeStepper, args...) =
    initialize!(clock, time_stepping.Δt_millisec, args...)

"""$(TYPEDSIGNATURES)
Initialize the clock with the time step `Δt` and `period` to integrate for."""
function initialize!(clock::Clock, Δt::Period, period::Period)
    clock.Δt = Δt
    clock.period = period
    clock.n_timesteps = ceil(Int, Millisecond(period).value/Millisecond(Δt).value)
    initialize!(clock)      # set start time, reset counter
end

"""$(TYPEDSIGNATURES)
Initialize the clock with the time step `Δt` and the number of time steps `n_timesteps`."""
function initialize!(clock::Clock, Δt::Period, n_timesteps::Integer)
    clock.Δt = Δt
    clock.n_timesteps = n_timesteps
    clock.period = Δt * n_timesteps
    initialize!(clock)
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
    initialize!(clock, time_stepping, clock.period)
end

"""
    Century <: Period

Convenience time type representing a 100-year period.
"""
struct Century <: Period
    value::Int64
end

# Constructor for Period types to resolve ambiguous dispatch
Century(p::Period) = convert(Century, p)

Dates._units(m::Century) = m.value == 1 ? " century" : " centuries"

# convert Century -> Year
Base.convert(::Type{Year}, c::Century) = Year(c.value*100)

# promotion rules
Base.promote_rule(::Type{Century}, ::Type{Year}) = Year
Base.promote_rule(::Type{Century}, ::Type{Month}) = Month
Base.promote_rule(::Type{Century}, ::Type{Day}) = Day
Base.promote_rule(::Type{Century}, ::Type{Hour}) = Hour
Base.promote_rule(::Type{Century}, ::Type{Second}) = Second

"""
    Millenium <: Period

Convenience time type representing a 1000-year period.
"""
struct Millenium <: Period
    value::Int64
end

# Constructor for Period types to resolve ambiguous dispatch
Millenium(p::Period) = convert(Millenium, p)

Dates._units(m::Millenium) = m.value == 1 ? " millenium" : " millenia"

# convert Millenium -> Century and Year
Base.convert(::Type{Century}, m::Millenium) = Century(m.value*10)
Base.convert(::Type{Year}, m::Millenium) = Year(Century(m))

# promotion rules for converting to common types, e.g. in collections
Base.promote_rule(::Type{Millenium}, ::Type{Century}) = Century
Base.promote_rule(::Type{Millenium}, ::Type{Year}) = Year
Base.promote_rule(::Type{Millenium}, ::Type{Month}) = Month
Base.promote_rule(::Type{Millenium}, ::Type{Day}) = Day
Base.promote_rule(::Type{Millenium}, ::Type{Hour}) = Hour
Base.promote_rule(::Type{Millenium}, ::Type{Second}) = Second

# add coarserperiod dispatches for Century and Millenium
Dates.coarserperiod(::Type{Year}) = (Century, 100)
Dates.coarserperiod(::Type{Century}) = (Millenium, 10)

# conversion rules for floating point -> time types
Dates.Second(       x::AbstractFloat) = convert(Second, x)
Dates.Minute(       x::AbstractFloat) = Second(60x)
Dates.Hour(         x::AbstractFloat) = Minute(60x)
Dates.Day(          x::AbstractFloat) = Hour(24x)
Dates.Week(         x::AbstractFloat) = Day(7x)
Dates.Month(        x::AbstractFloat) = Day(30x)  # approximate
Dates.Year(         x::AbstractFloat) = Day(365x) # approximate
Century(      x::AbstractFloat) = Year(100x)
Millenium(    x::AbstractFloat) = Century(10x)

# use Dates.second to round to integer seconds
Dates.second(x::Dates.Nanosecond) = round(Int, x.value*1e-9)
Dates.second(x::Dates.Microsecond) = round(Int, x.value*1e-6)
Dates.second(x::Dates.Millisecond) = round(Int, x.value*1e-3)

# defined to convert from floats to Dates.Second (which require ints by default) via rounding
function Base.convert(::Type{Second}, x::AbstractFloat)
    xr = round(Int64, x)
    x == xr || @warn "Rounding and converting $x to $xr for integer seconds."
    return Second(xr)
end

# conversion rule that allows integers to be autmoatically converted into Seconds
function Base.convert(::Type{Second}, x::Integer)
    @warn "Input '$x' assumed to have units of seconds. Use Minute($x), Hour($x), or Day($x) otherwise."
    return Second(round(Int64, x))
end

# month conversions
Base.convert(::Type{Dates.Day}, m::Dates.Month) = Day(Hour(m))
Base.convert(::Type{Dates.Hour}, m::Dates.Month) = Hour(Second(m))
Base.convert(::Type{Dates.Minute}, m::Dates.Month) = Minute(Second(m))
function Base.convert(::Type{Dates.Second}, m::Dates.Month)
    return Second(m.value * 30 * 24 * 60 * 60) # approximate
end
Base.convert(::Type{Dates.Millisecond}, m::Dates.Month) = Millisecond(Second(m))

# year conversions
Base.convert(::Type{Dates.Day}, y::Dates.Year) = Day(Hour(y))
Base.convert(::Type{Dates.Hour}, y::Dates.Year) = Hour(Second(y))
Base.convert(::Type{Dates.Minute}, y::Dates.Year) = Minute(Second(y))

function Base.convert(::Type{Dates.Second}, y::Dates.Year)
    return Second(y.value * 365 * 24 * 60 * 60) # approximate
end

Base.convert(::Type{Dates.Millisecond}, y::Dates.Year) = Millisecond(Second(y))

# additional century conversions
Base.convert(::Type{Millisecond}, c::Century) = Millisecond(Second(c))
Base.convert(::Type{Second}, c::Century) = Second(Year(c))
Base.convert(::Type{Hour}, c::Century) = Hour(Year(c))
Base.convert(::Type{Day}, c::Century) = Day(Year(c))
Base.convert(::Type{Month}, c::Century) = Month(Year(c))

# additional month conversions
Base.convert(::Type{Millisecond}, m::Millenium) = Millisecond(Second(m))
Base.convert(::Type{Second}, m::Millenium) = Second(Year(m))
Base.convert(::Type{Hour}, m::Millenium) = Hour(Year(m))
Base.convert(::Type{Day}, m::Millenium) = Day(Year(m))
Base.convert(::Type{Month}, m::Millenium) = Month(Year(m))