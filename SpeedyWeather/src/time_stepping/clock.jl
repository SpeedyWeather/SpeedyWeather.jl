const DEFAULT_DATE = DateTime(2000, 1, 1)

abstract type AbstractClock <: AbstractModelComponent end

export Clock
"""
Clock struct keeps track of the model time, how many days to integrate for
and how many time steps this takes.
$(TYPEDFIELDS)"""
@kwdef mutable struct Clock{DT, S, NF, I, MS} <: AbstractClock
    "current model time"
    time::DT = DEFAULT_DATE

    "start time of simulation"
    start::DT = DEFAULT_DATE

    "Rotation time for daily cycle"
    rotation_time::DT = time

    "Orbit time for seasonal cycle"
    orbit_time::DT = time

    "period to integrate for"
    period::S = Second(0)

    "Dilation (relative speed wrt time) of rotation time for daily cycle, used for solar zenith calculations"
    rotation_dilation::NF = 1f0

    "Dilation (relative speed wrt time) of orbit time for seasonal cycle, used for solar zenith calculations"
    orbit_dilation::NF = 1f0

    "Counting all steps the time stepper takes during simulation"
    step_counter::I = 0

    "Counting all time steps of size Δt, ignoring additional spin-up steps"
    time_step_counter::I = 0

    "Number of time stepper steps, regardless their step size"
    n_steps::I = 0
    
    "Number of time steps of size Δt"
    n_time_steps::I = 0

    "Time step"
    Δt::MS = Millisecond(0)
end

# default constructor (extended e.g. by ReactantExt)
Clock(architecture::AbstractArchitecture) = Clock()

# we don't want to adapt the clock, it has to stay mutable,
# we just explicitly transfer nothing to the kernels
Adapt.adapt_structure(to, ::Clock) = nothing

# Clock has no GPU arrays, pass through unchanged
on_architecture(::AbstractArchitecture, clock::Clock) = clock

time_step!(clock::Clock, time_stepping::AbstractTimeStepper) =
    time_step!(clock, time_stepping.Δt_millisec)

"""$(TYPEDSIGNATURES)
Time step `Δt` dilated by `dilation`, rounded to whole milliseconds. Used to advance
(or, on the leapfrog spin-up step, rewind) the rotation and orbit times."""
dilate(Δt, dilation) = Millisecond(round(Int, Millisecond(Δt).value * dilation))

function time_step!(clock::Clock, Δt; increase_counter::Bool = true)
    clock.time += Δt
    clock.step_counter += 1                     # always increased, counts time stepper steps
    clock.time_step_counter += increase_counter # spin up steps may not count for clock

    # step the rotation and orbit time separately for solar zenith calculations
    # using the dilation factors to speed up or slow down the rotation and orbit time relative to the model time
    clock.rotation_time += dilate(Δt, clock.rotation_dilation)
    clock.orbit_time += dilate(Δt, clock.orbit_dilation)

    return nothing
end

# copy! (converts on the fly for Reactant types to work as well)
function Base.copy!(clock::Clock{DT, S, NF, I, MS}, clock_old::Clock) where {DT, S, NF, I, MS}
    # explicitly convert to the new types too
    clock.time = convert(DT, clock_old.time)
    clock.start = convert(DT, clock_old.start)
    clock.rotation_time = convert(DT, clock_old.rotation_time)
    clock.orbit_time = convert(DT, clock_old.orbit_time)
    clock.period = convert(S, clock_old.period)
    clock.rotation_dilation = convert(NF, clock_old.rotation_dilation)
    clock.orbit_dilation = convert(NF, clock_old.orbit_dilation)
    clock.step_counter = convert(I, clock_old.step_counter)
    clock.time_step_counter = convert(I, clock_old.time_step_counter)
    clock.n_steps = convert(I, clock_old.n_steps)
    clock.n_time_steps = convert(I, clock_old.n_time_steps)
    clock.Δt = convert(MS, clock_old.Δt)
    return nothing
end

# for copy!(::Variables, ::Variables)
_copy_entry!(dest::AbstractClock, src::AbstractClock) = copy!(dest, src)

"""$(TYPEDSIGNATURES)
Initialize the clock with the time step `Δt` from `time_stepping` and
computes the time dilation given a `planet` with `length_of_day` and `length_of_year`."""
function initialize!(clock::Clock, time_stepping::AbstractTimeStepper, planet::AbstractPlanet, args...)

    # given a planet's length of day and year, we can compute the dilation factors for the rotation and orbit time
    @assert planet.length_of_day.value != 0 "planet.length_of_day is 0 (rotating infinitely fast)," *
        " choose Second(typemax(Int)) for a tidally locked planet or daily_cycle=false for average insulation." 
    @assert planet.length_of_year.value != 0 "planet.length_of_year is 0 (orbiting infinitely fast)," *
        " choose Second(typemax(Int)) for a seasonally locked planet or seasonal_cycle=false." 
    clock.rotation_dilation = Second(EARTH_DAY).value / Second(planet.length_of_day).value
    clock.orbit_dilation = Second(EARTH_YEAR).value / Second(planet.length_of_year).value

    # now initialize rest of the clock as normal
    return initialize!(clock, time_stepping, args...)
end

"""$(TYPEDSIGNATURES)
Initialize the clock with the time step `Δt` and `period` to integrate for.
`n_time_steps` is for the clock only, spin up steps (e.g. Leapfrog with 1 Euler to start)
are not counted for the clock."""
function initialize!(clock::Clock, time_stepping::AbstractTimeStepper, period::Period)
    Δt = time_stepping.Δt_millisec
    n_time_steps = ceil(Int, Millisecond(period).value / Millisecond(Δt).value)
    clock.Δt = Δt
    clock.period = period
    clock.n_time_steps = n_time_steps
    clock.n_steps = n_time_steps + spin_up_steps(time_stepping)
    return initialize!(clock)      # set start time, reset counter
end

"""$(TYPEDSIGNATURES)
Initialize the clock with the time step `Δt` and the number of time stepper steps `n_steps`."""
function initialize!(clock::Clock, time_stepping::AbstractTimeStepper, n_steps::Integer)
    Δt = time_stepping.Δt_millisec
    n_time_steps = max(0, n_steps - spin_up_steps(time_stepping)) # in case there is spin up the clock may not advance
    clock.Δt = Δt
    clock.n_steps = n_steps     # number of steps the time stepper should take, regardless step size
    clock.n_time_steps = n_time_steps
    clock.period = Δt * n_time_steps    # therefore calculate period from n_time_steps not n_steps
    return initialize!(clock)
end

"""$(TYPEDSIGNATURES)
Initialize the clock with setting the start time and resetting the (time) step counters."""
function initialize!(clock::Clock)
    clock.start = clock.time        # store the start time
    clock.step_counter = 0          # reset counter for time stepper steps, regardless step size
    clock.time_step_counter = 0     # reset counter for steps of size Δt
    return clock
end

function set!(clock::Clock; time::DateTime, start::DateTime)
    clock.time = time
    clock.start = start
    clock.rotation_time = time      # when time is explicitly set, synchronize orbit & rotation times too
    clock.orbit_time = time         # after this they generally run apart 
    return clock
end

"""
$(TYPEDSIGNATURES)
Create and initialize a clock from `time_stepping`."""
function Clock(time_stepping::AbstractTimeStepper; kwargs...)
    clock = Clock(; kwargs...)
    return initialize!(clock, time_stepping, clock.period)
end

"""
    Century{IntType} <: Period

Convenience time type representing a 100-year period.
Parametric on `IntType` for Reactant compatibility."""
struct Century{IntType} <: Period
    value::IntType
end

Century(v::Integer) = Century{typeof(v)}(v)

# Constructor for Period types to resolve ambiguous dispatch
Century(p::Period) = convert(Century, p)

Dates._units(m::Century) = abs(m.value) == 1 ? " century" : " centuries"

# convert Century -> Year
Base.convert(::Type{Year}, c::Century) = Year(c.value * 100)

# promotion rules
Base.promote_rule(::Type{<:Century}, ::Type{Year}) = Year
Base.promote_rule(::Type{<:Century}, ::Type{Month}) = Month
Base.promote_rule(::Type{<:Century}, ::Type{Day}) = Day
Base.promote_rule(::Type{<:Century}, ::Type{Hour}) = Hour
Base.promote_rule(::Type{<:Century}, ::Type{Second}) = Second

"""
    Millenium{IntType} <: Period

Convenience time type representing a 1000-year period.
Parametric on `IntType` for Reactant compatibility."""
struct Millenium{IntType} <: Period
    value::IntType
end

Millenium(v::Integer) = Millenium{typeof(v)}(v)

# Constructor for Period types to resolve ambiguous dispatch
Millenium(p::Period) = convert(Millenium, p)

Dates._units(m::Millenium) = abs(m.value) == 1 ? " millenium" : " millenia"

# convert Millenium -> Century and Year
Base.convert(::Type{Century}, m::Millenium) = Century(m.value * 10)
Base.convert(::Type{Year}, m::Millenium) = Year(Century(m))

# promotion rules for converting to common types, e.g. in collections
Base.promote_rule(::Type{<:Millenium}, ::Type{<:Century}) = Century
Base.promote_rule(::Type{<:Millenium}, ::Type{Year}) = Year
Base.promote_rule(::Type{<:Millenium}, ::Type{Month}) = Month
Base.promote_rule(::Type{<:Millenium}, ::Type{Day}) = Day
Base.promote_rule(::Type{<:Millenium}, ::Type{Hour}) = Hour
Base.promote_rule(::Type{<:Millenium}, ::Type{Second}) = Second

# add coarserperiod dispatches for Century and Millenium
Dates.coarserperiod(::Type{Year}) = (Century, 100)
Dates.coarserperiod(::Type{Century}) = (Millenium, 10)

# secondofday: seconds since midnight for a DateTime
secondofday(dt::DateTime) = div(Dates.value(Dates.Time(dt).instant), 1_000_000_000)

# conversion rules for floating point -> time types
Dates.Second(x::AbstractFloat) = convert(Second, x)
Dates.Minute(x::AbstractFloat) = Second(60x)
Dates.Hour(x::AbstractFloat) = Minute(60x)
Dates.Day(x::AbstractFloat) = Hour(24x)
Dates.Week(x::AbstractFloat) = Day(7x)
Dates.Month(x::AbstractFloat) = Day(30x)  # approximate
Dates.Year(x::AbstractFloat) = Day(365x) # approximate
Century(x::AbstractFloat) = Year(100x)
Millenium(x::AbstractFloat) = Century(10x)

Dates.Second(c::Dates.CompoundPeriod) = sum(Dates.Second.(c.periods))
Dates.Minute(c::Dates.CompoundPeriod) = sum(Dates.Minute.(c.periods))
Dates.Hour(c::Dates.CompoundPeriod) = sum(Dates.Hour.(c.periods))
Dates.Day(c::Dates.CompoundPeriod) = sum(Dates.Day.(c.periods))

# use Dates.second to round to integer seconds
Dates.second(x::Dates.Nanosecond) = round(Int, x.value * 1.0e-9)
Dates.second(x::Dates.Microsecond) = round(Int, x.value * 1.0e-6)
Dates.second(x::Dates.Millisecond) = round(Int, x.value * 1.0e-3)

# defined to convert from floats to Second (which require ints by default) via rounding
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
Base.convert(::Type{Day}, m::Month) = Day(Hour(m))
Base.convert(::Type{Hour}, m::Month) = Hour(Second(m))
Base.convert(::Type{Minute}, m::Month) = Minute(Second(m))
function Base.convert(::Type{Second}, m::Month)
    return Second(m.value * 30 * 24 * 60 * 60)  # approximate
end
Base.convert(::Type{Millisecond}, m::Month) = Millisecond(Second(m))

# year conversions
Base.convert(::Type{Day}, y::Year) = Day(Hour(y))
Base.convert(::Type{Hour}, y::Year) = Hour(Second(y))
Base.convert(::Type{Minute}, y::Year) = Minute(Second(y))

function Base.convert(::Type{Second}, y::Year)
    return Second(y.value * 365 * 24 * 60 * 60) # approximate
end

Base.convert(::Type{Millisecond}, y::Year) = Millisecond(Second(y))

# additional century conversions
Base.convert(::Type{Millisecond}, c::Century) = Millisecond(Second(c))
Base.convert(::Type{Second}, c::Century) = Second(Year(c))
Base.convert(::Type{Hour}, c::Century) = Hour(Year(c))
Base.convert(::Type{Day}, c::Century) = Day(Year(c))
Base.convert(::Type{Month}, c::Century) = Month(Year(c))

# additional millenium conversions
Base.convert(::Type{Millisecond}, m::Millenium) = Millisecond(Second(m))
Base.convert(::Type{Second}, m::Millenium) = Second(Year(m))
Base.convert(::Type{Hour}, m::Millenium) = Hour(Year(m))
Base.convert(::Type{Day}, m::Millenium) = Day(Year(m))
Base.convert(::Type{Month}, m::Millenium) = Month(Year(m))
