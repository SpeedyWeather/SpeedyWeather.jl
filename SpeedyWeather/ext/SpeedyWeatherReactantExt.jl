module SpeedyWeatherReactantExt

using SpeedyWeather
using Reactant
using DocStringExtensions
using Dates

using SpeedyWeather: ReactantDevice, scale!, get_step, unpack, timestep!, first_timesteps!, later_timestep!

const ReactantDatesExt = Base.get_extension(
    Reactant, :ReactantDatesExt
)

const ReactantSimulation = Union{
    Simulation{V, <:BarotropicModel{SG, <:ReactantDevice}},
    Simulation{V, <:ShallowWaterModel{SG, <:ReactantDevice}}, Simulation{V, <:PrimitiveDryModel{SG, <:ReactantDevice}}, Simulation{V, <:PrimitiveWetModel{SG, <:ReactantDevice}},
} where {V, SG}

# time stepping functions with Reactant, take compiled functions as optional arguments
# in case they are not provided, they are compiled on the fly

"""
$(TYPEDSIGNATURES)

Time-stepping function for a simulation with Reactant, take compiled functions as optional arguments.
In case they are not provided, they are compiled on the fly.

Example usage:

```julia
simulation = initialize!(model) 
initialize!(simulation; steps=10) # don't forget this! 
r_first! = @compile SpeedyWeather.first_timesteps!(simulation)
r_later! = @compile SpeedyWeather.later_timestep!(simulation)
SpeedyWeather.time_stepping!(simulation, r_first!, r_later!)
SpeedyWeather.finalize!(simulation)
```
"""
function SpeedyWeather.time_stepping!(simulation::ReactantSimulation, r_first_timesteps! = nothing, r_later_timestep! = nothing, enable_checkpointing = true)
    
    clock = simulation.variables.prognostic.clock

    if isnothing(r_first_timesteps!)
        @info "Reactant compiling first_timesteps!"
        r_first_timesteps! = @compile first_timesteps!(simulation)
    end

    #TODO: reenable @trace once Reactant issues fixed
    #@trace checkpointing = enable_checkpointing for _ in clock.timestep_counter:clock.n_timesteps
    #    r_later_timestep!(simulation)
    #end
    
    r_first_timesteps!(simulation)

    if isnothing(r_later_timestep!)
        @info "Reactant compiling later_timestep!"
        r_later_timestep! = @compile later_timestep!(simulation)
    end

    for i in (Int(clock.timestep_counter) + 1):Int(clock.n_timesteps)
        r_later_timestep!(simulation)
    end
    return
end

# that's for Reactant TracableDateTime
SpeedyWeather.secondofday(dt::ReactantDatesExt.ReactantDateTime) = Dates.second(convert(ReactantDatesExt.ReactantTime, dt).instant)

# For Reactant tracing, dispatch on the ReactantDateTime argument and forward
# to the internal `_year_angle`/`_solar_hour_angle` implementations. The
# `length_of_day`/`length_of_year` arguments may be regular `Second` (from the
# untraced SolarZenith struct fields) or Reactant period types, so we accept any.
# Both signatures are spelled out to avoid ambiguity with the SpeedyWeather methods
# that dispatch on `Dates.AbstractDateTime` + `Second`.
@inline SpeedyWeather.year_angle(::Type{T}, time::ReactantDatesExt.ReactantDateTime, length_of_day::Second, length_of_year::Second) where {T} = SpeedyWeather._year_angle(T, time, length_of_day, length_of_year)
@inline SpeedyWeather.year_angle(::Type{T}, time::ReactantDatesExt.ReactantDateTime, length_of_day::ReactantDatesExt.ReactantSecond, length_of_year::ReactantDatesExt.ReactantSecond) where {T} = SpeedyWeather._year_angle(T, time, length_of_day, length_of_year)
@inline SpeedyWeather.solar_hour_angle(::Type{T}, time::ReactantDatesExt.ReactantDateTime, λ, length_of_day::Second) where {T} = SpeedyWeather._solar_hour_angle(T, time, λ, length_of_day)
@inline SpeedyWeather.solar_hour_angle(::Type{T}, time::ReactantDatesExt.ReactantDateTime, λ, length_of_day::ReactantDatesExt.ReactantSecond) where {T} = SpeedyWeather._solar_hour_angle(T, time, λ, length_of_day)

#TODO: move those the the ReactantDatesExt once I am sure it's all additional functionality I need to add
Dates.firstdayofmonth(dt::ReactantDatesExt.ReactantDate{I}) where {I} =
    ReactantDatesExt.ReactantDate{I}(Date(Dates.year(dt), Dates.month(dt)))
Dates.firstdayofmonth(dt::ReactantDatesExt.ReactantDateTime{I}) where {I} =
    ReactantDatesExt.ReactantDateTime{I}(DateTime(Dates.year(dt), Dates.month(dt)))

# ReactantDatesExt versions of Base.convert methods (mirrors Dates.jl conversions)
Base.convert(::Type{ReactantDatesExt.ReactantDateTime}, dt::ReactantDatesExt.ReactantDate) =
    ReactantDatesExt.ReactantDateTime(Dates.UTInstant(ReactantDatesExt.ReactantMillisecond(Dates.value(dt) * 86400000)))
Base.convert(::Type{ReactantDatesExt.ReactantDate}, dt::ReactantDatesExt.ReactantDateTime) =
    ReactantDatesExt.ReactantDate(Dates.UTInstant(ReactantDatesExt.ReactantDay(div(Dates.value(dt), 86400000))))
Base.convert(::Type{ReactantDatesExt.ReactantTime}, dt::ReactantDatesExt.ReactantDateTime) =
    ReactantDatesExt.ReactantTime(ReactantDatesExt.ReactantNanosecond((Dates.value(dt) % 86400000) * 1000000))

Base.convert(::Type{ReactantDatesExt.ReactantDateTime}, x::Millisecond) =
    ReactantDatesExt.ReactantDateTime(Dates.UTInstant(ReactantDatesExt.ReactantMillisecond(Dates.value(x))))
Base.convert(::Type{Millisecond}, dt::ReactantDatesExt.ReactantDateTime) =
    Millisecond(Dates.value(dt))
Base.convert(::Type{ReactantDatesExt.ReactantDate}, x::Day) =
    ReactantDatesExt.ReactantDate(Dates.UTInstant(ReactantDatesExt.ReactantDay(Dates.value(x))))
Base.convert(::Type{Day}, dt::ReactantDatesExt.ReactantDate) =
    Day(Dates.value(dt))

# Typed outer constructors for Reactant period types so that Period arithmetic
# (e.g. -(x::P, y::P) = P(value(x)-value(y))) preserves the type parameter I
# when value arithmetic falls back to a plain Int via ConcretePJRTNumber.to_number.
for T in (
    :ReactantYear, :ReactantMonth, :ReactantDay,
    :ReactantHour, :ReactantMinute, :ReactantSecond,
    :ReactantMillisecond, :ReactantMicrosecond, :ReactantNanosecond,
)
    @eval (::Type{ReactantDatesExt.$T{I}})(v::Number) where {I} = ReactantDatesExt.$T(convert(I, v))
end

# Dates.days for ReactantMillisecond (mirrors days(c::Millisecond) = div(value(c), 86400000))
Dates.days(rm::ReactantDatesExt.ReactantMillisecond) = div(Dates.value(rm), 86400000)

# Reactant-compatible versions of Dates.yearmonthday, yearmonth, monthday
# These replace ternary ? : operators with if/else for Reactant tracing compatibility.
# The algorithm is identical to Dates/src/accessors.jl (proleptic Gregorian calendar).
function Dates.yearmonthday(days::Reactant.TracedRNumber)
    z = days + 306; h = 100z - 25; a = fld(h, 3652425); b = a - fld(a, 4)
    y = fld(100b + h, 36525); c = b + z - 365y - fld(y, 4); m = div(5c + 456, 153)
    d = c - div(153m - 457, 5)
    @trace if m > 12
        res = (y + 1, m - 12, d)
    else
        res = (y, m, d)
    end
    return res 
end

function Dates.yearmonth(days::Reactant.TracedRNumber)
    z = days + 306; h = 100z - 25; a = fld(h, 3652425); b = a - fld(a, 4)
    y = fld(100b + h, 36525); c = b + z - 365y - fld(y, 4); m = div(5c + 456, 153)
    @trace if m > 12
        res = (y + 1, m - 12)
    else
        res = (y, m)
    end
    return res
end

function Dates.monthday(days::Reactant.TracedRNumber)
    z = days + 306; h = 100z - 25; a = fld(h, 3652425); b = a - fld(a, 4)
    y = fld(100b + h, 36525); c = b + z - 365y - fld(y, 4); m = div(5c + 456, 153)
    d = c - div(153m - 457, 5)
    @trace if m > 12
        res = (m - 12, d)
    else
        res = (m, d)
    end
    return res
end

# Reactant-compatible version of Dates.isleapyear that replaces the && / ||
# short-circuiting (which requires concrete Bool) with bitwise operators that
# work on TracedRNumber{Bool}.
Dates.isleapyear(y::Reactant.TracedRNumber) =
    (y % 4 == 0) & ((y % 100 != 0) | (y % 400 == 0))

# Reactant-compatible version of Dates.dayofyear(y, m, d)
# Replaces the MONTHDAYS[m] tuple lookup (which fails on TracedRNumber) with a
# branching computation over the cumulative month-day offsets
# (0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334).
function Dates.dayofyear(y::Reactant.TracedRNumber, m::Reactant.TracedRNumber, d::Reactant.TracedRNumber)
    @trace if m == 1
        monthdays = zero(m)
    elseif m == 2
        monthdays = oftype(m, 31)
    elseif m == 3
        monthdays = oftype(m, 59)
    elseif m == 4
        monthdays = oftype(m, 90)
    elseif m == 5
        monthdays = oftype(m, 120)
    elseif m == 6
        monthdays = oftype(m, 151)
    elseif m == 7
        monthdays = oftype(m, 181)
    elseif m == 8
        monthdays = oftype(m, 212)
    elseif m == 9
        monthdays = oftype(m, 243)
    elseif m == 10
        monthdays = oftype(m, 273)
    elseif m == 11
        monthdays = oftype(m, 304)
    else
        monthdays = oftype(m, 334)
    end
    leap_correction = ifelse((m > 2) & Dates.isleapyear(y), one(m), zero(m))
    return monthdays + d + leap_correction
end

# dayofyear accessor for Reactant date/datetime types mirrors the pattern in
# ReactantDatesExt/accessors.jl (yearmonthday then forward to dayofyear(y, m, d))
function Dates.dayofyear(dt::Union{ReactantDatesExt.ReactantDate, ReactantDatesExt.ReactantDateTime})
    y, m, d = Dates.yearmonthday(dt)
    return Dates.dayofyear(y, m, d)
end

# These function extend those defined in SpeedyWeather/src/dynamics/clock.jl
# They will not move to ReactantDatesExt as they aren't part of stdlib Dates.jl
Dates.second(x::ReactantDatesExt.ReactantNanosecond) = round(Int, x.value * 1.0e-9)
Dates.second(x::ReactantDatesExt.ReactantMicrosecond) = round(Int, x.value * 1.0e-6)
Dates.second(x::ReactantDatesExt.ReactantMillisecond) = round(Int, x.value * 1.0e-3)

SpeedyWeather.Clock(architecture::ReactantDevice) = Reactant.to_rarray(SpeedyWeather.Clock(), track_numbers = true)

# _units is defined per concrete Period type in Dates with singular/plural logic.
# For Reactant period types the value is a traced number so we can't branch on it;
# always return the plural form which is the common case.
for (T, unit) in (
    (:ReactantYear, "years"),
    (:ReactantMonth, "months"),
    (:ReactantDay, "days"),
    (:ReactantHour, "hours"),
    (:ReactantMinute, "minutes"),
    (:ReactantSecond, "seconds"),
    (:ReactantMillisecond, "milliseconds"),
    (:ReactantMicrosecond, "microseconds"),
    (:ReactantNanosecond, "nanoseconds"),
)
    @eval Dates._units(::ReactantDatesExt.$T) = " " * $unit
end

end
