module TracableDates

# CUSTOM DATES IMPLEMENTATION WITH PARAMETRIC VALUE TYPES
# Standalone replacement for Julia's Dates module with `value::IntType` instead of `value::Int64`
# so that Reactant can trace through the integer operations.
# This code is mostly LLM generated, instructed to copy the Julia.Base implementation
# I really dislike that we have to do this here.
# Alternatively, we would just abandon having a Date like type anywhere in the model
# If we decide to keep this implementation, we could try to minimize it a bit further
# and really only keep the bits we use in anything that's traced by Reactant, so in a run!
# and convert between regular Dates and SpeedyDates when needed
# Type names mirror Dates directly: Second, Millisecond, DateTime, etc.

# types (matching Dates exports)
export Period, DatePeriod, TimePeriod
export Millisecond, Second, Minute, Hour, Day, Week, Month, Year
export DateTime
export CompoundPeriod
# additional types not in Dates
export Century, Millenium
# accessors (matching Dates exports)
export year, month, day, hour, minute, second, millisecond, yearmonthday
export dayofyear, daysinmonth, isleapyear, firstdayofmonth
# functions (matching Dates exports)
export canonicalize, now

# ============================================================================
# ABSTRACT TYPES
# ============================================================================
abstract type AbstractTime end
abstract type Period <: AbstractTime end
abstract type DatePeriod <: Period end
abstract type TimePeriod <: Period end
abstract type Instant <: AbstractTime end
abstract type AbstractDateTime <: AbstractTime end

# ============================================================================
# PERIOD TYPES — parametric on IntType
# ============================================================================
struct Millisecond{IntType} <: TimePeriod
    value::IntType
end

struct Second{IntType} <: TimePeriod
    value::IntType
end

struct Minute{IntType} <: TimePeriod
    value::IntType
end

struct Hour{IntType} <: TimePeriod
    value::IntType
end

struct Day{IntType} <: DatePeriod
    value::IntType
end

struct Week{IntType} <: DatePeriod
    value::IntType
end

struct Month{IntType} <: DatePeriod
    value::IntType
end

struct Year{IntType} <: DatePeriod
    value::IntType
end

struct Century{IntType} <: DatePeriod
    value::IntType
end

struct Millenium{IntType} <: DatePeriod
    value::IntType
end

# Convenience constructors: infer IntType from input
Millisecond(v::Integer) = Millisecond{typeof(v)}(v)
Second(v::Integer) = Second{typeof(v)}(v)
Minute(v::Integer) = Minute{typeof(v)}(v)
Hour(v::Integer) = Hour{typeof(v)}(v)
Day(v::Integer) = Day{typeof(v)}(v)
Week(v::Integer) = Week{typeof(v)}(v)
Month(v::Integer) = Month{typeof(v)}(v)
Year(v::Integer) = Year{typeof(v)}(v)
Century(v::Integer) = Century{typeof(v)}(v)
Millenium(v::Integer) = Millenium{typeof(v)}(v)

# ============================================================================
# VALUE ACCESSOR
# ============================================================================
value(x::Period) = x.value

# ============================================================================
# PRINT / SHOW
# ============================================================================
_units(x::Millisecond) = " millisecond" * (abs(value(x)) == 1 ? "" : "s")
_units(x::Second) = " second" * (abs(value(x)) == 1 ? "" : "s")
_units(x::Minute) = " minute" * (abs(value(x)) == 1 ? "" : "s")
_units(x::Hour) = " hour" * (abs(value(x)) == 1 ? "" : "s")
_units(x::Day) = " day" * (abs(value(x)) == 1 ? "" : "s")
_units(x::Week) = " week" * (abs(value(x)) == 1 ? "" : "s")
_units(x::Month) = " month" * (abs(value(x)) == 1 ? "" : "s")
_units(x::Year) = " year" * (abs(value(x)) == 1 ? "" : "s")
_units(x::Century) = abs(value(x)) == 1 ? " century" : " centuries"
_units(x::Millenium) = abs(value(x)) == 1 ? " millenium" : " millenia"

Base.print(io::IO, x::Period) = print(io, value(x), _units(x))
Base.show(io::IO, ::MIME"text/plain", x::Period) = print(io, x)
Base.show(io::IO, p::P) where {P <: Period} = print(io, nameof(P), '(', value(p), ')')

# ============================================================================
# ZERO / ONE / TYPEMIN / TYPEMAX / ISFINITE
# ============================================================================
Base.zero(::Type{P}) where {P <: Period} = P(0)
Base.zero(::P) where {P <: Period} = P(0)
Base.one(::Union{Type{P}, P}) where {P <: Period} = 1
Base.isfinite(::Union{Type{P}, P}) where {P <: Period} = true
Base.iszero(x::Period) = iszero(value(x))

# ============================================================================
# PERIOD ARITHMETIC
# ============================================================================
# Unary
Base.:(-)(x::P) where {P <: Period} = P(-value(x))
Base.:(+)(x::Period) = x
Base.abs(x::P) where {P <: Period} = P(abs(value(x)))
Base.sign(x::Period) = sign(value(x))
Base.signbit(x::Period) = signbit(value(x))

# Same-type binary arithmetic (P matches exact concrete type including IntType)
Base.:(==)(x::P, y::P) where {P <: Period} = value(x) == value(y)
Base.isless(x::P, y::P) where {P <: Period} = isless(value(x), value(y))

for op in (:+, :-, :lcm, :gcd)
    @eval Base.$op(x::P, y::P) where {P <: Period} = P($op(value(x), value(y)))
end

Base.:(/)(x::P, y::P) where {P <: Period} = value(x) / value(y)
Base.:(/)(x::P, y::Real) where {P <: Period} = P(value(x) / y)
Base.div(x::P, y::P, r::RoundingMode) where {P <: Period} = div(value(x), value(y), r)
Base.div(x::P, y::P) where {P <: Period} = div(value(x), value(y))

for op in (:rem, :mod)
    @eval begin
        Base.$op(x::P, y::P) where {P <: Period} = P($op(value(x), value(y)))
        Base.$op(x::P, y::Real) where {P <: Period} = P($op(value(x), y))
    end
end

Base.:(*)(x::P, y::Real) where {P <: Period} = P(value(x) * y)
Base.:(*)(y::Real, x::Period) = x * y
Base.:(*)(x::P, y::Integer) where {P <: Period} = P(value(x) * y)
Base.:(*)(y::Integer, x::Period) = x * y

# ============================================================================
# CONVERSION FACTORS (FixedPeriod chain)
# ============================================================================
# toms: convert any FixedPeriod to milliseconds
toms(x::Millisecond) = value(x)
toms(x::Second) = value(x) * 1000
toms(x::Minute) = value(x) * 60000
toms(x::Hour) = value(x) * 3600000
toms(x::Day) = value(x) * 86400000
toms(x::Week) = value(x) * 604800000
toms(x::Century) = value(x) * 100 * 365 * 86400000
toms(x::Millenium) = value(x) * 1000 * 365 * 86400000

# ============================================================================
# CONVERSIONS BETWEEN FIXED PERIOD TYPES
# ============================================================================
function _divexact(x, y)
    q, r = divrem(x, y)
    r == 0 || throw(InexactError(:_divexact, Int, x / y))
    return q
end

# Millisecond -> coarser
Base.convert(::Type{<:Second}, x::Millisecond) = Second(_divexact(value(x), 1000))
Base.convert(::Type{<:Minute}, x::Millisecond) = Minute(_divexact(value(x), 60000))
Base.convert(::Type{<:Hour}, x::Millisecond) = Hour(_divexact(value(x), 3600000))
Base.convert(::Type{<:Day}, x::Millisecond) = Day(_divexact(value(x), 86400000))
Base.convert(::Type{<:Week}, x::Millisecond) = Week(_divexact(value(x), 604800000))

# Second -> coarser / finer
Base.convert(::Type{<:Millisecond}, x::Second) = Millisecond(value(x) * 1000)
Base.convert(::Type{<:Minute}, x::Second) = Minute(_divexact(value(x), 60))
Base.convert(::Type{<:Hour}, x::Second) = Hour(_divexact(value(x), 3600))
Base.convert(::Type{<:Day}, x::Second) = Day(_divexact(value(x), 86400))
Base.convert(::Type{<:Week}, x::Second) = Week(_divexact(value(x), 604800))

# Minute -> coarser / finer
Base.convert(::Type{<:Millisecond}, x::Minute) = Millisecond(value(x) * 60000)
Base.convert(::Type{<:Second}, x::Minute) = Second(value(x) * 60)
Base.convert(::Type{<:Hour}, x::Minute) = Hour(_divexact(value(x), 60))
Base.convert(::Type{<:Day}, x::Minute) = Day(_divexact(value(x), 1440))
Base.convert(::Type{<:Week}, x::Minute) = Week(_divexact(value(x), 10080))

# Hour -> coarser / finer
Base.convert(::Type{<:Millisecond}, x::Hour) = Millisecond(value(x) * 3600000)
Base.convert(::Type{<:Second}, x::Hour) = Second(value(x) * 3600)
Base.convert(::Type{<:Minute}, x::Hour) = Minute(value(x) * 60)
Base.convert(::Type{<:Day}, x::Hour) = Day(_divexact(value(x), 24))
Base.convert(::Type{<:Week}, x::Hour) = Week(_divexact(value(x), 168))

# Day -> coarser / finer
Base.convert(::Type{<:Millisecond}, x::Day) = Millisecond(value(x) * 86400000)
Base.convert(::Type{<:Second}, x::Day) = Second(value(x) * 86400)
Base.convert(::Type{<:Minute}, x::Day) = Minute(value(x) * 1440)
Base.convert(::Type{<:Hour}, x::Day) = Hour(value(x) * 24)
Base.convert(::Type{<:Week}, x::Day) = Week(_divexact(value(x), 7))

# Week -> finer
Base.convert(::Type{<:Millisecond}, x::Week) = Millisecond(value(x) * 604800000)
Base.convert(::Type{<:Second}, x::Week) = Second(value(x) * 604800)
Base.convert(::Type{<:Minute}, x::Week) = Minute(value(x) * 10080)
Base.convert(::Type{<:Hour}, x::Week) = Hour(value(x) * 168)
Base.convert(::Type{<:Day}, x::Week) = Day(value(x) * 7)

# Year <-> Month conversions
Base.convert(::Type{<:Month}, x::Year) = Month(value(x) * 12)
Base.convert(::Type{<:Year}, x::Month) = Year(_divexact(value(x), 12))

# Century -> finer
Base.convert(::Type{<:Year}, c::Century) = Year(value(c) * 100)
Base.convert(::Type{<:Month}, c::Century) = Month(Year(c))
Base.convert(::Type{<:Day}, c::Century) = Day(Year(c))
Base.convert(::Type{<:Hour}, c::Century) = Hour(Year(c))
Base.convert(::Type{<:Second}, c::Century) = Second(Year(c))
Base.convert(::Type{<:Millisecond}, c::Century) = Millisecond(Second(c))

# Millenium -> finer
Base.convert(::Type{Century}, m::Millenium) = Century(value(m) * 10)
Base.convert(::Type{<:Year}, m::Millenium) = Year(Century(m))
Base.convert(::Type{<:Month}, m::Millenium) = Month(Year(m))
Base.convert(::Type{<:Day}, m::Millenium) = Day(Year(m))
Base.convert(::Type{<:Hour}, m::Millenium) = Hour(Year(m))
Base.convert(::Type{<:Second}, m::Millenium) = Second(Year(m))
Base.convert(::Type{<:Millisecond}, m::Millenium) = Millisecond(Second(m))

# Identity conversions
Base.convert(::Type{P}, x::P) where {P <: Period} = x

# Construct from Period (type conversion) — explicit per-type to avoid ambiguity
Millisecond(p::Period) = convert(Millisecond, p)
Second(p::Period) = convert(Second, p)
Minute(p::Period) = convert(Minute, p)
Hour(p::Period) = convert(Hour, p)
Day(p::Period) = convert(Day, p)
Week(p::Period) = convert(Week, p)
Month(p::Period) = convert(Month, p)
Year(p::Period) = convert(Year, p)
Century(p::Period) = convert(Century, p)
Millenium(p::Period) = convert(Millenium, p)

# ============================================================================
# PROMOTION RULES (FixedPeriod: promote to finer resolution)
# ============================================================================
const FixedPeriod = Union{Millenium, Century, Week, Day, Hour, Minute, Second, Millisecond}

Base.promote_rule(::Type{Week{I1}}, ::Type{Day{I2}}) where {I1, I2} = Day{promote_type(I1, I2)}
Base.promote_rule(::Type{Week{I1}}, ::Type{Hour{I2}}) where {I1, I2} = Hour{promote_type(I1, I2)}
Base.promote_rule(::Type{Week{I1}}, ::Type{Minute{I2}}) where {I1, I2} = Minute{promote_type(I1, I2)}
Base.promote_rule(::Type{Week{I1}}, ::Type{Second{I2}}) where {I1, I2} = Second{promote_type(I1, I2)}
Base.promote_rule(::Type{Week{I1}}, ::Type{Millisecond{I2}}) where {I1, I2} = Millisecond{promote_type(I1, I2)}

Base.promote_rule(::Type{Day{I1}}, ::Type{Hour{I2}}) where {I1, I2} = Hour{promote_type(I1, I2)}
Base.promote_rule(::Type{Day{I1}}, ::Type{Minute{I2}}) where {I1, I2} = Minute{promote_type(I1, I2)}
Base.promote_rule(::Type{Day{I1}}, ::Type{Second{I2}}) where {I1, I2} = Second{promote_type(I1, I2)}
Base.promote_rule(::Type{Day{I1}}, ::Type{Millisecond{I2}}) where {I1, I2} = Millisecond{promote_type(I1, I2)}

Base.promote_rule(::Type{Hour{I1}}, ::Type{Minute{I2}}) where {I1, I2} = Minute{promote_type(I1, I2)}
Base.promote_rule(::Type{Hour{I1}}, ::Type{Second{I2}}) where {I1, I2} = Second{promote_type(I1, I2)}
Base.promote_rule(::Type{Hour{I1}}, ::Type{Millisecond{I2}}) where {I1, I2} = Millisecond{promote_type(I1, I2)}

Base.promote_rule(::Type{Minute{I1}}, ::Type{Second{I2}}) where {I1, I2} = Second{promote_type(I1, I2)}
Base.promote_rule(::Type{Minute{I1}}, ::Type{Millisecond{I2}}) where {I1, I2} = Millisecond{promote_type(I1, I2)}

Base.promote_rule(::Type{Second{I1}}, ::Type{Millisecond{I2}}) where {I1, I2} = Millisecond{promote_type(I1, I2)}

# Millenium -> Century -> Year (approximate, via Day chain)
Base.promote_rule(::Type{Millenium{I1}}, ::Type{Century{I2}}) where {I1, I2} = Day{promote_type(I1, I2)}
Base.promote_rule(::Type{Millenium{I1}}, ::Type{Week{I2}}) where {I1, I2} = Day{promote_type(I1, I2)}
Base.promote_rule(::Type{Millenium{I1}}, ::Type{Day{I2}}) where {I1, I2} = Day{promote_type(I1, I2)}
Base.promote_rule(::Type{Millenium{I1}}, ::Type{Hour{I2}}) where {I1, I2} = Hour{promote_type(I1, I2)}
Base.promote_rule(::Type{Millenium{I1}}, ::Type{Second{I2}}) where {I1, I2} = Second{promote_type(I1, I2)}
Base.promote_rule(::Type{Millenium{I1}}, ::Type{Millisecond{I2}}) where {I1, I2} = Millisecond{promote_type(I1, I2)}

Base.promote_rule(::Type{Century{I1}}, ::Type{Week{I2}}) where {I1, I2} = Day{promote_type(I1, I2)}
Base.promote_rule(::Type{Century{I1}}, ::Type{Day{I2}}) where {I1, I2} = Day{promote_type(I1, I2)}
Base.promote_rule(::Type{Century{I1}}, ::Type{Hour{I2}}) where {I1, I2} = Hour{promote_type(I1, I2)}
Base.promote_rule(::Type{Century{I1}}, ::Type{Second{I2}}) where {I1, I2} = Second{promote_type(I1, I2)}
Base.promote_rule(::Type{Century{I1}}, ::Type{Millisecond{I2}}) where {I1, I2} = Millisecond{promote_type(I1, I2)}

# OtherPeriod (Month, Year) promote among themselves
Base.promote_rule(::Type{Month{I1}}, ::Type{Year{I2}}) where {I1, I2} = Month{promote_type(I1, I2)}
Base.promote_rule(::Type{Millenium{I1}}, ::Type{Year{I2}}) where {I1, I2} = Month{promote_type(I1, I2)}
Base.promote_rule(::Type{Millenium{I1}}, ::Type{Month{I2}}) where {I1, I2} = Month{promote_type(I1, I2)}
Base.promote_rule(::Type{Century{I1}}, ::Type{Year{I2}}) where {I1, I2} = Month{promote_type(I1, I2)}
Base.promote_rule(::Type{Century{I1}}, ::Type{Month{I2}}) where {I1, I2} = Month{promote_type(I1, I2)}

# Cross-type arithmetic via promotion (different period types)
for op in (:(==), :isless, :/, :rem, :mod, :lcm, :gcd)
    @eval Base.$op(x::Period, y::Period) = $op(promote(x, y)...)
end
Base.div(x::Period, y::Period, r::RoundingMode) = div(promote(x, y)..., r)
Base.div(x::Period, y::Period) = div(promote(x, y)...)

# Cross-type +/- via promotion
Base.:(+)(x::Period, y::Period) = +(promote(x, y)...)
Base.:(-)(x::Period, y::Period) = -(promote(x, y)...)

# ============================================================================
# UTINSTANT — parametric
# ============================================================================
struct UTInstant{P <: Period} <: Instant
    periods::P
end

Base.:(+)(x::UTInstant) = x
Base.:(-)(x::T, y::T) where {T <: UTInstant} = x.periods - y.periods

# ============================================================================
# DATETIME — parametric
# ============================================================================

# Calendar helpers (same algorithms as Julia Base Dates)
const _SHIFTEDMONTHDAYS = (306, 337, 0, 31, 61, 92, 122, 153, 184, 214, 245, 275)
const _DAYSINMONTH = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

_isleapyear(y) = (y % 4 == 0) && ((y % 100 != 0) || (y % 400 == 0))
isleapyear(y::Integer) = _isleapyear(y)
daysinmonth(y::Integer, m::Integer) = _DAYSINMONTH[m] + (m == 2 && _isleapyear(y))

function _totaldays(y, m, d)
    z = m < 3 ? y - 1 : y
    mdays = _SHIFTEDMONTHDAYS[m]
    return d + mdays + 365z + fld(z, 4) - fld(z, 100) + fld(z, 400) - 306
end

function _yearmonthday(days)
    z = days + 306; h = 100z - 25; a = fld(h, 3652425); b = a - fld(a, 4)
    y = fld(100b + h, 36525); c = b + z - 365y - fld(y, 4); m = div(5c + 456, 153)
    d = c - div(153m - 457, 5); return m > 12 ? (y + 1, m - 12, d) : (y, m, d)
end

function _year(days)
    z = days + 306; h = 100z - 25; a = fld(h, 3652425); b = a - fld(a, 4)
    y = fld(100b + h, 36525); c = b + z - 365y - fld(y, 4); m = div(5c + 456, 153)
    return m > 12 ? y + 1 : y
end

function _month(days)
    z = days + 306; h = 100z - 25; a = fld(h, 3652425); b = a - fld(a, 4)
    y = fld(100b + h, 36525); c = b + z - 365y - fld(y, 4); m = div(5c + 456, 153)
    return m > 12 ? m - 12 : m
end

function _day(days)
    z = days + 306; h = 100z - 25; a = fld(h, 3652425); b = a - fld(a, 4)
    y = fld(100b + h, 36525); c = b + z - 365y - fld(y, 4); m = div(5c + 456, 153)
    return c - div(153m - 457, 5)
end

function _dayofyear(days)
    y, m, d = _yearmonthday(days)
    return _totaldays(y, m, d) - _totaldays(y, 1, 1) + 1
end

"""
    DateTime{IntType}

A parametric DateTime type where the underlying millisecond counter is of type `IntType`.
Mirrors `Dates.DateTime` but allows non-Int64 integer types for Reactant compatibility."""
struct DateTime{IntType} <: AbstractDateTime
    instant::UTInstant{Millisecond{IntType}}
    DateTime{IntType}(instant::UTInstant{Millisecond{IntType}}) where {IntType} = new{IntType}(instant)
end

# Convenience: infer IntType from instant
DateTime(instant::UTInstant{Millisecond{IntType}}) where {IntType} = DateTime{IntType}(instant)

# UTM helper
_UTM(x) = UTInstant(Millisecond(x))

# value accessor for DateTime (returns milliseconds since epoch)
value(dt::DateTime) = value(dt.instant.periods)

# Core constructor from parts
function DateTime(
        y::Int64, m::Int64 = 1, d::Int64 = 1,
        h::Int64 = 0, mi::Int64 = 0, s::Int64 = 0, ms::Int64 = 0
    )
    rata = ms + 1000 * (s + 60mi + 3600h + 86400 * _totaldays(y, m, d))
    return DateTime(_UTM(rata))
end

# Fallback constructor from any integers
DateTime(y, m = 1, d = 1, h = 0, mi = 0, s = 0, ms = 0) =
    DateTime(Int64(y), Int64(m), Int64(d), Int64(h), Int64(mi), Int64(s), Int64(ms))

# Convenience constructors from Period types
function DateTime(
        y::Year, m::Month = Month(1), d::Day = Day(1),
        h::Hour = Hour(0), mi::Minute = Minute(0),
        s::Second = Second(0), ms::Millisecond = Millisecond(0)
    )
    return DateTime(
        value(y), value(m), value(d),
        value(h), value(mi), value(s), value(ms)
    )
end

# ============================================================================
# DATETIME ACCESSORS
# ============================================================================
function year(dt::DateTime)
    days = fld(value(dt), 86400000)
    return _year(days)
end

function month(dt::DateTime)
    days = fld(value(dt), 86400000)
    return _month(days)
end

function day(dt::DateTime)
    days = fld(value(dt), 86400000)
    return _day(days)
end

function yearmonthday(dt::DateTime)
    days = fld(value(dt), 86400000)
    return _yearmonthday(days)
end

function hour(dt::DateTime)
    ms = value(dt)
    return mod(fld(ms, 3600000), 24)
end

function minute(dt::DateTime)
    ms = value(dt)
    return mod(fld(ms, 60000), 60)
end

function second(dt::DateTime)
    ms = value(dt)
    return mod(fld(ms, 1000), 60)
end

function millisecond(dt::DateTime)
    ms = value(dt)
    return mod(ms, 1000)
end

function dayofyear(dt::DateTime)
    days = fld(value(dt), 86400000)
    return _dayofyear(days)
end

# seconds since midnight (used by zenith calculations)
function secondofday(dt::DateTime)
    ms = value(dt)
    return mod(fld(ms, 1000), 86400)
end

# ============================================================================
# DATETIME PRINTING
# ============================================================================
_pad2(x) = lpad(x, 2, '0')
_pad3(x) = lpad(x, 3, '0')
_pad4(x) = lpad(x, 4, '0')

function Base.show(io::IO, dt::DateTime)
    y, m, d = yearmonthday(dt)
    h = hour(dt)
    mi = minute(dt)
    s = second(dt)
    ms = millisecond(dt)
    print(io, _pad4(y), '-', _pad2(m), '-', _pad2(d), 'T', _pad2(h), ':', _pad2(mi), ':', _pad2(s))
    return ms != 0 && print(io, '.', _pad3(ms))
end

Base.show(io::IO, ::MIME"text/plain", dt::DateTime) = show(io, dt)

# ============================================================================
# DATETIME EQUALITY AND COMPARISON
# ============================================================================
Base.:(==)(x::DateTime, y::DateTime) = value(x) == value(y)
Base.isless(x::DateTime, y::DateTime) = isless(value(x), value(y))
Base.isfinite(::Union{Type{<:DateTime}, DateTime}) = true
Base.zero(::Type{DateTime}) = DateTime(_UTM(0))
Base.zero(::Type{DateTime{I}}) where {I} = DateTime{I}(UTInstant(Millisecond{I}(zero(I))))
Base.zero(::DateTime{I}) where {I} = DateTime{I}(UTInstant(Millisecond{I}(zero(I))))

# ============================================================================
# DATETIME ARITHMETIC
# ============================================================================
# DateTime - DateTime -> Millisecond
Base.:(-)(x::DateTime, y::DateTime) = Millisecond(value(x) - value(y))

# DateTime +/- FixedPeriod (via milliseconds)
Base.:(+)(x::DateTime, y::FixedPeriod) = DateTime(_UTM(value(x) + toms(y)))
Base.:(-)(x::DateTime, y::FixedPeriod) = DateTime(_UTM(value(x) - toms(y)))
Base.:(+)(y::FixedPeriod, x::DateTime) = x + y

# DateTime + Year
function Base.:(+)(dt::DateTime, y::Year)
    oy, m, d = yearmonthday(dt)
    ny = oy + value(y)
    ld = daysinmonth(ny, m)
    return DateTime(ny, m, d <= ld ? d : ld, hour(dt), minute(dt), second(dt), millisecond(dt))
end

function Base.:(-)(dt::DateTime, y::Year)
    oy, m, d = yearmonthday(dt)
    ny = oy - value(y)
    ld = daysinmonth(ny, m)
    return DateTime(ny, m, d <= ld ? d : ld, hour(dt), minute(dt), second(dt), millisecond(dt))
end

# DateTime + Month
_monthwrap(m1, m2) = (v = mod1(m1 + m2, 12); return v < 0 ? 12 + v : v)
_yearwrap(y, m1, m2) = y + fld(m1 + m2 - 1, 12)

function Base.:(+)(dt::DateTime, z::Month)
    y, m, d = yearmonthday(dt)
    ny = _yearwrap(y, m, value(z))
    mm = _monthwrap(m, value(z))
    ld = daysinmonth(ny, mm)
    return DateTime(ny, mm, d <= ld ? d : ld, hour(dt), minute(dt), second(dt), millisecond(dt))
end

function Base.:(-)(dt::DateTime, z::Month)
    y, m, d = yearmonthday(dt)
    ny = _yearwrap(y, m, -value(z))
    mm = _monthwrap(m, -value(z))
    ld = daysinmonth(ny, mm)
    return DateTime(ny, mm, d <= ld ? d : ld, hour(dt), minute(dt), second(dt), millisecond(dt))
end

Base.:(+)(y::Year, x::DateTime) = x + y
Base.:(+)(y::Month, x::DateTime) = x + y

# Multiple periods
Base.:(+)(a::DateTime, b::Period, c::Period) = a + b + c
Base.:(+)(a::DateTime, b::Period, c::Period, d::Period...) = (+)((a + b + c), d...)

# ============================================================================
# ADDITIONAL DATETIME HELPERS
# ============================================================================
# days since start of month (as integer)
function days(ms::Millisecond)
    return div(value(ms), 86400000)
end

# first day of the month for a given DateTime
function firstdayofmonth(dt::DateTime)
    y, m, _ = yearmonthday(dt)
    return DateTime(y, m, 1)
end

# second accessor for Millisecond (used in schedule.jl: Dates.second(every_n_timesteps * clock.Δt))
second(ms::Millisecond) = div(value(ms), 1000)

# ============================================================================
# FORMAT HELPER (minimal, for output/feedback)
# ============================================================================
function format(dt::DateTime, fmt::String = "yyyy-mm-ddTHH:MM:SS")
    y, m, d = yearmonthday(dt)
    h = hour(dt)
    mi = minute(dt)
    s = second(dt)
    result = replace(fmt, "yyyy" => _pad4(y))
    result = replace(result, "mm" => _pad2(m))
    result = replace(result, "dd" => _pad2(d))
    result = replace(result, "HH" => _pad2(h))
    result = replace(result, "MM" => _pad2(mi))
    result = replace(result, "SS" => _pad2(s))
    return result
end

# RFC1123-like format for now() replacement
function now()
    # Return a fixed "now" — in practice this is only used for logging
    # and will be replaced by the caller with the actual system time
    error("now() is not available without the Dates package. Use a fixed DateTime or pass time explicitly.")
end

# ============================================================================
# BROADCASTING SUPPORT
# ============================================================================
Base.Broadcast.broadcastable(x::Union{Period, UTInstant, DateTime}) = Ref(x)

# ============================================================================
# COMPOUND PERIOD AND CANONICALIZE
# ============================================================================
"""
    CompoundPeriod

A combination of `Period` types representing a canonicalized duration.
Mirrors `Dates.CompoundPeriod`."""
struct CompoundPeriod
    periods::Vector{Period}
end

CompoundPeriod(ps::Period...) = CompoundPeriod(collect(Period, ps))

function Base.:(==)(x::CompoundPeriod, y::CompoundPeriod)
    # compare by converting both to total milliseconds
    return _compound_toms(x) == _compound_toms(y)
end

function _compound_toms(cp::CompoundPeriod)
    total = 0
    for p in cp.periods
        total += toms(p)
    end
    return total
end

function Base.show(io::IO, cp::CompoundPeriod)
    isempty(cp.periods) && return print(io, "empty period")
    first = true
    for p in cp.periods
        first || print(io, ", ")
        print(io, p)
        first = false
    end
    return
end

Base.show(io::IO, ::MIME"text/plain", cp::CompoundPeriod) = show(io, cp)

"""
    canonicalize(cp::CompoundPeriod)

Reduce a `CompoundPeriod` to its canonical form by converting
milliseconds up through the period hierarchy.
Mirrors `Dates.canonicalize`."""
function canonicalize(cp::CompoundPeriod)
    ms = _compound_toms(cp)
    return _canonicalize_ms(ms)
end

function _canonicalize_ms(ms::Integer)
    neg = ms < 0
    ms = abs(ms)
    periods = Period[]

    d, ms = divrem(ms, 86400000)
    d != 0 && push!(periods, Day(neg ? -d : d))

    h, ms = divrem(ms, 3600000)
    h != 0 && push!(periods, Hour(neg ? -h : h))

    m, ms = divrem(ms, 60000)
    m != 0 && push!(periods, Minute(neg ? -m : m))

    s, ms = divrem(ms, 1000)
    s != 0 && push!(periods, Second(neg ? -s : s))

    ms != 0 && push!(periods, Millisecond(neg ? -ms : ms))

    isempty(periods) && push!(periods, Millisecond(0))
    return CompoundPeriod(periods)
end

"""
    canonicalize(p::Period)

Canonicalize a single period by converting to milliseconds first."""
canonicalize(p::FixedPeriod) = _canonicalize_ms(toms(p))

"""
    round(p::Millisecond, ::Type{T}) where T <: Period

Round a `Millisecond` duration to the nearest multiple of period type `T`.
Returns a `CompoundPeriod`. Mirrors `Base.round(::Dates.Millisecond, ::Type)`."""
function Base.round(p::Millisecond, ::Type{T}) where {T <: FixedPeriod}
    target_ms = toms(T(1))
    rounded = Millisecond(div(value(p) + div(target_ms, 2), target_ms) * target_ms)
    return canonicalize(rounded)
end

"""
    round(p::Millisecond, target::FixedPeriod)

Round a `Millisecond` duration to the nearest multiple of `target`.
Returns a `CompoundPeriod`. Mirrors `Base.round(::Dates.Millisecond, ::Dates.Period)`."""
function Base.round(p::Millisecond, target::FixedPeriod)
    target_ms = toms(target)
    rounded = Millisecond(div(value(p) + div(target_ms, 2), target_ms) * target_ms)
    return canonicalize(rounded)
end

end # module TracableDates
