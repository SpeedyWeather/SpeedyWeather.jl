# CUSTOM DATES IMPLEMENTATION WITH PARAMETRIC VALUE TYPES
# Mirrors Julia Base Dates but with `value::IntType` instead of `value::Int64`
# so that Reactant can trace through the integer operations.
# This code is mostly LLM generated, instructed to copy the Julia.Base implementation 
# I really dislike that we have to do this here. 

# ============================================================================
# PERIOD TYPES — parametric on IntType, subtyping Dates abstract types
# ============================================================================
struct SpeedyMillisecond{IntType} <: Dates.TimePeriod
    value::IntType
end

struct SpeedySecond{IntType} <: Dates.TimePeriod
    value::IntType
end

struct SpeedyMinute{IntType} <: Dates.TimePeriod
    value::IntType
end

struct SpeedyHour{IntType} <: Dates.TimePeriod
    value::IntType
end

struct SpeedyDay{IntType} <: Dates.DatePeriod
    value::IntType
end

struct SpeedyWeek{IntType} <: Dates.DatePeriod
    value::IntType
end

struct SpeedyMonth{IntType} <: Dates.DatePeriod
    value::IntType
end

struct SpeedyYear{IntType} <: Dates.DatePeriod
    value::IntType
end

# Union alias for dispatch on all Speedy period types
const SpeedyPeriod = Union{SpeedyMillisecond, SpeedySecond, SpeedyMinute, SpeedyHour, SpeedyDay, SpeedyWeek, SpeedyMonth, SpeedyYear}

# Convenience constructors: infer IntType from input
SpeedyMillisecond(v::Integer) = SpeedyMillisecond{typeof(v)}(v)
SpeedySecond(v::Integer) = SpeedySecond{typeof(v)}(v)
SpeedyMinute(v::Integer) = SpeedyMinute{typeof(v)}(v)
SpeedyHour(v::Integer) = SpeedyHour{typeof(v)}(v)
SpeedyDay(v::Integer) = SpeedyDay{typeof(v)}(v)
SpeedyWeek(v::Integer) = SpeedyWeek{typeof(v)}(v)
SpeedyMonth(v::Integer) = SpeedyMonth{typeof(v)}(v)
SpeedyYear(v::Integer) = SpeedyYear{typeof(v)}(v)

# ============================================================================
# VALUE ACCESSOR — extend Dates.value
# ============================================================================
Dates.value(x::SpeedyPeriod) = x.value

# ============================================================================
# PRINT / SHOW
# ============================================================================
_units(x::SpeedyMillisecond) = " millisecond" * (abs(Dates.value(x)) == 1 ? "" : "s")
_units(x::SpeedySecond) = " second" * (abs(Dates.value(x)) == 1 ? "" : "s")
_units(x::SpeedyMinute) = " minute" * (abs(Dates.value(x)) == 1 ? "" : "s")
_units(x::SpeedyHour) = " hour" * (abs(Dates.value(x)) == 1 ? "" : "s")
_units(x::SpeedyDay) = " day" * (abs(Dates.value(x)) == 1 ? "" : "s")
_units(x::SpeedyWeek) = " week" * (abs(Dates.value(x)) == 1 ? "" : "s")
_units(x::SpeedyMonth) = " month" * (abs(Dates.value(x)) == 1 ? "" : "s")
_units(x::SpeedyYear) = " year" * (abs(Dates.value(x)) == 1 ? "" : "s")

Base.print(io::IO, x::SpeedyPeriod) = print(io, Dates.value(x), _units(x))
Base.show(io::IO, ::MIME"text/plain", x::SpeedyPeriod) = print(io, x)
Base.show(io::IO, p::P) where {P <: SpeedyPeriod} = print(io, nameof(P), '(', Dates.value(p), ')')

# ============================================================================
# ZERO / ONE / TYPEMIN / TYPEMAX / ISFINITE
# ============================================================================
Base.zero(::Type{P}) where {P <: SpeedyPeriod} = P(0)
Base.zero(::P) where {P <: SpeedyPeriod} = P(0)
Base.one(::Union{Type{P}, P}) where {P <: SpeedyPeriod} = 1
Base.isfinite(::Union{Type{P}, P}) where {P <: SpeedyPeriod} = true
Base.iszero(x::SpeedyPeriod) = iszero(Dates.value(x))

# ============================================================================
# PERIOD ARITHMETIC
# ============================================================================
# Unary
Base.:(-)(x::P) where {P <: SpeedyPeriod} = P(-Dates.value(x))
Base.:(+)(x::SpeedyPeriod) = x
Base.abs(x::P) where {P <: SpeedyPeriod} = P(abs(Dates.value(x)))
Base.sign(x::SpeedyPeriod) = sign(Dates.value(x))
Base.signbit(x::SpeedyPeriod) = signbit(Dates.value(x))

# Same-type binary arithmetic (P matches exact concrete type including IntType)
Base.:(==)(x::P, y::P) where {P <: SpeedyPeriod} = Dates.value(x) == Dates.value(y)
Base.isless(x::P, y::P) where {P <: SpeedyPeriod} = isless(Dates.value(x), Dates.value(y))

for op in (:+, :-, :lcm, :gcd)
    @eval Base.$op(x::P, y::P) where {P <: SpeedyPeriod} = P($op(Dates.value(x), Dates.value(y)))
end

Base.:(/)(x::P, y::P) where {P <: SpeedyPeriod} = Dates.value(x) / Dates.value(y)
Base.:(/)(x::P, y::Real) where {P <: SpeedyPeriod} = P(Dates.value(x) / y)
Base.div(x::P, y::P, r::RoundingMode) where {P <: SpeedyPeriod} = div(Dates.value(x), Dates.value(y), r)
Base.div(x::P, y::P) where {P <: SpeedyPeriod} = div(Dates.value(x), Dates.value(y))

for op in (:rem, :mod)
    @eval begin
        Base.$op(x::P, y::P) where {P <: SpeedyPeriod} = P($op(Dates.value(x), Dates.value(y)))
        Base.$op(x::P, y::Real) where {P <: SpeedyPeriod} = P($op(Dates.value(x), y))
    end
end

Base.:(*)(x::P, y::Real) where {P <: SpeedyPeriod} = P(Dates.value(x) * y)
Base.:(*)(y::Real, x::SpeedyPeriod) = x * y
Base.:(*)(x::P, y::Integer) where {P <: SpeedyPeriod} = P(Dates.value(x) * y)
Base.:(*)(y::Integer, x::SpeedyPeriod) = x * y

# ============================================================================
# CONVERSION FACTORS (FixedPeriod chain)
# ============================================================================
# toms: convert any FixedPeriod to milliseconds (extends Dates.toms)
Dates.toms(x::SpeedyMillisecond) = Dates.value(x)
Dates.toms(x::SpeedySecond) = Dates.value(x) * 1000
Dates.toms(x::SpeedyMinute) = Dates.value(x) * 60000
Dates.toms(x::SpeedyHour) = Dates.value(x) * 3600000
Dates.toms(x::SpeedyDay) = Dates.value(x) * 86400000
Dates.toms(x::SpeedyWeek) = Dates.value(x) * 604800000

# ============================================================================
# CONVERSIONS BETWEEN FIXED PERIOD TYPES
# ============================================================================
function _divexact(x, y)
    q, r = divrem(x, y)
    r == 0 || throw(InexactError(:_divexact, Int, x / y))
    return q
end

# Millisecond -> coarser
Base.convert(::Type{<:SpeedySecond}, x::SpeedyMillisecond) = SpeedySecond(_divexact(Dates.value(x), 1000))
Base.convert(::Type{<:SpeedyMinute}, x::SpeedyMillisecond) = SpeedyMinute(_divexact(Dates.value(x), 60000))
Base.convert(::Type{<:SpeedyHour}, x::SpeedyMillisecond) = SpeedyHour(_divexact(Dates.value(x), 3600000))
Base.convert(::Type{<:SpeedyDay}, x::SpeedyMillisecond) = SpeedyDay(_divexact(Dates.value(x), 86400000))
Base.convert(::Type{<:SpeedyWeek}, x::SpeedyMillisecond) = SpeedyWeek(_divexact(Dates.value(x), 604800000))

# Second -> coarser / finer
Base.convert(::Type{<:SpeedyMillisecond}, x::SpeedySecond) = SpeedyMillisecond(Dates.value(x) * 1000)
Base.convert(::Type{<:SpeedyMinute}, x::SpeedySecond) = SpeedyMinute(_divexact(Dates.value(x), 60))
Base.convert(::Type{<:SpeedyHour}, x::SpeedySecond) = SpeedyHour(_divexact(Dates.value(x), 3600))
Base.convert(::Type{<:SpeedyDay}, x::SpeedySecond) = SpeedyDay(_divexact(Dates.value(x), 86400))
Base.convert(::Type{<:SpeedyWeek}, x::SpeedySecond) = SpeedyWeek(_divexact(Dates.value(x), 604800))

# Minute -> coarser / finer
Base.convert(::Type{<:SpeedyMillisecond}, x::SpeedyMinute) = SpeedyMillisecond(Dates.value(x) * 60000)
Base.convert(::Type{<:SpeedySecond}, x::SpeedyMinute) = SpeedySecond(Dates.value(x) * 60)
Base.convert(::Type{<:SpeedyHour}, x::SpeedyMinute) = SpeedyHour(_divexact(Dates.value(x), 60))
Base.convert(::Type{<:SpeedyDay}, x::SpeedyMinute) = SpeedyDay(_divexact(Dates.value(x), 1440))
Base.convert(::Type{<:SpeedyWeek}, x::SpeedyMinute) = SpeedyWeek(_divexact(Dates.value(x), 10080))

# Hour -> coarser / finer
Base.convert(::Type{<:SpeedyMillisecond}, x::SpeedyHour) = SpeedyMillisecond(Dates.value(x) * 3600000)
Base.convert(::Type{<:SpeedySecond}, x::SpeedyHour) = SpeedySecond(Dates.value(x) * 3600)
Base.convert(::Type{<:SpeedyMinute}, x::SpeedyHour) = SpeedyMinute(Dates.value(x) * 60)
Base.convert(::Type{<:SpeedyDay}, x::SpeedyHour) = SpeedyDay(_divexact(Dates.value(x), 24))
Base.convert(::Type{<:SpeedyWeek}, x::SpeedyHour) = SpeedyWeek(_divexact(Dates.value(x), 168))

# Day -> coarser / finer
Base.convert(::Type{<:SpeedyMillisecond}, x::SpeedyDay) = SpeedyMillisecond(Dates.value(x) * 86400000)
Base.convert(::Type{<:SpeedySecond}, x::SpeedyDay) = SpeedySecond(Dates.value(x) * 86400)
Base.convert(::Type{<:SpeedyMinute}, x::SpeedyDay) = SpeedyMinute(Dates.value(x) * 1440)
Base.convert(::Type{<:SpeedyHour}, x::SpeedyDay) = SpeedyHour(Dates.value(x) * 24)
Base.convert(::Type{<:SpeedyWeek}, x::SpeedyDay) = SpeedyWeek(_divexact(Dates.value(x), 7))

# Week -> finer
Base.convert(::Type{<:SpeedyMillisecond}, x::SpeedyWeek) = SpeedyMillisecond(Dates.value(x) * 604800000)
Base.convert(::Type{<:SpeedySecond}, x::SpeedyWeek) = SpeedySecond(Dates.value(x) * 604800)
Base.convert(::Type{<:SpeedyMinute}, x::SpeedyWeek) = SpeedyMinute(Dates.value(x) * 10080)
Base.convert(::Type{<:SpeedyHour}, x::SpeedyWeek) = SpeedyHour(Dates.value(x) * 168)
Base.convert(::Type{<:SpeedyDay}, x::SpeedyWeek) = SpeedyDay(Dates.value(x) * 7)

# Year <-> Month conversions
Base.convert(::Type{<:SpeedyMonth}, x::SpeedyYear) = SpeedyMonth(Dates.value(x) * 12)
Base.convert(::Type{<:SpeedyYear}, x::SpeedyMonth) = SpeedyYear(_divexact(Dates.value(x), 12))

# Identity conversions
Base.convert(::Type{P}, x::P) where {P <: SpeedyPeriod} = x

# Construct from Period (type conversion) — explicit per-type to avoid ambiguity
SpeedyMillisecond(p::SpeedyPeriod) = convert(SpeedyMillisecond, p)
SpeedySecond(p::SpeedyPeriod) = convert(SpeedySecond, p)
SpeedyMinute(p::SpeedyPeriod) = convert(SpeedyMinute, p)
SpeedyHour(p::SpeedyPeriod) = convert(SpeedyHour, p)
SpeedyDay(p::SpeedyPeriod) = convert(SpeedyDay, p)
SpeedyWeek(p::SpeedyPeriod) = convert(SpeedyWeek, p)
SpeedyMonth(p::SpeedyPeriod) = convert(SpeedyMonth, p)
SpeedyYear(p::SpeedyPeriod) = convert(SpeedyYear, p)

# ============================================================================
# PROMOTION RULES (FixedPeriod: promote to finer resolution)
# ============================================================================
const SpeedyFixedPeriod = Union{SpeedyWeek, SpeedyDay, SpeedyHour, SpeedyMinute, SpeedySecond, SpeedyMillisecond}

Base.promote_rule(::Type{SpeedyWeek{I1}}, ::Type{SpeedyDay{I2}}) where {I1, I2} = SpeedyDay{promote_type(I1, I2)}
Base.promote_rule(::Type{SpeedyWeek{I1}}, ::Type{SpeedyHour{I2}}) where {I1, I2} = SpeedyHour{promote_type(I1, I2)}
Base.promote_rule(::Type{SpeedyWeek{I1}}, ::Type{SpeedyMinute{I2}}) where {I1, I2} = SpeedyMinute{promote_type(I1, I2)}
Base.promote_rule(::Type{SpeedyWeek{I1}}, ::Type{SpeedySecond{I2}}) where {I1, I2} = SpeedySecond{promote_type(I1, I2)}
Base.promote_rule(::Type{SpeedyWeek{I1}}, ::Type{SpeedyMillisecond{I2}}) where {I1, I2} = SpeedyMillisecond{promote_type(I1, I2)}

Base.promote_rule(::Type{SpeedyDay{I1}}, ::Type{SpeedyHour{I2}}) where {I1, I2} = SpeedyHour{promote_type(I1, I2)}
Base.promote_rule(::Type{SpeedyDay{I1}}, ::Type{SpeedyMinute{I2}}) where {I1, I2} = SpeedyMinute{promote_type(I1, I2)}
Base.promote_rule(::Type{SpeedyDay{I1}}, ::Type{SpeedySecond{I2}}) where {I1, I2} = SpeedySecond{promote_type(I1, I2)}
Base.promote_rule(::Type{SpeedyDay{I1}}, ::Type{SpeedyMillisecond{I2}}) where {I1, I2} = SpeedyMillisecond{promote_type(I1, I2)}

Base.promote_rule(::Type{SpeedyHour{I1}}, ::Type{SpeedyMinute{I2}}) where {I1, I2} = SpeedyMinute{promote_type(I1, I2)}
Base.promote_rule(::Type{SpeedyHour{I1}}, ::Type{SpeedySecond{I2}}) where {I1, I2} = SpeedySecond{promote_type(I1, I2)}
Base.promote_rule(::Type{SpeedyHour{I1}}, ::Type{SpeedyMillisecond{I2}}) where {I1, I2} = SpeedyMillisecond{promote_type(I1, I2)}

Base.promote_rule(::Type{SpeedyMinute{I1}}, ::Type{SpeedySecond{I2}}) where {I1, I2} = SpeedySecond{promote_type(I1, I2)}
Base.promote_rule(::Type{SpeedyMinute{I1}}, ::Type{SpeedyMillisecond{I2}}) where {I1, I2} = SpeedyMillisecond{promote_type(I1, I2)}

Base.promote_rule(::Type{SpeedySecond{I1}}, ::Type{SpeedyMillisecond{I2}}) where {I1, I2} = SpeedyMillisecond{promote_type(I1, I2)}

# OtherPeriod (Month, Year) promote among themselves
Base.promote_rule(::Type{SpeedyMonth{I1}}, ::Type{SpeedyYear{I2}}) where {I1, I2} = SpeedyMonth{promote_type(I1, I2)}

# Cross-type arithmetic via promotion (different period types)
for op in (:(==), :isless, :/, :rem, :mod, :lcm, :gcd)
    @eval Base.$op(x::SpeedyPeriod, y::SpeedyPeriod) = $op(promote(x, y)...)
end
Base.div(x::SpeedyPeriod, y::SpeedyPeriod, r::RoundingMode) = div(promote(x, y)..., r)
Base.div(x::SpeedyPeriod, y::SpeedyPeriod) = div(promote(x, y)...)

# Cross-type +/- via promotion
Base.:(+)(x::SpeedyPeriod, y::SpeedyPeriod) = +(promote(x, y)...)
Base.:(-)(x::SpeedyPeriod, y::SpeedyPeriod) = -(promote(x, y)...)

# ============================================================================
# UTINSTANT — parametric
# ============================================================================
struct SpeedyUTInstant{P<:Dates.Period} <: Dates.Instant
    periods::P
end

Base.:(+)(x::SpeedyUTInstant) = x
Base.:(-)(x::T, y::T) where {T <: SpeedyUTInstant} = x.periods - y.periods

# ============================================================================
# DATETIME — parametric
# ============================================================================

# Calendar helpers (same algorithms as Julia Base Dates)
const _SHIFTEDMONTHDAYS = (306, 337, 0, 31, 61, 92, 122, 153, 184, 214, 245, 275)
const _DAYSINMONTH = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

# Dates.isleapyear and _daysinmonth already exist for Integer in Base Dates.
# We use a local _daysinmonth that works on plain integers for our internal helpers.
_isleapyear(y) = (y % 4 == 0) && ((y % 100 != 0) || (y % 400 == 0))
_daysinmonth(y, m) = _DAYSINMONTH[m] + (m == 2 && _isleapyear(y))

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
    SpeedyDateTime{IntType}

A parametric DateTime type where the underlying millisecond counter is of type `IntType`.
Mirrors `Dates.DateTime` but allows non-Int64 integer types for Reactant compatibility."""
struct SpeedyDateTime{IntType} <: Dates.AbstractDateTime
    instant::SpeedyUTInstant{SpeedyMillisecond{IntType}}
    SpeedyDateTime{IntType}(instant::SpeedyUTInstant{SpeedyMillisecond{IntType}}) where {IntType} = new{IntType}(instant)
end

# Convenience: infer IntType from instant
SpeedyDateTime(instant::SpeedyUTInstant{SpeedyMillisecond{IntType}}) where {IntType} = SpeedyDateTime{IntType}(instant)

# UTM helper
_SpeedyUTM(x) = SpeedyUTInstant(SpeedyMillisecond(x))

# value accessor for DateTime (returns milliseconds since epoch)
Dates.value(dt::SpeedyDateTime) = Dates.value(dt.instant.periods)

# Core constructor from parts
function SpeedyDateTime(y::Int64, m::Int64=1, d::Int64=1,
                    h::Int64=0, mi::Int64=0, s::Int64=0, ms::Int64=0)
    rata = ms + 1000 * (s + 60mi + 3600h + 86400 * _totaldays(y, m, d))
    return SpeedyDateTime(_SpeedyUTM(rata))
end

# Fallback constructor from any integers
SpeedyDateTime(y, m=1, d=1, h=0, mi=0, s=0, ms=0) =
    SpeedyDateTime(Int64(y), Int64(m), Int64(d), Int64(h), Int64(mi), Int64(s), Int64(ms))

# Convenience constructors from Period types
function SpeedyDateTime(y::SpeedyYear, m::SpeedyMonth=SpeedyMonth(1), d::SpeedyDay=SpeedyDay(1),
                    h::SpeedyHour=SpeedyHour(0), mi::SpeedyMinute=SpeedyMinute(0),
                    s::SpeedySecond=SpeedySecond(0), ms::SpeedyMillisecond=SpeedyMillisecond(0))
    return SpeedyDateTime(Dates.value(y), Dates.value(m), Dates.value(d),
                      Dates.value(h), Dates.value(mi), Dates.value(s), Dates.value(ms))
end

# ============================================================================
# DATETIME ACCESSORS — extend Dates.year, Dates.month, etc.
# ============================================================================
function Dates.year(dt::SpeedyDateTime)
    days = fld(Dates.value(dt), 86400000)
    return _year(days)
end

function Dates.month(dt::SpeedyDateTime)
    days = fld(Dates.value(dt), 86400000)
    return _month(days)
end

function Dates.day(dt::SpeedyDateTime)
    days = fld(Dates.value(dt), 86400000)
    return _day(days)
end

function Dates.yearmonthday(dt::SpeedyDateTime)
    days = fld(Dates.value(dt), 86400000)
    return _yearmonthday(days)
end

function Dates.hour(dt::SpeedyDateTime)
    ms = Dates.value(dt)
    return mod(fld(ms, 3600000), 24)
end

function Dates.minute(dt::SpeedyDateTime)
    ms = Dates.value(dt)
    return mod(fld(ms, 60000), 60)
end

function Dates.second(dt::SpeedyDateTime)
    ms = Dates.value(dt)
    return mod(fld(ms, 1000), 60)
end

function Dates.millisecond(dt::SpeedyDateTime)
    ms = Dates.value(dt)
    return mod(ms, 1000)
end

function Dates.dayofyear(dt::SpeedyDateTime)
    days = fld(Dates.value(dt), 86400000)
    return _dayofyear(days)
end

# seconds since midnight (used by zenith calculations)
function secondofday(dt::SpeedyDateTime)
    ms = Dates.value(dt)
    return mod(fld(ms, 1000), 86400)
end

# ============================================================================
# DATETIME PRINTING
# ============================================================================
_pad2(x) = lpad(x, 2, '0')
_pad3(x) = lpad(x, 3, '0')
_pad4(x) = lpad(x, 4, '0')

function Base.show(io::IO, dt::SpeedyDateTime)
    y, m, d = Dates.yearmonthday(dt)
    h = Dates.hour(dt)
    mi = Dates.minute(dt)
    s = Dates.second(dt)
    ms = Dates.millisecond(dt)
    print(io, _pad4(y), '-', _pad2(m), '-', _pad2(d), 'T', _pad2(h), ':', _pad2(mi), ':', _pad2(s))
    ms != 0 && print(io, '.', _pad3(ms))
end

Base.show(io::IO, ::MIME"text/plain", dt::SpeedyDateTime) = show(io, dt)

# ============================================================================
# DATETIME EQUALITY AND COMPARISON
# ============================================================================
Base.:(==)(x::SpeedyDateTime, y::SpeedyDateTime) = Dates.value(x) == Dates.value(y)
Base.isless(x::SpeedyDateTime, y::SpeedyDateTime) = isless(Dates.value(x), Dates.value(y))
Base.isfinite(::Union{Type{<:SpeedyDateTime}, SpeedyDateTime}) = true
Base.zero(::Type{SpeedyDateTime}) = SpeedyMillisecond(0)
Base.zero(::Type{SpeedyDateTime{I}}) where {I} = SpeedyMillisecond{I}(0)
Base.zero(::SpeedyDateTime{I}) where {I} = SpeedyMillisecond{I}(0)

# ============================================================================
# DATETIME ARITHMETIC
# ============================================================================
# DateTime - DateTime -> Millisecond
Base.:(-)(x::SpeedyDateTime, y::SpeedyDateTime) = SpeedyMillisecond(Dates.value(x) - Dates.value(y))

# DateTime +/- FixedPeriod (via milliseconds)
Base.:(+)(x::SpeedyDateTime, y::SpeedyFixedPeriod) = SpeedyDateTime(_SpeedyUTM(Dates.value(x) + Dates.toms(y)))
Base.:(-)(x::SpeedyDateTime, y::SpeedyFixedPeriod) = SpeedyDateTime(_SpeedyUTM(Dates.value(x) - Dates.toms(y)))
Base.:(+)(y::SpeedyFixedPeriod, x::SpeedyDateTime) = x + y

# DateTime + Year
function Base.:(+)(dt::SpeedyDateTime, y::SpeedyYear)
    oy, m, d = Dates.yearmonthday(dt)
    ny = oy + Dates.value(y)
    ld = _daysinmonth(ny, m)
    return SpeedyDateTime(ny, m, d <= ld ? d : ld, Dates.hour(dt), Dates.minute(dt), Dates.second(dt), Dates.millisecond(dt))
end

function Base.:(-)(dt::SpeedyDateTime, y::SpeedyYear)
    oy, m, d = Dates.yearmonthday(dt)
    ny = oy - Dates.value(y)
    ld = _daysinmonth(ny, m)
    return SpeedyDateTime(ny, m, d <= ld ? d : ld, Dates.hour(dt), Dates.minute(dt), Dates.second(dt), Dates.millisecond(dt))
end

# DateTime + Month
_monthwrap(m1, m2) = (v = mod1(m1 + m2, 12); return v < 0 ? 12 + v : v)
_yearwrap(y, m1, m2) = y + fld(m1 + m2 - 1, 12)

function Base.:(+)(dt::SpeedyDateTime, z::SpeedyMonth)
    y, m, d = Dates.yearmonthday(dt)
    ny = _yearwrap(y, m, Dates.value(z))
    mm = _monthwrap(m, Dates.value(z))
    ld = _daysinmonth(ny, mm)
    return SpeedyDateTime(ny, mm, d <= ld ? d : ld, Dates.hour(dt), Dates.minute(dt), Dates.second(dt), Dates.millisecond(dt))
end

function Base.:(-)(dt::SpeedyDateTime, z::SpeedyMonth)
    y, m, d = Dates.yearmonthday(dt)
    ny = _yearwrap(y, m, -Dates.value(z))
    mm = _monthwrap(m, -Dates.value(z))
    ld = _daysinmonth(ny, mm)
    return SpeedyDateTime(ny, mm, d <= ld ? d : ld, Dates.hour(dt), Dates.minute(dt), Dates.second(dt), Dates.millisecond(dt))
end

Base.:(+)(y::SpeedyYear, x::SpeedyDateTime) = x + y
Base.:(+)(y::SpeedyMonth, x::SpeedyDateTime) = x + y

# Multiple periods
Base.:(+)(a::SpeedyDateTime, b::SpeedyPeriod, c::SpeedyPeriod) = a + b + c
Base.:(+)(a::SpeedyDateTime, b::SpeedyPeriod, c::SpeedyPeriod, d::SpeedyPeriod...) = (+)((a + b + c), d...)

# ============================================================================
# CONVERSION FROM/TO Julia Base Dates
# ============================================================================
# Convert from Dates.DateTime to SpeedyDateTime
SpeedyDateTime(dt::Dates.DateTime) = SpeedyDateTime(_SpeedyUTM(Dates.value(dt)))

# Convert from SpeedyDateTime to Dates.DateTime
Dates.DateTime(dt::SpeedyDateTime) = Dates.DateTime(Dates.UTM(Int64(Dates.value(dt))))

# Convert from Dates.Period to SpeedyPeriod
SpeedyMillisecond(p::Dates.Millisecond) = SpeedyMillisecond(Dates.value(p))
SpeedySecond(p::Dates.Second) = SpeedySecond(Dates.value(p))
SpeedyMinute(p::Dates.Minute) = SpeedyMinute(Dates.value(p))
SpeedyHour(p::Dates.Hour) = SpeedyHour(Dates.value(p))
SpeedyDay(p::Dates.Day) = SpeedyDay(Dates.value(p))
SpeedyWeek(p::Dates.Week) = SpeedyWeek(Dates.value(p))
SpeedyMonth(p::Dates.Month) = SpeedyMonth(Dates.value(p))
SpeedyYear(p::Dates.Year) = SpeedyYear(Dates.value(p))

# Convert from SpeedyPeriod to Dates.Period
Dates.Millisecond(p::SpeedyMillisecond) = Dates.Millisecond(Int64(Dates.value(p)))
Dates.Second(p::SpeedySecond) = Dates.Second(Int64(Dates.value(p)))
Dates.Minute(p::SpeedyMinute) = Dates.Minute(Int64(Dates.value(p)))
Dates.Hour(p::SpeedyHour) = Dates.Hour(Int64(Dates.value(p)))
Dates.Day(p::SpeedyDay) = Dates.Day(Int64(Dates.value(p)))
Dates.Week(p::SpeedyWeek) = Dates.Week(Int64(Dates.value(p)))
Dates.Month(p::SpeedyMonth) = Dates.Month(Int64(Dates.value(p)))
Dates.Year(p::SpeedyYear) = Dates.Year(Int64(Dates.value(p)))

# ============================================================================
# BROADCASTING SUPPORT
# ============================================================================
Base.Broadcast.broadcastable(x::Union{SpeedyPeriod, SpeedyUTInstant, SpeedyDateTime}) = Ref(x)
