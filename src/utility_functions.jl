"""
$(TYPEDSIGNATURES)
Check whether elements of a vector `v` are strictly increasing."""
function isincreasing(v::AbstractVector)
    is_increasing = true
    for i in 2:length(v)
        is_increasing &= v[i-1] < v[i] ? true : false
    end
    return is_increasing
end

"""
$(TYPEDSIGNATURES)

Check whether elements of a vector `v` are strictly decreasing."""
function isdecreasing(v::AbstractVector)
    is_decreasing = true
    for i in 2:length(v)
        is_decreasing &= v[i-1] > v[i] ? true : false
    end
    return is_decreasing
end

"""
    clip_negatives!(A::AbstractArray)

Set all negative entries `a` in `A` to zero."""
function clip_negatives!(A::AbstractArray{T}) where T
    @inbounds for i in eachindex(A)
        A[i] = max(A[i], zero(T))
    end
end

"""
    underflow!(A::AbstractArray, ϵ::Real)

Underflows element `a` in `A` to zero if `abs(a) < ϵ`."""
function underflow!(A::AbstractArray{T}, ϵ::Real) where T
    ϵT = convert(T, abs(ϵ))
    @inbounds for i in eachindex(A)
        A[i] = abs(A[i]) < ϵT ? zero(T) : A[i]
    end
end

"""
    flipgsign!(A::AbstractArray)

Like `-A` but in-place."""
function flipsign!(A::AbstractArray)
    @inbounds for i in eachindex(A)
        A[i] = -A[i]
    end
    A
end

"""
    A = nans(T, dims...)

Allocate array A with NaNs of type T. Similar to `zeros(T, dims...)`."""
function nans(::Type{T}, dims...) where T
    return fill(convert(T, NaN), dims...)
end

"""
    A = nans(dims...)

Allocate `A::Array{Float64}` with NaNs."""
nans(dims...) = nans(Float64, dims...)

"""
$(TYPEDSIGNATURES)
Prints to `io` all fields of a struct `A` identified by their
`keys`."""
function print_fields(io::IO, A, keys;arrays::Bool=false)
    keys_filtered = arrays ? keys : filter(key -> ~(getfield(A, key) isa AbstractArray), keys)
    n = length(keys_filtered)
    filtered = n < length(keys)
    for (i, key) in enumerate(keys_filtered)
        last = (i == n) & ~filtered
        key = keys_filtered[i]
        val = getfield(A, key)
        ~last ? println(io, "├ $key::$(typeof(val)) = $val") :
                print(io,  "└ $key::$(typeof(val)) = $val")
    end
    if filtered                 # add the names of arrays
        s = "└── arrays: "
        for key in keys
            if ~(key in keys_filtered)
                s *= "$key, "
            end
        end

        # remove last ", "
        s_without_comma = s[1:prevind(s, findlast(==(','), s))]
        print(io, s_without_comma)
    end
end

"""
    Century <: Period

Convenience time type representing a 100-year period.
"""
struct Century <: Period
    value::Int64
end

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
Dates.Second(x::AbstractFloat) = convert(Second, x)
Dates.Minute(x::AbstractFloat) = Second(60x)
Dates.Hour(  x::AbstractFloat) = Minute(60x)
Dates.Day(   x::AbstractFloat) = Hour(24x)
Dates.Week(  x::AbstractFloat) = Day(7x)
Dates.Month( x::AbstractFloat) = Day(30x)  # approximate
Dates.Year(  x::AbstractFloat) = Day(365x) # approximate
Century(     x::AbstractFloat) = Year(100x)
Millenium(   x::AbstractFloat) = Century(10x)

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
function Base.convert(::Type{Dates.Second}, m::Dates.Month)
    return Second(m.value * 30 * 24 * 60 * 60) # approximate
end

# year conversions
Base.convert(::Type{Dates.Day}, y::Dates.Year) = Day(Hour(y))
Base.convert(::Type{Dates.Hour}, y::Dates.Year) = Hour(Second(y))
function Base.convert(::Type{Dates.Second}, y::Dates.Year)
    return Second(y.value * 365 * 24 * 60 * 60) # approximate
end

# additional century conversions
Base.convert(::Type{Second}, c::Century) = Second(Year(c))
Base.convert(::Type{Hour}, c::Century) = Hour(Year(c))
Base.convert(::Type{Day}, c::Century) = Day(Year(c))
Base.convert(::Type{Month}, c::Century) = Month(Year(c))

# additional month conversions
Base.convert(::Type{Second}, m::Millenium) = Second(Year(m))
Base.convert(::Type{Hour}, m::Millenium) = Hour(Year(m))
Base.convert(::Type{Day}, m::Millenium) = Day(Year(m))
Base.convert(::Type{Month}, m::Millenium) = Month(Year(m))
