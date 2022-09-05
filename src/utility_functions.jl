"""
    true/false = isincreasing(v::Vector)

Check whether elements of a vector `v` are increasing."""
function isincreasing(x::Vector)
    is_increasing = true
    for i in 2:length(x)
        is_increasing &= x[i] > x[i-1] ? true : false
    end
    return is_increasing
end

"""
    true/false = is_power_2(i::Integer)

Checks whether an integer `i` is a power of 2, i.e. i = 2^k, with k = 0,1,2,3,...."""
is_power_2(i::Integer) = i != 0 ? i & (i-1) == 0 : false
is_power_2_or_0(i::Integer) = i & (i-1) == 0

"""
    m = roundup_fft(n::Int;
                    small_primes::Vector{Int}=[2,3,5])

Returns an integer `m >= n` with only small prime factors 2, 3, 5 (default, others can be specified
with the keyword argument `small_primes`) to obtain an efficiently fourier-transformable number of
longitudes, m = 2^i * 3^j * 5^k >= n, with i,j,k >=0.
"""
function roundup_fft(n::Integer;small_primes::Vector{T}=[2,3,5]) where {T<:Integer}
    factors_not_in_small_primes = true      # starting condition for while loop
    n += isodd(n) ? 1 : 0                   # start with an even n
    while factors_not_in_small_primes
        
        factors = Primes.factor(n)          # prime factorization
        all_factors_small = true            # starting condition
        
        for i in 1:length(factors)          # loop over factors and check they are small
            factor = factors.pe[i].first    # extract factor from factors
            all_factors_small &= factor in small_primes
        end
        
        factors_not_in_small_primes = ~all_factors_small    # all factors small will abort while loop
        n += 2                                              # test for next larger even n
    
    end
    return n-2      # subtract unnecessary last += 2 addition
end

"""
    clip_negatives!(A::AbstractArray)

Set all negative entries `a` in `A` to zero."""
function clip_negatives!(A::AbstractArray{T}) where T
    @inbounds for i in eachindex(A)
        A[i] = max(A[i],zero(T))
    end
end

"""
    underflow!(A::AbstractArray,ϵ::Real)

Underflows element `a` in `A` to zero if `abs(a) < ϵ`."""
function underflow!(A::AbstractArray{T},ϵ::Real) where T
    ϵT = convert(T,abs(ϵ))
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
    readable_secs(secs::Real) -> Dates.CompoundPeriod

Returns `Dates.CompoundPeriod` rounding to either (days, hours), (hours, minutes), (minutes,
seconds), or seconds with 1 decimal place accuracy for >10s and two for less.
E.g.
```julia
julia> readable_secs(12345)
3 hours, 26 minutes
```
"""
function readable_secs(secs::Real)
    millisecs = Dates.Millisecond(round(secs * 10 ^ 3))
    if millisecs >= Dates.Day(1)
        return Dates.canonicalize(round(millisecs, Dates.Hour))
    elseif millisecs >= Dates.Hour(1)
        return Dates.canonicalize(round(millisecs, Dates.Minute))
    elseif millisecs >= Dates.Minute(1)
        return Dates.canonicalize(round(millisecs, Dates.Second))
    elseif millisecs >= Dates.Second(10)
        return Dates.canonicalize(round(millisecs, Dates.Millisecond(100)))
    end
    return Dates.canonicalize(round(millisecs, Dates.Millisecond(10)))
end

# define as will only become availale in Julia 1.9
pkgversion(m::Module) = VersionNumber(TOML.parsefile(joinpath(
    dirname(string(first(methods(m.eval)).file)), "..", "Project.toml"))["version"])