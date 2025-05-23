const DEFAULT_DEALIASING = 2.0

"""
$(TYPEDSIGNATURES)
For the spectral truncation `trunc` (e.g. 31 for T31) return the
grid resolution parameter `nlat_half` (number of latitude rings on one hemisphere
including the Equator) following a dealiasing parameter (default 2)
to match spectral and grid resolution."""
function RingGrids.get_nlat_half(   trunc::Integer,
                                    dealiasing::Real=DEFAULT_DEALIASING)
                    
    return roundup_fft(ceil(Int, ((1+dealiasing)*trunc+1)/4))
end

# inverse of get_nlat_half to reobtain dealiasing
get_dealiasing(trunc, nlat_half) = (4nlat_half)/(trunc+1) - 1

"""
$(TYPEDSIGNATURES)
For the grid resolution parameter `nlat_half` (e.g. 24 for a 48-ring FullGaussianGrid)
return the spectral truncation `trunc` (max degree of spherical harmonics)
following a dealiasing parameter (default 2) to match spectral and grid resolution."""
function get_truncation(nlat_half::Integer,
                        dealiasing::Real=DEFAULT_DEALIASING)

    return floor(Int, (4nlat_half-1)/(dealiasing+1))
end

# unpack nlat_half from provided map
get_truncation(field::AbstractField, dealiasing::Real=DEFAULT_DEALIASING) =
    get_truncation(field.nlat_half, dealiasing)

"""
    m = roundup_fft(n::Int;
                    small_primes::Vector{Int}=[2, 3, 5])

Returns an integer `m >= n` with only small prime factors 2, 3 (default, others can be specified
with the keyword argument `small_primes`) to obtain an efficiently fourier-transformable number of
longitudes, m = 2^i * 3^j * 5^k >= n, with i, j, k >=0.
"""
function roundup_fft(n::Integer; small_primes::Vector{T}=[2, 3, 5]) where {T<:Integer}
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
    true/false = is_power_2(i::Integer)

Checks whether an integer `i` is a power of 2, i.e. i = 2^k, with k = 0, 1, 2, 3, ...."""
is_power_2(i::Integer) = i != 0 ? i & (i-1) == 0 : false
is_power_2_or_0(i::Integer) = i & (i-1) == 0