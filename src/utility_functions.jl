"""Check whether elements of a vector are increasing."""
function isincreasing(x::Vector)
    is_increasing = true
    for i in 2:length(x)
        is_increasing &= x[i] > x[i-1] ? true : false
    end
    return is_increasing
end

"""Set all negative entries in an array to zero."""
function clip_negatives!(A::AbstractArray{T}) where T
    @inbounds for i in eachindex(A)
        A[i] = max(A[i],zero(T))
    end
end

"""Returns an integer `m >= n` of argument `n` with only prime factors 2 and 3
to obtain an easily fourier-transformable number of longitudes.
    m = 2^i * 3^j >= n, with i>=0 but j = 0,1
"""
function roundup_fft(n::Int)
    lz = leading_zeros(n)                   # determine scale of n

    # create a mask that's 1 for all but the two most significant figures of n
    # for finding the next largest integer of n with factors 2 and 3
    mask = (1 << (8*sizeof(n)-lz-2))-1    
    
    # round up by adding mask as offset and then masking all insignificant bits
    n_roundedup = (n + mask) & ~mask
    return n_roundedup
end