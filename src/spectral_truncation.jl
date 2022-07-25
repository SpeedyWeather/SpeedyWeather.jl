"""
    m = roundup_fft(n::Int;
                    small_primes::Vector{Int}=[2,3,5])

Returns an integer `m >= n` with only small prime factors 2, 3, 5 (default, others can be specified
with the keyword argument `small_primes`) to obtain an efficiently fourier-transformable number of
longitudes, m = 2^i * 3^j * 5^k >= n, with i,j,k >=0.
"""
function roundup_fft(n::Int;small_primes::Vector{Int}=[2,3,5])
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

struct TriangularTruncation
    trunc::Int      # spectral truncation (trunc=lmax=mmax)
    nlon::Int       # number of longitudes
    nlat::Int       # number of latitudes

    TriangularTruncation(trunc,nlon,nlat) = is_triangular_truncation(trunc,nlon,nlat) ?
        new(trunc,nlon,nlat) : error("Not a triangular truncation.")
end

"""
    is_triangular_truncation(trunc::Int,nlon::Int,nlat::Int)

Tests whether the inputs `trunc, nlon, nlat` satisfy the triangular truncation constraints.
`trunc` is the maximum degree and order of the Legendre polynomials (0-based), `nlon` is the
number of longitudes, `nlat` the number of latitudes on the spatial grid. The constraints are

    - nlon >= 3T+1
    - nlat >= (3T+1)/2
    - nlon = 2nlat."""
function is_triangular_truncation(trunc::Int,nlon::Int,nlat::Int)
    info = "$(nlon)x$(nlat) grid and T$trunc violate triangular truncation constraints, "
    constraints = "nlon >= 3T+1, nlat >= (3T+1)/2, nlon = 2nlat, nlat even."
    feedback = info*constraints
    @assert nlon >= 3*trunc+1 feedback
    @assert nlat >= (3*trunc+1)/2 feedback
    @assert nlon == 2nlat feedback
    @assert iseven(nlat) feedback
    return true     # return true if all @asserts pass the constraints
end

# unpack 
is_triangular_truncation(T::TriangularTruncation) = is_triangular_truncation(T.trunc,T.nlon,T.nlat)

"""
    tri_trunc = triangular_truncation(;trunc::Int=0,nlon::Int=0,nlat::Int=0)

Returns a `tri_trunc::TriangularTruncation` struct, that contains the spectral truncation `trunc`,
and the grid sizes `nlon,nlat` that satisfy the triangular truncation constraints.
Provide either `trunc`, `nlon`, or `nlat` but not a combination. If `trunc` is provided,
`nlon` will be easily Fast Fourier-transformable, as determined in `roundup_fft`.
If `nlon` or `nlat` are provided, `tri_trunc.trunc` is the largest spectral truncation
that satisfies the constraints."""
function triangular_truncation(;trunc::Int=0,nlon::Int=0,nlat::Int=0)
    # if all trunc, nlon, nlat are provided, test if triangular truncation
    if trunc > 0 && nlon > 0 && nlat > 0    
        is_triangular_truncation(trunc,nlon,nlat)

    # if only trunc is provided, get nlon,nlat
    elseif trunc > 0 && nlon == 0 && nlat == 0
        nlon = roundup_fft(3*trunc+1)
        nlat = nlon÷2

    # if only nlon is provided, get trunc,nlat
    elseif trunc == 0 && nlon > 0 && nlat == 0
        nlat = nlon ÷ 2
        trunc = floor(Int,(nlon-1)/3)

    # if only nlat is provided, get trunc,nlon
    elseif trunc == 0 && nlon == 0 && nlat > 0
        nlon = 2nlat
        trunc = floor(Int,(nlon-1)/3)

    else
        throw(error("Please use either trunc, nlon or nlat keywords, not a combination of them."))
    end

    return TriangularTruncation(trunc,nlon,nlat)
end

"""
    spectral_truncation!(alms,ltrunc,mtrunc)

Truncate spectral coefficients `alms` in-place by setting (a) the upper right triangle to zero and (b)
all coefficients for which the degree `l` is larger than the truncation `ltrunc` or order `m` larger than
the truncaction `mtrunc`."""
function spectral_truncation!(  alms::AbstractMatrix{NF},   # spectral field to be truncated
                                ltrunc::Int,                # truncate to max degree ltrunc
                                mtrunc::Int                 # truncate to max order mtrunc
                                ) where NF                  # number format NF (can be complex)
    
    lmax,mmax = size(alms) .- 1    # 0-based degree l, order m of the legendre polynomials

    @inbounds for m in 1:mmax+1         # order m = 0,mmax but 1-based
        for l in 1:lmax+1               # degree l = 0,lmax but 1-based
            if m > l ||                 # upper triangle (m>l)
                l > ltrunc+1 ||         # and degrees l>ltrunc
                m > mtrunc+1            # and orders m>mtrunc

                alms[l,m] = zero(NF)    # set that coefficient to zero
            end
        end
    end
    return alms
end

"""
    spectral_truncation!(alms,trunc)

Truncate spectral coefficients `alms` in-place by setting (a) the upper right triangle to zero and (b)
all coefficients for which the degree l is larger than the truncation `trunc`."""
function spectral_truncation!(  alms::AbstractMatrix{NF},       # spectral field to be truncated
                                trunc::Int                      # truncate to degree/order trunc
                                ) where NF                      # number format NF (can be complex)
    return spectral_truncation!(alms,trunc,trunc)               # use trunc=ltrunc=mtrunc
end

"""
    spectral_truncation!(alms)

Truncate spectral coefficients `alms` in-place by setting the upper right triangle to zero. This is
to enforce that all coefficients for which the degree l is larger than order m are zero."""
spectral_truncation!(alms::AbstractMatrix) = spectral_truncation!(alms,size(alms)...)

"""
    alms_trunc = spectral_truncation(alms,trunc)

Returns a spectral coefficient matrix `alms_trunc` that is truncated from `alms` to the size (`trunc`+1)^2.
`alms_trunc` only contains those coefficient of `alms` for which m,l <= trunc, and l>=m are zero anyway.
If `trunc` is larger than the implicit truncation in `alms` obtained from its size than `spectral_interpolation`
is automatically called instead, returning `alms_interp`, a coefficient matrix that is larger than `alms`
with padded zero coefficients. Also works with higher dimensional arrays, but truncation is only applied to
the first two dimensions."""
function spectral_truncation(   alms::AbstractArray{NF},    # spectral field to be truncated
                                ltrunc::Int,                # truncate to max degree ltrunc
                                mtrunc::Int                 # truncate to max order mtrunc
                                ) where NF                  # number format NF (can be complex)
    
    @boundscheck length(size(alms)) >= 2 || throw(BoundsError(alms,(ltrunc,mtrunc)))

    size_alms = size(alms)
    lmax,mmax = size_alms[1:2] .- 1     # 0-based degree l, order m of the spherical harmonics
    
    # interpolate to higher resolution if output larger than input
    (ltrunc > lmax || mtrunc > mmax) && return spectral_interpolation(alms,ltrunc,mtrunc)

    size_alms_new = collect(size_alms)              # convert to vector for mutability
    size_alms_new[1] = ltrunc+1                     # new (smaller) size
    size_alms_new[2] = mtrunc+1
    size_alms_new = tuple(size_alms_new...)
    alms_trunc = Array{NF}(undef,size_alms_new...)  # preallocate new (smaller) array

    # copy data over
    copyto!(alms_trunc,@view(alms[CartesianIndices(size_alms_new)]))
    spectral_truncation!(alms_trunc,ltrunc,mtrunc)  # make sure undef is replaced with zero
    return alms_trunc
end

spectral_truncation(alms::AbstractArray,trunc::Int) = spectral_truncation(alms,trunc,trunc)

"""
    alms_interp = spectral_interpolation(   alms::AbstractArray{NF},    # spectral field to be truncated
                                            ltrunc::Int,                # truncate to max degree ltrunc
                                            mtrunc::Int                 # truncate to max order mtrunc
                                            ) where NF                  # number format NF (can be complex)

Returns a spectral coefficient matrix `alms_interp` that is `alms` padded with zeros to interpolate in
spectral space. If `trunc` is smaller or equal to the implicit truncation in `alms` obtained from its size
than `spectral_truncation` is automatically called instead, returning `alms_trunc`, a coefficient matrix that
is smaller than `alms`, implicitly setting higher degrees and orders to zero. Also works with higher
dimensional arrays, but interpolation is only applied to the first two dimensions."""
function spectral_interpolation(alms::AbstractArray{NF},    # spectral field to be truncated
                                ltrunc::Int,                # truncate to max degree ltrunc
                                mtrunc::Int                 # truncate to max order mtrunc
                                ) where NF                  # number format NF (can be complex)
    
    @boundscheck length(size(alms)) >= 2 || throw(BoundsError(alms,(ltrunc,mtrunc)))

    size_alms = size(alms)
    lmax,mmax = size_alms[1:2] .- 1     # 0-based degree l, order m of the spherical harmonics
    
    # interpolate to higher resolution if output larger than input
    (ltrunc <= lmax && mtrunc <= mmax) && return spectral_interpolation(alms,ltrunc,mtrunc)

    size_alms_new = collect(size_alms)              # convert to vector for mutability
    size_alms_new[1] = ltrunc+1                     # new size
    size_alms_new[2] = mtrunc+1
    size_alms_new = tuple(size_alms_new...)         # convert back to tuple
    alms_trunc = zeros(NF,size_alms_new...)         # allocate new (larger) array

    # copy data over
    copyto!(@view(alms_trunc[CartesianIndices(size_alms)]),alms)
    spectral_truncation!(alms_trunc,ltrunc,mtrunc)  # make sure upper triangle is zero
    return alms_trunc
end

spectral_interpolation(alms::AbstractArray,trunc::Int) = spectral_interpolation(alms,trunc,trunc)