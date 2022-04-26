"""
    m = roundup_fft(n::Int)

Returns an integer `m >= n` with only prime factors 2 and 3 to obtain an efficiently
fourier-transformable number of longitudes, m = 2^i * 3^j >= n, with i>=0 but j = 0,1.
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

"""
    triangular_truncation(trunc::Int,nlon::Int,nlat::Int)

Tests whether the inputs `trunc, nlon, nlat` satisfy the triangular truncation constraints.
`trunc` is the maximum degree and order of the Legendre polynomials (0-based), `nlon` is the
number of longitudes, `nlat` the number of latitudes on the spatial grid. The constraints are

    - nlon >= 3T+1
    - nlat >= (3T+1)/2
    - nlon = 2nlat."""
function triangular_truncation(trunc::Int,nlon::Int,nlat::Int)
    info = "$(nlon)x$(nlat) grid and T$trunc violate triangular truncation constraints, "
    constraints = "nlon >= 3T+1, nlat >= (3T+1)/2, nlon = 2nlat."
    feedback = info*constraints
    @assert nlon >= 3*trunc+1 feedback
    @assert nlat >= (3*trunc+1)/2 feedback
    @assert nlon == 2nlat feedback
end

"""
    nlon, nlat = triangular_truncation(trunc::Int)

Returns the grid size `nlon,nlat` for a spectral truncation `trunc` that satisfies the
triangular truncation constraints. Returns the smallest pair `nlon,nlat` that is also
easily Fast Fourier-transformable, as determined in `roundup_fft`."""
function triangular_truncation(trunc::Int)
    nlon = roundup_fft(3*trunc+1)
    nlat = nlon÷2
    triangular_truncation(trunc,nlon,nlat)
    return nlon,nlat
end

"""
    trunc = triangular_truncation(nlon::Int,nlat::Int)

Returns the largest spectral truncation `trunc` that satisfies the triangular truncation constraints
based on the grid size `nlon,nlat`, which may or may not be easily Fast Fourier-transformable."""
function triangular_truncation(nlon::Int,nlat::Int)
    trunc = floor(Int,(nlon-1)/3)
    triangular_truncation(trunc,nlon,nlat)
    return trunc
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
`alms_trunc` only contains those coefficient of `alms` for which m,l <= trunc, and m>l are zero anyway. If
`trunc` is larger than the implicit truncation in `alms` obtained from its size than `spectral_interpolation`
is automatically called instead, returning `alms_interp`, a coefficient matrix that is larger than `alms`
with padded zero coefficients."""
function spectral_truncation(   alms::AbstractMatrix{NF},   # spectral field to be truncated
                                ltrunc::Int,                # truncate to max degree ltrunc
                                mtrunc::Int                 # truncate to max order mtrunc
                                ) where NF                  # number format NF (can be complex)
    
    lmax,mmax = size(alms) .- 1    # 0-based degree l, order m of the spherical harmonics
    (ltrunc > lmax || mtrunc > mmax) && return spectral_interpolation(alms,ltrunc,mtrunc)
    alms_trunc = Matrix{NF}(undef,ltrunc+1,mtrunc+1)
    copyto!(alms_trunc,@view(alms[1:ltrunc+1,1:mtrunc+1]))
    spectral_truncation!(alms_trunc,ltrunc,mtrunc)
    return alms_trunc
end

spectral_truncation(alms::AbstractMatrix,trunc::Int) = spectral_truncation(alms,trunc,trunc)

"""
    alms_interp = spectral_interp(alms,trunc)

Returns a spectral coefficient matrix `alms_interp` that is `alms` padded with zeros to interpolate in
spectral space. If `trunc` is smaller or equal to the implicit truncation in `alms` obtained from its size
than `spectral_truncation` is automatically called instead, returning `alms_trunc`, a coefficient matrix that
is smaller than `alms`, implicitly setting higher degrees and orders to zero."""
function spectral_interpolation(    alms::AbstractMatrix{NF},   # spectral field to be truncated
                                    ltrunc::Int,                # truncate to max degree ltrunc
                                    mtrunc::Int                 # truncate to max order mtrunc
                                    ) where NF                  # number format NF (can be complex)
    
    lmax,mmax = size(alms) .- 1    # 0-based degree l, order m of the spherical harmonics
    (ltrunc <= lmax && mtrunc <= mmax) && return spectral_truncation(alms,ltrunc,mtrunc)
    alms_interp = zeros(NF,ltrunc+1,mtrunc+1)               # allocate larger array
    copyto!(@view(alms_interp[1:lmax+1,1:mmax+1]),alms)     # copy over coeffs from alms
    spectral_truncation!(alms_interp,ltrunc,mtrunc)         # set other coeffs to zero
    return alms_interp
end

spectral_interpolation(alms::AbstractMatrix,trunc::Int) = spectral_interpolation(alms,trunc,trunc)