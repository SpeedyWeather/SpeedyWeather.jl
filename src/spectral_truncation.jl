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
    spectral_truncation!(alms::LowerTriangularMatrix,ltrunc::Integer,mtrunc::Integer)

Truncate spectral coefficients `alms` in-place by setting all coefficients for which the degree `l`
is larger than the truncation `ltrunc` or order `m` larger than the truncaction `mtrunc`.
Similar to `spectral_truncation!(::AbstractMatrix,...) but skips the upper triangle which is
zero by design for LowerTriangularMatrix."""
function spectral_truncation!(  alms::LowerTriangularMatrix{NF},    # spectral field to be truncated
                                ltrunc::Integer,                    # truncate to max degree ltrunc
                                mtrunc::Integer,                    # truncate to max order mtrunc
                                ) where NF                          # number format NF (can be complex)
    
    lmax,mmax = size(alms) .- 1         # 0-based degree l, order m of the legendre polynomials

    lm = 1
    @inbounds for m in 1:mmax+1         # order m = 0,mmax but 1-based
        for l in m:lmax+1               # degree l = 0,lmax but 1-based
            if  l > ltrunc+1 ||         # and degrees l>ltrunc
                m > mtrunc+1            # and orders m>mtrunc

                alms[lm] = zero(NF)     # set that coefficient to zero
                
            end
            lm += 1
        end
    end
    return alms
end

"""
    spectral_truncation!(alms,trunc)

Truncate spectral coefficients `alms` in-place by setting (a) the upper right triangle to zero and (b)
all coefficients for which the degree l is larger than the truncation `trunc`."""
function spectral_truncation!(  alms::AbstractMatrix,   # spectral field to be truncated
                                trunc::Int)             # truncate to degree/order trunc
    return spectral_truncation!(alms,trunc,trunc)       # use trunc=ltrunc=mtrunc
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
    alms_trunc = LowerTriangularMatrix{NF}(undef,size_alms_new...)  # preallocate new (smaller) array

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
    alms_trunc = zeros(LowerTriangularMatrix{NF},size_alms_new...)  # allocate new (larger) array

    # copy data over
    copyto!(@view(alms_trunc[CartesianIndices(size_alms)]),alms)
    spectral_truncation!(alms_trunc,ltrunc,mtrunc)  # make sure upper triangle is zero
    return alms_trunc
end

spectral_interpolation(alms::AbstractArray,trunc::Int) = spectral_interpolation(alms,trunc,trunc)