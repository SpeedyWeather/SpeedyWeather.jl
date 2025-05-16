"""
$(TYPEDSIGNATURES)
Triangular truncation to degree `ltrunc` and order `mtrunc` (both 0-based). Truncate spectral coefficients `alms` in-place
by setting all coefficients for which the degree `l` is larger than the truncation `ltrunc` or order `m` larger
than the truncaction `mtrunc`."""
function spectral_truncation!(
    alms::LowerTriangularArray,     # spectral field to be truncated
    ltrunc::Integer,                # truncate to max degree ltrunc (0-based)
    mtrunc::Integer,                # truncate to max order mtrunc (0-based)
)   
    lmax, mmax = size(alms, OneBased; as=Matrix)   # 1-based degree l, order m of the legendre polynomials

    ltrunc += 1     # 0-based to 1-based
    mtrunc += 1

    for k in eachmatrix(alms)
        lm = 1
        for m in 1:mmax                 # order m = 0, mmax but 1-based
            for l in m:lmax             # degree l = 0, lmax but 1-based
                if  l > ltrunc ||       # and degrees l>ltrunc
                    m > mtrunc          # and orders m>mtrunc

                    alms[lm, k] = 0     # set that coefficient to zero
                end
                lm += 1
            end
        end
    end
    return alms
end

"""
$(TYPEDSIGNATURES)
Sets the upper triangle of `A` to zero."""
function spectral_truncation!(A::AbstractMatrix)
    lmax, mmax = size(A)

    for m in 1:mmax
        for l in 1:lmax
            if  m > l
                A[l, m] = 0
            end
        end
    end
    return A
end

"""
$(TYPEDSIGNATURES)
Triangular truncation of `alms` to degree and order `trunc` in-place."""
spectral_truncation!(alms::LowerTriangularArray, trunc::Integer) = spectral_truncation!(alms, trunc, trunc)

"""
$(TYPEDSIGNATURES)
Triangular truncation of `alms` to the size of it, sets additional rows to zero."""
spectral_truncation!(alms::LowerTriangularArray) = spectral_truncation!(alms, size(alms, 2, ZeroBased, as=Matrix))


"""
$(TYPEDSIGNATURES)
Returns a LowerTriangularArray that is truncated from `alms` to the size (`ltrunc`+1) x (`mtrunc`+1),
both inputs are 0-based. If `ltrunc` or `mtrunc` is larger than the corresponding size of`alms` than
`spectral_interpolation` is automatically called instead, returning a LowerTriangularArray padded zero
coefficients for higher wavenumbers."""
function spectral_truncation(
    ::Type{NF},                 # number format NF (can be complex)
    alms::LowerTriangularArray{T, N, ArrayType, S}, # spectral field to be truncated
    ltrunc::Integer,            # truncate to max degree ltrunc (0-based)
    mtrunc::Integer,            # truncate to max order mtrunc (0-based)
) where {NF, T, N, S, ArrayType}
    
    lmax, mmax, k... = size(alms, ZeroBased, as=Matrix)
    
    # interpolate to higher resolution if output larger than input
    (ltrunc > lmax || mtrunc > mmax) && return spectral_interpolation(NF, alms, ltrunc, mtrunc)

    # preallocate new (smaller) array
    ArrayType_ = LowerTriangularArrays.nonparametric_type(ArrayType)
    alms_trunc = zeros(LowerTriangularArray{NF, N, ArrayType_{NF, N}, S}, Spectrum(ltrunc+1, mtrunc+1), k...)  

    # copy data over, copyto! copies the largest matching subset of harmonics
    copyto!(alms_trunc, alms)
    return alms_trunc
end

spectral_truncation(alms::LowerTriangularArray, ltrunc::Integer, mtrunc::Integer) = spectral_truncation(eltype(alms), alms, ltrunc, mtrunc)
spectral_truncation(alms::LowerTriangularArray, trunc::Integer) = spectral_truncation(alms, trunc, trunc)

"""
$(TYPEDSIGNATURES)
Returns a LowerTriangularArray that is interpolated from `alms` to the size (`ltrunc`+1) x (`mtrunc`+1),
both inputs are 0-based, by padding zeros for higher wavenumbers. If `ltrunc` or `mtrunc` are smaller than the
corresponding size of`alms` than `spectral_truncation` is automatically called instead, returning a smaller
LowerTriangularArray."""
function spectral_interpolation(
    ::Type{NF},                 # number format NF (can be complex)
    alms::LowerTriangularArray{T, N, ArrayType, S}, # spectral field to be truncated
    ltrunc::Integer,            # truncate to max degree ltrunc (0-based)
    mtrunc::Integer,            # truncate to max order mtrunc (0-based)
) where {NF, T, N, S, ArrayType}                
    
    lmax, mmax, k... = size(alms, ZeroBased, as=Matrix)
    
    # truncate to lower resolution if output smaller than input
    (ltrunc <= lmax && mtrunc <= mmax) && return spectral_truncation(NF, alms, ltrunc, mtrunc)

    # preallocate new (larger) array
    ArrayType_ = LowerTriangularArrays.nonparametric_type(ArrayType)
    alms_interp = zeros(LowerTriangularArray{NF, N, ArrayType_{NF, N}, S}, Spectrum(ltrunc+1, mtrunc+1), k...)  

    # copy data over, copyto! copies the largest matching subset of harmonics
    copyto!(alms_interp, alms)
    return alms_interp
end

spectral_interpolation(alms::LowerTriangularArray, trunc::Integer) = spectral_interpolation(alms, trunc, trunc)

"""
$(TYPEDSIGNATURES)
Set imaginary component of m=0 modes (the zonal modes in the first column) to 0."""
function zero_imaginary_zonal_modes!(
    alms::LowerTriangularArray,
)
    lmax, mmax = size(alms, OneBased, as=Matrix)

    for k in eachmatrix(alms)
        for l in 1:lmax
            alms[l, k] = real(alms[l, k])
        end
    end
    return alms
end

"""$(TYPEDSIGNATURES)
Smooth the spectral field `A` following A_smooth = (1-c*∇²ⁿ)A with power n of a normalised Laplacian
so that the highest degree lmax is dampened by multiplication with c. Anti-diffusion for c<0."""
function spectral_smoothing(A::LowerTriangularArray, c::Real; power::Real=1)
    A_smooth = copy(A)
    spectral_smoothing!(A_smooth, c; power)
    return A_smooth
end

"""$(TYPEDSIGNATURES)
Smooth the spectral field `A` following A *= (1-(1-c)*∇²ⁿ) with power n of a normalised Laplacian
so that the highest degree lmax is dampened by multiplication with c. Anti-diffusion for c>1."""
function spectral_smoothing!(   L::LowerTriangularArray,
                                c::Real;
                                power::Real=1,          # power of Laplacian used for smoothing
                                truncation::Int=-1)     # smoothing wrt wavenumber (0 = largest)
                        
    lmax, mmax = size(L; as=Matrix)
        
    # normalize by largest eigenvalue by default, or wrt to given truncation
    eigenvalue_norm = truncation == -1 ? -mmax*(mmax+1) : -truncation*(truncation+1)

    for k in eachmatrix(L)
        lm = 1
        for m in 1:mmax
            for l in m:lmax
                eigenvalue_normalised = -l*(l-1)/eigenvalue_norm
                # for eigenvalue_norm < largest eigenvalue the factor becomes negative
                # set to zero in that case
                L[lm, k] *= max(1 - (1-c)*eigenvalue_normalised^power, 0)
                lm += 1
            end
        end
    end
end