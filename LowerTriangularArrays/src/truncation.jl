"""
$(TYPEDSIGNATURES)
Triangular truncation to degree `ltrunc` and order `mtrunc` (both 0-based). Truncate spectral coefficients `alms` in-place
by setting all coefficients for which the degree `l` is larger than the truncation `ltrunc` or order `m` larger
than the truncaction `mtrunc`."""
function truncate!(
        alms::LowerTriangularArray,     # spectral field to be truncated
        ltrunc::Integer,                # truncate to max degree ltrunc (0-based)
        mtrunc::Integer,                # truncate to max order mtrunc (0-based)
    )
    (; l_indices, m_indices) = alms.spectrum

    # Convert to 1-based indexing
    ltrunc += 1     # 0-based to 1-based
    mtrunc += 1

    # Launch kernel for GPU/Reactant compatibility
    arch = architecture(alms)
    launch!(arch, SpectralWorkOrder, size(alms), _truncate_kernel!,
        alms.data, l_indices, m_indices, ltrunc, mtrunc)

    return alms
end

@kernel inbounds = true function _truncate_kernel!(data, @Const(l_indices), @Const(m_indices), ltrunc, mtrunc)
    I = @index(Global, NTuple)
    lm = I[1]

    l = l_indices[lm]
    m = m_indices[lm]

    if l > ltrunc || m > mtrunc
        data[I...] = 0
    end
end

"""
$(TYPEDSIGNATURES)
Sets the upper triangle of `A` to zero."""
function truncate!(A::AbstractMatrix)
    lmax, mmax = size(A)

    for m in 1:mmax
        for l in 1:lmax
            if m > l
                A[l, m] = 0
            end
        end
    end
    return A
end

"""
$(TYPEDSIGNATURES)
Triangular truncation of `alms` to degree and order `trunc` in-place."""
truncate!(alms::LowerTriangularArray, trunc::Integer) = truncate!(alms, trunc, trunc)

"""
$(TYPEDSIGNATURES)
Triangular truncation of `alms` to the size of it, sets additional rows to zero."""
truncate!(alms::LowerTriangularArray) = truncate!(alms, size(alms, 2, ZeroBased, as = Matrix))


"""
$(TYPEDSIGNATURES)
Returns a LowerTriangularArray that is truncated from `alms` to the size (`ltrunc`+1) x (`mtrunc`+1),
both inputs are 0-based. If `ltrunc` or `mtrunc` is larger than the corresponding size of`alms` than
`truncate` is automatically called instead, returning a LowerTriangularArray padded zero
coefficients for higher wavenumbers."""
function truncate(
        ::Type{NF},                 # number format NF (can be complex)
        alms::LowerTriangularArray{T, N, ArrayType, S}, # spectral field to be truncated
        ltrunc::Integer,            # truncate to max degree ltrunc (0-based)
        mtrunc::Integer,            # truncate to max order mtrunc (0-based)
    ) where {NF, T, N, S, ArrayType}

    lmax, mmax, k... = size(alms, ZeroBased, as = Matrix)

    # interpolate to higher resolution if output larger than input
    (ltrunc > lmax || mtrunc > mmax) && return interpolate(NF, alms, ltrunc, mtrunc)

    # preallocate new (smaller) array
    ArrayType_ = nonparametric_type(ArrayType)
    alms_trunc = zeros(LowerTriangularArray{NF, N, ArrayType_{NF, N}, S}, Spectrum(ltrunc + 1, mtrunc + 1, architecture = architecture(alms)), k...)

    # copy data over, copyto! copies the largest matching subset of harmonics
    copyto!(alms_trunc, alms)
    return alms_trunc
end

truncate(alms::LowerTriangularArray, ltrunc::Integer, mtrunc::Integer) = truncate(eltype(alms), alms, ltrunc, mtrunc)
truncate(alms::LowerTriangularArray, trunc::Integer) = truncate(alms, trunc, trunc)

"""
$(TYPEDSIGNATURES)
Returns a LowerTriangularArray that is interpolated from `alms` to the size (`ltrunc`+1) x (`mtrunc`+1),
both inputs are 0-based, by padding zeros for higher wavenumbers. If `ltrunc` or `mtrunc` are smaller than the
corresponding size of `alms` than `truncate` is automatically called instead, returning a smaller
LowerTriangularArray."""
function interpolate(
        ::Type{NF},                 # number format NF (can be complex)
        alms::LowerTriangularArray{T, N, ArrayType, S}, # spectral field to be truncated
        ltrunc::Integer,            # truncate to max degree ltrunc (0-based)
        mtrunc::Integer,            # truncate to max order mtrunc (0-based)
    ) where {NF, T, N, S, ArrayType}

    lmax, mmax, k... = size(alms, ZeroBased, as = Matrix)

    # truncate to lower resolution if output smaller than input
    (ltrunc <= lmax && mtrunc <= mmax) && return truncate(NF, alms, ltrunc, mtrunc)

    # preallocate new (larger) array
    ArrayType_ = nonparametric_type(ArrayType)
    alms_interp = zeros(LowerTriangularArray{NF, N, ArrayType_{NF, N}, S}, Spectrum(ltrunc + 1, mtrunc + 1, architecture = architecture(alms)), k...)

    # copy data over, copyto! copies the largest matching subset of harmonics
    copyto!(alms_interp, alms)
    return alms_interp
end

interpolate(alms::LowerTriangularArray, trunc::Integer) = interpolate(alms, trunc, trunc)
