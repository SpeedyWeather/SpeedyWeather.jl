"""$(TYPEDSIGNATURES)
Triangular truncation to degree `lmax` and order `mmax` (both 1-based). Truncate spectral coefficients `alms` in-place
by setting all coefficients for which the degree `l` is larger than the truncation `lmax` or order `m` larger
than the truncaction `mmax`."""
function truncate!(
        alms::LowerTriangularArray,   # spectral field to be truncated
        lmax::Integer,                # truncate to max degree (1-based)
        mmax::Integer,                # truncate to max order (1-based)
    )
    (; l_indices, m_indices) = alms.spectrum

    # Launch kernel for GPU/Reactant compatibility
    arch = architecture(alms)
    launch!(
        arch, SpectralWorkOrder, size(alms), _truncate_kernel!,
        alms.data, l_indices, m_indices, lmax, mmax
    )

    return alms
end

@kernel inbounds = true function _truncate_kernel!(data, l_indices, m_indices, lmax, mmax)
    I = @index(Global, NTuple)
    lm = I[1]

    l = l_indices[lm]
    m = m_indices[lm]

    if l > lmax || m > mmax
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
Triangular truncation of `alms` to degree and order `truncation` in-place."""
truncate!(alms::LowerTriangularArray, truncation::Integer) = truncate!(alms, truncation, truncation)

"""
$(TYPEDSIGNATURES)
Triangular truncation of `alms` to the size of it, sets additional rows to zero."""
truncate!(alms::LowerTriangularArray) = truncate!(alms, size(alms, 2, OneBased, as = Matrix))

"""$(TYPEDSIGNATURES)
Set the extra degree of `L` to zero, i.e. all harmonics with a degree ``l`` larger than the
truncation (``l > m_{max}``, 1-based), a no-op for square `LowerTriangularArray`s.
Spectral variables in SpeedyWeather carry one more degree ``l_{max} = m_{max}+1`` as vector
quantities (e.g. ``u, v`` or gradients) need it for the meridional derivative. Scalar quantities
must not use it. Same as `truncate!(L)` but branch-free (for Reactant): the kernel reads and
writes all elements, not just those of the extra degree."""
function truncate_extra_degree!(L::LowerTriangularArray)
    (; mmax, l_indices) = L.spectrum
    launch!(
        architecture(L), SpectralWorkOrder, size(L),
        truncate_extra_degree_kernel!, L, l_indices, mmax
    )
    return L
end

"""$(TYPEDSIGNATURES)
Apply `truncate_extra_degree!` to all elements of `vars`, recursing into nested `NamedTuple`s
(e.g. namespaces like `vars.prognostic.tracers`). Non-spectral elements are left untouched."""
truncate_extra_degree!(vars::NamedTuple) = (foreach(truncate_extra_degree!, values(vars)); vars)

"""$(TYPEDSIGNATURES)
No-op fallback for non-spectral `x` (e.g. grid variables or the clock), returns `x`."""
truncate_extra_degree!(x) = x

@kernel inbounds = true function truncate_extra_degree_kernel!(data, l_indices, mmax)
    I = @index(Global, Cartesian)
    l = l_indices[I[1]]
    data[I] = ifelse(l > mmax, zero(eltype(data)), data[I])
end

"""
$(TYPEDSIGNATURES)
Returns a LowerTriangularArray that is truncated from `alms` to the size `ltrunc` x `mtrunc`,
both inputs are 1-based. If `ltrunc` or `mtrunc` is larger than the corresponding size of`alms` than
`truncate` is automatically called instead, returning a LowerTriangularArray padded zero
coefficients for higher wavenumbers."""
function truncate(
        ::Type{NF},                 # number format NF (can be complex)
        alms::LowerTriangularArray{T, N, ArrayType, S}, # spectral field to be truncated
        ltrunc::Integer,            # truncate to max degree (1-based)
        mtrunc::Integer,            # truncate to max order (1-based)
    ) where {NF, T, N, S, ArrayType}

    lmax, mmax, k... = size(alms, OneBased, as = Matrix)

    # interpolate to higher resolution if output larger than input
    (ltrunc > lmax || mtrunc > mmax) && return interpolate(NF, alms, ltrunc, mtrunc)

    # preallocate new (smaller) array
    ArrayType_ = nonparametric_type(ArrayType)
    alms_trunc = zeros(LowerTriangularArray{NF, N, ArrayType_{NF, N}, S}, Spectrum(ltrunc, mtrunc, architecture = architecture(alms)), alms.dims, k...)

    # copy data over, copyto! copies the largest matching subset of harmonics
    copyto!(alms_trunc, alms)
    return alms_trunc
end

truncate(alms::LowerTriangularArray, ltrunc::Integer, mtrunc::Integer) = truncate(eltype(alms), alms, ltrunc, mtrunc)
truncate(alms::LowerTriangularArray, truncation::Integer) = truncate(alms, truncation, truncation)

"""
$(TYPEDSIGNATURES)
Returns a LowerTriangularArray that is interpolated from `alms` to the size `ltrunc` x `mtrunc`,
both inputs are 1-based, by padding zeros for higher wavenumbers. If `ltrunc` or `mtrunc` are smaller than the
corresponding size of `alms` than `truncate` is automatically called instead, returning a smaller
LowerTriangularArray."""
function interpolate(
        ::Type{NF},                 # number format NF (can be complex)
        alms::LowerTriangularArray{T, N, ArrayType, S}, # spectral field to be truncated
        ltrunc::Integer,            # truncate to max degree ltrunc (1-based)
        mtrunc::Integer,            # truncate to max order mtrunc (1-based)
    ) where {NF, T, N, S, ArrayType}

    lmax, mmax, k... = size(alms, OneBased, as = Matrix)

    # truncate to lower resolution if output smaller than input
    (ltrunc <= lmax && mtrunc <= mmax) && return truncate(NF, alms, ltrunc, mtrunc)

    # preallocate new (larger) array
    ArrayType_ = nonparametric_type(ArrayType)
    alms_interp = zeros(LowerTriangularArray{NF, N, ArrayType_{NF, N}, S}, Spectrum(ltrunc, mtrunc, architecture = architecture(alms)), alms.dims, k...)

    # copy data over, copyto! copies the largest matching subset of harmonics
    copyto!(alms_interp, alms)
    return alms_interp
end

interpolate(alms::LowerTriangularArray, ltrunc::Integer, mtrunc::Integer) = interpolate(eltype(alms), alms, ltrunc, mtrunc)
interpolate(alms::LowerTriangularArray, truncation::Integer) = interpolate(alms, truncation, truncation)
