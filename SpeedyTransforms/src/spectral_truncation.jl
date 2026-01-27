# These have been moved to LowerTriangularArrays
# TODO remove, rename? Change 0 to 1-based indexing?
const spectral_truncation = LowerTriangularArrays.truncate
const spectral_truncation! = LowerTriangularArrays.truncate!
const spectral_interpolation = LowerTriangularArrays.interpolate

"""$(TYPEDSIGNATURES)
Set imaginary component of m=0 modes (the zonal modes in the first column) to 0."""
function zero_imaginary_zonal_modes!(
        alms::LowerTriangularArray,
    )
    lmax, mmax = size(alms, OneBased, as = Matrix)

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
function spectral_smoothing(A::LowerTriangularArray, args...; kwargs...)
    return spectral_smoothing!(copy(A), args...; kwargs...)
end

"""$(TYPEDSIGNATURES)
Smooth the spectral field `A` following A *= (1-(1-c)*∇²ⁿ) with power n of a normalised Laplacian
so that the highest degree lmax is dampened by multiplication with c. Anti-diffusion for c>1."""
function spectral_smoothing!(
        L::LowerTriangularArray,
        c::Real;
        power::Real = 1,          # power of Laplacian used for smoothing
        truncation::Int = -1      # smoothing wrt wavenumber (0 = largest)
    )
    lmax, mmax = size(L; as = Matrix)

    # normalize by largest eigenvalue by default, or wrt to given truncation
    eigenvalue_norm = truncation == -1 ? -mmax * (mmax + 1) : -truncation * (truncation + 1)

    # Launch kernel
    launch!(
        architecture(L), SpectralWorkOrder, size(L), spectral_smoothing_kernel!,
        L, c, power, eigenvalue_norm, L.spectrum.l_indices,
    )

    return L
end

@kernel function spectral_smoothing_kernel!(
        L,
        @Const(c),
        @Const(power),
        @Const(eigenvalue_norm),
        @Const(l_indices)
    )
    I = @index(Global, Cartesian)  # I[1] == lm, I[2] == k

    # Get the degree l for this coefficient (1-based)
    l = l_indices[I[1]]

    # Calculate eigenvalue_normalised (1-based)
    eigenvalue_normalised = -l * (l - 1) / eigenvalue_norm

    # Apply smoothing: for eigenvalue_norm < largest eigenvalue the factor becomes negative
    # set to zero in that case
    L[I] *= max(1 - (1 - c) * eigenvalue_normalised^power, 0)
end
