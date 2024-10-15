"""$(TYPEDSIGNATURES)
A fused dot product of two vectors `a` and `b` but split into odd and even elements,
i.e. `odd = a1*b1 + a3*b3 + ...` and `even = a2*b2 + a4*b4 + ...`. Returns `(odd, even)`."""
@inline function fused_oddeven_dot(a::AbstractVector, b::AbstractVector)
    @boundscheck axes(a) == axes(b) || throw(DimensionMismatch)
    Base.require_one_based_indexing(a, b)

    odd  = zero(eltype(a))      # dot prodcut with elements 1, 3, 5, ... of a, b
    even = zero(eltype(a))      # dot product with elements 2, 4, 6, ... of a, b
    
    n = length(a)
    isoddn = isodd(n)
    n_even = n - isoddn         # if n is odd do last odd element after the loop
    
    @inbounds for i in 1:2:n_even
        odd = muladd(a[i], b[i], odd)
        even = muladd(a[i+1], b[i+1], even)
    end
    
    # now do the last element if n is odd
    odd = muladd(a[end], isoddn*b[end], odd)
    return odd, even
end

"""$(TYPEDSIGNATURES)
Inverse Legendre transform of the spherical harmonic coefficients `specs`,
to be stored in `g_north` and `g_south` for northern and southern latitudes respectively.
Not to be called directly, use `transform!` instead."""
function _legendre!(
    g_north::AbstractArray{<:Complex, 3},   # Legendre-transformed output, northern latitudes
    g_south::AbstractArray{<:Complex, 3},   # and southern latitudes
    specs::LowerTriangularArray,            # Legendre-transformed input
    S::SpectralTransform;                   # precomputed transform
    unscale_coslat::Bool = false,           # unscale by cosine of latitude on the fly?
)
    (; nlat_half, nlayers) = S              # dimensions    
    (; lmax, mmax ) = S                     # 0-based max degree l, order m of spherical harmonics  
    (; legendre_polynomials) = S            # precomputed Legendre polynomials    
    (; mmax_truncation) = S                 # Legendre shortcut, shortens loop over m, 0-based  
    (; coslat⁻¹, lon_offsets ) = S

    @boundscheck ismatching(S, specs) || throw(DimensionMismatch(S, specs))
    @boundscheck size(g_north) == size(g_south) == (S.nfreq_max, nlayers, nlat_half) || throw(DimensionMismatch(S, specs))

    @inbounds for k in eachmatrix(specs)            # loop over all specs/grids (e.g. vertical dimension)
        for j in 1:nlat_half                        # symmetry: loop over northern latitudes only
            g_north[:, k, j] .= 0                   # reset scratch memory
            g_south[:, k, j] .= 0                   # reset scratch memory

            # INVERSE LEGENDRE TRANSFORM by looping over wavenumbers l, m
            lm = 1                                  # single running index for non-zero l, m indices
            for m in 1:mmax_truncation[j] + 1       # Σ_{m=0}^{mmax}, but 1-based index, shortened to mmax_truncation
                lm_end = lm + lmax-m+1              # last index in column

                # view on lower triangular column
                spec_view = view(specs.data, lm:lm_end, k)
                legendre_view = view(legendre_polynomials.data, lm:lm_end, j)
                
                # dot product but split into even and odd harmonics on the fly for better performance
                # function is 1-based (odd, even, odd, ...) but here use 0-based indexing to name
                # the "even" and "odd" harmonics
                even, odd = fused_oddeven_dot(spec_view, legendre_view)

                # # CORRECT FOR LONGITUDE OFFSETTS (if grid points don't start at 0°E)
                o = lon_offsets[m, j]
                g_north[m, k, j] += o*(even+odd)    # accumulate in phase factors for northern
                g_south[m, k, j] += o*(even-odd)    # and southern hemisphere

                lm = lm_end + 1                         # first index of next m column
            end

            if unscale_coslat
                g_north[:, k, j] .*= coslat⁻¹[j]        # scale in place
                g_south[:, k, j] .*= coslat⁻¹[j]
            end
        end
    end
end

@inline function fused_oddeven_matvec!(
    odd::AbstractVector,
    even::AbstractVector,
    M::AbstractMatrix,
    v::AbstractVector,
)
    m, n = size(M)              # indexed with i, j
    isoddm = isodd(m)
    m_even = m - isoddm         # if m is odd do last odd element after the loop

    @boundscheck axes(odd) == axes(even) || throw(DimensionMismatch)
    @boundscheck axes(M, 1) == axes(v) || throw(DimensionMismatch)
    @boundscheck axes(M, 2) == axes(odd) || throw(DimensionMismatch)
    
    @inbounds for j in eachindex(odd, even)
        odd_j  = zero(eltype(odd))
        even_j = zero(eltype(even))

        for i in 1:2:m_even
            odd_j = muladd(M[i, j], v[i], odd_j)
            even_j = muladd(M[i+1, j], v[i+1], even_j)
        end

        # now do the last row if m is odd
        odd[j] = muladd(M[end, j], isoddm*v[end], odd_j)
        even[j] = even_j
    end

    return odd, even
end

"""$(TYPEDSIGNATURES)
Inverse Legendre transform, batched in the vertical."""
function _legendre_batched!(
    g_north::AbstractArray{<:Complex, 3},   # Legendre-transformed output, northern latitudes
    g_south::AbstractArray{<:Complex, 3},   # and southern latitudes
    specs::LowerTriangularArray,            # Legendre-transformed input
    S::SpectralTransform;                   # precomputed transform
    unscale_coslat::Bool = false,           # unscale by cosine of latitude on the fly?
)
    (; nlat_half, nlayers) = S              # dimensions    
    (; lmax, mmax ) = S                     # 0-based max degree l, order m of spherical harmonics  
    (; legendre_polynomials) = S            # precomputed Legendre polynomials    
    (; mmax_truncation) = S                 # Legendre shortcut, shortens loop over m, 0-based  
    (; coslat⁻¹, lon_offsets ) = S

    @boundscheck ismatching(S, specs) || throw(DimensionMismatch(S, specs))
    @boundscheck size(g_north) == size(g_south) == (S.nfreq_max, nlayers, nlat_half) || throw(DimensionMismatch(S, specs))

    odd  = S.scratch_memory_column_odd      # use scratch memory for vertically-batched dot product
    even = S.scratch_memory_column_even

    @inbounds for j in 1:nlat_half          # symmetry: loop over northern latitudes only
        g_north[:, :, j] .= 0               # reset scratch memory
        g_south[:, :, j] .= 0               # reset scratch memory

        # INVERSE LEGENDRE TRANSFORM by looping over wavenumbers l, m
        lm = 1                              # single running index for non-zero l, m indices
        for m in 1:mmax_truncation[j] + 1   # Σ_{m=0}^{mmax}, but 1-based index, shortened to mmax_truncation
            lm_end = lm + lmax-m+1          # last index in column

            # view on lower triangular column, but batched in vertical
            spec_view = view(specs.data, lm:lm_end, :)
            legendre_view = view(legendre_polynomials.data, lm:lm_end, j)

            # dot product but split into even and odd harmonics on the fly for better performance
            # function is 1-based (odd, even, odd, ...) but here use 0-based indexing to name
            # the "even" and "odd" harmonics, batched in the vertical so it's a mat vec multiplication
            even, odd = fused_oddeven_matvec!(even, odd, spec_view, legendre_view)

            # CORRECT FOR LONGITUDE OFFSETTS (if grid points don't start at 0°E)
            o = lon_offsets[m, j]           # rotation through multiplication with complex unit vector
            for k in eachindex(even, odd)
                g_north[m, k, j] += o*(even[k] + odd[k])
                g_south[m, k, j] += o*(even[k] - odd[k])
            end

            lm = lm_end + 1                         # first index of next m column
        end

        if unscale_coslat
            g_north[:, :, j] .*= coslat⁻¹[j]        # scale in place
            g_south[:, :, j] .*= coslat⁻¹[j]
        end
    end
end