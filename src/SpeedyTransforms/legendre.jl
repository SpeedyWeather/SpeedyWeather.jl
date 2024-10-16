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
function _legendre_unbatched!(
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

"""$(TYPEDSIGNATURES)
A matrix-vector multiplication betwen a matrix `M` and a vector `v` but split into odd and even elements
along the first dimension of `M` and `v`, i.e. it implicitly takes the transpose of `M` as argument,
so that `size(M, 1) == size(v))`. It is the generalisation of `fused_oddeven_dot` to matrices,
batching the dot product in the 2nd dimension of `M`. `odd` and `even` are the results of the dot product
and have to be provided as preallocated arrays. Returns `(odd, even)`."""
@inline function _fused_oddeven_matvec!(
    north::AbstractVector,
    south::AbstractVector,
    M::AbstractMatrix,
    v::AbstractVector,
)
    m, n = size(M)              # indexed with i, j
    isoddm = isodd(m)
    m_even = m - isoddm         # if m is odd do last odd element after the loop

    @boundscheck axes(north) == axes(south) || throw(DimensionMismatch)
    @boundscheck axes(M, 1) == axes(v) || throw(DimensionMismatch)
    @boundscheck axes(M, 2) <= axes(north) || throw(DimensionMismatch)
    
    @inbounds for j in 1:n
        odd_j  = zero(eltype(north))    # dot prodcut with elements 1, 3, 5, ... of M, v
        even_j = zero(eltype(south))    # dot product with elements 2, 4, 6, ... of M, v

        for i in 1:2:m_even
            odd_j = muladd(M[i, j], v[i], odd_j)
            even_j = muladd(M[i+1, j], v[i+1], even_j)
        end

        # now do the last row if m is odd, all written as muladds
        odd_j = muladd(M[end, j], isoddm*v[end], odd_j)
        north[j] = muladd( 1, even_j, odd_j)    # north = odd + even
        south[j] = muladd(-1, even_j, odd_j)    # sotuh = odd - even (note 1-based here)
    end

    return north, south
end

"""$(TYPEDSIGNATURES)
Inverse Legendre transform, batched in the vertical."""
function _legendre!(
    g_north::AbstractArray{<:Complex, 3},   # Legendre-transformed output, northern latitudes
    g_south::AbstractArray{<:Complex, 3},   # and southern latitudes
    specs::LowerTriangularArray,            # Legendre-transformed input
    S::SpectralTransform;                   # precomputed transform
    unscale_coslat::Bool = false,           # unscale by cosine of latitude on the fly?
)
    (; nlat_half) = S                       # dimensions    
    (; lmax, mmax ) = S                     # 0-based max degree l, order m of spherical harmonics  
    (; legendre_polynomials) = S            # precomputed Legendre polynomials    
    (; mmax_truncation) = S                 # Legendre shortcut, shortens loop over m, 0-based  
    (; coslat⁻¹, lon_offsets ) = S
    nlayers = size(specs, 2)                # get number of layers of specs for fewer layers than precomputed in S

    @boundscheck ismatching(S, specs) || throw(DimensionMismatch(S, specs))
    @boundscheck size(g_north) == size(g_south) == (S.nfreq_max, S.nlayers, nlat_half) || throw(DimensionMismatch(S, specs))

    north = S.scratch_memory_column_north   # use scratch memory for vertically-batched dot product
    south = S.scratch_memory_column_south

    @inbounds for j in 1:nlat_half          # symmetry: loop over northern latitudes only
        g_north[:, 1:nlayers, j] .= 0       # reset scratch memory
        g_south[:, 1:nlayers, j] .= 0       # reset scratch memory

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
            north, south = _fused_oddeven_matvec!(north, south, spec_view, legendre_view)

            # CORRECT FOR LONGITUDE OFFSETTS (if grid points don't start at 0°E)
            o = lon_offsets[m, j]           # rotation through multiplication with complex unit vector
            for k in 1:nlayers
                g_north[m, k, j] = muladd(o, north[k], g_north[m, k, j])
                g_south[m, k, j] = muladd(o, south[k], g_south[m, k, j])
            end

            lm = lm_end + 1                         # first index of next m column
        end

        if unscale_coslat
            g_north[:, 1:nlayers, j] .*= coslat⁻¹[j]        # scale in place
            g_south[:, 1:nlayers, j] .*= coslat⁻¹[j]
        end
    end
end

function _fused_oddeven_outer_product_accumulate!(
    specs::AbstractMatrix,
    legendre::AbstractVector,
    even::AbstractVector,
    odd::AbstractVector,
)
    lmax, nlayers = size(specs)
    isoddlmax = isodd(lmax)
    lmax_even = lmax - isoddlmax

    @inbounds for k in 1:nlayers
        even_k, odd_k = even[k], odd[k]
        for l in 1:2:lmax_even
            specs[l,   k] = muladd(legendre[l],  even_k, specs[l,   k])
            specs[l+1, k] = muladd(legendre[l+1], odd_k, specs[l+1, k])
        end
        specs[end, k] = muladd(legendre[end], isoddlmax*even_k, specs[end, k])
    end
end

"""$(TYPEDSIGNATURES)
Legendre transform, batched in the vertical."""
function _legendre!(                        # GRID TO SPECTRAL
    specs::LowerTriangularArray,            # Fourier and Legendre-transformed output
    f_north::AbstractArray{<:Complex, 3},   # Fourier-transformed input, northern latitudes
    f_south::AbstractArray{<:Complex, 3},   # and southern latitudes
    S::SpectralTransform,                   # precomputed transform
)
    (; nlat, nlat_half, nlayers) = S        # dimensions    
    (; lmax, mmax) = S                      # 0-based max degree l, order m of spherical harmonics  
    (; legendre_polynomials) = S            # precomputed Legendre polynomials    
    (; mmax_truncation) = S                 # Legendre shortcut, shortens loop over m, 0-based  
    (; solid_angles, lon_offsets) = S

    @boundscheck ismatching(S, specs) || throw(DimensionMismatch(S, specs))
    @boundscheck size(f_north) == size(f_south) == (S.nfreq_max, nlayers, nlat_half) || throw(DimensionMismatch(S, specs))

    even = S.scratch_memory_column_north     # use scratch memory for outer product
    odd = S.scratch_memory_column_south

    fill!(specs, 0)                         # reset as we accumulate into specs

    @inbounds for j_north in 1:nlat_half    # symmetry: loop over northern latitudes only
        j = j_north                         # symmetric index / ring-away from pole index

        # SOLID ANGLES including quadrature weights (sinθ Δθ) and azimuth (Δϕ) on ring j
        ΔΩ = solid_angles[j]                # = sinθ Δθ Δϕ, solid angle for a grid point

        lm = 1                              # single running index for spherical harmonics
        for m in 1:mmax_truncation[j] + 1   # Σ_{m=0}^{mmax}, but 1-based index, shortened to mmax_truncation

            # SOLID ANGLE QUADRATURE WEIGHTS and LONGITUDE OFFSET
            o = lon_offsets[m, j]           # longitude offset rotation by multiplication with complex unit vector
            ΔΩ_rotated = ΔΩ*conj(o)         # complex conjugate for rotation back to prime meridian
            
            # LEGENDRE TRANSFORM
            for k in 1:nlayers
                fn, fs  = f_north[m, k, j], f_south[m, k, j]
                @fastmath even[k] = ΔΩ_rotated*(fn + fs)
                @fastmath odd[k]  = ΔΩ_rotated*(fn - fs)
            end

            # integration over l = m:lmax+1
            lm_end = lm + lmax-m+1                      # last index in column m
            spec_view = view(specs.data, lm:lm_end, :)
            legendre_view = view(legendre_polynomials.data, lm:lm_end, j)

            _fused_oddeven_outer_product_accumulate!(spec_view, legendre_view, even, odd)

            lm = lm_end + 1                             # first index of next column m+1
        end
    end
end