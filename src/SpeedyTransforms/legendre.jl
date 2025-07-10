# (inverse) legendre transform kernel, called from _legendre!
@inline function _fused_oddeven_matvec!(
    north::AbstractVector,      # output, accumulator vector, northern latitudes
    south::AbstractVector,      # output, accumulator vector, southern latitudes
    specs::AbstractMatrix,      # input, spherical harmonic coefficients
    legendre::AbstractVector,   # input, Legendre polynomials
)
    lmax, nlayers = axes(specs)             # lmax is the number of degrees at order m, 
    isoddlmax = isodd(length(lmax))
    lmax_even = length(lmax) - isoddlmax    # if lmax is odd do last odd element after the loop

    @boundscheck size(north) == size(south) || throw(DimensionMismatch)
    @boundscheck size(specs, 1) == length(legendre) || throw(DimensionMismatch)
    @boundscheck size(specs, 2) <= length(north) || throw(DimensionMismatch)
    
    @inbounds for k in nlayers
        # "even" and "odd" coined with 0-based indexing, i.e. the even l=0 mode is 1st element
        even_k = zero(eltype(south))    # dot product with elements 1, 3, 5, ...
        odd_k  = zero(eltype(north))    # dot prodcut with elements 2, 4, 6, ...

        for l in 1:2:lmax_even          # dot product in pairs for contiguous memory access
            even_k = muladd(specs[l,  k], legendre[l],  even_k)
            odd_k = muladd(specs[l+1, k], legendre[l+1], odd_k)
        end

        # now do the last row if lmax is odd, all written as muladds
        even_k = muladd(specs[end, k], isoddlmax*legendre[end], even_k)
        north[k] = muladd( 1, odd_k, even_k)    # north = even + odd
        south[k] = muladd(-1, odd_k, even_k)    # south = even - odd
    end

    return north, south
end

"""$(TYPEDSIGNATURES)
Inverse Legendre transform, batched in the vertical. Not to be used
directly, but called from transform!."""
function _legendre!(
    g_north::AbstractArray{<:Complex, 3},   # Legendre-transformed output, northern latitudes
    g_south::AbstractArray{<:Complex, 3},   # and southern latitudes
    specs::LowerTriangularArray,            # input: spherical harmonic coefficients
    S::SpectralTransform;                   # precomputed transform
    unscale_coslat::Bool = false,           # unscale by cosine of latitude on the fly?
)
    (; nlat_half) = S.grid                  # dimensions    
    (; lmax, mmax ) = S.spectrum            # 1-based max degree l, order m of spherical harmonics  
    (; legendre_polynomials) = S            # precomputed Legendre polynomials    
    (; mmax_truncation) = S                 # Legendre shortcut, shortens loop over m, 1-based  
    (; coslat⁻¹, lon_offsets ) = S
    nlayers = axes(specs, 2)                # get number of layers of specs for fewer layers than precomputed in S

    lmax = lmax-1                           # 0-based max degree l of spherical harmonics
    mmax = mmax-1                           # 0-based max order m of spherical harmonics

    @boundscheck ismatching(S, specs) || throw(DimensionMismatch(S, specs))
    @boundscheck size(g_north) == size(g_south) == (S.nfreq_max, S.nlayers, nlat_half) || throw(DimensionMismatch(S, specs))

    north = S.scratch_memory_column_north   # use scratch memory for vertically-batched dot product
    south = S.scratch_memory_column_south

    @inbounds for j in 1:nlat_half          # symmetry: loop over northern latitudes only
        g_north[:, nlayers, j] .= 0       # reset scratch memory
        g_south[:, nlayers, j] .= 0       # reset scratch memory

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
            for k in nlayers
                g_north[m, k, j] = muladd(o, north[k], g_north[m, k, j])
                g_south[m, k, j] = muladd(o, south[k], g_south[m, k, j])
            end

            lm = lm_end + 1                         # first index of next m column
        end

        if unscale_coslat
            g_north[:, nlayers, j] .*= coslat⁻¹[j]        # scale in place
            g_south[:, nlayers, j] .*= coslat⁻¹[j]
        end
    end
end

# (forward) Legendre kernel, called from _legendre!
@inline function _fused_oddeven_outer_product_accumulate!(
    specs::AbstractMatrix,      # output, accumulated spherical harmonic coefficients
    legendre::AbstractVector,   # input, Legendre polynomials
    even::AbstractVector,       # input, even harmonics
    odd::AbstractVector,        # input, odd harmonics
)
    lmax, nlayers = size(specs)
    isoddlmax = isodd(lmax)
    lmax_even = lmax - isoddlmax

    @boundscheck size(odd) == size(even) || throw(DimensionMismatch)
    @boundscheck size(specs, 1) == length(legendre) || throw(DimensionMismatch)
    @boundscheck size(specs, 2) <= length(even) || throw(DimensionMismatch)

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
(Forward) Legendre transform, batched in the vertical. Not to be used
directly, but called from transform!."""
function _legendre!(                        # GRID TO SPECTRAL
    specs::LowerTriangularArray,            # Fourier and Legendre-transformed output
    f_north::AbstractArray{<:Complex, 3},   # Fourier-transformed input, northern latitudes
    f_south::AbstractArray{<:Complex, 3},   # and southern latitudes
    S::SpectralTransform,                   # precomputed transform
)
    (; nlat) = S                            # dimensions
    (; nlat_half) = S.grid
    (; lmax, mmax) = S.spectrum             # 1-based max degree l, order m of spherical harmonics  
    (; legendre_polynomials) = S            # precomputed Legendre polynomials    
    (; mmax_truncation) = S                 # Legendre shortcut, shortens loop over m, 1-based  
    (; solid_angles, lon_offsets) = S
    nlayers = axes(specs, 2)                # get number of layers of specs for fewer layers than precomputed in S

    lmax = lmax-1                           # 0-based max degree l of spherical harmonics
    mmax = mmax-1                           # 0-based max order m of spherical harmonics

    @boundscheck ismatching(S, specs) || throw(DimensionMismatch(S, specs))
    @boundscheck size(f_north) == size(f_south) == (S.nfreq_max, S.nlayers, nlat_half) || throw(DimensionMismatch(S, specs))

    even = S.scratch_memory_column_north    # use scratch memory for outer product
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
            for k in nlayers
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

"""
$(TYPEDSIGNATURES)
Unscale by cosine of latitude on the fly.
"""
function unscale_coslat!(
    g_north::AbstractArray{<:Complex, 3}, 
    g_south::AbstractArray{<:Complex, 3}, 
    coslat⁻¹::AbstractArray{<:Real, 1},
)
    launch!(architecture(g_north), :array_3d, size(g_north), unscale_coslat_kernel!, 
            g_north, g_south, coslat⁻¹)
    synchronize(architecture(g_north))
end 

@kernel inbounds=true function unscale_coslat_kernel!(
    g_north,
    g_south,
    @Const(coslat⁻¹),
)
    i, k, j = @index(Global, NTuple)
    g_north[i, k, j] *= coslat⁻¹[j]
    g_south[i, k, j] *= coslat⁻¹[j]
end
    