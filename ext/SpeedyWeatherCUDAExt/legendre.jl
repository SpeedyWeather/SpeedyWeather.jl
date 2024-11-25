# Convenience functions for calculating lm and lm_end from a given m within a 
# Kernel
function get_lm_end(m, lmax)
    return Int((m*lmax) - ((m-2)*m) + 0.5*m*(m-1))
end

function get_lm(m, lmax)
    return get_lm_end(m-1, lmax) + 1
end

# (inverse) legendre transform kernel, called from _legendre!
function phase_factor_kernel!(
    g_north,
    g_south,
    specs_data,
    legendre_polynomials_data, 
    lmax,
    lon_offsets,
    kjm_indices
)
    tid = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x

    if tid <= length(kjm_indices)
        k, j, m  = kjm_indices[tid]

        lm = get_lm(m, lmax)
        lm_end = lm + lmax-m+1          # last index in column

        # view on lower triangular column, but batched in vertical
        spec_view = view(specs_data, lm:lm_end, :)
        legendre_view = view(legendre_polynomials_data, lm:lm_end, j)


        # dot product but split into even and odd harmonics on the fly as this 
        # is how the previous implementation was enacted

        lmax_range, _ = axes(spec_view)             # lmax is the number of degrees at order m, 
        isoddlmax = isodd(length(lmax_range))
        lmax_even = length(lmax_range) - isoddlmax    # if lmax is odd do last odd element after the loop

        # Got rid of bounds check, potentially unsafe?
        # @boundscheck size(north) == size(south) || throw(DimensionMismatch)
        # @boundscheck size(spec_view, 1) == length(legendre_view) || throw(DimensionMismatch)
        # @boundscheck size(spec_view, 2) <= length(north) || throw(DimensionMismatch)

        # "even" and "odd" coined with 0-based indexing, i.e. the even l=0 mode is 1st element
        even_k = zero(eltype(g_south))    # dot product with elements 1, 3, 5, ...
        odd_k  = zero(eltype(g_north))    # dot prodcut with elements 2, 4, 6, ...

        # Switched to while loop as more performant from inside a Kernel     
        l = 1
        while l < lmax_even     # dot product in pairs for contiguous memory access
            even_k += spec_view[l, k] * legendre_view[l]
            odd_k += spec_view[l+1, k] * legendre_view[l+1]
            l += 2
        end

        # now do the last row if lmax is odd
        even_k += spec_view[end, k] * (isoddlmax * legendre_view[end])
        north = even_k + odd_k
        south = even_k - odd_k

        # CORRECT FOR LONGITUDE OFFSETTS (if grid points don't start at 0°E)
        o = lon_offsets[m, j]           # rotation through multiplication with complex unit vector

        g_north[m, k, j] += o * north
        g_south[m, k, j] += o * south
    end
    return 
end


"""$(TYPEDSIGNATURES)
Inverse Legendre transform, adapted for CUDA and batched across j (lattitude), 
k (vertical layers) and m (spherical harmonic order). Not to be used directly, 
but called from transform! with CuArrays."""
function SpeedyTransforms._legendre!(
    g_north::CuArray{<:Complex, 3},   # Legendre-transformed output, northern latitudes
    g_south::CuArray{<:Complex, 3},   # and southern latitudes
    specs::LowerTriangularArray,            # input: spherical harmonic coefficients
    S::SpectralTransform,                   # precomputed transform
    kjm_indices::CuArray;            # precomputed jm index map
    unscale_coslat::Bool = false,           # unscale by cosine of latitude on the fly?
)
    (; nlat_half) = S                       # dimensions    
    (; lmax, mmax ) = S                     # 0-based max degree l, order m of spherical harmonics  
    (; legendre_polynomials) = S            # precomputed Legendre polynomials    
    (; mmax_truncation) = S                 # Legendre shortcut, shortens loop over m, 0-based  
    (; coslat⁻¹, lon_offsets ) = S
    nlayers = axes(specs, 2)                # get number of layers of specs for fewer layers than precomputed in S

    @boundscheck SpeedyTransforms.ismatching(S, specs) || throw(DimensionMismatch(S, specs))
    @boundscheck size(g_north) == size(g_south) == (S.nfreq_max, S.nlayers, nlat_half) || throw(DimensionMismatch(S, specs))

    g_north .= 0
    g_south .= 0

    # INVERSE LEGENDRE TRANSFORM by looping over wavenumbers l, m and layer k
    k = CUDA.@cuda launch=false phase_factor_kernel!(
        g_north,
        g_south,
        specs.data,
        legendre_polynomials.data, 
        lmax,
        lon_offsets,
        kjm_indices
    )
    config = CUDA.launch_configuration(k.fun)
    threads = min(length(kjm_indices), config.threads)
    blocks = cld(length(kjm_indices), threads)
    k(
        g_north,
        g_south,
        specs.data,
        legendre_polynomials.data, 
        lmax,
        lon_offsets,
        kjm_indices; 
        threads, 
        blocks
    )
    
    if unscale_coslat
        @inbounds for j in 1:nlat_half          # symmetry: loop over northern latitudes only
            g_north[:, nlayers, j] .*= coslat⁻¹[j]        # scale in place
            g_south[:, nlayers, j] .*= coslat⁻¹[j]
        end
    end
end