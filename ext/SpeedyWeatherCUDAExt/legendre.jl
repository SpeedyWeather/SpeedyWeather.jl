# convert i, j indices of a matrix (here 0-based l,m though...) to a single 1-based running index
import SpeedyWeather.LowerTriangularMatrices: ij2k

# range of the running indices lm in a l-column (degrees of spherical harmonics)
# given the column index m (order of harmonics) 
get_lm_range(m, lmax) = ij2k(2*m - 1, m, lmax):ij2k(lmax+m, m, lmax)
get_2lm_range(m, lmax) = 2*ij2k(2*m - 1, m, lmax)-1:2*ij2k(lmax+m, m, lmax)
 
# (inverse) legendre transform kernel, called from _legendre!
function inverse_legendre_kernel!(
    g_north,                        # Scratch storage for legendre coefficients
    g_south,                        # before fft
    specs_data,                     # Data passed from spectral grid
    legendre_polynomials_data,      # Pre-calculated Legendre coefficients
    lmax,                           # Max l-value, from SpectralTransform struct
    lon_offsets,                    # Longitude 
    kjm_indices
)
    tid = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x

    if tid <= size(kjm_indices, 1)
        # Unpack indices from precomputed kjm_indices using single thread index
        k = kjm_indices[tid, 1]
        j = kjm_indices[tid, 2]
        m = kjm_indices[tid, 3]

        # are m, lmax 0-based here or 1-based? 
        lm_range = get_lm_range(m, lmax)    # assumes 1-based

        # view on lower triangular column, but batched in vertical
        spec_view = view(specs_data, lm_range, :)
        legendre_view = view(legendre_polynomials_data, lm_range, j)
        
        # dot product but split into even and odd harmonics on the fly as this 
        # is how the previous implementation was enacted
        lmax_range = length(lm_range)           # number of degrees at order m, lmax-m
        isoddlmax = isodd(lmax_range)
        lmax_even = lmax_range - isoddlmax     # if odd do last odd element after the loop

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
    g_north::CuArray{<:Complex, 3},     # Legendre-transformed output, northern latitudes
    g_south::CuArray{<:Complex, 3},     # and southern latitudes
    specs::LowerTriangularArray,        # input: spherical harmonic coefficients
    S::SpectralTransform,               # precomputed transform
    unscale_coslat::Bool = false,       # unscale by cosine of latitude on the fly?
)
    (; nlat_half) = S                   # dimensions    
    (; lmax, mmax ) = S                 # 0-based max degree l, order m of spherical harmonics  
    (; legendre_polynomials) = S        # precomputed Legendre polynomials    
    (; jm_index_size, kjm_indices ) = S # kjm loop indices precomputed for threads  
    (; coslat⁻¹, lon_offsets ) = S
    # NOTE: this comes out as a range, not an integer
    nlayers = axes(specs, 2)            # get number of layers of specs for fewer layers than precomputed in S

    @boundscheck SpeedyTransforms.ismatching(S, specs) || throw(DimensionMismatch(S, specs))
    @boundscheck size(g_north) == size(g_south) == (S.nfreq_max, S.nlayers, nlat_half) || throw(DimensionMismatch(S, specs))
    # reduced_kjm = kjm_indices[1:(nlayers.stop * jm_index_size), :]  # get the reduced kjm indices

    # @show reduced_kjm
    g_north .= 0
    g_south .= 0

    # INVERSE LEGENDRE TRANSFORM by looping over wavenumbers l, m and layer k
    kernel = CUDA.@cuda launch=false inverse_legendre_kernel!(
        g_north,
        g_south,
        specs.data,
        legendre_polynomials.data, 
        lmax,
        lon_offsets,
        kjm_indices
    )
    config = CUDA.launch_configuration(kernel.fun)
    threads = min(size(kjm_indices, 1), config.threads)
    blocks = cld(size(kjm_indices, 1), threads)

    # actually launch kernel!
    kernel(
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


# (forward) Legendre kernel, called from _legendre!
function forward_legendre_kernel!(
    specs_data,                 # output, accumulated spherical harmonic coefficients
    legendre_polynomials_data,  # input, Legendre polynomials
    f_north,                    # input, Fourier-transformed northern latitudes
    f_south,                    # input, southern latitudes
    lmax,                       # Max l-value, from SpectralTransform struct
    lon_offsets,                # Longitude 
    solid_angles,               # Solid angles for each latitude
    kjm_indices                 # precomputed indices for thread
)
    tid = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x

    if tid <= size(kjm_indices, 1)
        # Unpack indices from precomputed kjm_indices using single thread index
        k = kjm_indices[tid, 1]
        j = kjm_indices[tid, 2]
        m = kjm_indices[tid, 3]

        # j = j_north                         # symmetric index / ring-away from pole index
        lm_range = get_lm_range(m, lmax) 
        # This is the lm_range for the reinterpreted specs_data, taking into 
        # account real and imaginary parts
        lm2_range = get_2lm_range(m, lmax)

        # SOLID ANGLES including quadrature weights (sinθ Δθ) and azimuth (Δϕ) on ring j
        ΔΩ = solid_angles[j]                # = sinθ Δθ Δϕ, solid angle for a grid point

        # SOLID ANGLE QUADRATURE WEIGHTS and LONGITUDE OFFSET
        o = lon_offsets[m, j]           # longitude offset rotation by multiplication with complex unit vector
        ΔΩ_rotated = ΔΩ*conj(o)         # complex conjugate for rotation back to prime meridian
        
        even_k = zero(eltype(f_north))
        odd_k = zero(eltype(f_north))
        even_spec = zero(eltype(f_north))
        odd_spec = zero(eltype(f_north))

        # CUDA.@cuprintln("reached 1, ($k, $j, $m)")

        # LEGENDRE TRANSFORM
        fn = f_north[m, k, j]
        fs = f_south[m, k, j]
        @fastmath even_k = ΔΩ_rotated*(fn + fs)
        @fastmath odd_k  = ΔΩ_rotated*(fn - fs)

        # CUDA.@cuprintln("reached 2, ($k, $j, $m)")

        # integration over l = m:lmax+1
        # lm_end = lm + lmax-m+1                      # last index in column m
        # spec_view = view(scratch_memory, lm_range, :, j)
        # legendre_view = view(legendre_polynomials_data, lm_range, j)
        legendre_view = view(legendre_polynomials_data, lm_range, j)

        # spec_view_ = reinterpret(real(eltype(spec_view)), spec_view)

        lmax_range = length(lm_range)           # number of degrees at order m, lmax-m
        isoddlmax = isodd(lmax_range)
        lmax_even = lmax_range - isoddlmax

        # CUDA.@cushow (lm_range.start, lm_range.stop, lmax_range, lmax_even)
        # CUDA.@cuprintln("reached 3, ($k, $j, $m)")
    
        # @boundscheck size(spec_view, 1) == length(legendre_view) || throw(DimensionMismatch)
    
        # even_k, odd_k = even[k], odd[k]
        # for l in 1:2:lmax_even
        l = 1
        # f_north[m, k, j] = 0
        # f_south[m, k, j] = 0
        lm_index = Int(0)
        while l < lmax_even     # dot product in pairs for contiguous memory access
            even_spec = legendre_view[l] * even_k
            odd_spec = legendre_view[l+1] * odd_k
            # precalculate the index for each of the 4 elements, the real and 
            # imaginary parts of the even and odd harmonics
            lm_index = lm2_range[2l - 1]
            CUDA.@atomic specs_data[lm_index, k] += even_spec.re
            lm_index = lm2_range[2l]
            CUDA.@atomic specs_data[lm_index, k] += even_spec.im
            lm_index = lm2_range[2l + 1]
            CUDA.@atomic specs_data[lm_index, k] += odd_spec.re
            lm_index = lm2_range[2l + 2]
            CUDA.@atomic specs_data[lm_index, k] += odd_spec.im
            
            l += 2
        end

        if isoddlmax == 1
            even_spec = legendre_view[end] * even_k
            lm_index = lm2_range[end-1]
            CUDA.@atomic specs_data[lm_index, k] += legendre_view[end] * even_k.re
            lm_index = lm2_range[end]
            CUDA.@atomic specs_data[lm_index, k] += legendre_view[end] * even_k.im
            # CUDA.@atomic spec_view[end-1, k] += legendre_view[end] * even_k.re
            # CUDA.@atomic spec_view[end,   k] += legendre_view[end] * even_k.im
        end
        # spec_view[end, k] += 1 * isoddlmax
        # spec_view[end, k] += 1

        # CUDA.@cuprintln("reached 4, ($k, $j, $m)")
    end
    return
end    


"""$(TYPEDSIGNATURES)
(Forward) Legendre transform, batched in the vertical. Not to be used
directly, but called from transform!."""
function SpeedyTransforms._legendre!(                        # GRID TO SPECTRAL
    specs::LowerTriangularArray,            # Fourier and Legendre-transformed output
    f_north::CuArray{<:Complex, 3},         # Fourier-transformed input, northern latitudes
    f_south::CuArray{<:Complex, 3},         # and southern latitudes
    S::SpectralTransform,                   # precomputed transform
)
    (; nlat_half) = S                       # dimensions
    (; lmax) = S                            # 0-based max degree l, order m of spherical harmonics  
    (; legendre_polynomials) = S            # precomputed Legendre polynomials    
    (; kjm_indices, jm_index_size) = S      # Legendre shortcut, shortens loop over m, 0-based  
    (; solid_angles, lon_offsets) = S
    nlayers = axes(specs, 2)                # get number of layers of specs for fewer layers than precomputed in S

    @boundscheck SpeedyTransforms.ismatching(S, specs) || throw(DimensionMismatch(S, specs))
    @boundscheck size(f_north) == size(f_south) == (S.nfreq_max, S.nlayers, nlat_half) || throw(DimensionMismatch(S, specs))

    # even = S.scratch_memory_column_north    # use scratch memory for outer product
    # odd = S.scratch_memory_column_south

    # fill!(S.scratch_memory_legendre, 0)
    fill!(specs, 0)                         # reset as we accumulate into specs
    specs_reinterpret = reinterpret(real(eltype(specs.data)), specs.data)

    # INVERSE LEGENDRE TRANSFORM by looping over wavenumbers l, m and layer k
    kernel = CUDA.@cuda launch=false forward_legendre_kernel!(
        specs_reinterpret,
        legendre_polynomials.data, 
        f_north,
        f_south,
        lmax,
        lon_offsets,
        solid_angles,
        kjm_indices
    )
    config = CUDA.launch_configuration(kernel.fun)
    threads = min(size(kjm_indices, 1), config.threads)
    blocks = cld(size(kjm_indices, 1), threads)

    # actually launch kernel!
    kernel(
        specs_reinterpret,
        legendre_polynomials.data, 
        f_north,
        f_south,
        lmax,
        lon_offsets,
        solid_angles,
        kjm_indices; 
        threads, 
        blocks
    )
    CUDA.synchronize()

    # Reduce across the extra dimension to get the final result
    # Base.mapreducedim!(identity, +, specs.data, S.scratch_memory_legendre);
end

function forward_legendre_kernel_alt!(
    specs_data,                 # output, accumulated spherical harmonic coefficients
    legendre_polynomials_data,  # input, Legendre polynomials
    f_north,                    # input, Fourier-transformed northern latitudes
    f_south,                    # input, southern latitudes
    lmax,                       # Max l-value, from SpectralTransform struct
    lon_offsets,                # Longitude 
    solid_angles,               # Solid angles for each latitude
    km_indices,                 # precomputed indices for thread
    j                           # latitude index
)
    tid = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x

    if tid <= size(km_indices, 1)
        # Unpack indices from precomputed kjm_indices using single thread index
        k = km_indices[tid, 1]
        m = km_indices[tid, 2]
        
        # j = j_north                         # symmetric index / ring-away from pole index
        lm_range = get_lm_range(m, lmax) 

        # SOLID ANGLES including quadrature weights (sinθ Δθ) and azimuth (Δϕ) on ring j
        ΔΩ = solid_angles[m]                # = sinθ Δθ Δϕ, solid angle for a grid point

        # SOLID ANGLE QUADRATURE WEIGHTS and LONGITUDE OFFSET
        o = lon_offsets[m]           # longitude offset rotation by multiplication with complex unit vector
        ΔΩ_rotated = ΔΩ*conj(o)         # complex conjugate for rotation back to prime meridian
        
        even_k = zero(eltype(f_north))
        odd_k = zero(eltype(f_north))

        # CUDA.@cuprintln("reached 1, ($k, $j, $m)")

        # LEGENDRE TRANSFORM
        fn = f_north[m, k, j]
        fs = f_south[m, k, j]
        @fastmath even_k = ΔΩ_rotated*(fn + fs)
        @fastmath odd_k  = ΔΩ_rotated*(fn - fs)

        # CUDA.@cuprintln("reached 2, ($k, $j, $m)")

        # integration over l = m:lmax+1
        # lm_end = lm + lmax-m+1                      # last index in column m
        spec_view = view(specs_data, lm_range, k)
        legendre_view = view(legendre_polynomials_data, lm_range, j)

        lmax_range = length(lm_range)           # number of degrees at order m, lmax-m
        isoddlmax = isodd(lmax_range)
        lmax_even = lmax_range - isoddlmax

        l = 1
        while l < lmax_even     # dot product in pairs for contiguous memory access
            spec_view[l] += legendre_view[l] * even_k
            spec_view[l+1] += legendre_view[l+1] * odd_k
            l += 2
        end
        spec_view[end] += legendre_view[end] * isoddlmax*even_k
    end
    return
end

function SpeedyTransforms._legendre!(                        # GRID TO SPECTRAL
    specs::LowerTriangularArray,            # Fourier and Legendre-transformed output
    f_north::CuArray{<:Complex, 3},         # Fourier-transformed input, northern latitudes
    f_south::CuArray{<:Complex, 3},         # and southern latitudes
    S::SpectralTransform,                   # precomputed transform
    alt_fl::Bool
)
    (; nlat_half) = S                       # dimensions
    (; lmax) = S                            # 0-based max degree l, order m of spherical harmonics  
    (; mmax_truncation) = S                 # precomputed mmax truncation
    (; legendre_polynomials) = S            # precomputed Legendre polynomials    
    (; solid_angles, lon_offsets) = S
    nlayers = axes(specs, 2)                # get number of layers of specs for fewer layers than precomputed in S

    @boundscheck SpeedyTransforms.ismatching(S, specs) || throw(DimensionMismatch(S, specs))
    @boundscheck size(f_north) == size(f_south) == (S.nfreq_max, S.nlayers, nlat_half) || throw(DimensionMismatch(S, specs))

    fill!(specs, 0)                         # reset as we accumulate into specs
    # FORWARD LEGENDRE TRANSFORM by looping over wavenumbers l, m and layer k
    # loop manually over j to avoid race condition
    for j in 1:nlat_half  
        # Construct our km_indices array for each j
        km_indices = zeros(Int, length(nlayers) * (mmax_truncation[j]+1), 2)
        i = 0
        for k in nlayers
            for m in 1:mmax_truncation[j]+1
                i += 1
                km_indices[i, :] .= [k, m]
            end
        end
        
        km_indices = CUDA.cu(km_indices)
        kernel = CUDA.@cuda launch=false forward_legendre_kernel_alt!(
            specs.data,
            legendre_polynomials.data, 
            f_north,
            f_south,
            lmax,
            lon_offsets,
            solid_angles,
            km_indices,
            j
        )
        config = CUDA.launch_configuration(kernel.fun)
        threads = min(size(km_indices, 1), config.threads)
        blocks = cld(size(km_indices, 1), threads)

        # actually launch kernel!
        kernel(
            specs.data,
            legendre_polynomials.data, 
            f_north,
            f_south,
            lmax,
            lon_offsets,
            solid_angles,
            km_indices,
            j; 
            threads, 
            blocks
        )
        CUDA.synchronize()
    end
end