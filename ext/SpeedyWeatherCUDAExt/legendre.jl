# convert l,m indices of a matrix (here 0-based l,m though...) to a single 1-based running index
import SpeedyWeather.LowerTriangularArrays: lm2i, get_lm_range, get_2lm_range
import SpeedyWeather.Architectures: ismatching

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
    S::SpectralTransform;               # precomputed transform
    unscale_coslat::Bool = false,       # unscale by cosine of latitude on the fly?
)
    (; nlat_half) = S.grid              # dimensions    
    (; lmax, mmax ) = S.spectrum        # 1-based max degree l, order m of spherical harmonics  
    (; legendre_polynomials) = S        # precomputed Legendre polynomials    
    (; jm_index_size, kjm_indices ) = S # kjm loop indices precomputed for threads  
    (; coslat⁻¹, lon_offsets ) = S
    # NOTE: this comes out as a range, not an integer
    nlayers = axes(specs, 2)            # get number of layers of specs for fewer layers than precomputed in S

    lmax = lmax-1                       # 0-based max degree l of spherical harmonics
    mmax = mmax-1                       # 0-based max order m of spherical harmonics

    @boundscheck ismatching(S, specs) || throw(DimensionMismatch(S, specs))
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
    
    #TODO: suboptimal adjustment to GPU made to avoid scalar indexing here. 
    #TODO: either integrate it into the kernel or write a new kernel for it
    if unscale_coslat
        @inbounds for j in 1:nlat_half          # symmetry: loop over northern latitudes only
            g_north[:, nlayers, j] .*= coslat⁻¹[j:j]        # scale in place
            g_south[:, nlayers, j] .*= coslat⁻¹[j:j]
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

        # LEGENDRE TRANSFORM
        fn = f_north[m, k, j]
        fs = f_south[m, k, j]
        @fastmath even_k = ΔΩ_rotated*(fn + fs)
        @fastmath odd_k  = ΔΩ_rotated*(fn - fs)

        # Create view of the legendre polynomials data for the lm range of interest. 
        # We can't do the saem for the specs data as @atomic cannot be used on views
        legendre_view = view(legendre_polynomials_data, lm_range, j)

        lmax_range = length(lm_range)           # number of degrees at order m, lmax-m
        isoddlmax = isodd(lmax_range)
        lmax_even = lmax_range - isoddlmax

        l = 1
        # Create 
        lm_index = Int(0)
        while l < lmax_even     # dot product in pairs for contiguous memory access
            # Calculate the even and odd harmonics for this l value
            even_spec = legendre_view[l] * even_k
            odd_spec = legendre_view[l+1] * odd_k
            # Precalculate the index for each of the 4 elements, the real and 
            # imaginary parts of the even and odd harmonics, as it can't be done 
            # on same line as the @atomic call
            # Even real part
            lm_index = lm2_range[2l - 1]
            CUDA.@atomic specs_data[lm_index, k] += even_spec.re
            # Even imaginary part
            lm_index = lm2_range[2l]
            CUDA.@atomic specs_data[lm_index, k] += even_spec.im
            # Odd real part
            lm_index = lm2_range[2l + 1]
            CUDA.@atomic specs_data[lm_index, k] += odd_spec.re
            # Odd imaginary part
            lm_index = lm2_range[2l + 2]
            CUDA.@atomic specs_data[lm_index, k] += odd_spec.im
            
            # We still increment l by 2 as we are processing 2 elements at a time
            l += 2
        end

        # now do the last row if lmax is odd
        if isoddlmax == 1
            even_spec = legendre_view[end] * even_k
            lm_index = lm2_range[end-1]
            CUDA.@atomic specs_data[lm_index, k] += legendre_view[end] * even_k.re
            lm_index = lm2_range[end]
            CUDA.@atomic specs_data[lm_index, k] += legendre_view[end] * even_k.im
        end

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
    (; nlat_half) = S.grid                  # dimensions
    (; lmax) = S.spectrum                   # 1-based max degree l, order m of spherical harmonics  
    (; legendre_polynomials) = S            # precomputed Legendre polynomials    
    (; kjm_indices, jm_index_size) = S      # Legendre shortcut, shortens loop over m, 0-based  
    (; solid_angles, lon_offsets) = S
    nlayers = axes(specs, 2)                # get number of layers of specs for fewer layers than precomputed in S

    lmax = lmax - 1                         # 0-based max degree l of spherical harmonics

    @boundscheck ismatching(S, specs) || throw(DimensionMismatch(S, specs))
    @boundscheck size(f_north) == size(f_south) == (S.nfreq_max, S.nlayers, nlat_half) || throw(DimensionMismatch(S, specs))

    fill!(specs, 0)                         # reset as we accumulate into specs

    # Reinterpret the specs data as a real array for atomic operations (makes a 
    # view into the original array with real and imaginary parts interleaved)
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
end
