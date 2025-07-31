# KernelAbstractions implementation of Legendre transform used only on GPU 

import SpeedyWeather.LowerTriangularArrays: lm2i, get_lm_range, get_2lm_range
import Atomix 

# (inverse) legendre transform kernel, called from _legendre!
@kernel inbounds=true function inverse_legendre_kernel!(
    g_north,                        # Scratch storage for legendre coefficients
    g_south,                        # before fft
    specs_data,                     # Data passed from spectral grid
    legendre_polynomials_data,      # Pre-calculated Legendre coefficients
    lmax,                           # Max l-value, from SpectralTransform struct
    lon_offsets,                    # Longitude 
    kjm_indices                     # precomputed indices for thread    
)
    tid = @index(Global, Linear)

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


"""$(TYPEDSIGNATURES)
Inverse Legendre transform, adapted for KernelAbstractions (GPU usage) and batched across j (lattitude), 
k (vertical layers) and m (spherical harmonic order). Not to be used directly, 
but called from transform! with CuArrays."""
function _legendre!(
    g_north::AbstractArray{<:Complex, 3},       # Legendre-transformed output, northern latitudes
    g_south::AbstractArray{<:Complex, 3},       # and southern latitudes
    specs::LowerTriangularArray,                # input: spherical harmonic coefficients
    scratch_memory::ColumnScratchMemory,        # scratch memory (unused here, but used in CPU _legendre!)
    S::SpectralTransform{NF,<:Architectures.GPU};             # precomputed transform
    unscale_coslat::Bool = false,               # unscale by cosine of latitude on the fly?
) where NF
    (; nlat_half) = S.grid              # dimensions    
    (; lmax ) = S.spectrum              # 1-based max degree l, order m of spherical harmonics  
    (; legendre_polynomials) = S        # precomputed Legendre polynomials    
    (; kjm_indices ) = S                # kjm loop indices precomputed for threads  
    (; coslat⁻¹, lon_offsets ) = S
    # NOTE: this comes out as a range, not an integer
    nlayers = size(specs, 2)            # get number of layers of specs for fewer layers than precomputed in S

    lmax = lmax-1                       # 0-based max degree l of spherical harmonics

    @boundscheck SpeedyTransforms.ismatching(S, specs) || throw(DimensionMismatch(S, specs))
    @boundscheck size(g_north) == size(g_south) == (S.nfreq_max, S.nlayers, nlat_half) || throw(DimensionMismatch(S, specs))
    @boundscheck nlayers <= S.nlayers || throw(DimensionMismatch(S, specs))

    g_north .= 0
    g_south .= 0
   
    # Launch the kernel with the specified configuration
    launch!(
        S.architecture,
        LinearWorkOrder,
        (S.jm_index_size*nlayers,),
        inverse_legendre_kernel!,
        g_north,
        g_south,
        specs.data,
        legendre_polynomials.data, 
        lmax,
        lon_offsets,
        kjm_indices;
    )

    synchronize(S.architecture)
    
    # unscale by cosine of latitude on the fly if requested
    if unscale_coslat
        unscale_coslat!(g_north, g_south, coslat⁻¹, architecture=S.architecture)
    end
end


@kernel inbounds=true function forward_legendre_kernel!(
    specs_data,                 # output, accumulated spherical harmonic coefficients
    legendre_polynomials_data,  # input, Legendre polynomials
    f_north,                    # input, Fourier-transformed northern latitudes
    f_south,                    # input, southern latitudes
    lmax,                       # Max l-value, from SpectralTransform struct
    lon_offsets,                # Longitude 
    solid_angles,               # Solid angles for each latitude
    kjm_indices                 # precomputed indices for thread
)
    tid = @index(Global, Linear)

    # Unpack indices from precomputed kjm_indices using single thread index
    k = kjm_indices[tid, 1]
    j = kjm_indices[tid, 2]
    m = kjm_indices[tid, 3]

    lm_range = get_lm_range(m, lmax) 
    lm2_range = get_2lm_range(m, lmax)

    ΔΩ = solid_angles[j]                # Solid angle for a grid point
    o = lon_offsets[m, j]               # Longitude offset rotation
    ΔΩ_rotated = ΔΩ * conj(o)           # Rotation back to prime meridian
        
    even_k = zero(eltype(f_north))
    odd_k = zero(eltype(f_north))
    even_spec = zero(eltype(f_north))
    odd_spec = zero(eltype(f_north))

    fn = f_north[m, k, j]
    fs = f_south[m, k, j]
    even_k = ΔΩ_rotated * (fn + fs)
    odd_k  = ΔΩ_rotated * (fn - fs)

    legendre_view = view(legendre_polynomials_data, lm_range, j)

    lmax_range = length(lm_range)
    isoddlmax = isodd(lmax_range)
    lmax_even = lmax_range - isoddlmax

    l = 1
    while l < lmax_even
        even_spec = legendre_view[l] * even_k
        odd_spec = legendre_view[l+1] * odd_k

        Atomix.@atomic specs_data[lm2_range[2l - 1], k] += even_spec.re
        Atomix.@atomic specs_data[lm2_range[2l    ], k] += even_spec.im
        Atomix.@atomic specs_data[lm2_range[2l + 1], k] += odd_spec.re
        Atomix.@atomic specs_data[lm2_range[2l + 2], k] += odd_spec.im

        l += 2
    end

    if isoddlmax == 1
        even_spec = legendre_view[end] * even_k
        Atomix.@atomic specs_data[lm2_range[end - 1], k] += even_spec.re
        Atomix.@atomic specs_data[lm2_range[end    ], k] += even_spec.im
    end
end

function _legendre!(                        # GRID TO SPECTRAL
    specs::LowerTriangularArray,            # Fourier and Legendre-transformed output
    f_north::AbstractArray{<:Complex, 3},   # Fourier-transformed input, northern latitudes
    f_south::AbstractArray{<:Complex, 3},   # and southern latitudes
    scratch_memory::ColumnScratchMemory,    # scratch memory (not used here, but for CPU _legendre!)
    S::SpectralTransform{NF,<:Architectures.GPU},        # precomputed transform
) where NF
    (; lmax) = S.spectrum                   # 1-based max degree l, order m of spherical harmonics  
    (; legendre_polynomials) = S            # precomputed Legendre polynomials    
    (; kjm_indices) = S                     # Legendre shortcut, shortens loop over m, 0-based  
    (; solid_angles, lon_offsets) = S

    nlayers = size(specs, 2)                # get number of layers of specs for fewer layers than precomputed in S

    @boundscheck SpeedyTransforms.ismatching(S, specs) || throw(DimensionMismatch(S, specs))
    @boundscheck size(f_north) == size(f_south) == (S.nfreq_max, S.nlayers, S.grid.nlat_half) || throw(DimensionMismatch(S, specs))
    @boundscheck nlayers <= S.nlayers || throw(DimensionMismatch(S, specs))

    fill!(specs, 0)                         # reset as we accumulate into specs
    lmax = lmax-1                           # 0-based max degree l of spherical harmonics

    specs_reinterpret = reinterpret(real(eltype(specs.data)), specs.data)

    launch!(
        S.architecture,
        LinearWorkOrder,
        (S.jm_index_size*nlayers,),
        forward_legendre_kernel!,
        specs_reinterpret,
        legendre_polynomials.data, 
        f_north,
        f_south,
        lmax,
        lon_offsets,
        solid_angles,
        kjm_indices;
    )

    # NOTE: synchronize here?
    synchronize(S.architecture)
end