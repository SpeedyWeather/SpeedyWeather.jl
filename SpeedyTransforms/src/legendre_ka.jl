# KernelAbstractions implementation of Legendre transform used only on GPU

"""$(TYPEDSIGNATURES)
`muladd` for a complex number scaled by a real one, i.e. `a + x*z` with the real and imaginary
parts fused separately. Julia does not contract `a + x*z` by itself, and calling `muladd` on the
complex number directly promotes the real factor to a complex one, turning this into a full
complex multiply-add (measurably slower). `muladd` rather than `fma` on the components: it
contracts into a fused multiply-add wherever the hardware has one and falls back to a plain
multiply and add where it does not, whereas `fma` would force a slow software emulation there."""
@inline muladd_complex(x::Real, z::Complex, a::Complex) =
    Complex(muladd(x, real(z), real(a)), muladd(x, imag(z), imag(a)))

# (inverse) legendre transform kernel, called from _legendre!
# One thread per (latitude ring j, order m, layer k); the (j, m) pair and the offset/length of
# the lower-triangular column at order m come from the precomputed `jm_indices` table, whose rows
# are ordered by order m first so that neighbouring threads hold neighbouring rings. That way the
# transposed Legendre polynomials are read coalesced and the coefficients, identical across the
# rings of one order, are a broadcast within the warp.
@kernel inbounds = true function inverse_legendre_kernel!(
        g_north,                        # Scratch storage for legendre coefficients
        g_south,                        # before fft
        specs_data,                     # Data passed from spectral grid
        legendre_polynomials_data,      # input, pre-calculated Legendre polynomials, transposed
        lon_offsets,                    # input, longitude offset rotation per (m, j)
        jm_indices,                     # precomputed (j, m, lm_offset, column length) per row
        nlat_half,                      # stride between harmonics in the transposed polynomials
    )
    i, k = @index(Global, NTuple)

    j = jm_indices[i, 1]
    m = jm_indices[i, 2]
    lm_offset = jm_indices[i, 3]        # offset to index the lower triangular column directly
    lmax_range = jm_indices[i, 4]       # number of degrees at order m, lmax-m

    # dot product but split into even and odd harmonics on the fly as this
    # is how the previous implementation was enacted
    isoddlmax = isodd(lmax_range)
    lmax_even = lmax_range - isoddlmax  # if odd do last odd element after the loop

    # "even" and "odd" coined with 0-based indexing, i.e. the even l=0 mode is 1st element
    even_k = zero(eltype(g_south))      # dot product with elements 1, 3, 5, ...
    odd_k = zero(eltype(g_north))       # dot product with elements 2, 4, 6, ...

    # a while loop here because the degree l and the index p into the transposed polynomials
    # advance together, and the last degree is handled after the loop when the column is odd
    p = lm_offset * nlat_half + j       # transposed polynomials at (j, lm_offset + 1)
    l = 1
    while l < lmax_even                 # dot product in pairs for contiguous memory access
        even_k = muladd_complex(legendre_polynomials_data[p], specs_data[lm_offset + l, k], even_k)
        odd_k = muladd_complex(legendre_polynomials_data[p + nlat_half], specs_data[lm_offset + l + 1, k], odd_k)
        p += 2nlat_half
        l += 2
    end

    if isoddlmax                        # now do the last row if the column length is odd
        even_k = muladd_complex(legendre_polynomials_data[p], specs_data[lm_offset + lmax_range, k], even_k)
    end

    # CORRECT FOR LONGITUDE OFFSETTS (if grid points don't start at 0°E)
    # rotation through multiplication with complex unit vector, zero for the rows beyond the
    # Legendre truncation: those write the zeros the Fourier transform reads for the frequencies
    # the Legendre transform never touches (`lon_offsets` isn't defined for them either)
    o = lmax_range > 0 ? lon_offsets[m, j] : zero(eltype(lon_offsets))

    g_north[m, k, j] = o * (even_k + odd_k)
    g_south[m, k, j] = o * (even_k - odd_k)
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
        S::SpectralTransform{NF, <:Architectures.GPU};             # precomputed transform
        unscale_coslat::Bool = false,               # unscale by cosine of latitude on the fly?
    ) where {NF}

    (; nlat_half) = S.grid              # dimensions
    (; legendre_polynomials_transposed) = S     # precomputed Legendre polynomials, (j, lm)
    (; jm_indices) = S                  # (j, m) loop indices precomputed for threads
    (; coslat⁻¹, lon_offsets) = S
    nlayers = size(specs, 2)            # get number of layers of specs for fewer layers than precomputed in S

    @boundscheck SpeedyTransforms.ismatching(S, specs) || throw(DimensionMismatch(S, specs))
    # Scratch dim 2 is the per-call capacity (= max(planned_K) on CPU, S.nlayers elsewhere);
    # allow it to exceed nlayers so chunked/sliced callers (e.g. test scratches copied
    # from CPU) still pass.
    @boundscheck (size(g_north) == size(g_south) && size(g_north, 1) == S.nfreq_max &&
                  size(g_north, 3) == nlat_half && size(g_north, 2) >= nlayers) ||
                 throw(DimensionMismatch(S, specs))
    @boundscheck nlayers <= S.nlayers || throw(DimensionMismatch(S, specs))

    # no need to reset g_north/g_south here: every frequency the Fourier transform reads is
    # written by the kernel below, the rows beyond the Legendre truncation with zeros

    # Launch the kernel with the specified configuration
    launch!(
        S.architecture,
        ArrayWorkOrder,
        (S.jm_index_size, nlayers),
        inverse_legendre_kernel!,
        g_north,
        g_south,
        specs.data,
        legendre_polynomials_transposed,
        lon_offsets,
        jm_indices,
        nlat_half
    )

    # unscale by cosine of latitude on the fly if requested
    if unscale_coslat
        unscale_coslat!(g_north, g_south, coslat⁻¹, nlayers, architecture = S.architecture)
    end
    return nothing
end


# (forward) legendre transform kernel, called from _legendre!
# Parallelised over the *output* coefficients: one thread per (spherical harmonic lm, layer k),
# each summing over the latitude rings that contribute at its order m. That way every coefficient
# is written by exactly one thread, so no atomic accumulation (and no reset of specs) is needed,
# and neighbouring threads read/write neighbouring coefficients.
@kernel inbounds = true function forward_legendre_kernel!(
        specs_data,                 # output, spherical harmonic coefficients
        legendre_polynomials_data,  # input, Legendre polynomials
        f_north,                    # input, Fourier-transformed northern latitudes
        f_south,                    # input, southern latitudes
        lon_offsets,                # input, longitude offset rotation per (m, j)
        solid_angles,               # input, solid angles for each latitude
        lm_indices,                 # precomputed (m, hemisphere sign, row range) per coefficient
        jm_indices,                 # precomputed (j, m, lm_offset, column length) per row
    )
    lm, k = @index(Global, NTuple)

    m = lm_indices[lm, 1]           # order m of this coefficient
    # north + south for even, north - south for odd degrees l-m, the symmetry split
    hemisphere_sign = convert(real(eltype(f_north)), lm_indices[lm, 2])

    spec = zero(eltype(specs_data))
    # rows of jm_indices holding the rings that contribute at order m, empty if none do
    for r in lm_indices[lm, 3]:lm_indices[lm, 4]
        j = jm_indices[r, 1]

        # SOLID ANGLE QUADRATURE WEIGHT (sinθ Δθ Δϕ) and LONGITUDE OFFSET, the complex
        # conjugate rotating back to the prime meridian
        ΔΩ_rotated = solid_angles[j] * conj(lon_offsets[m, j])

        f = muladd_complex(hemisphere_sign, f_south[m, k, j], f_north[m, k, j])
        spec = muladd_complex(legendre_polynomials_data[lm, j], ΔΩ_rotated * f, spec)
    end

    specs_data[lm, k] = spec
end

function _legendre!(                        # GRID TO SPECTRAL
        specs::LowerTriangularArray,            # Fourier and Legendre-transformed output
        f_north::AbstractArray{<:Complex, 3},   # Fourier-transformed input, northern latitudes
        f_south::AbstractArray{<:Complex, 3},   # and southern latitudes
        scratch_memory::ColumnScratchMemory,    # scratch memory (not used here, but for CPU _legendre!)
        S::SpectralTransform{NF, <:Architectures.GPU},        # precomputed transform
    ) where {NF}
    (; legendre_polynomials) = S            # precomputed Legendre polynomials
    (; lm_indices, jm_indices) = S          # coefficient indices precomputed for threads
    (; solid_angles, lon_offsets) = S

    nlayers = size(specs, 2)                # get number of layers of specs for fewer layers than precomputed in S

    @boundscheck SpeedyTransforms.ismatching(S, specs) || throw(DimensionMismatch(S, specs))
    @boundscheck (size(f_north) == size(f_south) && size(f_north, 1) == S.nfreq_max &&
                  size(f_north, 3) == S.grid.nlat_half && size(f_north, 2) >= nlayers) ||
                 throw(DimensionMismatch(S, specs))
    @boundscheck nlayers <= S.nlayers || throw(DimensionMismatch(S, specs))


    launch!(
        S.architecture,
        ArrayWorkOrder,
        (size(specs.data, 1), nlayers),
        forward_legendre_kernel!,
        specs.data,
        legendre_polynomials.data,
        f_north,
        f_south,
        lon_offsets,
        solid_angles,
        lm_indices,
        jm_indices
    )
    return nothing
end
