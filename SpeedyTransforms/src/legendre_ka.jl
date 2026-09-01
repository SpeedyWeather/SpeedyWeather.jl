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


# (inverse) legendre transform kernel for the RECOMPUTED polynomials, called from _legendre!
# Same thread layout as `inverse_legendre_kernel!` above — one thread per (ring j, order m, layer
# k), taken from `jm_indices` — but instead of loading the precomputed polynomials it runs the
# ScaledLegendre column recursion in registers as it walks down the column, so nothing is read from
# global memory per harmonic except the α/β recursion coefficients. Those are indexed by lm alone,
# which is identical for every thread of a warp (a warp holds neighbouring rings at the *same*
# order), so they broadcast out of cache instead of being fetched per thread.
@kernel inbounds = true function inverse_legendre_recompute_kernel!(
        g_north,                        # Scratch storage for legendre coefficients
        g_south,                        # before fft
        specs_data,                     # Data passed from spectral grid
        αhi, αlo, βhi, βlo,             # input, recursion coefficients (double-single), by lm
        sectoral_hi, sectoral_lo,       # input, starting values λ_m^m (double-single), (j, m)
        sectoral_scale,                 # input, their extended-exponent scale, (j, m)
        xhi, xlo,                       # input, cos(colat) (double-single), by j
        lon_offsets,                    # input, longitude offset rotation per (m, j)
        jm_indices,                     # precomputed (j, m, lm_offset, column length) per row
    )
    i, k = @index(Global, NTuple)

    j = jm_indices[i, 1]
    m = jm_indices[i, 2]
    lm_offset = jm_indices[i, 3]        # offset to index the lower triangular column directly
    ncolumn = jm_indices[i, 4]          # number of degrees at order m, lmax-m

    NF = eltype(αhi)
    # "even" and "odd" coined with 0-based indexing, i.e. the even l=0 mode is 1st element
    even_k = zero(eltype(g_north))      # dot product with elements 1, 3, 5, ...
    odd_k = zero(eltype(g_south))       # dot product with elements 2, 4, 6, ...

    # ncolumn == 0 marks the rows beyond the Legendre truncation, which only write zeros; skipping
    # here also keeps the sectoral lookup below in bounds, as those rows carry an order m that can
    # exceed mmax (they exist to cover the Fourier frequencies, see `jm_indices`' construction)
    if ncolumn > 0
        xh, xl = xhi[j], xlo[j]
        p1h, p1l = sectoral_hi[j, m], sectoral_lo[j, m]  # λ_{l-1}^m, starting at the sectoral mode
        p2h, p2l = zero(NF), zero(NF)                    # λ_{l-2}^m ≡ 0 above the sectoral mode
        s = Int32(sectoral_scale[j, m])

        # PHASE 1: climb out of the underflow region (λ_m^m ~ sinθ^m is far below floatmin near the
        # poles at high m). Every value here is zero to working precision by the `scale_shift`
        # invariant, so nothing is accumulated and the recursion just steps. If the column ends
        # before s reaches 0 the whole column is negligible and this thread writes zeros, which is
        # also what the precomputed path produces there.
        q = 1
        while q <= ncolumn && s > 0
            p1h, p1l, p2h, p2l, s = ScaledLegendre.legendre_advance(
                αhi[lm_offset + q + 1], αlo[lm_offset + q + 1],
                βhi[lm_offset + q + 1], βlo[lm_offset + q + 1],
                xh, xl, p1h, p1l, p2h, p2l, s,
            )
            q += 1
        end

        # PHASE 2: the hot loop, two recursion steps per iteration so the even/odd (north/south
        # symmetry) split of the dot product survives without a branch or a parity test inside.
        # `a` collects positions q_first, q_first+2, ..., `b` the ones in between; which of the two
        # is the "even" sum depends only on the parity of where phase 1 left off.
        q_first = q
        a = zero(eltype(g_north))
        b = zero(eltype(g_south))
        while q < ncolumn                   # a full pair still fits
            a = muladd_complex(ScaledLegendre.legendre_value(p1h, p1l), specs_data[lm_offset + q, k], a)
            p1h, p1l, p2h, p2l, s = ScaledLegendre.legendre_advance(
                αhi[lm_offset + q + 1], αlo[lm_offset + q + 1],
                βhi[lm_offset + q + 1], βlo[lm_offset + q + 1],
                xh, xl, p1h, p1l, p2h, p2l, s,
            )
            b = muladd_complex(ScaledLegendre.legendre_value(p1h, p1l), specs_data[lm_offset + q + 1, k], b)
            p1h, p1l, p2h, p2l, s = ScaledLegendre.legendre_advance(
                αhi[lm_offset + q + 2], αlo[lm_offset + q + 2],
                βhi[lm_offset + q + 2], βlo[lm_offset + q + 2],
                xh, xl, p1h, p1l, p2h, p2l, s,
            )
            q += 2
        end

        if q == ncolumn                     # odd number of remaining degrees, do the last one
            a = muladd_complex(ScaledLegendre.legendre_value(p1h, p1l), specs_data[lm_offset + ncolumn, k], a)
        end

        # 1-based positions: an odd position is an even degree l-m, hence the parity test
        even_k, odd_k = isodd(q_first) ? (a, b) : (b, a)
    end

    # CORRECT FOR LONGITUDE OFFSETTS (if grid points don't start at 0°E), see the kernel above
    o = ncolumn > 0 ? lon_offsets[m, j] : zero(eltype(lon_offsets))

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
    legendre = S.legendre               # precomputed or recomputed Legendre polynomials
    (; jm_indices) = S                  # (j, m) loop indices precomputed for threads
    (; coslat⁻¹, lon_offsets) = S
    nlayers = size(specs, 2)            # get number of layers of specs for fewer layers than precomputed in S

    @boundscheck SpeedyTransforms.ismatching(S, specs) || throw(DimensionMismatch(S, specs))
    # Scratch dim 2 is the per-call capacity (= max(planned_K) on CPU, S.nlayers elsewhere);
    # allow it to exceed nlayers so chunked/sliced callers (e.g. test scratches copied
    # from CPU) still pass.
    @boundscheck (
        size(g_north) == size(g_south) && size(g_north, 1) == S.nfreq_max &&
            size(g_north, 3) == nlat_half && size(g_north, 2) >= nlayers
    ) ||
        throw(DimensionMismatch(S, specs))
    @boundscheck nlayers <= S.nlayers || throw(DimensionMismatch(S, specs))

    # no need to reset g_north/g_south here: every frequency the Fourier transform reads is
    # written by the kernel below, the rows beyond the Legendre truncation with zeros

    # Launch the kernel with the specified configuration. The two modes differ only in where the
    # polynomials come from — a global-memory load from the transposed table, or the recursion run
    # in registers — so the thread layout and the work size are identical.
    if legendre isa RecomputedLegendre
        launch!(
            S.architecture,
            ArrayWorkOrder,
            (S.jm_index_size, nlayers),
            inverse_legendre_recompute_kernel!,
            g_north,
            g_south,
            specs.data,
            legendre.αhi, legendre.αlo, legendre.βhi, legendre.βlo,
            legendre.sectoral_hi, legendre.sectoral_lo, legendre.sectoral_scale,
            legendre.xhi, legendre.xlo,
            lon_offsets,
            jm_indices
        )
    else
        launch!(
            S.architecture,
            ArrayWorkOrder,
            (S.jm_index_size, nlayers),
            inverse_legendre_kernel!,
            g_north,
            g_south,
            specs.data,
            legendre.polynomials_transposed,   # precomputed Legendre polynomials, (j, lm)
            lon_offsets,
            jm_indices,
            nlat_half
        )
    end

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
# `lm_offset` shifts the thread index onto a contiguous *block* of coefficients: it is 0 for the
# precomputed path (one launch covering every coefficient, `legendre_polynomials_data` the whole
# table) and the base of the current tile for the recomputed path, where the launch is repeated per
# tile and `legendre_polynomials_data` holds only that tile's rows, see the recompute kernel below.
@kernel inbounds = true function forward_legendre_kernel!(
        specs_data,                 # output, spherical harmonic coefficients
        legendre_polynomials_data,  # input, Legendre polynomials, (lm - lm_offset, j)
        f_north,                    # input, Fourier-transformed northern latitudes
        f_south,                    # input, southern latitudes
        lon_offsets,                # input, longitude offset rotation per (m, j)
        solid_angles,               # input, solid angles for each latitude
        lm_indices,                 # precomputed (m, hemisphere sign, row range) per coefficient
        jm_indices,                 # precomputed (j, m, lm_offset, column length) per row
        lm_offset,                  # index of the coefficient just before this launch's block
    )
    lm_local, k = @index(Global, NTuple)
    lm = lm_local + lm_offset       # global coefficient index

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
        spec = muladd_complex(legendre_polynomials_data[lm_local, j], ΔΩ_rotated * f, spec)
    end

    specs_data[lm, k] = spec
end

# Recompute one tile of the Legendre polynomials for the forward transform, called from _legendre!
# One thread per (ring j, order m) walking the whole column at that order down into `tile`, which
# holds a contiguous block of orders in the same (lm, j) layout the precomputed table uses — so
# `forward_legendre_kernel!` above reads it with exactly the same indexing and the same coalescing.
# The ring j is the *fast* thread index so that neighbouring threads of a warp share the order m
# and therefore the α/β lookups (broadcast from cache) and the column length (no divergence); the
# writes down a column are strided, which is the right trade because the tile is written once per
# transform but read once per layer.
@kernel inbounds = true function recompute_legendre_tile_kernel!(
        tile,                       # output, one tile of Legendre polynomials, (lm - lm_base, j)
        αhi, αlo, βhi, βlo,         # input, recursion coefficients (double-single), by lm
        sectoral_hi, sectoral_lo,   # input, starting values λ_m^m (double-single), (j, m)
        sectoral_scale,             # input, their extended-exponent scale, (j, m)
        xhi, xlo,                   # input, cos(colat) (double-single), by j
        lmax,                       # 1-based number of degrees, to get each column's length
        m_first,                    # first order m of this tile
        lm_base,                    # index of the coefficient just before this tile's first column
    )
    j, mi = @index(Global, NTuple)
    m = m_first + mi - 1
    ncolumn = lmax - m + 1                              # degrees at order m

    # index just before this column's first entry, LowerTriangularArrays.lm2i(m, m, lmax) - 1
    lm_column = m + (m - 1) * lmax - m * (m - 1) ÷ 2 - 1
    row = lm_column - lm_base                           # same, but relative to this tile

    NF = eltype(tile)
    xh, xl = xhi[j], xlo[j]
    p1h, p1l = sectoral_hi[j, m], sectoral_lo[j, m]     # λ_{l-1}^m, starting at the sectoral mode
    p2h, p2l = zero(NF), zero(NF)                       # λ_{l-2}^m ≡ 0 above the sectoral mode
    s = Int32(sectoral_scale[j, m])

    # identical to ScaledLegendre.legendre_column!, written out here because the output is a
    # strided column of a matrix rather than a contiguous vector
    for q in 1:ncolumn
        tile[row + q, j] = s == 0 ? ScaledLegendre.legendre_value(p1h, p1l) : zero(NF)
        i = lm_column + q + 1
        p1h, p1l, p2h, p2l, s = ScaledLegendre.legendre_advance(
            αhi[i], αlo[i], βhi[i], βlo[i], xh, xl, p1h, p1l, p2h, p2l, s,
        )
    end
end

function _legendre!(                        # GRID TO SPECTRAL
        specs::LowerTriangularArray,            # Fourier and Legendre-transformed output
        f_north::AbstractArray{<:Complex, 3},   # Fourier-transformed input, northern latitudes
        f_south::AbstractArray{<:Complex, 3},   # and southern latitudes
        scratch_memory::ColumnScratchMemory,    # scratch memory (not used here, but for CPU _legendre!)
        S::SpectralTransform{NF, <:Architectures.GPU},        # precomputed transform
    ) where {NF}
    legendre = S.legendre                   # precomputed or recomputed Legendre polynomials
    (; lm_indices, jm_indices) = S          # coefficient indices precomputed for threads
    (; solid_angles, lon_offsets) = S
    (; lmax) = S.spectrum                   # 1-based max degree l
    (; nlat_half) = S.grid

    nlayers = size(specs, 2)                # get number of layers of specs for fewer layers than precomputed in S

    @boundscheck SpeedyTransforms.ismatching(S, specs) || throw(DimensionMismatch(S, specs))
    @boundscheck (
        size(f_north) == size(f_south) && size(f_north, 1) == S.nfreq_max &&
            size(f_north, 3) == S.grid.nlat_half && size(f_north, 2) >= nlayers
    ) ||
        throw(DimensionMismatch(S, specs))
    @boundscheck nlayers <= S.nlayers || throw(DimensionMismatch(S, specs))


    if legendre isa RecomputedLegendre
        # The forward transform is parallelised over its output coefficients, summing over rings,
        # which runs orthogonally to the recursion (down a column at fixed order) — so unlike the
        # inverse it cannot recompute in registers. Instead recompute into a `tile` buffer, one
        # contiguous block of orders m at a time. A block of orders is a contiguous range of lm in
        # the lower-triangular layout, so every coefficient belongs to exactly one tile and is
        # written exactly once: no accumulation across tiles and the result does not depend on the
        # tile size. The recompute cost is per transform call, not per layer, so it amortises over
        # nlayers. See `tile_order_blocks` for how the blocks are sized.
        (; tile, tile_orders) = legendre
        for block in tile_orders
            m_first, m_last = first(block), last(block)
            lm_base = first(LowerTriangularArrays.get_lm_range(m_first, lmax - 1)) - 1
            n_lm = last(LowerTriangularArrays.get_lm_range(m_last, lmax - 1)) - lm_base

            launch!(
                S.architecture,
                ArrayWorkOrder,
                (nlat_half, length(block)),
                recompute_legendre_tile_kernel!,
                tile,
                legendre.αhi, legendre.αlo, legendre.βhi, legendre.βlo,
                legendre.sectoral_hi, legendre.sectoral_lo, legendre.sectoral_scale,
                legendre.xhi, legendre.xlo,
                lmax,
                m_first,
                lm_base
            )

            launch!(
                S.architecture,
                ArrayWorkOrder,
                (n_lm, nlayers),
                forward_legendre_kernel!,
                specs.data,
                tile,
                f_north,
                f_south,
                lon_offsets,
                solid_angles,
                lm_indices,
                jm_indices,
                lm_base
            )
        end
    else
        launch!(
            S.architecture,
            ArrayWorkOrder,
            (size(specs.data, 1), nlayers),
            forward_legendre_kernel!,
            specs.data,
            legendre.polynomials.data,      # precomputed Legendre polynomials, (lm, j)
            f_north,
            f_south,
            lon_offsets,
            solid_angles,
            lm_indices,
            jm_indices,
            0                               # no offset: one launch covers every coefficient
        )
    end
    return nothing
end
