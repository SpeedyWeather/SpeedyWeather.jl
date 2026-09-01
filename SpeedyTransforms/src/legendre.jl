"""$(TYPEDSIGNATURES)
Return the Legendre polynomials for latitude ring `j`, as a vector indexed by the running `lm`
index (all orders `m`, all degrees `l`, one ring). For `PrecomputedLegendre` this is a `view` into
the precomputed table; for `RecomputedLegendre` it fills the struct's `ring` scratch vector by
recomputing every column (looping over orders `m`) via `ScaledLegendre.legendre_column!` and
returns it. Meant to be hoisted to the top of the `j` loop in the CPU `_legendre!` methods below,
so the recursion for `RecomputedLegendre` is amortised over every order and every layer of one
ring — one recursion pass per ring per transform, run in `Float64` when `NF == Float64` and in
double-single arithmetic otherwise (see `ScaledLegendre`), so CPU and GPU agree."""
legendre_ring!(legendre::PrecomputedLegendre, j::Integer) = view(legendre.polynomials.data, :, j)

function legendre_ring!(legendre::RecomputedLegendre, j::Integer)
    (; lmax, ring, x, αhi, αlo, βhi, βlo, sectoral_hi, sectoral_lo, sectoral_scale) = legendre
    mmax = size(sectoral_hi, 2)
    lm_offset = 0
    @inbounds for m in 1:mmax
        ncolumn = lmax - m + 1
        out = view(ring, (lm_offset + 1):(lm_offset + ncolumn))
        ScaledLegendre.legendre_column!(
            out, x[j],
            αhi, αlo, βhi, βlo,
            lm_offset, ncolumn,
            sectoral_hi[j, m], sectoral_lo[j, m], sectoral_scale[j, m],
        )
        lm_offset += ncolumn
    end
    return ring
end

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
        odd_k = zero(eltype(north))    # dot prodcut with elements 2, 4, 6, ...

        for l in 1:2:lmax_even          # dot product in pairs for contiguous memory access
            even_k = muladd(specs[l, k], legendre[l], even_k)
            odd_k = muladd(specs[l + 1, k], legendre[l + 1], odd_k)
        end

        # now do the last row if lmax is odd, all written as muladds
        even_k = muladd(specs[end, k], isoddlmax * legendre[end], even_k)
        north[k] = muladd(1, odd_k, even_k)    # north = even + odd
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
        scratch_memory::ColumnScratchMemory,    # scratch memory for vertically batched Legendre transform
        S::SpectralTransform;                   # precomputed transform
        unscale_coslat::Bool = false,           # unscale by cosine of latitude on the fly?
    )
    (; nlat_half) = S.grid                  # dimensions
    (; lmax, mmax) = S.spectrum            # 1-based max degree l, order m of spherical harmonics
    legendre = S.legendre                   # precomputed or recomputed Legendre polynomials
    (; mmax_truncation) = S                 # Legendre shortcut, shortens loop over m, 1-based
    (; coslat⁻¹, lon_offsets) = S
    nlayers = axes(specs, 2)                # get number of layers of specs for fewer layers than precomputed in S

    lmax = lmax - 1                           # 0-based max degree l of spherical harmonics
    mmax = mmax - 1                           # 0-based max order m of spherical harmonics

    @boundscheck ismatching(S, specs) || throw(DimensionMismatch(S, specs))
    # scratch dim 2 is the per-call capacity (= max(planned_K) on CPU, S.nlayers elsewhere);
    # allow it to exceed length(nlayers) so chunked CPU calls pass the bound.
    @boundscheck (
        size(g_north) == size(g_south) && size(g_north, 1) == S.nfreq_max &&
            size(g_north, 3) == nlat_half && size(g_north, 2) >= length(nlayers)
    ) ||
        throw(DimensionMismatch(S, specs))

    north = scratch_memory.north     # use scratch memory for vertically-batched dot product
    south = scratch_memory.south

    return @inbounds for j in 1:nlat_half          # symmetry: loop over northern latitudes only
        g_north[:, nlayers, j] .= 0       # reset scratch memory
        g_south[:, nlayers, j] .= 0       # reset scratch memory

        legendre_j = legendre_ring!(legendre, j)  # precomputed view or recomputed ring, amortised over m, k

        # INVERSE LEGENDRE TRANSFORM by looping over wavenumbers l, m
        lm = 1                              # single running index for non-zero l, m indices
        for m in 1:(mmax_truncation[j] + 1)   # Σ_{m=0}^{mmax}, but 1-based index, shortened to mmax_truncation
            lm_end = lm + lmax - m + 1          # last index in column

            # view on lower triangular column, but batched in vertical
            spec_view = view(specs.data, lm:lm_end, :)
            legendre_view = view(legendre_j, lm:lm_end)

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
            @views g_north[:, nlayers, j] .*= coslat⁻¹[j]        # scale in place
            @views g_south[:, nlayers, j] .*= coslat⁻¹[j]
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

    return @inbounds for k in 1:nlayers
        even_k, odd_k = even[k], odd[k]
        for l in 1:2:lmax_even
            specs[l, k] = muladd(legendre[l], even_k, specs[l, k])
            specs[l + 1, k] = muladd(legendre[l + 1], odd_k, specs[l + 1, k])
        end
        specs[end, k] = muladd(legendre[end], isoddlmax * even_k, specs[end, k])
    end
end

"""$(TYPEDSIGNATURES)
(Forward) Legendre transform, batched in the vertical. Not to be used
directly, but called from transform!."""
function _legendre!(                        # GRID TO SPECTRAL
        specs::LowerTriangularArray,            # Fourier and Legendre-transformed output
        f_north::AbstractArray{<:Complex, 3},   # Fourier-transformed input, northern latitudes
        f_south::AbstractArray{<:Complex, 3},   # and southern latitudes
        scratch_memory::ColumnScratchMemory,    # scratch memory for vertically batched Legendre transform
        S::SpectralTransform;                   # precomputed transform
        add::Bool = false,                      # accumulate onto `specs` instead of overwriting it?
    )
    (; nlat) = S                            # dimensions
    (; nlat_half) = S.grid
    (; lmax, mmax) = S.spectrum             # 1-based max degree l, order m of spherical harmonics
    legendre = S.legendre                   # precomputed or recomputed Legendre polynomials
    (; mmax_truncation) = S                 # Legendre shortcut, shortens loop over m, 1-based
    (; solid_angles, lon_offsets) = S
    nlayers = axes(specs, 2)                # get number of layers of specs for fewer layers than precomputed in S

    lmax = lmax - 1                           # 0-based max degree l of spherical harmonics
    mmax = mmax - 1                           # 0-based max order m of spherical harmonics

    @boundscheck ismatching(S, specs) || throw(DimensionMismatch(S, specs))
    @boundscheck (
        size(f_north) == size(f_south) && size(f_north, 1) == S.nfreq_max &&
            size(f_north, 3) == nlat_half && size(f_north, 2) >= length(nlayers)
    ) ||
        throw(DimensionMismatch(S, specs))

    even = scratch_memory.north      # use scratch memory for outer product
    odd = scratch_memory.south

    # by default reset specs to 0 (this transform accumulates into it internally); `add=true`
    # (used by the Enzyme transform! adjoint rule, matching the inverse `_fourier!` `add` kwarg)
    # accumulates onto the existing contents instead, so no scratch spectral array is needed
    add || fill!(specs, 0)

    return @inbounds for j_north in 1:nlat_half    # symmetry: loop over northern latitudes only
        j = j_north                         # symmetric index / ring-away from pole index

        # SOLID ANGLES including quadrature weights (sinθ Δθ) and azimuth (Δϕ) on ring j
        ΔΩ = solid_angles[j]                # = sinθ Δθ Δϕ, solid angle for a grid point

        legendre_j = legendre_ring!(legendre, j)  # precomputed view or recomputed ring, amortised over m, k

        lm = 1                              # single running index for spherical harmonics
        for m in 1:(mmax_truncation[j] + 1)   # Σ_{m=0}^{mmax}, but 1-based index, shortened to mmax_truncation

            # SOLID ANGLE QUADRATURE WEIGHTS and LONGITUDE OFFSET
            o = lon_offsets[m, j]           # longitude offset rotation by multiplication with complex unit vector
            ΔΩ_rotated = ΔΩ * conj(o)         # complex conjugate for rotation back to prime meridian

            # LEGENDRE TRANSFORM
            for k in nlayers
                fn, fs = f_north[m, k, j], f_south[m, k, j]
                @fastmath even[k] = ΔΩ_rotated * (fn + fs)
                @fastmath odd[k] = ΔΩ_rotated * (fn - fs)
            end

            # integration over l = m:lmax+1
            lm_end = lm + lmax - m + 1                      # last index in column m
            spec_view = view(specs.data, lm:lm_end, :)
            legendre_view = view(legendre_j, lm:lm_end)

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
        nlayers::Integer = size(g_north, 2);     # scale only the layers this transform wrote
        architecture::AbstractArchitecture = DEFAULT_ARCHITECTURE
    )

    launch!(
        architecture, ArrayWorkOrder, (size(g_north, 1), nlayers, size(g_north, 3)),
        unscale_coslat_kernel!, g_north, g_south, coslat⁻¹
    )
    return nothing
end

@kernel inbounds = true function unscale_coslat_kernel!(
        g_north,
        g_south,
        coslat⁻¹,
    )
    i, k, j = @index(Global, NTuple)
    g_north[i, k, j] *= coslat⁻¹[j]
    g_south[i, k, j] *= coslat⁻¹[j]
end
