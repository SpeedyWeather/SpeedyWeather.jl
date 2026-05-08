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
    (; legendre_polynomials) = S            # precomputed Legendre polynomials
    (; mmax_truncation) = S                 # Legendre shortcut, shortens loop over m, 1-based
    (; coslat⁻¹, lon_offsets) = S
    nlayers = axes(specs, 2)                # get number of layers of specs for fewer layers than precomputed in S

    lmax = lmax - 1                           # 0-based max degree l of spherical harmonics
    mmax = mmax - 1                           # 0-based max order m of spherical harmonics

    @boundscheck ismatching(S, specs) || throw(DimensionMismatch(S, specs))
    @boundscheck size(g_north) == size(g_south) == (S.nfreq_max, S.nlayers, nlat_half) || throw(DimensionMismatch(S, specs))

    # Multithreaded over latitudes. Each task gets its own column-scratch slice
    # (view(north_threads, :, c) / view(south_threads, :, c)) so that the inner
    # `_fused_oddeven_matvec!` writes don't race. Latitudes are partitioned into
    # `nchunks` contiguous ranges so that the chunk index `c` is a stable per-task
    # buffer index — avoiding the task-migration pitfalls of `Threads.threadid()`.
    north_threads = scratch_memory.north_threads        # (nlayers, nthreads)
    south_threads = scratch_memory.south_threads        # (nlayers, nthreads)
    nchunks = min(size(north_threads, 2), nlat_half)
    chunks = _chunk_ranges(nlat_half, nchunks)

    @sync for c in 1:nchunks
        Threads.@spawn begin
            north_c = view(north_threads, :, c)
            south_c = view(south_threads, :, c)
            @inbounds for j in chunks[c]    # symmetry: loop over northern latitudes only
                g_north[:, nlayers, j] .= 0
                g_south[:, nlayers, j] .= 0

                # INVERSE LEGENDRE TRANSFORM by looping over wavenumbers l, m
                lm = 1                          # single running index for non-zero l, m indices
                for m in 1:(mmax_truncation[j] + 1) # Σ_{m=0}^{mmax}, 1-based, shortened
                    lm_end = lm + lmax - m + 1  # last index in column

                    # view on lower triangular column, but batched in vertical
                    spec_view = view(specs.data, lm:lm_end, :)
                    legendre_view = view(legendre_polynomials.data, lm:lm_end, j)

                    # dot product but split into even and odd harmonics on the fly
                    _fused_oddeven_matvec!(north_c, south_c, spec_view, legendre_view)

                    # CORRECT FOR LONGITUDE OFFSETTS (if grid points don't start at 0°E)
                    o = lon_offsets[m, j]
                    for k in nlayers
                        g_north[m, k, j] = muladd(o, north_c[k], g_north[m, k, j])
                        g_south[m, k, j] = muladd(o, south_c[k], g_south[m, k, j])
                    end

                    lm = lm_end + 1             # first index of next m column
                end

                if unscale_coslat
                    g_north[:, nlayers, j] .*= coslat⁻¹[j]
                    g_south[:, nlayers, j] .*= coslat⁻¹[j]
                end
            end
        end
    end
    return nothing
end

# Partition 1:n into `nchunks` contiguous UnitRanges. Used to give each task a
# stable index 1..nchunks for indexing per-thread scratch buffers, which avoids
# the task-migration pitfalls of Threads.threadid().
function _chunk_ranges(n::Integer, nchunks::Integer)
    nchunks = max(1, min(nchunks, n))
    base, rem = divrem(n, nchunks)
    ranges = Vector{UnitRange{Int}}(undef, nchunks)
    start = 1
    for c in 1:nchunks
        len = base + (c <= rem ? 1 : 0)
        ranges[c] = start:(start + len - 1)
        start += len
    end
    return ranges
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
        S::SpectralTransform,                   # precomputed transform
    )
    (; nlat) = S                            # dimensions
    (; nlat_half) = S.grid
    (; lmax, mmax) = S.spectrum             # 1-based max degree l, order m of spherical harmonics
    (; legendre_polynomials) = S            # precomputed Legendre polynomials
    (; mmax_truncation) = S                 # Legendre shortcut, shortens loop over m, 1-based
    (; solid_angles, lon_offsets) = S
    nlayers = axes(specs, 2)                # get number of layers of specs for fewer layers than precomputed in S

    lmax = lmax - 1                           # 0-based max degree l of spherical harmonics
    mmax = mmax - 1                           # 0-based max order m of spherical harmonics

    @boundscheck ismatching(S, specs) || throw(DimensionMismatch(S, specs))
    @boundscheck size(f_north) == size(f_south) == (S.nfreq_max, S.nlayers, nlat_half) || throw(DimensionMismatch(S, specs))

    fill!(specs, 0)                         # reset as we accumulate into specs

    # Multithreaded over latitudes. Each task accumulates into its own slice of
    # specs_threads (`view(specs_threads, :, :, c)`) and uses its own column
    # scratch slices from north_threads/south_threads (repurposed here as
    # `even`/`odd`). After all tasks finish we sum each thread-local accumulator
    # into the output `specs.data` to avoid races. As in the inverse transform
    # we partition latitudes into `nchunks` contiguous ranges so the chunk
    # index `c` is a stable per-task buffer index.
    even_threads = scratch_memory.north_threads         # (nlayers, nthreads), repurposed as `even`
    odd_threads = scratch_memory.south_threads          # (nlayers, nthreads), repurposed as `odd`
    specs_threads = scratch_memory.specs_threads        # (nspec, nlayers, nthreads)
    nchunks = min(size(specs_threads, 3), nlat_half)
    chunks = _chunk_ranges(nlat_half, nchunks)
    nlayers_used = length(nlayers)

    @sync for c in 1:nchunks
        Threads.@spawn begin
            even = view(even_threads, :, c)
            odd = view(odd_threads, :, c)
            specs_local = view(specs_threads, :, :, c)
            fill!(specs_local, 0)
            @inbounds for j_north in chunks[c]      # symmetry: only northern latitudes
                j = j_north                         # symmetric / ring-away-from-pole index

                # SOLID ANGLES including quadrature weights (sinθ Δθ) and azimuth (Δϕ)
                ΔΩ = solid_angles[j]

                lm = 1                              # single running index for spherical harmonics
                for m in 1:(mmax_truncation[j] + 1) # Σ_{m=0}^{mmax}, 1-based, shortened
                    o = lon_offsets[m, j]
                    ΔΩ_rotated = ΔΩ * conj(o)

                    for k in nlayers
                        fn, fs = f_north[m, k, j], f_south[m, k, j]
                        @fastmath even[k] = ΔΩ_rotated * (fn + fs)
                        @fastmath odd[k] = ΔΩ_rotated * (fn - fs)
                    end

                    lm_end = lm + lmax - m + 1
                    spec_view = view(specs_local, lm:lm_end, :)
                    legendre_view = view(legendre_polynomials.data, lm:lm_end, j)

                    _fused_oddeven_outer_product_accumulate!(spec_view, legendre_view, even, odd)

                    lm = lm_end + 1
                end
            end
        end
    end

    # REDUCE: sum each thread-local accumulator into the output specs.data.
    # specs.data may be 1D (LowerTriangularMatrix, single layer) or 2D
    # (LowerTriangularArray, multiple layers); restrict the accumulator to the
    # active nlayers and reshape if needed.
    @inbounds for c in 1:nchunks
        src = view(specs_threads, :, 1:nlayers_used, c)
        if ndims(specs.data) == 1
            specs.data .+= vec(src)
        else
            specs.data .+= src
        end
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)
Unscale by cosine of latitude on the fly.
"""
function unscale_coslat!(
        g_north::AbstractArray{<:Complex, 3},
        g_south::AbstractArray{<:Complex, 3},
        coslat⁻¹::AbstractArray{<:Real, 1};
        architecture::AbstractArchitecture = DEFAULT_ARCHITECTURE
    )

    launch!(
        architecture, Array3DWorkOrder, size(g_north), unscale_coslat_kernel!,
        g_north, g_south, coslat⁻¹
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
