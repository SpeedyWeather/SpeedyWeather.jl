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

    fill!(g_north, 0)
    fill!(g_south, 0)

    north_threads = scratch_memory.north_threads        # (nlayers, nthreads)
    south_threads = scratch_memory.south_threads        # (nlayers, nthreads)

    if Threads.nthreads() > 1
        # Multithreaded over latitudes using @threads :static, which pins each
        # iteration to a fixed OS thread so Threads.threadid() is stable and
        # safe to use as a buffer index into north_threads/south_threads.
        Threads.@threads :static for j in 1:nlat_half
            t = Threads.threadid()
            _legendre_inverse_lat!(
                g_north, g_south, specs,
                view(north_threads, :, t), view(south_threads, :, t),
                legendre_polynomials, mmax_truncation, lon_offsets,
                lmax, nlayers, j,
            )
        end
    else
        # Skip the @threads runtime entirely on a single thread (it allocates
        # the threading_run machinery on every call and dominates the runtime
        # for cheap iterations).
        north = view(north_threads, :, 1)
        south = view(south_threads, :, 1)
        for j in 1:nlat_half
            _legendre_inverse_lat!(
                g_north, g_south, specs, north, south,
                legendre_polynomials, mmax_truncation, lon_offsets,
                lmax, nlayers, j,
            )
        end
    end

    if unscale_coslat
        unscale_coslat!(g_north, g_south, coslat⁻¹; architecture = architecture(specs))
    end
    return nothing
end

# Per-latitude body of the inverse Legendre transform, shared between the
# threaded and serial paths.
@inline function _legendre_inverse_lat!(
        g_north, g_south, specs, north, south,
        legendre_polynomials, mmax_truncation, lon_offsets,
        lmax, nlayers, j::Integer,
    )
    @inbounds begin
        lm = 1
        for m in 1:(mmax_truncation[j] + 1)
            lm_end = lm + lmax - m + 1

            spec_view = view(specs.data, lm:lm_end, :)
            legendre_view = view(legendre_polynomials.data, lm:lm_end, j)

            _fused_oddeven_matvec!(north, south, spec_view, legendre_view)

            o = lon_offsets[m, j]
            for k in nlayers
                g_north[m, k, j] = muladd(o, north[k], g_north[m, k, j])
                g_south[m, k, j] = muladd(o, south[k], g_south[m, k, j])
            end

            lm = lm_end + 1
        end
    end
    return nothing
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

    even_threads = scratch_memory.north_threads         # (nlayers, nthreads), repurposed as `even`
    odd_threads = scratch_memory.south_threads          # (nlayers, nthreads), repurposed as `odd`

    if Threads.nthreads() > 1
        # Multithreaded: each thread accumulates into its own specs_threads
        # slice; we then reduce by summing into specs.data. specs_threads
        # views avoid races on the shared output.
        specs_threads = scratch_memory.specs_threads    # (nspec, nlayers, nthreads)
        fill!(specs_threads, 0)
        nlayers_used = length(nlayers)

        Threads.@threads :static for j in 1:nlat_half
            t = Threads.threadid()
            _legendre_forward_lat!(
                view(specs_threads, :, :, t), f_north, f_south,
                view(even_threads, :, t), view(odd_threads, :, t),
                legendre_polynomials, mmax_truncation, lon_offsets, solid_angles,
                lmax, nlayers, j,
            )
        end

        # REDUCE: sum each thread-local accumulator into specs.data.
        # specs.data may be 1D (LowerTriangularMatrix, single layer) or 2D.
        nthreads = size(specs_threads, 3)
        @inbounds for t in 1:nthreads
            src = view(specs_threads, :, 1:nlayers_used, t)
            if ndims(specs.data) == 1
                specs.data .+= vec(src)
            else
                specs.data .+= src
            end
        end
    else
        # Single-thread fast path: skip the @threads runtime *and* the
        # per-thread accumulator + reduction. Write straight into specs.data,
        # which matches what the old (pre-threading) implementation did.
        even = view(even_threads, :, 1)
        odd = view(odd_threads, :, 1)
        specs_data = ndims(specs.data) == 1 ?
            reshape(specs.data, length(specs.data), 1) : specs.data
        for j in 1:nlat_half
            _legendre_forward_lat!(
                specs_data, f_north, f_south, even, odd,
                legendre_polynomials, mmax_truncation, lon_offsets, solid_angles,
                lmax, nlayers, j,
            )
        end
    end
    return nothing
end

# Per-latitude body of the forward Legendre transform, accumulating into
# `specs_local`. Shared between threaded and serial paths.
@inline function _legendre_forward_lat!(
        specs_local, f_north, f_south, even, odd,
        legendre_polynomials, mmax_truncation, lon_offsets, solid_angles,
        lmax, nlayers, j::Integer,
    )
    @inbounds begin
        ΔΩ = solid_angles[j]

        lm = 1
        for m in 1:(mmax_truncation[j] + 1)
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
