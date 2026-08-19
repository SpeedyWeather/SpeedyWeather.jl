module SpeedyTransformsMetalExt

import Metal: Metal, MtlArray
import AbstractFFTs
import LinearAlgebra
using DocStringExtensions

using SpeedyTransforms
using SpeedyTransforms.RingGrids
using SpeedyTransforms.LowerTriangularArrays

import SpeedyTransforms: SpectralTransform, _fourier_batched!
import SpeedyTransforms.RingGrids: AbstractField
import SpeedyWeatherInternals.Architectures: architecture

# =====================================================================================
# METAL COMMAND-BUFFER-BATCHED FOURIER TRANSFORM
#
# The batched Fourier transform on a (reduced) grid applies one FFT per latitude ring
# (≈ 2*nlat_half tiny FFTs). Metal's native FFT (Metal.jl, backed by MPSGraph) no longer
# blocks on a per-call synchronize as of Metal.jl 1.10.1 (JuliaGPU/Metal.jl#872), but each
# call to `_fft!` still builds and commits its OWN dedicated `MPSCommandBuffer`
# (`Metal.MPS.MPSCommandBuffer(Metal.global_queue(Metal.device()))`), and constructing that
# from Metal.jl's task-local *batched* queue forces a full flush of any pending batched
# work first (`Base.cconvert` on `BatchedCommandQueue` calls `flush!`). So a `transform!`
# call that fires ~2*nlat_half tiny per-ring FFTs still pays that command-buffer-creation
# and flush cost once per ring — benchmarking confirms this: bumping to Metal.jl 1.10.2
# alone (which removed the blocking synchronize) did not measurably improve `transform!`
# performance, because the bottleneck was never the wait, it's the per-call submission
# overhead.
#
# This override removes that overhead the direct way: every ring's FFT is `encode!`d onto
# a shared `MPSCommandBuffer`, then `commit!`ed. Gather (grid -> per-ring buffer) and
# scatter (per-ring buffer -> scratch) stay as plain GPU array copies — Metal.jl's
# `BatchedCommandQueue` already auto-batches those (kernel launches and blit ops
# accumulate into one command buffer and commit lazily), so they don't need manual
# batching.
#
# Vertical layers are already batched for free within each ring's single FFT call: the FFT
# plan's transform axis is the longitude dimension only ("region=1"), so MPSGraph treats
# every other axis of the (nlon × nlayers) input, i.e. `nlayers`, as a free batch dimension
# — one graph execution already covers every vertical level for that ring, nothing to gain
# there.
#
# TWO STRATEGIES, both kept and selectable via `S.metal_merge_hemispheres`:
#
#   - `true` (default): north AND south are encoded onto ONE shared command buffer and
#     committed once per `transform!` call ("north-south merge").
#   - `false`: north and south are each encoded onto their OWN command buffer, committed
#     separately ("unmerged") — the original, simpler design.
#
# Dense A/B benchmarking (nlayers=8) showed a genuine crossover, not noise:
#   - trunc=32..197 (step 5, 34 points): merged wins 31/33 points, often by 10-20%.
#   - trunc=200..300 (step 10, 11 points): unmerged wins 10/11 points, several by
#     15-30% (e.g. trunc=230: 55.6ms merged vs 43.4ms unmerged; trunc=300: 71.8ms vs
#     63.7ms).
# Neither dominates everywhere, hence the toggle rather than a single hardcoded choice.
# Default is `true` (merged) since most usage is expected at the lower end of that range.
# A plausible mechanism for the crossover: at low nlat_half, one fewer command-buffer
# flush (merged) wins; at high nlat_half, the larger single command buffer's construction/
# encode cost (or lost host/device overlap between hemispheres — the GPU executing north's
# buffer while the host builds south's gather/encode) starts to dominate instead. Not
# fully root-caused; revisit if profiling ever explains it precisely.
#
# Explicit ring-length fusion (combining same-`nlon` rings, e.g. on HEALPix-type grids,
# into a single larger MPSGraph call) was also tried but reverted: it requires FFT-ing
# views into a larger shared buffer, and Metal.jl 1.10.2's FFT plans silently return zeros
# when applied to a buffer view at a nonzero offset (confirmed in isolation; likely an
# `MPSGraphTensorData` offset-handling bug, worth reporting upstream). It also wouldn't
# have helped SpeedyWeather's default `OctahedralGaussianGrid`, where every ring already
# has a distinct length. Each ring here therefore gets its own, independently-allocated
# (zero-offset) buffer, reused call-to-call via a small per-`SpectralTransform` cache.
# =====================================================================================

"""$(TYPEDSIGNATURES)
Per-`SpectralTransform` cache of reusable, independently-allocated per-ring buffers for
the Metal command-buffer-batched Fourier transform. One real (`nlon_j × nlayers`) and one
complex (`nfreq_j × nlayers`) buffer per ring, separately for north and south (needed for
the north-south merge strategy; the unmerged strategy just uses the north-half
of these sequentially for both passes). The real buffers double as the forward input /
inverse output, the complex buffers as the forward output / inverse input."""
struct MetalFourierCommandBatchCache{RT, CT}
    ring_real_north::Vector{RT}
    ring_real_south::Vector{RT}
    ring_complex_north::Vector{CT}
    ring_complex_south::Vector{CT}
    nlat_half::Int
    nlat::Int
    has_equator::Bool
    j_equator::Int
end

# One cache per SpectralTransform, keyed by object identity (stable for the transform's
# lifetime; saved as a global and not a field of `S` so `S` stays a concrete, stable type)
const FOURIER_CMDBUF_CACHES = IdDict{SpectralTransform, MetalFourierCommandBatchCache}()

function build_cache(S::SpectralTransform)
    NF = eltype(S)
    nlat_half = S.grid.nlat_half
    nlat = S.nlat
    nlayers = S.nlayers
    nlons = S.nlons

    ring_real_north = [Metal.zeros(NF, nlons[j], nlayers) for j in 1:nlat_half]
    ring_real_south = [Metal.zeros(NF, nlons[j], nlayers) for j in 1:nlat_half]
    ring_complex_north = [Metal.zeros(Complex{NF}, nlons[j] ÷ 2 + 1, nlayers) for j in 1:nlat_half]
    ring_complex_south = [Metal.zeros(Complex{NF}, nlons[j] ÷ 2 + 1, nlayers) for j in 1:nlat_half]

    has_equator = isodd(nlat)
    j_equator = (nlat + 1) ÷ 2

    return MetalFourierCommandBatchCache(
        ring_real_north, ring_real_south, ring_complex_north, ring_complex_south,
        nlat_half, nlat, has_equator, j_equator
    )
end

get_cache(S::SpectralTransform) = get!(() -> build_cache(S), FOURIER_CMDBUF_CACHES, S)::MetalFourierCommandBatchCache

"""$(TYPEDSIGNATURES)
Clear all cached Metal command-batched Fourier buffers (frees the associated GPU memory).
Mainly useful for tests/benchmarks."""
clear_metal_fourier_cache!() = (empty!(FOURIER_CMDBUF_CACHES); nothing)

# Encode one ring's FFT graph execution onto a shared, not-yet-committed command buffer.
# Reuses Metal.jl's own FFT graph cache (`_fft_graph_cache`) so the graph itself is built
# once per (ring size, direction) and shared with the plan's normal `mul!`/`*` usage
# elsewhere; only the command-buffer submission is batched here, everything else is
# identical to what `Metal.jl`'s own `_fft!` does internally.
@inline function encode_fft!(cmdbuf, plan, x, y)
    key = Metal.FFTGraphKey(plan)
    cached = Base.@lock Metal._fft_graph_cache_lock get!(() -> Metal.CachedFFTGraph(key), Metal._fft_graph_cache, key)
    feeds = Dict(cached.placeholder => Metal.MPSGraphTensorData(x))
    results = Dict(cached.result => Metal.MPSGraphTensorData(y))
    Metal.MPS.encode!(
        cmdbuf, cached.graph, Metal.NSDictionary(feeds), Metal.NSDictionary(results),
        Metal.ObjectiveC.nil, Metal.MPSGraphs.default_exec_desc()
    )
    return nothing
end

# =====================================================================================
# North-south merge strategy: north and south share a single command buffer.
# =====================================================================================

function forward_loop_merged!(cache::MetalFourierCommandBatchCache, f_north, f_south, field::AbstractField, S::SpectralTransform)
    (; rfft_plans, rings, nlons) = S
    (; ring_real_north, ring_real_south, ring_complex_north, ring_complex_south) = cache
    (; nlat_half, nlat, has_equator, j_equator) = cache
    queue = Metal.global_queue(Metal.device())

    @inbounds for j in 1:nlat_half
        j_south = nlat - j + 1
        copyto!(ring_real_north[j], view(field.data, rings[j], :))
        copyto!(ring_real_south[j], view(field.data, rings[j_south], :))
    end

    cmdbuf = Metal.MPS.MPSCommandBuffer(queue)
    @inbounds for j in 1:nlat_half
        encode_fft!(cmdbuf, rfft_plans[j], ring_real_north[j], ring_complex_north[j])
        # southern equator ring skipped entirely: same physical latitude the northern
        # pass already transformed, transforming it again would double-count it downstream
        (has_equator && j == j_equator) && continue
        encode_fft!(cmdbuf, rfft_plans[j], ring_real_south[j], ring_complex_south[j])
    end
    Metal.commit!(cmdbuf)

    @inbounds for j in 1:nlat_half
        nfreq = nlons[j] ÷ 2 + 1
        copyto!(view(f_north, 1:nfreq, :, j), ring_complex_north[j])
        if has_equator && j == j_equator
            fill!(view(f_south, 1:nfreq, :, j), 0)
        else
            copyto!(view(f_south, 1:nfreq, :, j), ring_complex_south[j])
        end
    end
    return nothing
end

function inverse_loop_merged!(cache::MetalFourierCommandBatchCache, field::AbstractField, g_north, g_south, S::SpectralTransform)
    (; brfft_plans, rings, nlons) = S
    (; ring_real_north, ring_real_south, ring_complex_north, ring_complex_south) = cache
    (; nlat_half, nlat, has_equator, j_equator) = cache
    queue = Metal.global_queue(Metal.device())

    @inbounds for j in 1:nlat_half
        nfreq = nlons[j] ÷ 2 + 1
        copyto!(ring_complex_north[j], view(g_north, 1:nfreq, :, j))
        (has_equator && j == j_equator) && continue   # south equator ring: north already covers it
        copyto!(ring_complex_south[j], view(g_south, 1:nfreq, :, j))
    end

    cmdbuf = Metal.MPS.MPSCommandBuffer(queue)
    @inbounds for j in 1:nlat_half
        encode_fft!(cmdbuf, brfft_plans[j], ring_complex_north[j], ring_real_north[j])
        (has_equator && j == j_equator) && continue
        encode_fft!(cmdbuf, brfft_plans[j], ring_complex_south[j], ring_real_south[j])
    end
    Metal.commit!(cmdbuf)

    @inbounds for j in 1:nlat_half
        copyto!(view(field.data, rings[j], :), ring_real_north[j])
        (has_equator && j == j_equator) && continue   # north already wrote this grid row
        j_south = nlat - j + 1
        copyto!(view(field.data, rings[j_south], :), ring_real_south[j])
    end
    return nothing
end

# =====================================================================================
# Unmerged strategy: north and south each get their own command buffer.
# Reuses the `_north` buffers for both passes sequentially (no need for the `_south`
# buffers here, but the cache always allocates them for the merged strategy above).
# =====================================================================================

function forward_loop_unmerged!(cache::MetalFourierCommandBatchCache, f_north, f_south, field::AbstractField, S::SpectralTransform)
    (; rfft_plans, rings, nlons) = S
    (; ring_real_north, ring_complex_north, nlat_half, nlat, has_equator, j_equator) = cache
    queue = Metal.global_queue(Metal.device())

    # northern rings: always transformed
    @inbounds for j in 1:nlat_half
        copyto!(ring_real_north[j], view(field.data, rings[j], :))
    end
    cmdbuf = Metal.MPS.MPSCommandBuffer(queue)
    @inbounds for j in 1:nlat_half
        encode_fft!(cmdbuf, rfft_plans[j], ring_real_north[j], ring_complex_north[j])
    end
    Metal.commit!(cmdbuf)
    @inbounds for j in 1:nlat_half
        nfreq = nlons[j] ÷ 2 + 1
        copyto!(view(f_north, 1:nfreq, :, j), ring_complex_north[j])
    end

    # southern rings: the equator ring (if any) is skipped entirely and zeroed instead —
    # it's the same physical latitude the northern pass already transformed, so
    # transforming it again here would double-count it downstream
    @inbounds for j in 1:nlat_half
        j_south = nlat - j + 1
        copyto!(ring_real_north[j], view(field.data, rings[j_south], :))
    end
    cmdbuf = Metal.MPS.MPSCommandBuffer(queue)
    @inbounds for j in 1:nlat_half
        (has_equator && j == j_equator) && continue
        encode_fft!(cmdbuf, rfft_plans[j], ring_real_north[j], ring_complex_north[j])
    end
    Metal.commit!(cmdbuf)
    @inbounds for j in 1:nlat_half
        nfreq = nlons[j] ÷ 2 + 1
        if has_equator && j == j_equator
            fill!(view(f_south, 1:nfreq, :, j), 0)
        else
            copyto!(view(f_south, 1:nfreq, :, j), ring_complex_north[j])
        end
    end
    return nothing
end

function inverse_loop_unmerged!(cache::MetalFourierCommandBatchCache, field::AbstractField, g_north, g_south, S::SpectralTransform)
    (; brfft_plans, rings, nlons) = S
    (; ring_real_north, ring_complex_north, nlat_half, nlat, has_equator, j_equator) = cache
    queue = Metal.global_queue(Metal.device())

    # northern rings: always transformed
    @inbounds for j in 1:nlat_half
        nfreq = nlons[j] ÷ 2 + 1
        copyto!(ring_complex_north[j], view(g_north, 1:nfreq, :, j))
    end
    cmdbuf = Metal.MPS.MPSCommandBuffer(queue)
    @inbounds for j in 1:nlat_half
        encode_fft!(cmdbuf, brfft_plans[j], ring_complex_north[j], ring_real_north[j])
    end
    Metal.commit!(cmdbuf)
    @inbounds for j in 1:nlat_half
        copyto!(view(field.data, rings[j], :), ring_real_north[j])
    end

    # southern rings: the equator ring, if any, is skipped entirely (both FFT and
    # scatter) — north already wrote those grid rows
    @inbounds for j in 1:nlat_half
        (has_equator && j == j_equator) && continue
        nfreq = nlons[j] ÷ 2 + 1
        copyto!(ring_complex_north[j], view(g_south, 1:nfreq, :, j))
    end
    cmdbuf = Metal.MPS.MPSCommandBuffer(queue)
    @inbounds for j in 1:nlat_half
        (has_equator && j == j_equator) && continue
        encode_fft!(cmdbuf, brfft_plans[j], ring_complex_north[j], ring_real_north[j])
    end
    Metal.commit!(cmdbuf)
    @inbounds for j in 1:nlat_half
        (has_equator && j == j_equator) && continue
        j_south = nlat - j + 1
        copyto!(view(field.data, rings[j_south], :), ring_real_north[j])
    end
    return nothing
end

# =====================================================================================
# Method overrides: dispatch on MtlArray scratch (more specific than the generic
# AbstractArray{<:Complex,3} methods in fourier.jl), then dispatch on
# `S.metal_merge_hemispheres` to pick the strategy.
# =====================================================================================

"""$(TYPEDSIGNATURES)
Metal command-buffer-batched forward (grid → spectral) batched Fourier transform. Picks
between the merged (one commit for north+south) and unmerged (one commit each) strategies
based on `S.metal_merge_hemispheres`; see the module-level comment for the benchmarked
crossover between them."""
function SpeedyTransforms._fourier_batched!(
        f_north::MtlArray{<:Complex, 3},
        f_south::MtlArray{<:Complex, 3},
        field::AbstractField,
        S::SpectralTransform,
    )
    @assert eltype(field) == eltype(S) "Number format of grid $(eltype(field)) and SpectralTransform $(eltype(S)) need to match."
    if !S.cuda_graphs   # generic toggle for "use accelerated batched Fourier path"; also gates Metal batching
        return Base.@invoke _fourier_batched!(
            f_north::AbstractArray{<:Complex, 3}, f_south::AbstractArray{<:Complex, 3},
            field::AbstractField, S::SpectralTransform,
        )
    end
    cache = get_cache(S)
    if S.metal_merge_hemispheres
        forward_loop_merged!(cache, f_north, f_south, field, S)
    else
        forward_loop_unmerged!(cache, f_north, f_south, field, S)
    end
    return nothing
end

"""$(TYPEDSIGNATURES)
Metal command-buffer-batched inverse (spectral → grid) batched Fourier transform. Picks
between the merged and unmerged strategies based on `S.metal_merge_hemispheres`; see the
module-level comment for the benchmarked crossover between them."""
function SpeedyTransforms._fourier_batched!(
        field::AbstractField,
        g_north::MtlArray{<:Complex, 3},
        g_south::MtlArray{<:Complex, 3},
        S::SpectralTransform,
    )
    if !S.cuda_graphs
        return Base.@invoke _fourier_batched!(
            field::AbstractField, g_north::AbstractArray{<:Complex, 3},
            g_south::AbstractArray{<:Complex, 3}, S::SpectralTransform,
        )
    end
    cache = get_cache(S)
    if S.metal_merge_hemispheres
        inverse_loop_merged!(cache, field, g_north, g_south, S)
    else
        inverse_loop_unmerged!(cache, field, g_north, g_south, S)
    end
    return nothing
end

end
