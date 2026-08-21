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
# There are two strategies kept for now, due to differences in performance at different resolution ranges.
#
#   - `true` (default): north AND south are encoded onto ONE shared command buffer and
#     committed once per `transform!` call ("north-south merge").
#   - `false`: north and south are each encoded onto their OWN command buffer, committed
#     separately ("unmerged") — the original, simpler design.
# =====================================================================================

"""$(TYPEDSIGNATURES)
Per-`(SpectralTransform, nlayers)` cache of reusable, independently-allocated per-ring
buffers for the Metal command-buffer-batched Fourier transform. One real (`nlon_j ×
nlayers`) and one complex (`nfreq_j × nlayers`) buffer per ring, separately for north and
south (needed for the north-south merge strategy; the unmerged strategy just uses the
north-half of these sequentially for both passes). The real buffers double as the forward
input / inverse output, the complex buffers as the forward output / inverse input.
Buffers are sized to the actual per-call `nlayers` (the batch `K`), not `S.nlayers` (the
scratch-memory upper bound) — a `SpectralTransform` built via `SpectralGrid` plans several
distinct `K`s (e.g. `1, 2, nlayers, 2nlayers, ...`) and each needs its own correctly-sized
buffers to match its own `rfft_plans_batched[K]`/`brfft_plans_batched[K]` plans."""
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

# One cache per (SpectralTransform, nlayers), the outer IdDict keyed by object identity
# (stable for the transform's lifetime; saved as a global and not a field of `S` so `S`
# stays a concrete, stable type), the inner Dict keyed by the batch size `K` actually used
# — a single `S` serves multiple `K`s (see `MetalFourierCommandBatchCache` docstring).
const FOURIER_CMDBUF_CACHES = IdDict{SpectralTransform, Dict{Int, MetalFourierCommandBatchCache}}()

function build_cache(S::SpectralTransform, nlayers::Int)
    NF = eltype(S)
    nlat_half = S.grid.nlat_half
    nlat = S.nlat
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

function get_cache(S::SpectralTransform, nlayers::Int)
    caches_by_K = get!(() -> Dict{Int, MetalFourierCommandBatchCache}(), FOURIER_CMDBUF_CACHES, S)
    return get!(() -> build_cache(S, nlayers), caches_by_K, nlayers)::MetalFourierCommandBatchCache
end

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
    (; rings, nlons) = S
    rfft_plans = S.rfft_plans_batched[size(field, 2)]
    (; ring_real_north, ring_real_south, ring_complex_north, ring_complex_south) = cache
    (; nlat_half, nlat, has_equator, j_equator) = cache
    queue = Metal.global_queue(Metal.device())

    # south hemisphere gathered in forward (ascending) physical row order through
    # field.data — j_south increases 1:1 with the loop, unlike `nlat - j + 1`, since GPU
    # global-memory traffic benefits from ascending access even though each individual
    # `rings[j]` is already a contiguous block regardless of visitation order
    @inbounds for j in 1:nlat_half
        copyto!(ring_real_north[j], view(field.data, rings[j], :))
    end
    @inbounds for j_south in (nlat_half + 1):nlat
        j = nlat - j_south + 1
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

function inverse_loop_merged!(cache::MetalFourierCommandBatchCache, field::AbstractField, g_north, g_south, S::SpectralTransform; add::Bool = false)
    (; rings, nlons) = S
    nlayers = size(field, 2)
    brfft_plans = S.brfft_plans_batched[nlayers]
    (; ring_real_north, ring_real_south, ring_complex_north, ring_complex_south) = cache
    (; nlat_half, nlat, has_equator, j_equator) = cache
    queue = Metal.global_queue(Metal.device())

    # g_north/g_south are scratch, oversized to the transform's max batch width — slice to
    # the actual per-call `nlayers`, not `:`
    @inbounds for j in 1:nlat_half
        nfreq = nlons[j] ÷ 2 + 1
        copyto!(ring_complex_north[j], view(g_north, 1:nfreq, 1:nlayers, j))
        (has_equator && j == j_equator) && continue   # south equator ring: north already covers it
        copyto!(ring_complex_south[j], view(g_south, 1:nfreq, 1:nlayers, j))
    end

    cmdbuf = Metal.MPS.MPSCommandBuffer(queue)
    @inbounds for j in 1:nlat_half
        encode_fft!(cmdbuf, brfft_plans[j], ring_complex_north[j], ring_real_north[j])
        (has_equator && j == j_equator) && continue
        encode_fft!(cmdbuf, brfft_plans[j], ring_complex_south[j], ring_real_south[j])
    end
    Metal.commit!(cmdbuf)

    @inbounds for j in 1:nlat_half
        dest = view(field.data, rings[j], :)
        add ? (dest .+= ring_real_north[j]) : copyto!(dest, ring_real_north[j])
    end
    # south hemisphere scattered in forward (ascending) physical row order, mirroring the
    # gather above
    @inbounds for j_south in (nlat_half + 1):nlat
        j = nlat - j_south + 1
        (has_equator && j == j_equator) && continue   # north already wrote this grid row
        dest = view(field.data, rings[j_south], :)
        add ? (dest .+= ring_real_south[j]) : copyto!(dest, ring_real_south[j])
    end
    return nothing
end

# =====================================================================================
# Unmerged strategy: north and south each get their own command buffer.
# Reuses the `_north` buffers for both passes sequentially (no need for the `_south`
# buffers here, but the cache always allocates them for the merged strategy above).
# =====================================================================================

function forward_loop_unmerged!(cache::MetalFourierCommandBatchCache, f_north, f_south, field::AbstractField, S::SpectralTransform)
    (; rings, nlons) = S
    rfft_plans = S.rfft_plans_batched[size(field, 2)]
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
    # transforming it again here would double-count it downstream.
    # Gathered in forward (ascending) physical row order through field.data, mirroring
    # the north loop above (see forward_loop_merged! for the same reasoning)
    @inbounds for j_south in (nlat_half + 1):nlat
        j = nlat - j_south + 1
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

function inverse_loop_unmerged!(cache::MetalFourierCommandBatchCache, field::AbstractField, g_north, g_south, S::SpectralTransform; add::Bool = false)
    (; rings, nlons) = S
    nlayers = size(field, 2)
    brfft_plans = S.brfft_plans_batched[nlayers]
    (; ring_real_north, ring_complex_north, nlat_half, nlat, has_equator, j_equator) = cache
    queue = Metal.global_queue(Metal.device())

    # northern rings: always transformed. g_north/g_south are scratch, oversized to the
    # transform's max batch width — slice to the actual per-call `nlayers`, not `:`
    @inbounds for j in 1:nlat_half
        nfreq = nlons[j] ÷ 2 + 1
        copyto!(ring_complex_north[j], view(g_north, 1:nfreq, 1:nlayers, j))
    end
    cmdbuf = Metal.MPS.MPSCommandBuffer(queue)
    @inbounds for j in 1:nlat_half
        encode_fft!(cmdbuf, brfft_plans[j], ring_complex_north[j], ring_real_north[j])
    end
    Metal.commit!(cmdbuf)
    @inbounds for j in 1:nlat_half
        dest = view(field.data, rings[j], :)
        add ? (dest .+= ring_real_north[j]) : copyto!(dest, ring_real_north[j])
    end

    # southern rings: the equator ring, if any, is skipped entirely (both FFT and
    # scatter) — north already wrote those grid rows
    @inbounds for j in 1:nlat_half
        (has_equator && j == j_equator) && continue
        nfreq = nlons[j] ÷ 2 + 1
        copyto!(ring_complex_north[j], view(g_south, 1:nfreq, 1:nlayers, j))
    end
    cmdbuf = Metal.MPS.MPSCommandBuffer(queue)
    @inbounds for j in 1:nlat_half
        (has_equator && j == j_equator) && continue
        encode_fft!(cmdbuf, brfft_plans[j], ring_complex_north[j], ring_real_north[j])
    end
    Metal.commit!(cmdbuf)
    # scattered in forward (ascending) physical row order through field.data, mirroring
    # the gather loops above
    @inbounds for j_south in (nlat_half + 1):nlat
        j = nlat - j_south + 1
        (has_equator && j == j_equator) && continue
        dest = view(field.data, rings[j_south], :)
        add ? (dest .+= ring_real_north[j]) : copyto!(dest, ring_real_north[j])
    end
    return nothing
end

# =====================================================================================
# Fused single-graph strategy: every ring's FFT (both hemispheres) is a separate,
# independent op inside ONE shared MPSGraph, built once and cached. A `transform!` call
# then issues exactly ONE `encode!` dispatch per direction (covering every ring at once)
# instead of one per ring — this is what actually eliminates the fixed ~0.07-0.08ms
# per-`encode!` Objective-C dispatch overhead identified as the dominant cost in
# `metal_profiling.md`, rather than just batching dispatches onto a shared command buffer
# (which the merged/unmerged strategies above already do, but still issue nlat_half or
# 2*nlat_half separate `encode!` calls). Reuses the same per-ring gather/scatter `MtlArray`
# buffers as the merged/unmerged strategies (`MetalFourierCommandBatchCache`, via
# `get_cache`) — only the graph construction and dispatch differ.
# =====================================================================================

"""$(TYPEDSIGNATURES)
Per-`(SpectralTransform, nlayers)` cache holding one shared `MPSGraph` for the forward
(real→complex) direction and one for the inverse (complex→real) direction, each containing
one independent FFT op per ring (north `1:nlat_half`, south `1:nlat_half` minus the
equator-duplicate ring). Built once; every `transform!` call reuses the same graphs and
just re-binds fresh `MPSGraphTensorData` (wrapping the current per-ring buffers) to the
cached placeholder/result tensors."""
struct MetalFusedFourierGraphCache
    forward_graph::Metal.MPSGraph
    forward_ph_north::Vector{Metal.MPSGraphTensor}
    forward_res_north::Vector{Metal.MPSGraphTensor}
    forward_ph_south::Vector{Metal.MPSGraphTensor}
    forward_res_south::Vector{Metal.MPSGraphTensor}
    inverse_graph::Metal.MPSGraph
    inverse_ph_north::Vector{Metal.MPSGraphTensor}
    inverse_res_north::Vector{Metal.MPSGraphTensor}
    inverse_ph_south::Vector{Metal.MPSGraphTensor}
    inverse_res_south::Vector{Metal.MPSGraphTensor}
    nlat_half::Int
    nlat::Int
    has_equator::Bool
    j_equator::Int
end

const FUSED_FOURIER_GRAPH_CACHES = IdDict{SpectralTransform, Dict{Int, MetalFusedFourierGraphCache}}()

# axis for a (nlon_j, nlayers) or (nfreq_j, nlayers) ring array with FFT along dim 1:
# Julia axis 1 of a 2-D array -> Metal axis (ndims - axis) = 2 - 1 = 1 (see the identical
# comment/derivation in Metal.jl's own `CachedFFTGraph`, fft.jl:246-249)
fft_axis() = Metal.NSArray([Metal.NSNumber(1)])

function build_fused_forward_graph(NF, nlayers::Int, nlat_half::Int, has_equator::Bool, j_equator::Int, nlons::Vector{Int})
    graph = Metal.MPSGraph()
    fft_desc = Metal.MPSGraphFFTDescriptor(; inverse = false)
    axes = fft_axis()
    ph_north = Vector{Metal.MPSGraphTensor}(undef, nlat_half)
    res_north = Vector{Metal.MPSGraphTensor}(undef, nlat_half)
    ph_south = Vector{Metal.MPSGraphTensor}(undef, nlat_half)
    res_south = Vector{Metal.MPSGraphTensor}(undef, nlat_half)
    @inbounds for j in 1:nlat_half
        nlon = nlons[j]
        ph_north[j] = Metal.placeholderTensor(graph, (nlon, nlayers), NF)
        res_north[j] = Metal.realToHermiteanFFTWithTensor(graph, ph_north[j], axes, fft_desc)
        (has_equator && j == j_equator) && continue   # south equator ring: north already covers it
        ph_south[j] = Metal.placeholderTensor(graph, (nlon, nlayers), NF)
        res_south[j] = Metal.realToHermiteanFFTWithTensor(graph, ph_south[j], axes, fft_desc)
    end
    return graph, ph_north, res_north, ph_south, res_south
end

function build_fused_inverse_graph(NF, nlayers::Int, nlat_half::Int, has_equator::Bool, j_equator::Int, nlons::Vector{Int})
    graph = Metal.MPSGraph()
    fft_desc = Metal.MPSGraphFFTDescriptor(; inverse = true)
    axes = fft_axis()
    ph_north = Vector{Metal.MPSGraphTensor}(undef, nlat_half)
    res_north = Vector{Metal.MPSGraphTensor}(undef, nlat_half)
    ph_south = Vector{Metal.MPSGraphTensor}(undef, nlat_half)
    res_south = Vector{Metal.MPSGraphTensor}(undef, nlat_half)
    @inbounds for j in 1:nlat_half
        nfreq = nlons[j] ÷ 2 + 1
        ph_north[j] = Metal.placeholderTensor(graph, (nfreq, nlayers), Complex{NF})
        res_north[j] = Metal.HermiteanToRealFFTWithTensor(graph, ph_north[j], axes, fft_desc)
        (has_equator && j == j_equator) && continue
        ph_south[j] = Metal.placeholderTensor(graph, (nfreq, nlayers), Complex{NF})
        res_south[j] = Metal.HermiteanToRealFFTWithTensor(graph, ph_south[j], axes, fft_desc)
    end
    return graph, ph_north, res_north, ph_south, res_south
end

function build_fused_cache(S::SpectralTransform, nlayers::Int)
    NF = eltype(S)
    nlat_half = S.grid.nlat_half
    nlat = S.nlat
    nlons = S.nlons
    has_equator = isodd(nlat)
    j_equator = (nlat + 1) ÷ 2

    fgraph, fph_n, fres_n, fph_s, fres_s = build_fused_forward_graph(NF, nlayers, nlat_half, has_equator, j_equator, nlons)
    igraph, iph_n, ires_n, iph_s, ires_s = build_fused_inverse_graph(NF, nlayers, nlat_half, has_equator, j_equator, nlons)

    return MetalFusedFourierGraphCache(
        fgraph, fph_n, fres_n, fph_s, fres_s,
        igraph, iph_n, ires_n, iph_s, ires_s,
        nlat_half, nlat, has_equator, j_equator,
    )
end

function get_fused_cache(S::SpectralTransform, nlayers::Int)
    caches_by_K = get!(() -> Dict{Int, MetalFusedFourierGraphCache}(), FUSED_FOURIER_GRAPH_CACHES, S)
    return get!(() -> build_fused_cache(S, nlayers), caches_by_K, nlayers)::MetalFusedFourierGraphCache
end

"""$(TYPEDSIGNATURES)
Clear all cached Metal fused-graph Fourier caches (frees the associated GPU memory).
Mainly useful for tests/benchmarks."""
clear_metal_fused_fourier_cache!() = (empty!(FUSED_FOURIER_GRAPH_CACHES); nothing)

function forward_loop_fused!(buf::MetalFourierCommandBatchCache, fused::MetalFusedFourierGraphCache, f_north, f_south, field::AbstractField, S::SpectralTransform)
    (; rings, nlons) = S
    (; ring_real_north, ring_real_south, ring_complex_north, ring_complex_south) = buf
    (; nlat_half, nlat, has_equator, j_equator) = fused
    queue = Metal.global_queue(Metal.device())

    @inbounds for j in 1:nlat_half
        copyto!(ring_real_north[j], view(field.data, rings[j], :))
    end
    @inbounds for j_south in (nlat_half + 1):nlat
        j = nlat - j_south + 1
        copyto!(ring_real_south[j], view(field.data, rings[j_south], :))
    end

    feeds = Dict{Metal.MPSGraphTensor, Metal.MPSGraphTensorData}()
    resultsdict = Dict{Metal.MPSGraphTensor, Metal.MPSGraphTensorData}()
    @inbounds for j in 1:nlat_half
        feeds[fused.forward_ph_north[j]] = Metal.MPSGraphTensorData(ring_real_north[j])
        resultsdict[fused.forward_res_north[j]] = Metal.MPSGraphTensorData(ring_complex_north[j])
        (has_equator && j == j_equator) && continue
        feeds[fused.forward_ph_south[j]] = Metal.MPSGraphTensorData(ring_real_south[j])
        resultsdict[fused.forward_res_south[j]] = Metal.MPSGraphTensorData(ring_complex_south[j])
    end

    cmdbuf = Metal.MPS.MPSCommandBuffer(queue)
    Metal.MPS.encode!(
        cmdbuf, fused.forward_graph, Metal.NSDictionary(feeds), Metal.NSDictionary(resultsdict),
        Metal.ObjectiveC.nil, Metal.MPSGraphs.default_exec_desc()
    )
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

function inverse_loop_fused!(buf::MetalFourierCommandBatchCache, fused::MetalFusedFourierGraphCache, field::AbstractField, g_north, g_south, S::SpectralTransform; add::Bool = false)
    (; rings, nlons) = S
    nlayers = size(field, 2)
    (; ring_real_north, ring_real_south, ring_complex_north, ring_complex_south) = buf
    (; nlat_half, nlat, has_equator, j_equator) = fused
    queue = Metal.global_queue(Metal.device())

    @inbounds for j in 1:nlat_half
        nfreq = nlons[j] ÷ 2 + 1
        copyto!(ring_complex_north[j], view(g_north, 1:nfreq, 1:nlayers, j))
        (has_equator && j == j_equator) && continue
        copyto!(ring_complex_south[j], view(g_south, 1:nfreq, 1:nlayers, j))
    end

    feeds = Dict{Metal.MPSGraphTensor, Metal.MPSGraphTensorData}()
    resultsdict = Dict{Metal.MPSGraphTensor, Metal.MPSGraphTensorData}()
    @inbounds for j in 1:nlat_half
        feeds[fused.inverse_ph_north[j]] = Metal.MPSGraphTensorData(ring_complex_north[j])
        resultsdict[fused.inverse_res_north[j]] = Metal.MPSGraphTensorData(ring_real_north[j])
        (has_equator && j == j_equator) && continue
        feeds[fused.inverse_ph_south[j]] = Metal.MPSGraphTensorData(ring_complex_south[j])
        resultsdict[fused.inverse_res_south[j]] = Metal.MPSGraphTensorData(ring_real_south[j])
    end

    cmdbuf = Metal.MPS.MPSCommandBuffer(queue)
    Metal.MPS.encode!(
        cmdbuf, fused.inverse_graph, Metal.NSDictionary(feeds), Metal.NSDictionary(resultsdict),
        Metal.ObjectiveC.nil, Metal.MPSGraphs.default_exec_desc()
    )
    Metal.commit!(cmdbuf)

    @inbounds for j in 1:nlat_half
        dest = view(field.data, rings[j], :)
        add ? (dest .+= ring_real_north[j]) : copyto!(dest, ring_real_north[j])
    end
    @inbounds for j_south in (nlat_half + 1):nlat
        j = nlat - j_south + 1
        (has_equator && j == j_equator) && continue
        dest = view(field.data, rings[j_south], :)
        add ? (dest .+= ring_real_south[j]) : copyto!(dest, ring_real_south[j])
    end
    return nothing
end

# =====================================================================================
# Method overrides: dispatch on MtlArray scratch (more specific than the generic
# AbstractArray{<:Complex,3} methods in fourier.jl), then dispatch on
# `S.metal_fused_fourier`/`S.metal_merge_hemispheres` to pick the strategy.
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
    cache = get_cache(S, size(field, 2))
    if S.metal_fused_fourier
        forward_loop_fused!(cache, get_fused_cache(S, size(field, 2)), f_north, f_south, field, S)
    elseif S.metal_merge_hemispheres
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
        S::SpectralTransform;
        add::Bool = false,
    )
    if !S.cuda_graphs
        return Base.@invoke _fourier_batched!(
            field::AbstractField, g_north::AbstractArray{<:Complex, 3},
            g_south::AbstractArray{<:Complex, 3}, S::SpectralTransform;
            add,
        )
    end
    cache = get_cache(S, size(field, 2))
    if S.metal_fused_fourier
        inverse_loop_fused!(cache, get_fused_cache(S, size(field, 2)), field, g_north, g_south, S; add)
    elseif S.metal_merge_hemispheres
        inverse_loop_merged!(cache, field, g_north, g_south, S; add)
    else
        inverse_loop_unmerged!(cache, field, g_north, g_south, S; add)
    end
    return nothing
end

end
