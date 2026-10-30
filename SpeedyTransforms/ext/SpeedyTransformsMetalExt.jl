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

# Dissable gpu_graphs for debugging with Metal GPU CI pipeline
import SpeedyTransforms: SpectralTransform, _fourier_batched!, default_gpu_graphs
default_gpu_graphs(::Metal.MetalBackend) = false

# =====================================================================================
# Single-graph Fourier, north+south combined per ring: north ring j and its mirrored
# south ring's data are combined into ONE tensor (`2*nlayers` columns: north in
# 1:nlayers, south in nlayers+1:2*nlayers) and hit with ONE FFT op, instead of two
# separate ops. Halves the ring-op count (~nlat_half ops instead of ~2*nlat_half) and the
# per-ring Dict/NSDictionary construction overhead.
# =====================================================================================

"""$(TYPEDSIGNATURES)
Per-`(SpectralTransform, nlayers)` cache of reusable per-ring buffers for the Metal
Fourier transform: one real (`nlon_j × 2nlayers`) and one complex (`nfreq_j × 2nlayers`)
buffer per ring `j in 1:nlat_half`, columns `1:nlayers` holding northern data and
`nlayers+1:2nlayers` the mirrored southern ring's data (unused/stale at the equator ring,
which has no distinct southern counterpart — never read back out, see
`forward_loop_fused!`). The real buffers double as the forward input / inverse output, the
complex buffers as the forward output / inverse input. Buffers are sized to the actual
per-call `nlayers` (the batch `K`), not `S.nlayers` (the scratch-memory upper bound) — a
`SpectralTransform` built via `SpectralGrid` plans several distinct `K`s (e.g. `1, 2,
nlayers, 2nlayers, ...`) and each needs its own correctly-sized buffers to match its own
`rfft_plans_batched[K]`/`brfft_plans_batched[K]` plans."""
struct MetalFourierBufferCache{RT, CT}
    ring_real_both::Vector{RT}
    ring_complex_both::Vector{CT}
    nlat_half::Int
    nlat::Int
    has_equator::Bool
    j_equator::Int
end

# One cache per (SpectralTransform, nlayers), the outer IdDict keyed by object identity
# (stable for the transform's lifetime; saved as a global and not a field of `S` so `S`
# stays a concrete, stable type), the inner Dict keyed by the batch size `K` actually used
# — a single `S` serves multiple `K`s (see `MetalFourierBufferCache` docstring).
const FOURIER_BUFFER_CACHES = IdDict{SpectralTransform, Dict{Int, MetalFourierBufferCache}}()

function build_buffer_cache(S::SpectralTransform, nlayers::Int)
    NF = eltype(S)
    nlat_half = S.grid.nlat_half
    nlat = S.nlat
    nlons = S.nlons

    ring_real_both = [Metal.zeros(NF, nlons[j], 2 * nlayers) for j in 1:nlat_half]
    ring_complex_both = [Metal.zeros(Complex{NF}, nlons[j] ÷ 2 + 1, 2 * nlayers) for j in 1:nlat_half]

    has_equator = isodd(nlat)
    j_equator = (nlat + 1) ÷ 2

    return MetalFourierBufferCache(ring_real_both, ring_complex_both, nlat_half, nlat, has_equator, j_equator)
end

function get_buffer_cache(S::SpectralTransform, nlayers::Int)
    caches_by_K = get!(() -> Dict{Int, MetalFourierBufferCache}(), FOURIER_BUFFER_CACHES, S)
    return get!(() -> build_buffer_cache(S, nlayers), caches_by_K, nlayers)::MetalFourierBufferCache
end

"""$(TYPEDSIGNATURES)
Per-`(SpectralTransform, nlayers)` cache holding one shared `MPSGraph` for the forward
direction and one for the inverse direction, each containing ONE FFT op per ring `j in
1:nlat_half` — operating on the combined `2*nlayers`-wide north+south buffer for that
ring. Built once; every `transform!` call reuses the same graphs and just re-binds fresh
`MPSGraphTensorData` (wrapping the current per-ring buffers) to the cached
placeholder/result tensors."""
struct MetalFourierGraphCache
    forward_graph::Metal.MPSGraph
    forward_ph::Vector{Metal.MPSGraphTensor}
    forward_res::Vector{Metal.MPSGraphTensor}
    inverse_graph::Metal.MPSGraph
    inverse_ph::Vector{Metal.MPSGraphTensor}
    inverse_res::Vector{Metal.MPSGraphTensor}
    nlat_half::Int
end

const FOURIER_GRAPH_CACHES = IdDict{SpectralTransform, Dict{Int, MetalFourierGraphCache}}()

# axis for a (nlon_j, 2*nlayers) or (nfreq_j, 2*nlayers) ring array with FFT along dim 1:
# Julia axis 1 of a 2-D array -> Metal axis (ndims - axis) = 2 - 1 = 1 (see the identical
# comment/derivation in Metal.jl's own `CachedFFTGraph`, fft.jl:246-249)
fft_axis() = Metal.NSArray([Metal.NSNumber(1)])

function build_fourier_forward_graph(NF, nlayers::Int, nlat_half::Int, nlons::Vector{Int})
    graph = Metal.MPSGraph()
    fft_desc = Metal.MPSGraphFFTDescriptor(; inverse = false)
    axes = fft_axis()
    ph = Vector{Metal.MPSGraphTensor}(undef, nlat_half)
    res = Vector{Metal.MPSGraphTensor}(undef, nlat_half)
    @inbounds for j in 1:nlat_half
        nlon = nlons[j]
        ph[j] = Metal.placeholderTensor(graph, (nlon, 2 * nlayers), NF)
        res[j] = Metal.realToHermiteanFFTWithTensor(graph, ph[j], axes, fft_desc)
    end
    return graph, ph, res
end

function build_fourier_inverse_graph(NF, nlayers::Int, nlat_half::Int, nlons::Vector{Int})
    graph = Metal.MPSGraph()
    fft_desc = Metal.MPSGraphFFTDescriptor(; inverse = true)
    axes = fft_axis()
    ph = Vector{Metal.MPSGraphTensor}(undef, nlat_half)
    res = Vector{Metal.MPSGraphTensor}(undef, nlat_half)
    @inbounds for j in 1:nlat_half
        nfreq = nlons[j] ÷ 2 + 1
        ph[j] = Metal.placeholderTensor(graph, (nfreq, 2 * nlayers), Complex{NF})
        res[j] = Metal.HermiteanToRealFFTWithTensor(graph, ph[j], axes, fft_desc)
    end
    return graph, ph, res
end

function build_fourier_graph_cache(S::SpectralTransform, nlayers::Int)
    NF = eltype(S)
    nlat_half = S.grid.nlat_half
    nlons = S.nlons

    fgraph, fph, fres = build_fourier_forward_graph(NF, nlayers, nlat_half, nlons)
    igraph, iph, ires = build_fourier_inverse_graph(NF, nlayers, nlat_half, nlons)

    return MetalFourierGraphCache(fgraph, fph, fres, igraph, iph, ires, nlat_half)
end

function get_fourier_graph_cache(S::SpectralTransform, nlayers::Int)
    caches_by_K = get!(() -> Dict{Int, MetalFourierGraphCache}(), FOURIER_GRAPH_CACHES, S)
    return get!(() -> build_fourier_graph_cache(S, nlayers), caches_by_K, nlayers)::MetalFourierGraphCache
end

"""$(TYPEDSIGNATURES)
Clear all cached Metal Fourier buffer/graph caches (frees the associated GPU memory).
Mainly useful for tests/benchmarks."""
clear_metal_fourier_cache!() = (empty!(FOURIER_BUFFER_CACHES); empty!(FOURIER_GRAPH_CACHES); nothing)

function forward_loop_fused!(buf::MetalFourierBufferCache, fused::MetalFourierGraphCache, f_north, f_south, field::AbstractField, S::SpectralTransform)
    (; rings, nlons) = S
    (; ring_real_both, ring_complex_both) = buf
    (; nlat_half, nlat, has_equator, j_equator) = buf
    nlayers = size(f_north, 2)
    queue = Metal.global_queue(Metal.device())

    @inbounds for j in 1:nlat_half
        copyto!(view(ring_real_both[j], :, 1:nlayers), view(field.data, rings[j], :))
    end
    @inbounds for j_south in (nlat_half + 1):nlat
        j = nlat - j_south + 1
        (has_equator && j == j_equator) && continue
        copyto!(view(ring_real_both[j], :, (nlayers + 1):(2 * nlayers)), view(field.data, rings[j_south], :))
    end

    feeds = Dict{Metal.MPSGraphTensor, Metal.MPSGraphTensorData}()
    results = Dict{Metal.MPSGraphTensor, Metal.MPSGraphTensorData}()
    @inbounds for j in 1:nlat_half
        feeds[fused.forward_ph[j]] = Metal.MPSGraphTensorData(ring_real_both[j])
        results[fused.forward_res[j]] = Metal.MPSGraphTensorData(ring_complex_both[j])
    end
    cmdbuf = Metal.MPS.MPSCommandBuffer(queue)
    Metal.MPS.encode!(
        cmdbuf, fused.forward_graph, Metal.NSDictionary(feeds), Metal.NSDictionary(results),
        Metal.ObjectiveC.nil, Metal.MPSGraphs.default_exec_desc()
    )
    Metal.commit!(cmdbuf)

    @inbounds for j in 1:nlat_half
        nfreq = nlons[j] ÷ 2 + 1
        copyto!(view(f_north, 1:nfreq, :, j), view(ring_complex_both[j], 1:nfreq, 1:nlayers))
        if has_equator && j == j_equator
            fill!(view(f_south, 1:nfreq, :, j), 0)
        else
            copyto!(view(f_south, 1:nfreq, :, j), view(ring_complex_both[j], 1:nfreq, (nlayers + 1):(2 * nlayers)))
        end
    end
    return nothing
end

function inverse_loop_fused!(buf::MetalFourierBufferCache, fused::MetalFourierGraphCache, field::AbstractField, g_north, g_south, S::SpectralTransform; add::Bool = false)
    (; rings, nlons) = S
    nlayers = size(field, 2)
    (; ring_real_both, ring_complex_both) = buf
    (; nlat_half, nlat, has_equator, j_equator) = buf
    queue = Metal.global_queue(Metal.device())

    @inbounds for j in 1:nlat_half
        nfreq = nlons[j] ÷ 2 + 1
        copyto!(view(ring_complex_both[j], 1:nfreq, 1:nlayers), view(g_north, 1:nfreq, 1:nlayers, j))
        (has_equator && j == j_equator) && continue
        copyto!(view(ring_complex_both[j], 1:nfreq, (nlayers + 1):(2 * nlayers)), view(g_south, 1:nfreq, 1:nlayers, j))
    end

    feeds = Dict{Metal.MPSGraphTensor, Metal.MPSGraphTensorData}()
    results = Dict{Metal.MPSGraphTensor, Metal.MPSGraphTensorData}()
    @inbounds for j in 1:nlat_half
        feeds[fused.inverse_ph[j]] = Metal.MPSGraphTensorData(ring_complex_both[j])
        results[fused.inverse_res[j]] = Metal.MPSGraphTensorData(ring_real_both[j])
    end
    cmdbuf = Metal.MPS.MPSCommandBuffer(queue)
    Metal.MPS.encode!(
        cmdbuf, fused.inverse_graph, Metal.NSDictionary(feeds), Metal.NSDictionary(results),
        Metal.ObjectiveC.nil, Metal.MPSGraphs.default_exec_desc()
    )
    Metal.commit!(cmdbuf)

    @inbounds for j in 1:nlat_half
        dest = view(field.data, rings[j], :)
        src = view(ring_real_both[j], :, 1:nlayers)
        add ? (dest .+= src) : copyto!(dest, src)
    end
    @inbounds for j_south in (nlat_half + 1):nlat
        j = nlat - j_south + 1
        (has_equator && j == j_equator) && continue
        dest = view(field.data, rings[j_south], :)
        src = view(ring_real_both[j], :, (nlayers + 1):(2 * nlayers))
        add ? (dest .+= src) : copyto!(dest, src)
    end
    return nothing
end

# =====================================================================================
# Method overrides: dispatch on MtlArray scratch (more specific than the generic
# AbstractArray{<:Complex,3} methods in fourier.jl) onto the fused single-graph strategy.
# =====================================================================================

"""$(TYPEDSIGNATURES)
Metal fused single-graph forward (grid → spectral) batched Fourier transform — see the
module-level comment."""
function SpeedyTransforms._fourier_batched!(
        f_north::MtlArray{<:Complex, 3},
        f_south::MtlArray{<:Complex, 3},
        field::AbstractField,
        S::SpectralTransform,
    )
    @assert eltype(field) == eltype(S) "Number format of grid $(eltype(field)) and SpectralTransform $(eltype(S)) need to match."
    if !S.gpu_graphs   # generic toggle for "use accelerated batched Fourier path"; also gates Metal fused Fourier
        return Base.@invoke _fourier_batched!(
            f_north::AbstractArray{<:Complex, 3}, f_south::AbstractArray{<:Complex, 3},
            field::AbstractField, S::SpectralTransform,
        )
    end
    forward_loop_fused!(get_buffer_cache(S, size(field, 2)), get_fourier_graph_cache(S, size(field, 2)), f_north, f_south, field, S)
    return nothing
end

"""$(TYPEDSIGNATURES)
Metal fused single-graph inverse (spectral → grid) batched Fourier transform — see the
module-level comment."""
function SpeedyTransforms._fourier_batched!(
        field::AbstractField,
        g_north::MtlArray{<:Complex, 3},
        g_south::MtlArray{<:Complex, 3},
        S::SpectralTransform;
        add::Bool = false,
    )
    if !S.gpu_graphs
        return Base.@invoke _fourier_batched!(
            field::AbstractField, g_north::AbstractArray{<:Complex, 3},
            g_south::AbstractArray{<:Complex, 3}, S::SpectralTransform;
            add,
        )
    end
    inverse_loop_fused!(get_buffer_cache(S, size(field, 2)), get_fourier_graph_cache(S, size(field, 2)), field, g_north, g_south, S; add)
    return nothing
end

end
