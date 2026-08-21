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
# Fused single-graph Fourier: every ring's FFT (both hemispheres, one direction) is a
# separate op inside ONE shared MPSGraph, built once and cached. Forward and inverse each
# issue two `encode!`/`commit!` calls per `transform!` (one for north, one for south, each
# on its own command buffer) — found to consistently beat sharing a single command buffer
# for both hemispheres, see metal_fourier_transform_status.md.
# =====================================================================================

"""$(TYPEDSIGNATURES)
Per-`(SpectralTransform, nlayers)` cache of reusable, independently-allocated per-ring
buffers for the Metal command-buffer-batched Fourier transform. One real (`nlon_j ×
nlayers`) and one complex (`nfreq_j × nlayers`) buffer per ring, separately for north and
south (both hemispheres are gathered before either is committed). The real buffers double as the forward
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

# =====================================================================================
# Fused single-graph strategy: every ring's FFT (both hemispheres) is a separate,
# independent op inside ONE shared MPSGraph, built once and cached. A `transform!` call
# then issues exactly TWO `encode!` dispatches per direction (one for north, one for
# south, covering every ring in that hemisphere at once) instead of one per ring.
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

# North and south are each encoded onto their own command buffer and committed
# separately (rather than sharing one commit)
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

    feeds_n = Dict{Metal.MPSGraphTensor, Metal.MPSGraphTensorData}()
    results_n = Dict{Metal.MPSGraphTensor, Metal.MPSGraphTensorData}()
    @inbounds for j in 1:nlat_half
        feeds_n[fused.forward_ph_north[j]] = Metal.MPSGraphTensorData(ring_real_north[j])
        results_n[fused.forward_res_north[j]] = Metal.MPSGraphTensorData(ring_complex_north[j])
    end
    cmdbuf_n = Metal.MPS.MPSCommandBuffer(queue)
    Metal.MPS.encode!(
        cmdbuf_n, fused.forward_graph, Metal.NSDictionary(feeds_n), Metal.NSDictionary(results_n),
        Metal.ObjectiveC.nil, Metal.MPSGraphs.default_exec_desc()
    )
    Metal.commit!(cmdbuf_n)

    feeds_s = Dict{Metal.MPSGraphTensor, Metal.MPSGraphTensorData}()
    results_s = Dict{Metal.MPSGraphTensor, Metal.MPSGraphTensorData}()
    @inbounds for j in 1:nlat_half
        (has_equator && j == j_equator) && continue
        feeds_s[fused.forward_ph_south[j]] = Metal.MPSGraphTensorData(ring_real_south[j])
        results_s[fused.forward_res_south[j]] = Metal.MPSGraphTensorData(ring_complex_south[j])
    end
    cmdbuf_s = Metal.MPS.MPSCommandBuffer(queue)
    Metal.MPS.encode!(
        cmdbuf_s, fused.forward_graph, Metal.NSDictionary(feeds_s), Metal.NSDictionary(results_s),
        Metal.ObjectiveC.nil, Metal.MPSGraphs.default_exec_desc()
    )
    Metal.commit!(cmdbuf_s)

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

    feeds_n = Dict{Metal.MPSGraphTensor, Metal.MPSGraphTensorData}()
    results_n = Dict{Metal.MPSGraphTensor, Metal.MPSGraphTensorData}()
    @inbounds for j in 1:nlat_half
        feeds_n[fused.inverse_ph_north[j]] = Metal.MPSGraphTensorData(ring_complex_north[j])
        results_n[fused.inverse_res_north[j]] = Metal.MPSGraphTensorData(ring_real_north[j])
    end
    cmdbuf_n = Metal.MPS.MPSCommandBuffer(queue)
    Metal.MPS.encode!(
        cmdbuf_n, fused.inverse_graph, Metal.NSDictionary(feeds_n), Metal.NSDictionary(results_n),
        Metal.ObjectiveC.nil, Metal.MPSGraphs.default_exec_desc()
    )
    Metal.commit!(cmdbuf_n)

    feeds_s = Dict{Metal.MPSGraphTensor, Metal.MPSGraphTensorData}()
    results_s = Dict{Metal.MPSGraphTensor, Metal.MPSGraphTensorData}()
    @inbounds for j in 1:nlat_half
        (has_equator && j == j_equator) && continue
        feeds_s[fused.inverse_ph_south[j]] = Metal.MPSGraphTensorData(ring_complex_south[j])
        results_s[fused.inverse_res_south[j]] = Metal.MPSGraphTensorData(ring_real_south[j])
    end
    cmdbuf_s = Metal.MPS.MPSCommandBuffer(queue)
    Metal.MPS.encode!(
        cmdbuf_s, fused.inverse_graph, Metal.NSDictionary(feeds_s), Metal.NSDictionary(results_s),
        Metal.ObjectiveC.nil, Metal.MPSGraphs.default_exec_desc()
    )
    Metal.commit!(cmdbuf_s)

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
    if !S.cuda_graphs   # generic toggle for "use accelerated batched Fourier path"; also gates Metal fused Fourier
        return Base.@invoke _fourier_batched!(
            f_north::AbstractArray{<:Complex, 3}, f_south::AbstractArray{<:Complex, 3},
            field::AbstractField, S::SpectralTransform,
        )
    end
    forward_loop_fused!(get_cache(S, size(field, 2)), get_fused_cache(S, size(field, 2)), f_north, f_south, field, S)
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
    if !S.cuda_graphs
        return Base.@invoke _fourier_batched!(
            field::AbstractField, g_north::AbstractArray{<:Complex, 3},
            g_south::AbstractArray{<:Complex, 3}, S::SpectralTransform;
            add,
        )
    end
    inverse_loop_fused!(get_cache(S, size(field, 2)), get_fused_cache(S, size(field, 2)), field, g_north, g_south, S; add)
    return nothing
end

end
