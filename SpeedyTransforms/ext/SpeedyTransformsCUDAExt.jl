module SpeedyTransformsCUDAExt

import CUDA: CUDA, CuArray, CuVector, CuGraphExec, capture, instantiate, launch
using cuFFT
import AbstractFFTs
import LinearAlgebra
import LinearAlgebra: mul!
using KernelAbstractions
using DocStringExtensions

using SpeedyTransforms
using SpeedyTransforms.RingGrids
using SpeedyTransforms.LowerTriangularArrays

import SpeedyTransforms: SpectralTransform, _fourier_batched!
import SpeedyTransforms.RingGrids: AbstractField

import SpeedyWeatherInternals.KernelLaunching: launch!, ArrayWorkOrder
import SpeedyWeatherInternals.Architectures: on_architecture, synchronize, CUDAGPU

# =====================================================================================
# CUDA GRAPHS ACCELERATION OF THE BATCHED FOURIER TRANSFORM
#
# The batched Fourier transform on a (reduced) grid applies one cuFFT per latitude ring
# (≈ 2*nlat_half tiny FFTs). On the GPU this is heavily launch-bound: the CPU spends far
# more time enqueuing the many small kernels than the GPU spends computing them.
#
# CUDA Graphs let us record this whole sequence ONCE and replay it with a single `launch`,
# eliminating the per-operation CPU launch overhead. Two requirements shape the design:
#
#  1. No device allocations inside the captured region (the generic GPU path allocates a
#     temporary per ring via `field.data[ilons, :]` and `plan * x`). We therefore pre-
#     allocate one contiguous packed buffer that holds every ring's dense block, and use
#     in-place `mul!` reading/writing reshaped views into it.
#  2. Few graph nodes. A single KernelAbstractions gather/scatter kernel moves ALL rings
#     between the strided grid/scratch layout and the packed buffer in one launch.
#     This collapses the graph from ~6*nlat_half nodes to ~2*nlat_half + 4.
#
# The captured graph bakes in the device pointers of the input `field.data`,  and the
# scratch buffers (`scratch_memory.north/.south`) and the packed work
# buffer `GPUFourierGraphCache`. In the SpeedyWeather time loop the same variable buffers are reused every
# timestep, so a graph captured for a given `field` is replayed on all subsequent steps.
# Caches are held per transform SIZE (keyed by the FFT plan, so an `S` that batches
# several layer counts gets one cache each); within a cache, graphs are keyed by the device
# pointer of `field.data` (see `graph_key`), keying on the wrapper object identity wouldn't
# work as a view on a CuArray is again a CuArray and not a SubArray which causes problems here.
#
# The backend-agnostic parts of this (kernels, cache struct, allocation-free loops,
# capture/replay control flow) live in gpu_graphs_common.jl, included below and shared with
# SpeedyTransformsAMDGPUExt.jl; only the capture/instantiate/launch primitives differ per
# backend (see `GraphBackend`).
# =====================================================================================

include("gpu_graphs_common.jl")

# The CUDA capture/instantiate/launch primitives `run_graph!` needs, plus the architecture
# used to `synchronize` before capture; see `GraphBackend` in gpu_graphs_common.jl. `capture`
# must not throw on an uncapturable region, hence `throw_error = false`.
const CUDA_GRAPH_BACKEND = GraphBackend(
    loop! -> capture(loop!; throw_error = false),
    instantiate,
    launch,
    CUDAGPU(),
)

# =====================================================================================
# Method overrides: dispatch on CuArray scratch (more specific than the generic
# AbstractArray{<:Complex,3} methods in fourier.jl)
# =====================================================================================

"""$(TYPEDSIGNATURES)
CUDA-Graphs accelerated forward (grid → spectral) batched Fourier transform.
Replays a cached CUDA graph of the fused gather + per-ring cuFFTs + scatter; see
[`run_graph!`](@ref)."""
function _fourier_batched!(
        f_north::CuArray{<:Complex, 3},
        f_south::CuArray{<:Complex, 3},
        field::AbstractField,
        S::SpectralTransform,
    )
    @assert eltype(field) == eltype(S) "Number format of grid $(eltype(field)) and SpectralTransform $(eltype(S)) need to match."
    if !S.gpu_graphs
        return Base.@invoke _fourier_batched!(
            f_north::AbstractArray{<:Complex, 3}, f_south::AbstractArray{<:Complex, 3},
            field::AbstractField, S::SpectralTransform,
        )
    end
    # the cache is selected by (and sized for) this transform size (= its FFT plan set); within
    # it the only thing that varies between calls is the field buffer, so that is the graph key
    cache = get_cache(CuGraphExec, S, size(field, 2))
    run_graph!(cache.forward_execs, graph_key(field.data), () -> forward_loop!(cache, f_north, f_south, field, S), CUDA_GRAPH_BACKEND)
    return nothing
end

"""$(TYPEDSIGNATURES)
CUDA-Graphs accelerated inverse (spectral → grid) batched Fourier transform.
Replays a cached CUDA graph of the fused gather + per-ring inverse cuFFTs + scatter; see
[`run_graph!`](@ref)."""
function _fourier_batched!(
        field::AbstractField,
        g_north::CuArray{<:Complex, 3},
        g_south::CuArray{<:Complex, 3},
        S::SpectralTransform;
        add::Bool = false,          # accumulate onto `field` instead of overwriting? (Enzyme adjoint rule)
    )
    if !S.gpu_graphs
        return Base.@invoke _fourier_batched!(
            field::AbstractField, g_north::AbstractArray{<:Complex, 3},
            g_south::AbstractArray{<:Complex, 3}, S::SpectralTransform; add,
        )
    end
    cache = get_cache(CuGraphExec, S, size(field, 2))
    # `add` is part of the graph key: a captured graph bakes in overwrite-vs-accumulate, so replaying
    # the wrong one would silently produce incorrect results.
    run_graph!(cache.inverse_execs, (graph_key(field.data), add), () -> inverse_loop!(cache, field, g_north, g_south, S, add), CUDA_GRAPH_BACKEND)
    return nothing
end

end
