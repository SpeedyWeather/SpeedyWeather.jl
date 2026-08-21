module SpeedyTransformsAMDGPUExt

import AMDGPU: AMDGPU, ROCArray
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
import SpeedyWeatherInternals.Architectures: on_architecture, GPU

# =====================================================================================
# HIP GRAPHS ACCELERATION OF THE BATCHED FOURIER TRANSFORM
#
# HIP-graphs equivalent of SpeedyTransformsCUDAExt.jl: captures the fused gather/rocFFT/scatter
# loop of the batched Fourier transform once and replays it with a single `launch`, eliminating
# per-operation CPU launch overhead on the GPU-launch-bound Fourier transform. See
# SpeedyTransformsCUDAExt.jl for the full rationale (packed buffer, single gather/scatter
# kernel per direction, etc.) — identical here, only the capture/instantiate/launch primitives
# differ (`AMDGPU.HIP.capture`/`instantiate`/`launch` instead of `CUDA.capture`/`instantiate`/
# `launch`).
#
# The backend-agnostic parts of this (kernels, cache struct, allocation-free loops,
# capture/replay control flow) live in gpu_graphs_common.jl, included below and shared with
# SpeedyTransformsCUDAExt.jl; only the capture/instantiate/launch primitives differ per backend
# (see `GraphBackend`).
# =====================================================================================

include("gpu_graphs_common.jl")

# Probe for the high-level HIP graph API. AMDGPU.HIP exports these in newer versions; on
# older installs only the raw C bindings (hipGraph_t, hipGraphExec_t, …) are present, and
# `_fourier_batched!` falls back to the unaccelerated path with a one-time warning.
const _HIP_GRAPHS_AVAILABLE =
    isdefined(AMDGPU, :HIP) &&
    isdefined(AMDGPU.HIP, :HIPGraphExec) &&
    isdefined(AMDGPU.HIP, :capture)

# The HIP capture/instantiate/launch primitives `run_graph!` needs, plus the architecture used
# to `synchronize` before capture; see `GraphBackend` in gpu_graphs_common.jl. `capture` must
# not throw on an uncapturable region, hence `throw_error = false`. `nothing` on older AMDGPU.jl
# installs where the high-level HIP graph API isn't available at all.
const HIP_GRAPH_BACKEND = if _HIP_GRAPHS_AVAILABLE
    GraphBackend(
        loop! -> AMDGPU.HIP.capture(loop!; throw_error = false),
        AMDGPU.HIP.instantiate,
        AMDGPU.HIP.launch,
        GPU(AMDGPU.ROCBackend()),
    )
else
    nothing
end

const GPU_GRAPHS_WARNED = Ref(false)

function warn_gpu_graphs_unavailable()
    if !GPU_GRAPHS_WARNED[]
        @warn "GPU-graphs acceleration of the Fourier transform (`gpu_graphs = true`) requires a newer AMDGPU.jl exposing the high-level HIP graph API (`AMDGPU.HIP.capture`/`instantiate`/`HIPGraphExec`); falling back to the unaccelerated path. Set `gpu_graphs = false` to silence this warning."
        GPU_GRAPHS_WARNED[] = true
    end
    return nothing
end

# =====================================================================================
# Method overrides: dispatch on ROCArray scratch (more specific than the generic
# AbstractArray{<:Complex,3} methods in fourier.jl)
# =====================================================================================

"""$(TYPEDSIGNATURES)
HIP-Graphs accelerated forward (grid → spectral) batched Fourier transform.
Replays a cached HIP graph of the fused gather + per-ring rocFFTs + scatter; see
[`run_graph!`](@ref)."""
function _fourier_batched!(
        f_north::ROCArray{<:Complex, 3},
        f_south::ROCArray{<:Complex, 3},
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
    if !_HIP_GRAPHS_AVAILABLE
        warn_gpu_graphs_unavailable()
        return Base.@invoke _fourier_batched!(
            f_north::AbstractArray{<:Complex, 3}, f_south::AbstractArray{<:Complex, 3},
            field::AbstractField, S::SpectralTransform,
        )
    end
    # the cache is selected by (and sized for) this transform size (= its FFT plan set); within
    # it the only thing that varies between calls is the field buffer, so that is the graph key
    cache = get_cache(AMDGPU.HIP.HIPGraphExec, S, size(field, 2))
    run_graph!(cache.forward_execs, graph_key(field.data), () -> forward_loop!(cache, f_north, f_south, field, S), HIP_GRAPH_BACKEND)
    return nothing
end

"""$(TYPEDSIGNATURES)
HIP-Graphs accelerated inverse (spectral → grid) batched Fourier transform.
Replays a cached HIP graph of the fused gather + per-ring inverse rocFFTs + scatter; see
[`run_graph!`](@ref)."""
function _fourier_batched!(
        field::AbstractField,
        g_north::ROCArray{<:Complex, 3},
        g_south::ROCArray{<:Complex, 3},
        S::SpectralTransform;
        add::Bool = false,          # accumulate onto `field` instead of overwriting? (Enzyme adjoint rule)
    )
    if !S.gpu_graphs
        return Base.@invoke _fourier_batched!(
            field::AbstractField, g_north::AbstractArray{<:Complex, 3},
            g_south::AbstractArray{<:Complex, 3}, S::SpectralTransform; add,
        )
    end
    if !_HIP_GRAPHS_AVAILABLE
        warn_gpu_graphs_unavailable()
        return Base.@invoke _fourier_batched!(
            field::AbstractField, g_north::AbstractArray{<:Complex, 3},
            g_south::AbstractArray{<:Complex, 3}, S::SpectralTransform; add,
        )
    end
    cache = get_cache(AMDGPU.HIP.HIPGraphExec, S, size(field, 2))
    # `add` is part of the graph key: a captured graph bakes in overwrite-vs-accumulate, so replaying
    # the wrong one would silently produce incorrect results.
    run_graph!(cache.inverse_execs, (graph_key(field.data), add), () -> inverse_loop!(cache, field, g_north, g_south, S, add), HIP_GRAPH_BACKEND)
    return nothing
end

end
