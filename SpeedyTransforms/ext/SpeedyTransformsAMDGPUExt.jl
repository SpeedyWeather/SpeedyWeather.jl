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
import SpeedyWeatherInternals.Architectures: on_architecture, synchronize, ROCGPU

# =====================================================================================
# HIP GRAPHS ACCELERATION OF THE BATCHED FOURIER TRANSFORM — CAPTURE NOT ENABLED
#
# CUDA graphs give the batched Fourier transform a large speedup on NVIDIA GPUs (see
# SpeedyTransformsCUDAExt.jl) by capturing the fused gather/FFT/scatter loop once and
# replaying it on later calls. The HIP-graphs equivalent for AMDGPU is not enabled: ROCm's
# stream-capture validator does not reliably reject operations that are illegal to capture —
# some raise a catchable HIPError, but others are silently accepted at capture time and only
# surface later as a GPU memory access fault when the (corrupt) graph is replayed, confirmed
# on real hardware (LUMI). Since capture failures aren't reliably catchable, "try capture and
# fall back on error" isn't safe, so `_fourier_batched!` below always falls back
# unconditionally to the generic (unaccelerated but correct) per-ring Fourier path in
# SpeedyTransforms/src/fourier.jl, with a one-time warning if the user asked for
# `gpu_graphs = true`.
#
# The backend-agnostic machinery (kernels, cache struct, allocation-free loops,
# capture/replay control flow) is shared with SpeedyTransformsCUDAExt.jl via
# gpu_graphs_common.jl. `HIP_GRAPH_BACKEND` below is wired up against that shared interface
# with real AMDGPU.HIP bindings, guarded by `_HIP_GRAPHS_AVAILABLE` so loading this extension
# never fails on older AMDGPU.jl installs that lack them — but it is not referenced by
# `_fourier_batched!`, so none of it runs today. Re-enabling capture, once ROCm/AMDGPU.jl
# provides reliable stream-capture validation, is then: swap the two `_fourier_batched!`
# methods below for `get_cache`/`run_graph!`-based ones mirroring SpeedyTransformsCUDAExt.jl,
# and test on real hardware before trusting it.
# =====================================================================================

include("gpu_graphs_common.jl")

# Probe for the high-level HIP graph API. AMDGPU.HIP exports these in newer versions; on
# older installs only the raw C bindings (hipGraph_t, hipGraphExec_t, …) are present.
const _HIP_GRAPHS_AVAILABLE =
    isdefined(AMDGPU, :HIP) &&
    isdefined(AMDGPU.HIP, :HIPGraphExec) &&
    isdefined(AMDGPU.HIP, :capture)

# Ready for when capture is re-enabled (see module comment above); NOT used by
# `_fourier_batched!` today. `nothing` on older AMDGPU.jl installs where the HIP graph API
# isn't available at all.
const HIP_GRAPH_BACKEND = if _HIP_GRAPHS_AVAILABLE
    GraphBackend(
        loop! -> AMDGPU.HIP.capture(loop!; throw_error = false),
        AMDGPU.HIP.instantiate,
        AMDGPU.HIP.launch,
        ROCGPU(),
    )
else
    nothing
end

const GPU_GRAPHS_WARNED = Ref(false)

function warn_gpu_graphs_unavailable()
    if !GPU_GRAPHS_WARNED[]
        @warn "GPU-graphs acceleration of the Fourier transform (`gpu_graphs = true`) is not implemented for AMDGPU/HIP; falling back to the unaccelerated path. Set `gpu_graphs = false` to silence this warning."
        GPU_GRAPHS_WARNED[] = true
    end
    return nothing
end

function _fourier_batched!(
        f_north::ROCArray{<:Complex, 3},
        f_south::ROCArray{<:Complex, 3},
        field::AbstractField,
        S::SpectralTransform,
    )
    S.gpu_graphs && warn_gpu_graphs_unavailable()
    return Base.@invoke _fourier_batched!(
        f_north::AbstractArray{<:Complex, 3}, f_south::AbstractArray{<:Complex, 3},
        field::AbstractField, S::SpectralTransform,
    )
end

function _fourier_batched!(
        field::AbstractField,
        g_north::ROCArray{<:Complex, 3},
        g_south::ROCArray{<:Complex, 3},
        S::SpectralTransform;
        add::Bool = false,          # accumulate onto `field` instead of overwriting? (Enzyme adjoint rule)
    )
    S.gpu_graphs && warn_gpu_graphs_unavailable()
    return Base.@invoke _fourier_batched!(
        field::AbstractField, g_north::AbstractArray{<:Complex, 3},
        g_south::AbstractArray{<:Complex, 3}, S::SpectralTransform; add,
    )
end

end
