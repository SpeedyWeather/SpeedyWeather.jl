module SpeedyTransformsAMDGPUExt

import AMDGPU: ROCArray

using SpeedyTransforms
using SpeedyTransforms.RingGrids

import SpeedyTransforms: SpectralTransform, _fourier_batched!
import SpeedyTransforms.RingGrids: AbstractField

# =====================================================================================
# HIP GRAPHS ACCELERATION OF THE BATCHED FOURIER TRANSFORM — NOT IMPLEMENTED
#
# CUDA graphs give the batched Fourier transform a large speedup on NVIDIA GPUs (see
# SpeedyTransformsCUDAExt.jl) by capturing the fused gather/FFT/scatter loop once and
# replaying it on later calls. An equivalent HIP-graphs implementation for AMDGPU was
# attempted (git history: 0c114dfc..e4a19ac7 on gd/hip-graphs) but had to be abandoned:
# ROCm's stream-capture validator does not reliably reject operations that are illegal to
# capture — some raise a catchable HIPError, but others are silently accepted at capture
# time and only surface later as a GPU memory access fault when the (corrupt) graph is
# replayed. This was confirmed on real hardware (LUMI) twice, including a 2026-08-05
# session that re-enabled capture on raw hip{Stream,Graph}* C bindings specifically to
# check whether the bug still reproduces — it does. Since capture failures aren't reliably
# catchable, "try capture and fall back on error" isn't safe, so this optimisation is not
# offered for AMDGPU: the generic (unaccelerated but correct) per-ring Fourier path in
# SpeedyTransforms/src/fourier.jl is used unconditionally, with a one-time warning if the
# user asked for `gpu_graphs = true`.
# =====================================================================================

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
