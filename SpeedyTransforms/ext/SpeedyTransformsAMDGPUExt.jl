module SpeedyTransformsAMDGPUExt

import AMDGPU: AMDGPU, ROCArray, ROCBackend

using SpeedyTransforms
using SpeedyTransforms.RingGrids
using SpeedyTransforms.LowerTriangularArrays

import SpeedyTransforms: SpectralTransform, GPUFourierGraphCache, build_cache, run_graph!, fft_plans

import SpeedyWeatherInternals.Architectures: GPU, architecture

# =====================================================================================
# build_cache: allocate the packed work buffers and per-ring views on AMDGPU device memory.
# =====================================================================================

function build_cache(S::SpectralTransform, nlayers::Integer, ::GPU{<:ROCBackend})
    NF = eltype(S)
    nlat = S.nlat
    nlat_half = S.grid.nlat_half
    rfft_plans, brfft_plans = fft_plans(S, nlayers)
    rings = S.rings
    nlons = collect(Int, S.nlons[1:nlat_half])
    nfreqs = nlons .÷ 2 .+ 1

    block_real = nlons .* nlayers
    block_complex = nfreqs .* nlayers
    real_offset = [0; cumsum(block_real)[1:(end - 1)]]
    complex_offset = [0; cumsum(block_complex)[1:(end - 1)]]

    istart_n = [rings[j].start for j in 1:nlat_half]
    istart_s = [rings[nlat - j + 1].start for j in 1:nlat_half]

    has_equator = isodd(nlat)
    j_equator = (nlat + 1) ÷ 2
    nlons_s = copy(nlons)
    has_equator && (nlons_s[j_equator] = 0)

    packed_real = AMDGPU.zeros(NF, sum(block_real))
    packed_complex = AMDGPU.zeros(Complex{NF}, sum(block_complex))
    real_view = [reshape(view(packed_real, real_offset[j] + 1:real_offset[j] + block_real[j]), nlons[j], nlayers) for j in 1:nlat_half]
    complex_view = [reshape(view(packed_complex, complex_offset[j] + 1:complex_offset[j] + block_complex[j]), nfreqs[j], nlayers) for j in 1:nlat_half]

    dev(x) = ROCArray(x)
    return GPUFourierGraphCache(
        packed_real, packed_complex, real_view, complex_view,
        rfft_plans, brfft_plans,
        dev(real_offset), dev(complex_offset), dev(nlons), dev(nfreqs),
        dev(istart_n), dev(istart_s), dev(nlons_s),
        S.nlon_max, S.nfreq_max, nlat_half, nlayers, has_equator, j_equator,
        architecture(packed_real),
        Dict{UInt, Nothing}(),  # unused: run_graph! never captures on AMDGPU, see below
        Dict{UInt, Nothing}(),
    )
end

# =====================================================================================
# run_graph!: HIP graph capture is currently disabled for AMDGPU. ROCm's stream-capture
# validator does not reliably reject operations that are illegal to capture — some raise a
# catchable HIPError (e.g. hipErrorStreamCaptureUnsupported, seen from rocFFT's `mul!`), but
# others are silently accepted at capture time and only surface later as a GPU memory access
# fault when the (corrupt) graph is replayed. That failure mode was confirmed on real hardware
# (LUMI) even after removing the one identifiable offending call, and matches a known upstream
# gap: HIP Graph capture on AMD GPUs does not raise `operation not permitted` for illegal
# operations the way CUDA Graph capture does (see pytorch/pytorch#155684, #155720). Since
# capture failures aren't reliably catchable, "try capture and fall back on error" isn't safe
# here — always run the allocation-free direct loop instead. Revisit once ROCm/AMDGPU.jl
# provides reliable stream-capture validation for library calls like rocFFT execution.
# =====================================================================================

function run_graph!(
        ::GPUFourierGraphCache{<:Any, <:Any, <:Any, <:Any, <:Any, <:GPU{<:ROCBackend}, E},
        ::Dict{UInt, E}, key, loop!::F,
    ) where {E, F}
    loop!()
    return nothing
end

end
