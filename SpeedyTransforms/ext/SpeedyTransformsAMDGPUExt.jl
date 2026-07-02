module SpeedyTransformsAMDGPUExt

import AMDGPU: AMDGPU, ROCArray, ROCBackend

using KernelAbstractions

using SpeedyTransforms
using SpeedyTransforms.RingGrids
using SpeedyTransforms.LowerTriangularArrays

import SpeedyTransforms: SpectralTransform, GPUFourierGraphCache, build_cache, run_graph!, MAX_GRAPHS, fft_plans

import SpeedyWeatherInternals.Architectures: GPU, architecture

# Probe for the high-level HIP graph API. AMDGPU.HIP exports these in newer versions;
# on older installs only the raw C bindings (hipGraph_t, hipGraphExec_t, …) are present.
const _HIP_GRAPHS_AVAILABLE =
    isdefined(AMDGPU, :HIP) &&
    isdefined(AMDGPU.HIP, :HIPGraphExec) &&
    isdefined(AMDGPU.HIP, :capture)

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
        Dict{UInt, Any}(),  # HIPGraphExec may not exist at compile time, so Any is unavoidable
        Dict{UInt, Any}(),
    )
end

# =====================================================================================
# run_graph!: HIP-specific graph capture and replay.
# Dispatches on the A=GPU{<:ROCBackend} type parameter of GPUFourierGraphCache.
# Falls back to running the allocation-free loop directly when the HIP graph API is
# unavailable (older AMDGPU.jl) or when capture fails.
# =====================================================================================

function run_graph!(
        ::GPUFourierGraphCache{<:Any, <:Any, <:Any, <:Any, <:Any, <:GPU{<:ROCBackend}, E},
        execs::Dict{UInt, E}, key, loop!::F,
    ) where {E, F}
    if !_HIP_GRAPHS_AVAILABLE
        loop!()                              # no graph API: run the fused loop directly
        return nothing
    end

    HIPGraphExecT = AMDGPU.HIP.HIPGraphExec  # only evaluated when _HIP_GRAPHS_AVAILABLE

    exec = get(execs, key, missing)
    if exec isa HIPGraphExecT
        AMDGPU.HIP.launch(exec)              # hot path: pure replay on AMDGPU.stream()
        return nothing
    elseif exec === nothing                  # capture previously failed; run directly
        loop!()
        return nothing
    end

    # first time we see this buffer (exec === missing)
    if length(execs) >= MAX_GRAPHS
        loop!()                              # cache full: don't capture, just run
        return nothing
    end

    # warm up so that one-time work (rocFFT init, kernel JIT, memory-pool growth) happens
    # OUTSIDE the capture region where it is not allowed
    loop!()
    KernelAbstractions.synchronize(ROCBackend())

    # `throw_error = false` only swallows `hipErrorStreamCaptureInvalidated` — capture can
    # also fail with a *different* HIPError (e.g. hipErrorStreamCaptureUnsupported, seen when
    # rocFFT's `mul!` performs an operation that isn't graph-capturable) which `capture` still
    # rethrows. Catch that too: the warm-up run above already produced the correct result for
    # this call, so it's safe to just record capture as failed and fall back for future calls.
    graph = try
        AMDGPU.HIP.capture(throw_error = false) do
            loop!()
        end
    catch err
        err isa AMDGPU.HIP.HIPError || rethrow()
        nothing
    end

    if graph === nothing
        # capture invalidated or failed with a hard HIP error; warmup already produced
        # the correct result
        execs[key] = nothing
        return nothing
    end

    exec = AMDGPU.HIP.instantiate(graph)
    execs[key] = exec
    AMDGPU.HIP.launch(exec)
    return nothing
end

end
