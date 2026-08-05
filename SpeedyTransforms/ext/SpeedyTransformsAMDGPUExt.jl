module SpeedyTransformsAMDGPUExt

import AMDGPU: AMDGPU, ROCArray, ROCBackend
import AMDGPU.HIP:
    hipStreamBeginCapture, hipStreamEndCapture,
    hipGraphInstantiateWithFlags, hipGraphLaunch,
    hipGraphExecDestroy, hipGraphDestroy,
    hipStreamCaptureModeGlobal, hipGraph_t, hipGraphExec_t,
    hipSuccess, hipGetErrorString, HIPError

using SpeedyTransforms
using SpeedyTransforms.RingGrids
using SpeedyTransforms.LowerTriangularArrays

import SpeedyTransforms: SpectralTransform, GPUFourierGraphCache, build_cache, run_graph!, MAX_GRAPHS, fft_plans

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
    exec_dict() = Dict{UInt, Union{Nothing, hipGraphExec_t}}()
    return GPUFourierGraphCache(
        packed_real, packed_complex, real_view, complex_view,
        rfft_plans, brfft_plans,
        dev(real_offset), dev(complex_offset), dev(nlons), dev(nfreqs),
        dev(istart_n), dev(istart_s), dev(nlons_s),
        S.nlon_max, S.nfreq_max, nlat_half, nlayers, has_equator, j_equator,
        architecture(packed_real),
        exec_dict(),
        exec_dict(),
    )
end

# =====================================================================================
# EXPERIMENTAL, re-enabled 2026-08-05 for testing on LUMI (session investigating whether the
# original crash still reproduces). Previously HIP graph capture was disabled outright: ROCm's
# stream-capture validator was found not to reliably reject illegal-to-capture operations —
# some raised a catchable HIPError (hipErrorStreamCaptureUnsupported, from rocFFT's `mul!`),
# but others were silently accepted at capture time and only surfaced as a GPU memory access
# fault on replay, confirmed on real hardware (LUMI) even after removing the one identifiable
# offending call (see git history: 0c114dfc..c7dc49b1 on gd/hip-graphs, and matching upstream
# gap pytorch/pytorch#155684, #155720).
#
# Standalone repro scripts (scratch_hip_graph_repro*.jl) run against the CURRENT AMDGPU.jl on
# this LUMI environment found NO reproduction: isolated single rocFFT calls (both directions)
# and full transform! calls (Legendre + all-ring Fourier, both directions) all captured cleanly
# and replayed correctly across 5 iterations with fresh data each time. This may mean the
# upstream bug is fixed; it may also just mean these tests aren't the right shape/scale to
# trigger it. This re-enablement is to find out under the ACTUAL model test suite / a real
# run!() -- do not consider this validated until that has been run and reviewed.
#
# Uses raw hip{Stream,Graph}* C bindings rather than a high-level capture/instantiate/launch
# wrapper because AMDGPU.jl 2.1.2 (as resolved on LUMI) does not expose one: `AMDGPU.HIP.capture`
# and `HIPGraphExec` are undefined on this version; only the raw C bindings are present.
# =====================================================================================

raw_stream() = AMDGPU.stream().stream

function hip_ok!(err, label)
    err == hipSuccess && return err
    error("$label failed: $(unsafe_string(hipGetErrorString(err))) (code $err)")
end

# Low-level equivalent of AMDGPU.HIP.capture(f; throw_error=false), built on the raw
# hipStream*Capture* bindings (see module docstring above for why). Returns the captured
# hipGraph_t, or `nothing` if capture was invalidated. A HIPError raised inside f() propagates
# to the caller, which is expected to catch it (mirrors capture(throw_error=false)'s contract
# of only swallowing invalidation, not arbitrary capture errors).
function raw_capture(f)
    stream = raw_stream()
    hip_ok!(hipStreamBeginCapture(stream, hipStreamCaptureModeGlobal), "hipStreamBeginCapture")

    thrown = nothing
    try
        f()
    catch err
        thrown = err
    end

    graph_ref = Ref{hipGraph_t}()
    end_err = hipStreamEndCapture(stream, graph_ref)

    if thrown !== nothing
        rethrow(thrown)
    end

    end_err == hipSuccess || return nothing   # capture invalidated
    return graph_ref[]
end

function raw_instantiate(graph)
    exec_ref = Ref{hipGraphExec_t}()
    hip_ok!(hipGraphInstantiateWithFlags(exec_ref, graph, UInt64(0)), "hipGraphInstantiateWithFlags")
    return exec_ref[]
end

function raw_launch!(exec)
    hip_ok!(hipGraphLaunch(exec, raw_stream()), "hipGraphLaunch")
    return nothing
end

# =====================================================================================
# run_graph!: HIP-specific graph capture and replay.
# Dispatches on the A=GPU{<:ROCBackend} type parameter of GPUFourierGraphCache.
#
# DEBUG INSTRUMENTATION (temporary, 2026-08-05 LUMI crash isolation session): every phase is
# announced with a flushed println tagged by direction/nlayers/buffer key, and AMDGPU.synchronize()
# is inserted after every risky GPU operation (warm-up, replay) to force async faults to surface
# at the point that actually caused them rather than at some later unrelated sync. This makes
# replay slower than the real design intends -- REMOVE before considering this production code.
# =====================================================================================

function run_graph!(
        cache::GPUFourierGraphCache{<:Any, <:Any, <:Any, <:Any, <:Any, <:GPU{<:ROCBackend}, E},
        execs::Dict{UInt, E}, key, loop!::F,
    ) where {E, F}
    direction = execs === cache.forward_execs ? "forward" : (execs === cache.inverse_execs ? "inverse" : "unknown")
    tag = "[$direction nlayers=$(cache.nlayers) key=0x$(string(key, base = 16))]"

    exec = get(execs, key, missing)
    if exec isa hipGraphExec_t
        println("$tag replaying cached graph"); flush(stdout)
        raw_launch!(exec)
        AMDGPU.synchronize()
        println("$tag replay OK"); flush(stdout)
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

    println("$tag first time seen -- warming up outside capture"); flush(stdout)
    # warm up so that one-time work (rocFFT init, kernel JIT, memory-pool growth) happens
    # OUTSIDE the capture region where it is not allowed
    loop!()
    AMDGPU.synchronize()
    println("$tag warm-up OK"); flush(stdout)

    println("$tag beginning capture"); flush(stdout)
    graph = try
        raw_capture() do
            loop!()
        end
    catch err
        err isa HIPError || rethrow()
        println("$tag capture threw HIPError: $err"); flush(stdout)
        nothing
    end

    if graph === nothing
        # capture invalidated or failed with a hard HIP error; warmup already produced
        # the correct result
        println("$tag capture invalidated/failed -- falling back to direct loop for this key"); flush(stdout)
        execs[key] = nothing
        return nothing
    end
    println("$tag capture OK -- instantiating"); flush(stdout)

    exec = raw_instantiate(graph)
    hipGraphDestroy(graph)   # exec is independently allocated; the graph template is not needed after this
    execs[key] = exec
    println("$tag instantiated -- launching first replay"); flush(stdout)
    raw_launch!(exec)
    AMDGPU.synchronize()
    println("$tag first replay OK"); flush(stdout)
    return nothing
end

end
