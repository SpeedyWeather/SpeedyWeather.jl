module SpeedyTransformsAMDGPUExt

import AMDGPU: AMDGPU, ROCArray, ROCBackend
import AMDGPU.HIP:
    hipStreamBeginCapture, hipStreamEndCapture,
    hipGraphInstantiateWithFlags, hipGraphLaunch,
    hipGraphExecDestroy, hipGraphDestroy,
    hipStreamCaptureModeGlobal, hipGraph_t, hipGraphExec_t,
    hipSuccess, hipGetErrorString, HIPError
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
import SpeedyWeatherInternals.Architectures: GPU, architecture

# =====================================================================================
# HIP GRAPHS ACCELERATION OF THE BATCHED FOURIER TRANSFORM
#
# This is the AMDGPU equivalent of SpeedyTransformsCUDAExt.jl — mirrors it structurally
# (kernels, cache struct, allocation-free loops, capture/replay bookkeeping) but is a
# self-contained, independent copy rather than sharing code via the main package. That
# duplication is deliberate: an earlier version of this file shared the gather/scatter
# kernels, GPUFourierGraphCache struct, and forward/inverse loops with the CUDA extension via
# SpeedyTransforms/src/fourier_gpu_graphs.jl. That worked, but because that file lived in the
# main package (not an extension), it defined extra methods of _fourier_batched!/fft_plans
# that exist in EVERY Julia session regardless of which (if any) GPU extension is loaded —
# even a plain CPU session. Upstream's extensive Enzyme differentiability test suite turned
# out to be sensitive to the total method table of these generic functions: the extra
# GPU-only methods (never dynamically reachable for CPU array types) still confused Enzyme's
# static activity analysis, causing a `roots_activep != activep` internal assertion and a
# hard crash in CPU-only AD tests that never touch a GPU. Keeping each backend's code
# entirely inside its own extension (as upstream's CUDA extension already does) means none of
# this exists unless that specific GPU package is actually loaded, sidestepping the issue.
#
# See SpeedyTransformsCUDAExt.jl for the full explanation of why graphs are needed at all; the
# design here is otherwise identical modulo the backend-specific capture/instantiate/replay
# API (HIP vs CUDA).
# =====================================================================================

"""Maximum number of cached graphs per direction per `SpectralTransform`. Prevents
unbounded growth (and host-side capture cost) when the transform is called with a stream
of freshly-allocated `field` buffers (e.g. the allocating `transform(field, S)`). Beyond
this many distinct buffers the allocation-free loop is run directly without capturing."""
const MAX_GRAPHS = 64

# =====================================================================================
# Fused gather/scatter kernels — move all rings at once between the strided grid/scratch
# layout and the contiguous packed work buffer (one launch instead of one copy per ring).
# Each kernel is launched over (nlon_max | nfreq_max, nlat_half, nlayers) and masks the
# threads beyond a ring's actual length.
# =====================================================================================

# grid field (npoints × nlayers) ring rows  ->  packed real buffer  (forward gather)
@kernel inbounds = true function gather_real_kernel!(packed, src, real_offset, nlons, istart)
    i, j, k = @index(Global, NTuple)
    nlon = nlons[j]
    if i <= nlon
        packed[real_offset[j] + (k - 1) * nlon + i] = src[istart[j] + i - 1, k]
    end
end

# packed complex buffer  ->  scratch (nfreq_max × nlayers × nlat_half)  (forward scatter)
@kernel inbounds = true function scatter_complex_kernel!(dst, packed, complex_offset, nfreqs)
    m, j, k = @index(Global, NTuple)
    nfreq = nfreqs[j]
    if m <= nfreq
        dst[m, k, j] = packed[complex_offset[j] + (k - 1) * nfreq + m]
    end
end

# scratch  ->  packed complex buffer  (inverse gather)
@kernel inbounds = true function gather_complex_kernel!(packed, src, complex_offset, nfreqs)
    m, j, k = @index(Global, NTuple)
    nfreq = nfreqs[j]
    if m <= nfreq
        packed[complex_offset[j] + (k - 1) * nfreq + m] = src[m, k, j]
    end
end

# packed real buffer  ->  grid field ring rows  (inverse scatter)
# `add` selects overwrite (`=`, the plain transform) vs accumulate (`+=`, used by the Enzyme adjoint
# rules). It is a plain `Bool` kernel argument rather than a type parameter: the branch is uniform
# across all threads (no divergence) and graphs are captured per `add` mode anyway (see
# `inverse_execs`), so a second specialisation would buy nothing.
@kernel inbounds = true function scatter_real_kernel!(dst, packed, real_offset, nlons, istart, add)
    i, j, k = @index(Global, NTuple)
    nlon = nlons[j]
    if i <= nlon
        i_dst = istart[j] + i - 1
        val = packed[real_offset[j] + (k - 1) * nlon + i]
        if add
            dst[i_dst, k] += val
        else
            dst[i_dst, k] = val
        end
    end
end

# =====================================================================================
# Per-SpectralTransform cache: contiguous packed work buffers, reshaped per-ring views
# for the FFTs, device-side gather/scatter metadata, and the instantiated HIP graphs.
# =====================================================================================

"""$(TYPEDSIGNATURES)
Cache (one per transform SIZE, i.e. per FFT plan) holding the pre-allocated packed
work buffers, the per-ring reshaped views and FFT plans used by the transforms, the device
gather/scatter index metadata, and the instantiated HIP graphs (one per distinct `field`
buffer and direction). A graph value of `nothing` marks a buffer for which capture failed
(fall back to direct loop)."""
struct GPUFourierGraphCache{PR, PC, RV, CV, RFP, BFP, IV, A}
    packed_real::PR             # ROCVector{NF}          — all rings' dense real blocks
    packed_complex::PC          # ROCVector{Complex{NF}} — all rings' dense complex blocks
    real_view::RV               # Vector of per-ring reshaped (nlon_j × nlayers) views into packed_real
    complex_view::CV            # Vector of per-ring reshaped (nfreq_j × nlayers) views into packed_complex
    rfft_plans::RFP             # forward FFT plans for THIS size (nlayers batches), per ring
    brfft_plans::BFP            # inverse FFT plans for THIS size, per ring
    real_offset::IV             # 0-based real-block offset per ring
    complex_offset::IV          # 0-based complex-block offset per ring
    nlons::IV                   # longitudes per ring
    nfreqs::IV                  # Fourier frequencies per ring
    istart_n::IV                # first grid row of each northern ring
    istart_s::IV                # first grid row of each southern ring
    nlons_s::IV                 # like nlons but 0 at the equator ring (south skip)
    nlon_max::Int               # maximum number of longitudes across all rings (for launch size)
    nfreq_max::Int              # maximum number of frequencies across all rings (for launch size)
    nlat_half::Int              # number of rings on one hemisphere, equator included
    nlayers::Int                # number of vertical layers (for launch size)
    has_equator::Bool           # whether the grid has a ring on the equator that needs special handling
    j_equator::Int              # latitude index of the equator ring (if any)
    arch::A                     # SpeedyWeather GPU architecture (for launch!)
    forward_execs::Dict{UInt, Union{Nothing, hipGraphExec_t}} # forward graphs, keyed by field buffer
    # inverse graphs, keyed by (field buffer, `add`): a captured graph bakes in whether the scatter
    # overwrites or accumulates, so the two modes must never share a graph.
    inverse_execs::Dict{Tuple{UInt, Bool}, Union{Nothing, hipGraphExec_t}}
end

# One cache per (SpectralTransform, transform size), keyed by the *forward FFT plan set*.
const GRAPH_CACHES = IdDict{Any, GPUFourierGraphCache}()

# Return the (forward, inverse) per-ring FFT plan sets for a transform of `nlayers` layers.
fft_plans(S::SpectralTransform, nlayers::Integer) = (S.rfft_plans_batched[nlayers], S.brfft_plans_batched[nlayers])

# build/allocate the cache for a transform of `nlayers` layers
function build_cache(S::SpectralTransform, nlayers::Integer)
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
        Dict{UInt, Union{Nothing, hipGraphExec_t}}(),
        Dict{Tuple{UInt, Bool}, Union{Nothing, hipGraphExec_t}}(),
    )
end

"""$(TYPEDSIGNATURES)
The per-size resource that identifies a graph cache. The forward FFT plan set is unique per
transform size and stable for the lifetime of `S`. Must return *the stored*
plan object (not a fresh allocation) so the identity is stable across calls."""
cache_key(S::SpectralTransform, nlayers::Integer) = fft_plans(S, nlayers)[1]

get_cache(S::SpectralTransform, nlayers::Integer) =
    get!(() -> build_cache(S, nlayers), GRAPH_CACHES, cache_key(S, nlayers))::GPUFourierGraphCache

"""$(TYPEDSIGNATURES)
Clear all cached HIP-graphs Fourier buffers and graphs (frees the associated GPU memory).
Mainly useful for tests/benchmarks."""
clear_fourier_graph_cache!() = (empty!(GRAPH_CACHES); nothing)

# =====================================================================================
# Allocation-free fused loops (capturable). Identical in structure to the CUDA extension's —
# they only use launch!, mul!, and fill! which are all device-agnostic.
# =====================================================================================

"""$(TYPEDSIGNATURES)
Allocation-free, fused forward (grid → spectral) batched Fourier loop: one gather kernel
packs all rings, per-ring in-place FFTs run on reshaped views, one scatter kernel writes
the scratch. Suitable for HIP-graph capture."""
function forward_loop!(cache::GPUFourierGraphCache, f_north, f_south, field::AbstractField, S::SpectralTransform)
    (; arch) = cache
    (; nlat_half, nlayers) = cache
    (; nlon_max, nfreq_max, nlons, nfreqs, j_equator) = cache
    (; real_view, complex_view, packed_real, packed_complex, rfft_plans) = cache
    (; real_offset, complex_offset, istart_n, istart_s) = cache

    real_size = (nlon_max, nlat_half, nlayers)
    complex_size = (nfreq_max, nlat_half, nlayers)

    # northern rings
    launch!(arch, ArrayWorkOrder, real_size, gather_real_kernel!, packed_real, field.data, real_offset, nlons, istart_n)
    @inbounds for j in 1:nlat_half
        mul!(complex_view[j], rfft_plans[j], real_view[j])
    end
    launch!(arch, ArrayWorkOrder, complex_size, scatter_complex_kernel!, f_north, packed_complex, complex_offset, nfreqs)

    # southern rings (the equator ring, if any, is zeroed rather than transformed)
    launch!(arch, ArrayWorkOrder, real_size, gather_real_kernel!, packed_real, field.data, real_offset, nlons, istart_s)
    @inbounds for j in 1:nlat_half
        if cache.has_equator && j == j_equator
            fill!(complex_view[j], 0)
        else
            mul!(complex_view[j], rfft_plans[j], real_view[j])
        end
    end
    launch!(arch, ArrayWorkOrder, complex_size, scatter_complex_kernel!, f_south, packed_complex, complex_offset, nfreqs)
    return nothing
end

"""$(TYPEDSIGNATURES)
Allocation-free, fused inverse (spectral → grid) batched Fourier loop: one gather kernel
packs all rings from the scratch, per-ring in-place inverse FFTs run on reshaped views, one
scatter kernel writes the grid field. `add` accumulates onto `field` instead of overwriting it
(Enzyme adjoint rule); it is baked into the captured graph, hence the per-`add` graph cache
key. Suitable for HIP-graph capture."""
function inverse_loop!(cache::GPUFourierGraphCache, field::AbstractField, g_north, g_south, S::SpectralTransform, add::Bool = false)
    (; arch) = cache
    (; nlat_half, nlayers) = cache
    (; nlon_max, nfreq_max, nlons, nlons_s, nfreqs, j_equator) = cache
    (; real_view, complex_view, packed_real, packed_complex, brfft_plans) = cache
    (; real_offset, complex_offset, istart_n, istart_s) = cache

    real_size = (nlon_max, nlat_half, nlayers)
    complex_size = (nfreq_max, nlat_half, nlayers)

    # northern rings
    launch!(arch, ArrayWorkOrder, complex_size, gather_complex_kernel!, packed_complex, g_north, complex_offset, nfreqs)
    @inbounds for j in 1:nlat_half
        mul!(real_view[j], brfft_plans[j], complex_view[j])
    end
    launch!(arch, ArrayWorkOrder, real_size, scatter_real_kernel!, field.data, packed_real, real_offset, nlons, istart_n, add)

    # southern rings (the equator ring, if any, is skipped: north already wrote those rows)
    launch!(arch, ArrayWorkOrder, complex_size, gather_complex_kernel!, packed_complex, g_south, complex_offset, nfreqs)
    @inbounds for j in 1:nlat_half
        (cache.has_equator && j == j_equator) && continue
        mul!(real_view[j], brfft_plans[j], complex_view[j])
    end
    launch!(arch, ArrayWorkOrder, real_size, scatter_real_kernel!, field.data, packed_real, real_offset, nlons_s, istart_s, add)
    return nothing
end

# =====================================================================================
# Capture / replay management
#
# EXPERIMENTAL, re-enabled 2026-08-05 for testing on LUMI (session investigating whether the
# original crash — a GPU memory access fault confirmed on real hardware, see git history:
# 0c114dfc..c7dc49b1 on gd/hip-graphs — still reproduces). Uses raw hip{Stream,Graph}* C
# bindings rather than a high-level capture/instantiate/launch wrapper because AMDGPU.jl 2.1.2
# (as resolved on LUMI) does not expose one: `AMDGPU.HIP.capture` and `HIPGraphExec` are
# undefined on this version; only the raw C bindings are present.
#
# DEBUG INSTRUMENTATION (temporary, same session): every phase is announced with a flushed
# println tagged by direction/nlayers/buffer key, and AMDGPU.synchronize() is inserted after
# every risky GPU operation (warm-up, replay) to force async faults to surface at the point
# that actually caused them rather than at some later unrelated sync. This makes replay slower
# than the real design intends -- REMOVE before considering this production code.
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

function run_graph!(cache::GPUFourierGraphCache, execs::AbstractDict, key, loop!::F) where {F}
    direction = execs === cache.forward_execs ? "forward" : (execs === cache.inverse_execs ? "inverse" : "unknown")
    # key is a UInt (forward) or a (UInt, Bool) tuple (inverse, `add` baked into the key)
    ptr, add = key isa Tuple ? key : (key, nothing)
    add_tag = add === nothing ? "" : " add=$add"
    tag = "[$direction nlayers=$(cache.nlayers) key=0x$(string(ptr, base = 16))$add_tag]"

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

# Stable per-buffer graph-cache key: the device address of `field.data`. See CUDA extension's
# graph_key for the full explanation (per-step view wrappers churn but the device pointer they
# alias is stable for the lifetime of the model buffer).
@inline graph_key(data) = reinterpret(UInt, pointer(data))

# =====================================================================================
# Method overrides: dispatch on ROCArray scratch (more specific than the generic
# AbstractArray{<:Complex,3} methods in fourier.jl)
# =====================================================================================

"""$(TYPEDSIGNATURES)
HIP-Graphs accelerated forward (grid → spectral) batched Fourier transform.
Replays a cached HIP graph of the fused gather + per-ring FFTs + scatter; see
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
    cache = get_cache(S, size(field, 2))
    run_graph!(cache, cache.forward_execs, graph_key(field.data), () -> forward_loop!(cache, f_north, f_south, field, S))
    return nothing
end

"""$(TYPEDSIGNATURES)
HIP-Graphs accelerated inverse (spectral → grid) batched Fourier transform.
Replays a cached HIP graph of the fused gather + per-ring inverse FFTs + scatter; see
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
    cache = get_cache(S, size(field, 2))
    # `add` is part of the graph key: a captured graph bakes in overwrite-vs-accumulate, so replaying
    # the wrong one would silently produce incorrect results.
    run_graph!(cache, cache.inverse_execs, (graph_key(field.data), add), () -> inverse_loop!(cache, field, g_north, g_south, S, add))
    return nothing
end

end
