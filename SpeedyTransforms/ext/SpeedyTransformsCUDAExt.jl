module SpeedyTransformsCUDAExt

import CUDA: CUDA, CUFFT, CuArray, CuGraphExec, capture, instantiate, launch
import AbstractFFTs
import LinearAlgebra
import LinearAlgebra: mul!
using DocStringExtensions

using SpeedyTransforms
using SpeedyTransforms.RingGrids
using SpeedyTransforms.LowerTriangularArrays

import SpeedyTransforms: SpectralTransform, _fourier_batched!
import SpeedyTransforms.RingGrids: AbstractField

# =====================================================================================
# CUDA GRAPHS ACCELERATION OF THE BATCHED FOURIER TRANSFORM
#
# The batched Fourier transform on a (reduced) grid applies one cuFFT per latitude ring
# (≈ 2*nlat_half tiny FFTs). On the GPU this is heavily launch-bound: the CPU spends far
# more time enqueuing the many small kernels than the GPU spends computing them.
#
# CUDA Graphs let us record this whole sequence ONCE and replay it with a single
# `launch`, eliminating the per-operation CPU launch overhead. To make the captured
# region replay correctly we must avoid device allocations inside it (the generic GPU
# path allocates a temporary per ring via `field.data[ilons, :]` and `plan * x`), so this
# implementation pre-allocates dense per-ring buffers once and uses in-place `mul!`.
#
# The captured graph bakes in the device pointers of `field.data`, the scratch buffers
# (`S.scratch_memory.north/.south`, stable for the lifetime of `S`) and the per-ring
# work buffers. In the SpeedyWeather time loop the same variable buffers are reused every
# timestep, so a graph captured for a given `field` is replayed on all subsequent steps.
# Graphs are cached per `SpectralTransform` and keyed by the `field.data` object.
# =====================================================================================

"""Toggle for the CUDA-Graphs accelerated batched Fourier transform. Set to `false` to
fall back to the generic (allocating) per-ring GPU path, e.g. for benchmarking."""
const FOURIER_GRAPHS_ENABLED = Ref(true)

"""Maximum number of cached graphs per direction per `SpectralTransform`. Prevents
unbounded growth (and host-side capture cost) when the transform is called with a stream
of freshly-allocated `field` buffers (e.g. the allocating `transform(field, S)`). Beyond
this many distinct buffers the allocation-free loop is run directly without capturing."""
const MAX_GRAPHS = 128

"""$(TYPEDSIGNATURES)
Per-`SpectralTransform` cache holding the pre-allocated dense per-ring work buffers and
the instantiated CUDA graphs (one per distinct `field` buffer and direction). A graph
value of `nothing` marks a buffer for which capture failed (fall back to direct loop)."""
mutable struct GPUFourierGraphCache{RealBuffer <: CuArray, ComplexBuffer <: CuArray}
    real_buf::Vector{RealBuffer}        # per ring: (nlon_j  × nlayers) real, FFT in/out
    cplx_buf::Vector{ComplexBuffer}     # per ring: (nfreq_j × nlayers) complex, FFT out/in
    forward_execs::IdDict{Any, Union{Nothing, CuGraphExec}}
    inverse_execs::IdDict{Any, Union{Nothing, CuGraphExec}}
end

# one cache per SpectralTransform, keyed by the (stable, unique) north scratch array
const GRAPH_CACHES = IdDict{Any, GPUFourierGraphCache}()

function build_cache(S::SpectralTransform)
    NF = eltype(S)
    nlat_half = S.grid.nlat_half
    nlayers = S.nlayers
    real_buf = [CUDA.zeros(NF, S.nlons[j], nlayers) for j in 1:nlat_half]
    cplx_buf = [CUDA.zeros(Complex{NF}, S.nlons[j] ÷ 2 + 1, nlayers) for j in 1:nlat_half]
    return GPUFourierGraphCache(
        real_buf, cplx_buf,
        IdDict{Any, Union{Nothing, CuGraphExec}}(),
        IdDict{Any, Union{Nothing, CuGraphExec}}(),
    )
end

get_cache(S::SpectralTransform) = get!(() -> build_cache(S), GRAPH_CACHES, S.scratch_memory.north)::GPUFourierGraphCache

"""$(TYPEDSIGNATURES)
Clear all cached CUDA-Graphs Fourier buffers and graphs (frees the associated GPU memory).
Mainly useful for tests/benchmarks."""
clear_fourier_graph_cache!() = (empty!(GRAPH_CACHES); nothing)

# =====================================================================================
# Allocation-free per-ring loops (capturable). These reproduce exactly the regions that
# the generic `_fourier_batched!` writes (rows 1:nfreq of each ring's scratch slice;
# the full ring rows of `field` for the inverse).
# =====================================================================================

"""$(TYPEDSIGNATURES)
Allocation-free forward (grid → spectral) batched Fourier loop using pre-allocated dense
per-ring buffers and in-place `mul!`. Suitable for CUDA-graph capture."""
function forward_loop!(cache::GPUFourierGraphCache, f_north, f_south, field::AbstractField, S::SpectralTransform)
    (; nlat, nlons, rings) = S
    nlat_half = S.grid.nlat_half
    real_buf = cache.real_buf
    cplx_buf = cache.cplx_buf
    @inbounds for j_north in 1:nlat_half
        j = j_north
        j_south = nlat - j_north + 1
        nlon = nlons[j]
        nfreq = nlon ÷ 2 + 1
        not_equator = j_north != j_south

        # northern ring (always)
        ilons = rings[j_north]
        copyto!(real_buf[j], view(field.data, ilons, :))
        mul!(cplx_buf[j], S.rfft_plans[j], real_buf[j])
        copyto!(view(f_north, 1:nfreq, :, j), cplx_buf[j])

        # southern ring (skip redundant FFT on the equator)
        if not_equator
            ilons = rings[j_south]
            copyto!(real_buf[j], view(field.data, ilons, :))
            mul!(cplx_buf[j], S.rfft_plans[j], real_buf[j])
            copyto!(view(f_south, 1:nfreq, :, j), cplx_buf[j])
        else
            fill!(view(f_south, 1:nfreq, :, j), 0)
        end
    end
    return nothing
end

"""$(TYPEDSIGNATURES)
Allocation-free inverse (spectral → grid) batched Fourier loop using pre-allocated dense
per-ring buffers and in-place `mul!`. Suitable for CUDA-graph capture."""
function inverse_loop!(cache::GPUFourierGraphCache, field::AbstractField, g_north, g_south, S::SpectralTransform)
    (; nlat, nlons, rings) = S
    nlat_half = S.grid.nlat_half
    real_buf = cache.real_buf
    cplx_buf = cache.cplx_buf
    @inbounds for j_north in 1:nlat_half
        j = j_north
        j_south = nlat - j_north + 1
        nlon = nlons[j]
        nfreq = nlon ÷ 2 + 1
        not_equator = j_north != j_south

        # northern ring (always)
        ilons = rings[j_north]
        copyto!(cplx_buf[j], view(g_north, 1:nfreq, :, j))
        mul!(real_buf[j], S.brfft_plans[j], cplx_buf[j])
        copyto!(view(field.data, ilons, :), real_buf[j])

        # southern ring (skip redundant FFT on the equator)
        if not_equator
            ilons = rings[j_south]
            copyto!(cplx_buf[j], view(g_south, 1:nfreq, :, j))
            mul!(real_buf[j], S.brfft_plans[j], cplx_buf[j])
            copyto!(view(field.data, ilons, :), real_buf[j])
        end
    end
    return nothing
end

# =====================================================================================
# Capture / replay management
# =====================================================================================

"""$(TYPEDSIGNATURES)
Run `loop!` (a `() -> ...` closure over the allocation-free batched FFT) either by
replaying a cached CUDA graph keyed by `key`, or — on first use — by warming up, capturing,
instantiating and caching a graph. Falls back to running `loop!` directly if capture fails
or the per-direction cache is full."""
function run_graph!(execs::IdDict, key, loop!::F) where {F}
    exec = get(execs, key, missing)
    if exec isa CuGraphExec
        launch(exec)                         # hot path: pure replay
        return nothing
    elseif exec === nothing
        loop!()                              # known non-capturable buffer
        return nothing
    end

    # first time we see this buffer (exec === missing)
    if length(execs) >= MAX_GRAPHS
        loop!()                              # cache full: don't capture, just run
        return nothing
    end

    # warm up so that one-time work (cuFFT init, kernel JIT, memory-pool growth) happens
    # OUTSIDE the capture, where it is allowed
    loop!()
    CUDA.synchronize()

    graph = capture(throw_error = false) do
        loop!()
    end

    if graph === nothing
        # capture invalidated (e.g. an unexpected allocation/sync); the warmup already
        # produced the correct result. Remember not to retry capture for this buffer.
        execs[key] = nothing
        return nothing
    end

    exec = instantiate(graph)
    execs[key] = exec
    launch(exec)                             # produce the result via the graph
    return nothing
end

# =====================================================================================
# Method overrides: dispatch on CuArray scratch (more specific than the generic
# AbstractArray{<:Complex,3} methods in fourier.jl)
# =====================================================================================

"""$(TYPEDSIGNATURES)
CUDA-Graphs accelerated forward (grid → spectral) batched Fourier transform.
Replays a cached CUDA graph of the per-ring cuFFTs; see [`run_graph!`](@ref)."""
function _fourier_batched!(
        f_north::CuArray{<:Complex, 3},
        f_south::CuArray{<:Complex, 3},
        field::AbstractField,
        S::SpectralTransform,
    )
    @assert eltype(field) == eltype(S) "Number format of grid $(eltype(field)) and SpectralTransform $(eltype(S)) need to match."
    if !FOURIER_GRAPHS_ENABLED[]
        return Base.@invoke _fourier_batched!(
            f_north::AbstractArray{<:Complex, 3}, f_south::AbstractArray{<:Complex, 3},
            field::AbstractField, S::SpectralTransform,
        )
    end
    cache = get_cache(S)
    run_graph!(cache.forward_execs, field.data, () -> forward_loop!(cache, f_north, f_south, field, S))
    return nothing
end

"""$(TYPEDSIGNATURES)
CUDA-Graphs accelerated inverse (spectral → grid) batched Fourier transform.
Replays a cached CUDA graph of the per-ring inverse cuFFTs; see [`run_graph!`](@ref)."""
function _fourier_batched!(
        field::AbstractField,
        g_north::CuArray{<:Complex, 3},
        g_south::CuArray{<:Complex, 3},
        S::SpectralTransform,
    )
    if !FOURIER_GRAPHS_ENABLED[]
        return Base.@invoke _fourier_batched!(
            field::AbstractField, g_north::AbstractArray{<:Complex, 3},
            g_south::AbstractArray{<:Complex, 3}, S::SpectralTransform,
        )
    end
    cache = get_cache(S)
    run_graph!(cache.inverse_execs, field.data, () -> inverse_loop!(cache, field, g_north, g_south, S))
    return nothing
end

end
