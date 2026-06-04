module SpeedyTransformsCUDAExt

import CUDA: CUDA, CUFFT, CuArray, CuVector, CuGraphExec, capture, instantiate, launch
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

import SpeedyWeatherInternals.KernelLaunching: launch!, Array3DWorkOrder
import SpeedyWeatherInternals.Architectures: architecture

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
# several layer counts gets one cache each); within a cache, graphs are keyed by `field.data`.
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
@kernel inbounds = true function gather_real_kernel!(packed, src, roff, nlon, rstart)
    r, j, k = @index(Global, NTuple)
    nl = nlon[j]
    if r <= nl
        packed[roff[j] + (k - 1) * nl + r] = src[rstart[j] + r - 1, k]
    end
end

# packed complex buffer  ->  scratch (nfreq_max × nlayers × nlat_half)  (forward scatter)
@kernel inbounds = true function scatter_cplx_kernel!(dst, packed, coff, nfreq)
    m, j, k = @index(Global, NTuple)
    nf = nfreq[j]
    if m <= nf
        dst[m, k, j] = packed[coff[j] + (k - 1) * nf + m]
    end
end

# scratch  ->  packed complex buffer  (inverse gather)
@kernel inbounds = true function gather_cplx_kernel!(packed, src, coff, nfreq)
    m, j, k = @index(Global, NTuple)
    nf = nfreq[j]
    if m <= nf
        packed[coff[j] + (k - 1) * nf + m] = src[m, k, j]
    end
end

# packed real buffer  ->  grid field ring rows  (inverse scatter)
@kernel inbounds = true function scatter_real_kernel!(dst, packed, roff, nlon, rstart)
    r, j, k = @index(Global, NTuple)
    nl = nlon[j]
    if r <= nl
        dst[rstart[j] + r - 1, k] = packed[roff[j] + (k - 1) * nl + r]
    end
end

# =====================================================================================
# Per-SpectralTransform cache: contiguous packed work buffers, reshaped per-ring views
# for the FFTs, device-side gather/scatter metadata, and the instantiated CUDA graphs.
# =====================================================================================

"""$(TYPEDSIGNATURES)
Cache (one per transform SIZE, i.e. per FFT plan) holding the pre-allocated packed
work buffers, the per-ring reshaped views and FFT plans used by the transforms, the device
gather/scatter index metadata, and the instantiated CUDA graphs (one per distinct `field`
buffer and direction). A graph value of `nothing` marks a buffer for which capture failed
(fall back to direct loop)."""
struct GPUFourierGraphCache{PR, PC, RV, CV, IV, A}
    packed_real::PR             # CuVector{NF}          — all rings' dense real blocks
    packed_cplx::PC             # CuVector{Complex{NF}} — all rings' dense complex blocks
    real_view::RV               # Vector of per-ring reshaped (nlon_j × nlayers) views into packed_real
    cplx_view::CV               # Vector of per-ring reshaped (nfreq_j × nlayers) views into packed_cplx
    rfft_plans::Vector{AbstractFFTs.Plan}   # forward FFT plans for THIS size (nlayers batches), per ring
    brfft_plans::Vector{AbstractFFTs.Plan}  # inverse FFT plans for THIS size, per ring
    roff::IV                    # 0-based real-block offset per ring
    coff::IV                    # 0-based complex-block offset per ring
    nlons::IV                   # longitudes per ring
    nfreq::IV                   # Fourier frequencies per ring
    rstart_n::IV                # first grid row of each northern ring
    rstart_s::IV                # first grid row of each southern ring
    nlons_s::IV                 # like nlons but 0 at the equator ring (south skip)
    nlon_max::Int
    nfreq_max::Int
    nlat_half::Int
    nlayers::Int
    has_equator::Bool
    jeq::Int
    arch::A                     # SpeedyWeather GPU architecture (for launch!)
    forward_execs::IdDict{Any, Union{Nothing, CuGraphExec}} # the actual graphs forward transforms
    inverse_execs::IdDict{Any, Union{Nothing, CuGraphExec}} # the actual graphs inverse transforms
end

# One cache per (SpectralTransform, transform size), keyed by the *forward FFT plan set*
# (`fft_plans(S, nlayers)[1]`)
# this is saved here directly as a global variable and not as a field of `S` because 
# it would be very hard to make `S` a concrete type again otherwise, and this seems to work
# without problems and performance mali
const GRAPH_CACHES = IdDict{Any, GPUFourierGraphCache}()

# EXTENSION POINT for batching: return the
# (forward, inverse) per-ring FFT plan sets for a transform of `nlayers` layers. The default
# assumes a single size; override this for an `S` that stores several plan sets (e.g. keyed
# by layer count) and the rest — caches, packing, capture, replay — follows automatically.
# gets called by `cache_key` when accessing the different caches and graphs
fft_plans(S::SpectralTransform, nlayers::Integer) = (S.rfft_plans, S.brfft_plans)

# build/allocate the cache for a transform of `nlayers` layers
function build_cache(S::SpectralTransform, nlayers::Integer)
    NF = eltype(S)
    nlat = S.nlat
    nlat_half = S.grid.nlat_half
    rfft_plans, brfft_plans = fft_plans(S, nlayers)
    rings = S.rings
    nlons = collect(Int, S.nlons[1:nlat_half])
    nfreqs = nlons .÷ 2 .+ 1

    block_r = nlons .* nlayers
    block_c = nfreqs .* nlayers
    roff = [0; cumsum(block_r)[1:(end - 1)]]
    coff = [0; cumsum(block_c)[1:(end - 1)]]

    rstart_n = [rings[j].start for j in 1:nlat_half]
    rstart_s = [rings[nlat - j + 1].start for j in 1:nlat_half]

    has_equator = isodd(nlat)
    jeq = (nlat + 1) ÷ 2
    nlons_s = copy(nlons)
    has_equator && (nlons_s[jeq] = 0)        # south pass skips the equator ring

    packed_real = CUDA.zeros(NF, sum(block_r))
    packed_cplx = CUDA.zeros(Complex{NF}, sum(block_c))
    real_view = [reshape(view(packed_real, roff[j] + 1:roff[j] + block_r[j]), nlons[j], nlayers) for j in 1:nlat_half]
    cplx_view = [reshape(view(packed_cplx, coff[j] + 1:coff[j] + block_c[j]), nfreqs[j], nlayers) for j in 1:nlat_half]

    dev(x) = CuArray(x)
    return GPUFourierGraphCache(
        packed_real, packed_cplx, real_view, cplx_view,
        rfft_plans, brfft_plans,
        dev(roff), dev(coff), dev(nlons), dev(nfreqs),
        dev(rstart_n), dev(rstart_s), dev(nlons_s),
        S.nlon_max, S.nfreq_max, nlat_half, nlayers, has_equator, jeq,
        architecture(packed_real),
        IdDict{Any, Union{Nothing, CuGraphExec}}(),
        IdDict{Any, Union{Nothing, CuGraphExec}}(),
    )
end

"""$(TYPEDSIGNATURES)
The per-size resource that identifies a graph cache. The forward FFT plan set is unique per
transform size and stable for the lifetime of `S`. Must return *the stored*
plan object (not a fresh allocation) so the identity is stable across calls."""
cache_key(S::SpectralTransform, nlayers::Integer) = fft_plans(S, nlayers)[1]

# keyed by the forward FFT plan set (= the per-size resource); `nlayers` sizes the cache
get_cache(S::SpectralTransform, nlayers::Integer) =
    get!(() -> build_cache(S, nlayers), GRAPH_CACHES, cache_key(S, nlayers))::GPUFourierGraphCache

"""$(TYPEDSIGNATURES)
Clear all cached CUDA-Graphs Fourier buffers and graphs (frees the associated GPU memory).
Mainly useful for tests/benchmarks."""
clear_fourier_graph_cache!() = (empty!(GRAPH_CACHES); nothing)

# =====================================================================================
# Allocation-free fused loops (capturable). They write exactly the regions the generic
# `_fourier_batched!` writes (rows 1:nfreq of each ring's scratch slice for the forward;
# the full ring rows of `field` for the inverse), so the result is identical.
# =====================================================================================

"""$(TYPEDSIGNATURES)
Allocation-free, fused forward (grid → spectral) batched Fourier loop: one gather kernel
packs all rings, per-ring in-place cuFFTs run on reshaped views, one scatter kernel writes
the scratch. Suitable for CUDA-graph capture."""
function forward_loop!(cache::GPUFourierGraphCache, f_north, f_south, field::AbstractField, S::SpectralTransform)
    arch = cache.arch
    nlat_half = cache.nlat_half
    L = cache.nlayers
    real_view = cache.real_view
    cplx_view = cache.cplx_view

    rfft_plans = cache.rfft_plans

    # northern rings
    launch!(arch, Array3DWorkOrder, (cache.nlon_max, nlat_half, L), gather_real_kernel!, cache.packed_real, field.data, cache.roff, cache.nlons, cache.rstart_n)
    @inbounds for j in 1:nlat_half
        mul!(cplx_view[j], rfft_plans[j], real_view[j])
    end
    launch!(arch, Array3DWorkOrder, (cache.nfreq_max, nlat_half, L), scatter_cplx_kernel!, f_north, cache.packed_cplx, cache.coff, cache.nfreq)

    # southern rings (the equator ring, if any, is zeroed rather than transformed)
    launch!(arch, Array3DWorkOrder, (cache.nlon_max, nlat_half, L), gather_real_kernel!, cache.packed_real, field.data, cache.roff, cache.nlons, cache.rstart_s)
    @inbounds for j in 1:nlat_half
        if cache.has_equator && j == cache.jeq
            fill!(cplx_view[j], 0)
        else
            mul!(cplx_view[j], rfft_plans[j], real_view[j])
        end
    end
    launch!(arch, Array3DWorkOrder, (cache.nfreq_max, nlat_half, L), scatter_cplx_kernel!, f_south, cache.packed_cplx, cache.coff, cache.nfreq)
    return nothing
end

"""$(TYPEDSIGNATURES)
Allocation-free, fused inverse (spectral → grid) batched Fourier loop: one gather kernel
packs all rings from the scratch, per-ring in-place inverse cuFFTs run on reshaped views,
one scatter kernel writes the grid field. Suitable for CUDA-graph capture."""
function inverse_loop!(cache::GPUFourierGraphCache, field::AbstractField, g_north, g_south, S::SpectralTransform)
    arch = cache.arch
    nlat_half = cache.nlat_half
    L = cache.nlayers
    real_view = cache.real_view
    cplx_view = cache.cplx_view

    brfft_plans = cache.brfft_plans

    # northern rings
    launch!(arch, Array3DWorkOrder, (cache.nfreq_max, nlat_half, L), gather_cplx_kernel!, cache.packed_cplx, g_north, cache.coff, cache.nfreq)
    @inbounds for j in 1:nlat_half
        mul!(real_view[j], brfft_plans[j], cplx_view[j])
    end
    launch!(arch, Array3DWorkOrder, (cache.nlon_max, nlat_half, L), scatter_real_kernel!, field.data, cache.packed_real, cache.roff, cache.nlons, cache.rstart_n)

    # southern rings (the equator ring, if any, is skipped: north already wrote those rows)
    launch!(arch, Array3DWorkOrder, (cache.nfreq_max, nlat_half, L), gather_cplx_kernel!, cache.packed_cplx, g_south, cache.coff, cache.nfreq)
    @inbounds for j in 1:nlat_half
        (cache.has_equator && j == cache.jeq) && continue
        mul!(real_view[j], brfft_plans[j], cplx_view[j])
    end
    launch!(arch, Array3DWorkOrder, (cache.nlon_max, nlat_half, L), scatter_real_kernel!, field.data, cache.packed_real, cache.roff, cache.nlons_s, cache.rstart_s)
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
    exec = get(execs, key, missing)          # `missing` is fallback value that gets return when there's no cached exec yet
    if exec isa CuGraphExec
        launch(exec)                         # hot path: pure replay
        return nothing
    elseif exec === nothing                  # fallback (`nothing` -> tried to capture, but failed (e.g. because MAX_GRAPHS reached))
        loop!()                              # fallback to non-graph direct loop for unknown non-capturable buffer
        return nothing
    end

    # first time we see this buffer (exec === missing)
    if length(execs) >= MAX_GRAPHS
        loop!()                              # cache full: don't capture, just run
        return nothing
    end

    # if we haven't exited yet, so far we know we have to capture

    # warm up so that one-time work (cuFFT init, kernel JIT, memory-pool growth) happens
    # OUTSIDE the capture, where it is allowed
    loop!()
    CUDA.synchronize()

    # do the capture
    graph = capture(throw_error = false) do
        loop!()
    end

    if graph === nothing
        # capture invalidated (e.g. an unexpected allocation/sync); the warmup already
        # produced the correct result. Remember not to retry capture for this buffer.
        execs[key] = nothing
        return nothing
    end

    # save the graph 
    exec = instantiate(graph)
    execs[key] = exec

    # run the graph
    launch(exec)                             # produce the result via the graph
    return nothing
end

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
    if !S.cuda_graphs
        return Base.@invoke _fourier_batched!(
            f_north::AbstractArray{<:Complex, 3}, f_south::AbstractArray{<:Complex, 3},
            field::AbstractField, S::SpectralTransform,
        )
    end
    # the cache is selected by (and sized for) this transform size (= its FFT plan set); within
    # it the only thing that varies between calls is the field buffer, so that is the graph key
    cache = get_cache(S, size(field, 2))
    run_graph!(cache.forward_execs, field.data, () -> forward_loop!(cache, f_north, f_south, field, S))
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
        S::SpectralTransform,
    )
    if !S.cuda_graphs
        return Base.@invoke _fourier_batched!(
            field::AbstractField, g_north::AbstractArray{<:Complex, 3},
            g_south::AbstractArray{<:Complex, 3}, S::SpectralTransform,
        )
    end
    cache = get_cache(S, size(field, 2))
    run_graph!(cache.inverse_execs, field.data, () -> inverse_loop!(cache, field, g_north, g_south, S))
    return nothing
end

end
