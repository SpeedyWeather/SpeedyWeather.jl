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

import SpeedyWeatherInternals.KernelLaunching: launch!, ArrayWorkOrder
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
# several layer counts gets one cache each); within a cache, graphs are keyed by the device
# pointer of `field.data` (see `graph_key`), keying on the wrapper object identity wouldn't
# work as a view on a CuArray is again a CuArray and not a SubArray which causes problems here.
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
@kernel inbounds = true function scatter_real_kernel!(dst, packed, real_offset, nlons, istart)
    i, j, k = @index(Global, NTuple)
    nlon = nlons[j]
    if i <= nlon
        dst[istart[j] + i - 1, k] = packed[real_offset[j] + (k - 1) * nlon + i]
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
    packed_complex::PC          # CuVector{Complex{NF}} — all rings' dense complex blocks
    real_view::RV               # Vector of per-ring reshaped (nlon_j × nlayers) views into packed_real
    complex_view::CV            # Vector of per-ring reshaped (nfreq_j × nlayers) views into packed_complex
    rfft_plans::Vector{AbstractFFTs.Plan}   # forward FFT plans for THIS size (nlayers batches), per ring
    brfft_plans::Vector{AbstractFFTs.Plan}  # inverse FFT plans for THIS size, per ring
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
    forward_execs::Dict{UInt, Union{Nothing, CuGraphExec}} # forward graphs
    inverse_execs::Dict{UInt, Union{Nothing, CuGraphExec}} # inverse graphs
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
fft_plans(S::SpectralTransform, nlayers::Integer) = (S.rfft_plans[nlayers], S.brfft_plans[nlayers])

# build/allocate the cache for a transform of `nlayers` layers
function build_cache(S::SpectralTransform, nlayers::Integer)
    NF = eltype(S)
    nlat = S.nlat
    nlat_half = S.grid.nlat_half
    rfft_plans, brfft_plans = fft_plans(S, nlayers)
    rings = S.rings
    nlons = collect(Int, S.nlons[1:nlat_half])
    nfreqs = nlons .÷ 2 .+ 1

    # packed buffer block offset
    block_real = nlons .* nlayers
    block_complex = nfreqs .* nlayers
    real_offset = [0; cumsum(block_real)[1:(end - 1)]]
    complex_offset = [0; cumsum(block_complex)[1:(end - 1)]]

    istart_n = [rings[j].start for j in 1:nlat_half]
    istart_s = [rings[nlat - j + 1].start for j in 1:nlat_half]

    has_equator = isodd(nlat)
    j_equator = (nlat + 1) ÷ 2
    nlons_s = copy(nlons)
    has_equator && (nlons_s[j_equator] = 0)        # south pass skips the equator ring

    packed_real = CUDA.zeros(NF, sum(block_real))
    packed_complex = CUDA.zeros(Complex{NF}, sum(block_complex))
    real_view = [reshape(view(packed_real, real_offset[j] + 1:real_offset[j] + block_real[j]), nlons[j], nlayers) for j in 1:nlat_half]
    complex_view = [reshape(view(packed_complex, complex_offset[j] + 1:complex_offset[j] + block_complex[j]), nfreqs[j], nlayers) for j in 1:nlat_half]

    dev(x) = CuArray(x)
    return GPUFourierGraphCache(
        packed_real, packed_complex, real_view, complex_view,
        rfft_plans, brfft_plans,
        dev(real_offset), dev(complex_offset), dev(nlons), dev(nfreqs),
        dev(istart_n), dev(istart_s), dev(nlons_s),
        S.nlon_max, S.nfreq_max, nlat_half, nlayers, has_equator, j_equator,
        architecture(packed_real),
        Dict{UInt, Union{Nothing, CuGraphExec}}(),
        Dict{UInt, Union{Nothing, CuGraphExec}}(),
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
    (; arch ) = cache
    (; nlat_half, nlayers) = cache
    (; nlon_max, nfreq_max, nlons, nfreqs, j_equator) = cache
    (; real_view, complex_view, packed_real, packed_complex, rfft_plans) = cache
    (; real_offset, complex_offset, istart_n, istart_s) = cache

    real_size = (nlon_max, nlat_half, nlayers)          # launch scatter/gather kernels over 
    complex_size = (nfreq_max, nlat_half, nlayers)

    # northern rings, 
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
packs all rings from the scratch, per-ring in-place inverse cuFFTs run on reshaped views,
one scatter kernel writes the grid field. Suitable for CUDA-graph capture."""
function inverse_loop!(cache::GPUFourierGraphCache, field::AbstractField, g_north, g_south, S::SpectralTransform)
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
    launch!(arch, ArrayWorkOrder, real_size, scatter_real_kernel!, field.data, packed_real, real_offset, nlons, istart_n)

    # southern rings (the equator ring, if any, is skipped: north already wrote those rows)
    launch!(arch, ArrayWorkOrder, complex_size, gather_complex_kernel!, packed_complex, g_south, complex_offset, nfreqs)
    @inbounds for j in 1:nlat_half
        (cache.has_equator && j == j_equator) && continue
        mul!(real_view[j], brfft_plans[j], complex_view[j])
    end
    launch!(arch, ArrayWorkOrder, real_size, scatter_real_kernel!, field.data, packed_real, real_offset, nlons_s, istart_s)
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
function run_graph!(execs::AbstractDict, key, loop!::F) where {F}
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

# Stable per-buffer graph-cache key: the device address of `field.data`.
#
# The time stepping fetches the grid field for each transform via `get_prognostic_step` /
# `get_tendency_step`, i.e. a `field_view`/`get_step` on a step-dimensioned grid variable. On
# the GPU `view(::CuArray, :, :, step)` returns again a fresh `CuArray` wrapper each call
# , so keying the cache on the wrapper's object identity would capture a
# new graph every timestep and the cache would grow without bound. The wrapper churns but it
# always aliases the same device memory at the same offset, so the device pointer — exactly what
# the captured graph bakes in — is stable; we key on that instead. (Safe because model buffers live
# for the whole run, so an address is never freed and recycled under a stale graph.)
@inline graph_key(data) = reinterpret(UInt, pointer(data))

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
    run_graph!(cache.forward_execs, graph_key(field.data), () -> forward_loop!(cache, f_north, f_south, field, S))
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
    run_graph!(cache.inverse_execs, graph_key(field.data), () -> inverse_loop!(cache, field, g_north, g_south, S))
    return nothing
end

end
