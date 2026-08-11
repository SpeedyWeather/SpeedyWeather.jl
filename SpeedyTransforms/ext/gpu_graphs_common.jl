# =====================================================================================
# BACKEND-AGNOSTIC GPU-GRAPHS MACHINERY FOR THE BATCHED FOURIER TRANSFORM
#
# `include()`-d by each graph-capturing backend extension (currently SpeedyTransformsCUDAExt.jl
# and SpeedyTransformsAMDGPUExt.jl). See SpeedyTransformsCUDAExt.jl for the rationale of
# GPU-graphs acceleration itself.
#
# This file must stay under SpeedyTransforms/ext/ and be `include()`-d from within an
# extension module, NOT moved to SpeedyTransforms/src/: the `_fourier_batched!` methods that
# actually use this machinery live in each extension, and the main package loads
# unconditionally in every Julia session — so if any of this were in src/ instead, those
# methods would exist even in CPU-only sessions with no GPU backend loaded, which breaks
# Enzyme's static activity analysis (see docs/dev/2026-08/gpu-graphs-common-interface.md).
#
# The only thing that differs between backends is the graph capture/instantiate/launch API
# itself (e.g. `CUDA.capture`/`instantiate`/`launch` vs a HIP equivalent) and the graph-exec
# type it produces (e.g. `CuGraphExec`). That's bundled into a `GraphBackend`, supplied by the
# including extension. Everything else — kernels, the per-transform-size cache, the
# allocation-free fused loops, and capture/replay control flow — is identical across backends
# and lives here once.
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
# for the FFTs, device-side gather/scatter metadata, and the instantiated GPU graphs.
# =====================================================================================

"""$(TYPEDSIGNATURES)
Cache (one per transform SIZE, i.e. per FFT plan) holding the pre-allocated packed
work buffers, the per-ring reshaped views and FFT plans used by the transforms, the device
gather/scatter index metadata, and the instantiated GPU graphs (one per distinct `field`
buffer and direction). A graph value of `nothing` marks a buffer for which capture failed
(fall back to direct loop). `E` is the backend's graph-exec type (e.g. `CuGraphExec`)."""
struct GPUFourierGraphCache{PR, PC, RV, CV, RFP, BFP, IV, A, E}
    packed_real::PR             # device Vector{NF}          — all rings' dense real blocks
    packed_complex::PC          # device Vector{Complex{NF}} — all rings' dense complex blocks
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
    forward_execs::Dict{UInt, Union{Nothing, E}} # forward graphs, keyed by field buffer
    # inverse graphs, keyed by (field buffer, `add`): a captured graph bakes in whether the scatter
    # overwrites or accumulates, so the two modes must never share a graph.
    inverse_execs::Dict{Tuple{UInt, Bool}, Union{Nothing, E}}
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
fft_plans(S::SpectralTransform, nlayers::Integer) = (S.rfft_plans_batched[nlayers], S.brfft_plans_batched[nlayers])

# build/allocate the cache for a transform of `nlayers` layers. `E` is the backend's
# graph-exec type, supplied by the including extension (e.g. `CuGraphExec`).
function build_cache(::Type{E}, S::SpectralTransform, nlayers::Integer) where {E}
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

    arch = S.architecture
    dev(x) = on_architecture(arch, x)              # device-agnostic transfer (works for any GPU backend)
    packed_real = dev(zeros(NF, sum(block_real)))
    packed_complex = dev(zeros(Complex{NF}, sum(block_complex)))
    real_view = [reshape(view(packed_real, real_offset[j] + 1:real_offset[j] + block_real[j]), nlons[j], nlayers) for j in 1:nlat_half]
    complex_view = [reshape(view(packed_complex, complex_offset[j] + 1:complex_offset[j] + block_complex[j]), nfreqs[j], nlayers) for j in 1:nlat_half]

    return GPUFourierGraphCache(
        packed_real, packed_complex, real_view, complex_view,
        rfft_plans, brfft_plans,
        dev(real_offset), dev(complex_offset), dev(nlons), dev(nfreqs),
        dev(istart_n), dev(istart_s), dev(nlons_s),
        S.nlon_max, S.nfreq_max, nlat_half, nlayers, has_equator, j_equator,
        arch,
        Dict{UInt, Union{Nothing, E}}(),
        Dict{Tuple{UInt, Bool}, Union{Nothing, E}}(),
    )
end

"""$(TYPEDSIGNATURES)
The per-size resource that identifies a graph cache. The forward FFT plan set is unique per
transform size and stable for the lifetime of `S`. Must return *the stored*
plan object (not a fresh allocation) so the identity is stable across calls."""
cache_key(S::SpectralTransform, nlayers::Integer) = fft_plans(S, nlayers)[1]

# keyed by the forward FFT plan set (= the per-size resource); `nlayers` sizes the cache
get_cache(::Type{E}, S::SpectralTransform, nlayers::Integer) where {E} =
    get!(() -> build_cache(E, S, nlayers), GRAPH_CACHES, cache_key(S, nlayers))::GPUFourierGraphCache

"""$(TYPEDSIGNATURES)
Clear all cached GPU-graphs Fourier buffers and graphs (frees the associated GPU memory).
Mainly useful for tests/benchmarks."""
clear_fourier_graph_cache!() = (empty!(GRAPH_CACHES); nothing)

# =====================================================================================
# Allocation-free fused loops (capturable). They write exactly the regions the generic
# `_fourier_batched!` writes (rows 1:nfreq of each ring's scratch slice for the forward;
# the full ring rows of `field` for the inverse), so the result is identical.
# =====================================================================================

"""$(TYPEDSIGNATURES)
Allocation-free, fused forward (grid → spectral) batched Fourier loop: one gather kernel
packs all rings, per-ring in-place FFTs run on reshaped views, one scatter kernel writes
the scratch. Suitable for GPU-graph capture."""
function forward_loop!(cache::GPUFourierGraphCache, f_north, f_south, field::AbstractField, S::SpectralTransform)
    (; arch) = cache
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
packs all rings from the scratch, per-ring in-place inverse FFTs run on reshaped views,
one scatter kernel writes the grid field. `add` accumulates onto `field` instead of overwriting it
(Enzyme adjoint rule); it is baked into the captured graph, hence the per-`add` graph cache key.
Suitable for GPU-graph capture."""
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
# =====================================================================================

"""$(TYPEDSIGNATURES)
Bundles the four backend-specific primitives `run_graph!` needs to capture and replay a GPU
graph, so the control flow in `run_graph!` itself can stay backend-agnostic:

- `capture(loop!)`: run `loop!` under stream capture and return the recorded graph, or
  `nothing` if capture failed/was invalidated (must not throw).
- `instantiate(graph)`: turn a captured graph into an executable.
- `launch(exec)`: replay an instantiated graph.
- `synchronize()`: block until the device is idle (used once, before capture, so that
  one-time warmup work happens outside the captured region).

e.g. for CUDA: `GraphBackend(loop! -> CUDA.capture(loop!; throw_error = false), CUDA.instantiate, CUDA.launch, CUDA.synchronize)`."""
struct GraphBackend{C, I, L, Y}
    capture::C
    instantiate::I
    launch::L
    synchronize::Y
end

"""$(TYPEDSIGNATURES)
Run `loop!` (a `() -> ...` closure over the allocation-free batched FFT) either by
replaying a cached GPU graph keyed by `key`, or — on first use — by warming up, capturing,
instantiating and caching a graph via the given `backend`. Falls back to running `loop!`
directly if capture fails or the per-direction cache is full."""
function run_graph!(execs::AbstractDict, key, loop!::F, backend::GraphBackend) where {F}
    exec = get(execs, key, missing)          # `missing` is fallback value that gets return when there's no cached exec yet
    if exec !== missing && exec !== nothing
        backend.launch(exec)                 # hot path: pure replay
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

    # warm up so that one-time work (FFT plan init, kernel JIT, memory-pool growth) happens
    # OUTSIDE the capture, where it is allowed
    loop!()
    backend.synchronize()

    # do the capture
    graph = backend.capture(loop!)

    if graph === nothing
        # capture invalidated (e.g. an unexpected allocation/sync); the warmup already
        # produced the correct result. Remember not to retry capture for this buffer.
        execs[key] = nothing
        return nothing
    end

    # save the graph
    exec = backend.instantiate(graph)
    execs[key] = exec

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
