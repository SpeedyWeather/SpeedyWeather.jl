# GPU-graph-accelerated batched Fourier transform — shared infrastructure for CUDA and AMDGPU.
#
# Both backends share: the gather/scatter KernelAbstractions kernels, the GPUFourierGraphCache
# struct, the allocation-free forward/inverse loops, graph-cache bookkeeping, and the
# _fourier_batched! entry points (which dispatch on AbstractGPUArray).
#
# What lives in each backend extension: build_cache (backend-specific allocation) and
# run_graph! (backend-specific graph capture/replay API), dispatched via the A type param.
#
# See SpeedyTransformsCUDAExt.jl for a full explanation of why graphs are needed.

# =====================================================================================
# Fused gather/scatter kernels — device-agnostic (KernelAbstractions), shared between
# CUDA and AMDGPU graph paths.
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
# for the FFTs, device-side gather/scatter metadata, and the instantiated GPU graphs.
# =====================================================================================

"""Maximum number of cached graphs per direction per `SpectralTransform`. Prevents
unbounded growth (and host-side capture cost) when the transform is called with a stream
of freshly-allocated `field` buffers (e.g. the allocating `transform(field, S)`). Beyond
this many distinct buffers the allocation-free loop is run directly without capturing."""
const MAX_GRAPHS = 64

"""$(TYPEDSIGNATURES)
Cache (one per transform SIZE, i.e. per FFT plan) holding the pre-allocated packed work
buffers, the per-ring reshaped views and FFT plans used by the transforms, the device
gather/scatter index metadata, and the instantiated GPU graphs keyed by field buffer.
A graph value of `nothing` marks a buffer for which capture failed (fall back to direct
loop). The `A` type parameter holds the architecture and is used to dispatch `run_graph!`
and `build_cache` to the correct backend extension."""
struct GPUFourierGraphCache{PR, PC, RV, CV, IV, A, E}
    packed_real::PR             # GPU vector of NF — all rings' dense real blocks
    packed_complex::PC          # GPU vector of Complex{NF} — all rings' dense complex blocks
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
    has_equator::Bool           # whether the grid has a ring on the equator
    j_equator::Int              # latitude index of the equator ring (if any)
    arch::A                     # SpeedyWeather GPU architecture (for launch! and run_graph! dispatch)
    forward_execs::Dict{UInt, E}  # per-buffer forward graph: missing=unseen, nothing=failed, else exec
    inverse_execs::Dict{UInt, E}  # per-buffer inverse graph
end

# One cache per (SpectralTransform, transform size), keyed by the *forward FFT plan set*.
const GRAPH_CACHES = IdDict{Any, GPUFourierGraphCache}()

# Return the (forward, inverse) per-ring FFT plan sets for a transform of `nlayers` layers.
fft_plans(S::SpectralTransform, nlayers::Integer) = (S.rfft_plans, S.brfft_plans)

cache_key(S::SpectralTransform, nlayers::Integer) = fft_plans(S, nlayers)[1]

get_cache(S::SpectralTransform, nlayers::Integer) =
    get!(() -> build_cache(S, nlayers, S.architecture), GRAPH_CACHES, cache_key(S, nlayers))::GPUFourierGraphCache

"""$(TYPEDSIGNATURES)
Clear all cached GPU-graph Fourier buffers and graphs (frees the associated GPU memory).
Mainly useful for tests/benchmarks."""
clear_fourier_graph_cache!() = (empty!(GRAPH_CACHES); nothing)

# Stable per-buffer graph-cache key: the device address of `field.data`.
# The time stepping fetches the grid field for each transform via a per-step view; on the GPU
# each view is a FRESH array wrapper aliasing the same device memory, so keying the cache on
# wrapper identity would capture a new graph every timestep. The device pointer is stable.
@inline graph_key(data) = reinterpret(UInt, pointer(data))

# =====================================================================================
# Allocation-free fused loops (capturable). Identical for CUDA and AMDGPU — they only
# use launch!, LinearAlgebra.mul!, and fill! which are all device-agnostic.
# =====================================================================================

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
        LinearAlgebra.mul!(complex_view[j], rfft_plans[j], real_view[j])
    end
    launch!(arch, ArrayWorkOrder, complex_size, scatter_complex_kernel!, f_north, packed_complex, complex_offset, nfreqs)

    # southern rings (the equator ring, if any, is zeroed rather than transformed)
    launch!(arch, ArrayWorkOrder, real_size, gather_real_kernel!, packed_real, field.data, real_offset, nlons, istart_s)
    @inbounds for j in 1:nlat_half
        if cache.has_equator && j == j_equator
            fill!(complex_view[j], 0)
        else
            LinearAlgebra.mul!(complex_view[j], rfft_plans[j], real_view[j])
        end
    end
    launch!(arch, ArrayWorkOrder, complex_size, scatter_complex_kernel!, f_south, packed_complex, complex_offset, nfreqs)
    return nothing
end

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
        fft_inverse_mul!(real_view[j], brfft_plans[j], complex_view[j])
    end
    launch!(arch, ArrayWorkOrder, real_size, scatter_real_kernel!, field.data, packed_real, real_offset, nlons, istart_n)

    # southern rings (the equator ring, if any, is skipped: north already wrote those rows)
    launch!(arch, ArrayWorkOrder, complex_size, gather_complex_kernel!, packed_complex, g_south, complex_offset, nfreqs)
    @inbounds for j in 1:nlat_half
        (cache.has_equator && j == j_equator) && continue
        fft_inverse_mul!(real_view[j], brfft_plans[j], complex_view[j])
    end
    launch!(arch, ArrayWorkOrder, real_size, scatter_real_kernel!, field.data, packed_real, real_offset, nlons_s, istart_s)
    return nothing
end

# Stubs extended by backend extensions with methods dispatching on the A type param.
function build_cache end
function run_graph! end

"""$(TYPEDSIGNATURES)
Execute the inverse (complex → real) FFT plan `plan` on `x`, writing into `y`. Defaults to the
generic AbstractFFTs `mul!`. Backend extensions may override this for plan types whose generic
`mul!` performs extra work that isn't HIP/CUDA-graph-capture-safe — e.g. AMDGPU's rocFFT wrapper
defensively copies `x` into a scratch buffer before executing (rocFFT's C2R transform destroys
its input), and that copy is not capturable on some ROCm versions (`hipErrorStreamCaptureUnsupported`).
That copy is unnecessary here: `complex_view[j]` is never read again after this call within a
timestep, so backends may bypass it and let the FFT destroy its input directly."""
fft_inverse_mul!(y, plan, x) = LinearAlgebra.mul!(y, plan, x)

# =====================================================================================
# _fourier_batched! for any GPU backend (CuArray, ROCArray, MtlArray all <: AbstractGPUArray).
# More specific than the AbstractArray methods in fourier.jl, so GPU calls land here.
# =====================================================================================

function _fourier_batched!(
        f_north::AbstractGPUArray{<:Complex, 3},
        f_south::AbstractGPUArray{<:Complex, 3},
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

function _fourier_batched!(
        field::AbstractField,
        g_north::AbstractGPUArray{<:Complex, 3},
        g_south::AbstractGPUArray{<:Complex, 3},
        S::SpectralTransform,
    )
    if !S.gpu_graphs
        return Base.@invoke _fourier_batched!(
            field::AbstractField, g_north::AbstractArray{<:Complex, 3},
            g_south::AbstractArray{<:Complex, 3}, S::SpectralTransform,
        )
    end
    cache = get_cache(S, size(field, 2))
    run_graph!(cache, cache.inverse_execs, graph_key(field.data), () -> inverse_loop!(cache, field, g_north, g_south, S))
    return nothing
end
