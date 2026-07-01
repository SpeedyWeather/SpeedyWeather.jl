module SpeedyTransformsCUDAExt

import CUDA: CUDA, CuArray, CuGraphExec, capture, instantiate, launch

using SpeedyTransforms
using SpeedyTransforms.RingGrids
using SpeedyTransforms.LowerTriangularArrays

import SpeedyTransforms: SpectralTransform, GPUFourierGraphCache, build_cache, run_graph!, MAX_GRAPHS, fft_plans

import SpeedyWeatherInternals.Architectures: GPU, architecture

# =====================================================================================
# build_cache: allocate the packed work buffers and per-ring views on CUDA device memory.
# =====================================================================================

function build_cache(S::SpectralTransform, nlayers::Integer, ::GPU{<:CUDA.CUDABackend})
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

    packed_real = CUDA.zeros(NF, sum(block_real))
    packed_complex = CUDA.zeros(Complex{NF}, sum(block_complex))
    real_view = [reshape(view(packed_real, real_offset[j] + 1:real_offset[j] + block_real[j]), nlons[j], nlayers) for j in 1:nlat_half]
    complex_view = [reshape(view(packed_complex, complex_offset[j] + 1:complex_offset[j] + block_complex[j]), nfreqs[j], nlayers) for j in 1:nlat_half]

    dev(x) = CuArray(x)
    exec_dict() = Dict{UInt, Union{Nothing, CuGraphExec}}()
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
# run_graph!: CUDA-specific graph capture and replay.
# Dispatches on the A=GPU{<:CUDA.CUDABackend} type parameter of GPUFourierGraphCache.
# The E type param on the dict ensures exec is Union{Nothing, CuGraphExec, Missing}
# — a small union Julia handles without boxing.
# =====================================================================================

function run_graph!(
        ::GPUFourierGraphCache{<:Any, <:Any, <:Any, <:Any, <:Any, <:GPU{<:CUDA.CUDABackend}, E},
        execs::Dict{UInt, E}, key, loop!::F,
    ) where {E, F}
    exec = get(execs, key, missing)
    if exec isa CuGraphExec
        launch(exec)                         # hot path: pure replay
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

    # warm up so that one-time work (cuFFT init, kernel JIT, memory-pool growth) happens
    # OUTSIDE the capture region where it is not allowed
    loop!()
    CUDA.synchronize()

    graph = capture(throw_error = false) do
        loop!()
    end

    if graph === nothing
        # capture invalidated (e.g. an unexpected allocation); warmup already produced
        # the correct result — remember not to retry for this buffer
        execs[key] = nothing
        return nothing
    end

    exec = instantiate(graph)
    execs[key] = exec
    launch(exec)
    return nothing
end

end
