# Phase 0b of docs/dev/2026-11/2026-11-17-fused-tracers.md: are N tracers faster when processed together?
#
# Times the spectral → grid transform and the horizontal advection of N tracers as the model does it now (one
# tracer after another, `foreach_active_tracer`) against a prototype that processes all tracers together: the
# tracers stacked along the layer axis of one array, one batched transform instead of N, one kernel instead of 2N.
# These are 4 of the transforms and 2 of the kernels per tracer and time step. The forward transforms of the
# advection are timed as three batches of N·L layers, and as one batch of 3N·L layers. Batches are split into
# chunks of whole tracers that fit the transform's largest batch (`S.nlayers`, set by `transform_batch`).
# Also times one full time step with N tracers, to relate the saving to it, and checks that the prototype
# gives the same results as the model.
#
# Run from the repository root with the benchmark environment, the architecture as argument (cpu, gpu for
# NVIDIA, metal: Metal is not in the benchmark environment, add it first):
#
#     julia --project=SpeedyWeather/benchmark SpeedyWeather/benchmark/fused_tracers_prototype.jl gpu
#     julia --project=SpeedyWeather/benchmark SpeedyWeather/benchmark/fused_tracers_prototype.jl gpu --output results.md
#
# The results are printed as a markdown table (and written to the file given with `--output`). `--quick` runs only
# T32 with 1 and 2 tracers, to check that the script works.

const ARCHITECTURE_ARG = isempty(ARGS) ? "cpu" : lowercase(ARGS[1])
const OUTPUT_INDEX = findfirst(==("--output"), ARGS)
const OUTPUT_PATH = isnothing(OUTPUT_INDEX) ? nothing : abspath(ARGS[OUTPUT_INDEX + 1])
const QUICK = "--quick" in ARGS

ARCHITECTURE_ARG == "gpu" && using CUDA
ARCHITECTURE_ARG == "metal" && using Metal

using SpeedyWeather, Printf, Statistics, InteractiveUtils
import SpeedyWeather.KernelAbstractions: @kernel, @index
using SpeedyWeather.KernelLaunching: launch!, RingGridWorkOrder
import SpeedyWeather: synchronize, foreach_active_tracer, transform_tracer!, tracer_advection!, divergence!,
    get_prognostic_step, get_tendency_step, DynamicalCore, time_step!
import SpeedyWeather.SpeedyTransforms: wrapped_view

const ARCHITECTURE = ARCHITECTURE_ARG == "cpu" ? SpeedyWeather.CPU() : SpeedyWeather.GPU()

# (model, truncation, nlayers) and the numbers of tracers, as in phase 0a without the case of no tracers
const CONFIGURATIONS = QUICK ? [(ShallowWaterModel, 32, 1), (PrimitiveWetModel, 32, 8)] : [
        (ShallowWaterModel, 32, 1), (ShallowWaterModel, 64, 1), (ShallowWaterModel, 128, 1),
        (PrimitiveWetModel, 32, 8), (PrimitiveWetModel, 64, 8), (PrimitiveWetModel, 128, 8), (PrimitiveWetModel, 256, 8),
    ]
const NTRACERS = QUICK ? (1, 2) : (1, 2, 4, 8)

# timing: median over samples of the time per call, each sample calls `f` repeatedly for at least MIN_SAMPLE_TIME
const NSAMPLES = 7
const MIN_SAMPLE_TIME = 0.2     # [s]

"""All tracers at once: `+A*div` (if `with_divergence`) into the tendencies `tend[:, 1:NL]`, and the fluxes
`u*A/coslat`, `v*A/coslat` into `tend[:, NL+1:2NL]`, `tend[:, 2NL+1:3NL]`. `A` holds N tracers of L layers stacked
along the layer axis, the shared fields `div`, `u`, `v` have L layers. The kernels `_horizontal_advection_kernel!`
and `_flux_divergence_kernel!` of the model in one, over all tracers."""
@kernel inbounds = true function fused_tracer_advection_kernel!(
        tend, A, div, u, v, coslat⁻¹, whichring, L, NL, with_divergence
    )
    ij, kn = @index(Global, NTuple)     # grid point ij, layer kn of the stacked tracers
    k = mod1(kn, L)                     # layer of the shared fields
    a = A[ij, kn]
    if with_divergence
        tend[ij, kn] = muladd(a, div[ij, k], tend[ij, kn])
    end
    a_coslat⁻¹ = a * coslat⁻¹[whichring[ij]]
    tend[ij, NL + kn] = u[ij, k] * a_coslat⁻¹
    tend[ij, 2NL + kn] = v[ij, k] * a_coslat⁻¹
end

"""Transform `input` to `output`, both holding `nslabs` slabs of `L` layers, in chunks of whole slabs that fit
the transform's largest batch `S.nlayers`. One call if all slabs fit."""
function transform_slabs!(output, input, nslabs, L, scratch_memory, S)
    slabs_per_call = max(1, S.nlayers ÷ L)
    slabs_per_call >= nslabs && return transform!(output, input, scratch_memory, S)
    for first_slab in 1:slabs_per_call:nslabs
        last_slab = min(first_slab + slabs_per_call - 1, nslabs)
        layers = ((first_slab - 1) * L + 1):(last_slab * L)
        transform!(wrapped_view(output, :, layers), wrapped_view(input, :, layers), scratch_memory, S)
    end
    return output
end

"""The arrays of the fused prototype: N tracers stacked along the layer axis of one array each."""
function fused_arrays(spectral_grid, N, L)
    (; spectrum, grid, NF) = spectral_grid
    return (
        spectral = zeros(Complex{NF}, spectrum, N * L),         # prognostic tracers, current step
        grid = zeros(NF, grid, N * L),                          # grid tracers, current step
        tend_grid = zeros(NF, grid, 3 * N * L),                 # [tendencies; u*A; v*A] on the grid
        tend_spectral = zeros(Complex{NF}, spectrum, 3 * N * L),    # and in spectral space
    )
end

"""Copy the current step of the model's prognostic tracers into the stacked array `fused.spectral`."""
function copy_tracers!(fused, vars, model, names, L)
    for (n, name) in enumerate(names)
        tracer = get_prognostic_step(getfield(vars.prognostic.tracers, name), model.time_stepping, model.spectral_transform, model)
        fused.spectral.data[:, ((n - 1) * L + 1):(n * L)] .= tracer.data
    end
    return fused
end

"""The prototype: spectral → grid of all tracers, then their horizontal advection, with the three forward
transforms as three batches of N·L layers (`single_batch = false`) or as one batch of 3N·L layers."""
function fused!(fused, vars, model, N, L; single_batch::Bool)
    S = model.spectral_transform
    scratch_memory = vars.scratch.transform_memory
    time_stepping = model.time_stepping
    NL = N * L

    # spectral → grid, all tracers
    transform_slabs!(fused.grid, fused.spectral, N, L, scratch_memory, S)

    # +A*div and the fluxes u*A, v*A of all tracers, one kernel
    u = get_prognostic_step(vars.grid.u, time_stepping, DynamicalCore(), model)
    v = get_prognostic_step(vars.grid.v, time_stepping, DynamicalCore(), model)
    with_divergence = haskey(vars.grid, :divergence)
    div = with_divergence ? get_prognostic_step(vars.grid.divergence, time_stepping, DynamicalCore(), model) : u
    (; coslat⁻¹) = model.geometry
    launch!(
        ARCHITECTURE, RingGridWorkOrder, (size(fused.grid, 1), NL), fused_tracer_advection_kernel!,
        fused.tend_grid, fused.grid, div, u, v, coslat⁻¹, fused.grid.grid.whichring, L, NL, with_divergence
    )

    # grid → spectral of the tendencies and both fluxes
    if single_batch
        transform_slabs!(fused.tend_spectral, fused.tend_grid, 3N, L, scratch_memory, S)
    else
        for part in 0:2
            layers = (part * NL + 1):((part + 1) * NL)
            transform_slabs!(
                wrapped_view(fused.tend_spectral, :, layers), wrapped_view(fused.tend_grid, :, layers),
                N, L, scratch_memory, S
            )
        end
    end

    # tendencies += -∇⋅(u*A, v*A), all tracers
    tend = wrapped_view(fused.tend_spectral, :, 1:NL)
    uA = wrapped_view(fused.tend_spectral, :, (NL + 1):(2NL))
    vA = wrapped_view(fused.tend_spectral, :, (2NL + 1):(3NL))
    divergence!(tend, uA, vA, S; add = true, flipsign = true)
    return nothing
end

"""As the model does it now: spectral → grid and the horizontal advection, one tracer after another."""
function per_tracer!(vars, model)
    S = model.spectral_transform
    foreach_active_tracer(transform_tracer!, vars, model, model.time_stepping, vars.scratch.transform_memory, S)
    tracer_advection!(vars, model)
    return nothing
end

"""Median time per call of `f` in seconds, see `NSAMPLES` and `MIN_SAMPLE_TIME`."""
function time_per_call(f)
    for _ in 1:3            # warm-up: compilation, FFT plans of new batch sizes, GPU graph capture
        f()
    end
    synchronize(ARCHITECTURE)
    t0 = time()             # number of calls per sample from one more call
    f()
    synchronize(ARCHITECTURE)
    ncalls = clamp(ceil(Int, MIN_SAMPLE_TIME / max(time() - t0, 1.0e-6)), 1, 10_000)
    samples = map(1:NSAMPLES) do _
        synchronize(ARCHITECTURE)
        t0 = time()
        for _ in 1:ncalls
            f()
        end
        synchronize(ARCHITECTURE)
        return (time() - t0) / ncalls
    end
    return median(samples)
end

"""Largest relative difference between the prototype and the model: of the grid tracers, and of the spectral
tendencies of the horizontal advection (both from zero grid tendencies, one call each from the same state)."""
function max_relative_difference(fused, vars, model, names, N, L)
    time_stepping = model.time_stepping
    S = model.spectral_transform
    copy_tracers!(fused, vars, model, names, L)
    fill!(fused.tend_grid.data, 0)
    for name in names
        fill!(get_tendency_step(getfield(vars.tendencies.grid_tracers, name), time_stepping, DynamicalCore()).data, 0)
    end
    per_tracer!(vars, model)
    fused!(fused, vars, model, N, L; single_batch = false)

    # Inf if the model's array is all zeros: then the comparison would say nothing
    function relative_difference(a, b)
        scale = maximum(abs, Array(b))
        return iszero(scale) ? Inf : maximum(abs, Array(a) .- Array(b)) / scale
    end
    differences = map(enumerate(names)) do (n, name)
        layers = ((n - 1) * L + 1):(n * L)
        grid = get_prognostic_step(getfield(vars.grid.tracers, name), time_stepping, S, model)
        tend = get_tendency_step(getfield(vars.tendencies.tracers, name), time_stepping, DynamicalCore())
        return max(
            relative_difference(fused.grid.data[:, layers], grid.data),
            relative_difference(fused.tend_spectral.data[:, layers], tend.data),
        )
    end
    return maximum(differences)
end

function benchmark(Model, truncation, nlayers, N)
    spectral_grid = SpectralGrid(; truncation, nlayers, architecture = ARCHITECTURE)
    names = ntuple(i -> Symbol(:tracer, i), N)
    tracers = NamedTuple{names}(ntuple(_ -> Tracer(), N))
    model = Model(spectral_grid; tracers)
    model.feedback.verbose = false
    simulation = initialize!(model)
    set!(simulation; namespace = :tracers, (name => ((λ, φ, σ) -> exp(-(λ - 180)^2 / 30^2) * (1 + σ)) for name in names)...)
    run!(simulation, steps = 4)     # a realistic state: winds, divergence and tracers on the grid
    vars, model = simulation.variables, simulation.model
    L = nlayers

    fused = fused_arrays(spectral_grid, N, L)
    difference = max_relative_difference(fused, vars, model, names, N, L)

    time_per_tracer = time_per_call(() -> per_tracer!(vars, model))
    time_fused = time_per_call(() -> fused!(fused, vars, model, N, L; single_batch = false))
    time_fused_single = time_per_call(() -> fused!(fused, vars, model, N, L; single_batch = true))
    initialize!(simulation)         # the clock etc. for time steps outside of `run!`
    time_step = time_per_call(() -> time_step!(vars, model.time_stepping, model))

    return (; Model, truncation, nlayers, N, time_per_tracer, time_fused, time_fused_single, time_step, difference)
end

function results_table(results)
    io = IOBuffer()
    println(io, "| Model | T | L | Tracers | per tracer [ms] | fused, 3 batches [ms] | fused, 1 batch [ms] | speedup | time step [ms] | saved of time step | max. rel. difference |")
    println(io, "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for r in results
        time_fused = min(r.time_fused, r.time_fused_single)
        @printf(
            io, "| %s | %d | %d | %d | %.3f | %.3f | %.3f | %.2f× | %.3f | %.0f%% | %.1e |\n",
            nameof(r.Model), r.truncation, r.nlayers, r.N, 1.0e3 * r.time_per_tracer, 1.0e3 * r.time_fused,
            1.0e3 * r.time_fused_single, r.time_per_tracer / time_fused, 1.0e3 * r.time_step,
            100 * (r.time_per_tracer - time_fused) / r.time_step, r.difference
        )
    end
    return String(take!(io))
end

function machine_summary()
    io = IOBuffer()
    println(
        io, "Architecture: `$ARCHITECTURE_ARG`, SpeedyWeather v$(pkgversion(SpeedyWeather)), Julia v$VERSION, ",
        "$(Sys.cpu_info()[1].model), $(Sys.CPU_THREADS) threads."
    )
    if ARCHITECTURE_ARG == "gpu"
        println(io, "GPU: ", CUDA.name(CUDA.device()), ", CUDA.jl v$(pkgversion(CUDA)).")
    elseif ARCHITECTURE_ARG == "metal"
        println(io, "GPU: ", Metal.current_device().name, ", Metal.jl v$(pkgversion(Metal)).")
    end
    return String(take!(io))
end

results = []
for (Model, truncation, nlayers) in CONFIGURATIONS, N in NTRACERS
    @info "$Model T$truncation L$nlayers with $N tracers"
    try
        push!(results, benchmark(Model, truncation, nlayers, N))
    catch error
        @warn "$Model T$truncation L$nlayers with $N tracers failed" exception = (error, catch_backtrace())
    end
end

report = machine_summary() * "\n" * results_table(results) * """

    "per tracer" is the spectral → grid transform and the horizontal advection of all N tracers as the model does it
    now, "fused" the prototype (best of the two variants for the speedup), "time step" one full time step with N
    tracers, "saved of time step" what the fused prototype would save of it. "max. rel. difference" between the
    prototype and the model (grid tracers and advection tendencies).
    """
print(report)
isnothing(OUTPUT_PATH) || (write(OUTPUT_PATH, report); @info "Results written to $OUTPUT_PATH")
