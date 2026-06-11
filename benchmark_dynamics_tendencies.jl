# Benchmark for SpeedyWeather.dynamics_tendencies!(vars, lf, model::PrimitiveEquation)
#
# Usage:
#   julia --project=SpeedyWeather benchmark_dynamics_tendencies.jl              # CPU
#   julia --project=SpeedyWeather benchmark_dynamics_tendencies.jl gpu          # GPU (loads CUDA)
#
# Sweeps over both PrimitiveDryModel and PrimitiveWetModel across a range of spectral
# truncations to surface where the mega-batch pays off most.

const ARCH_ARG = lowercase(get(ARGS, 1, "cpu"))
const USE_GPU = ARCH_ARG in ("gpu", "cuda")

using SpeedyWeather
using BenchmarkTools
using Printf

if USE_GPU
    using CUDA
    arch = SpeedyWeather.GPU()
    sync() = CUDA.synchronize()
else
    arch = SpeedyWeather.CPU()
    sync() = nothing
end

# resolutions to sweep, (trunc, nlayers)
const RESOLUTIONS = [
    (31, 8),
    (63, 8),
    (127, 8),
    (255, 8),
]

const MODEL_KINDS = [:PrimitiveWetModel]

"""Build a model + initialised variables, primed by one full timestep so all dynamics
variables (vars.dynamics.{u_mean_grid, dpres_dx, …}) are populated."""
function setup(arch, trunc, nlayers, model_kind::Symbol)
    spectral_grid = SpectralGrid(; trunc, nlayers, architecture = arch)
    ModelType = model_kind === :PrimitiveDryModel ? PrimitiveDryModel : PrimitiveWetModel
    model = ModelType(spectral_grid)
    sim = initialize!(model)
    # Run one timestep so all PHASE-0 inputs to dynamics_tendencies! are populated.
    run!(sim, period = Hour(1))
    return sim.model, sim.variables
end

"""Time a single configuration; returns (median_time_seconds, mem_bytes, allocs)."""
function bench_one(arch, trunc, nlayers, model_kind)
    model, vars = setup(arch, trunc, nlayers, model_kind)
    lf = 2   # later leapfrog index (post-init)

    # warm up + ensure correctness path
    SpeedyWeather.dynamics_tendencies!(vars, lf, model)
    sync()
    SpeedyWeather.dynamics_tendencies!(vars, lf, model)
    sync()

    if USE_GPU
        b = @benchmark begin
            SpeedyWeather.dynamics_tendencies!($vars, $lf, $model)
            CUDA.synchronize()
        end
    else
        b = @benchmark SpeedyWeather.dynamics_tendencies!($vars, $lf, $model)
    end

    return (median(b).time, median(b).memory, median(b).allocs)
end

function main()
    println("Benchmark: SpeedyWeather.dynamics_tendencies!(vars, lf, model::PrimitiveEquation)")
    println("Architecture: ", arch)
    println()
    @printf("%-22s %8s %8s %14s %14s %10s\n",
        "model", "trunc", "nlayers", "median (ms)", "memory (KiB)", "allocs")
    println(repeat('-', 82))

    for model_kind in MODEL_KINDS, (trunc, nlayers) in RESOLUTIONS
        try
            t_ns, mem, allocs = bench_one(arch, trunc, nlayers, model_kind)
            @printf("%-22s %8d %8d %14.3f %14.1f %10d\n",
                String(model_kind), trunc, nlayers,
                t_ns / 1.0e6, mem / 1024, allocs)
        catch err
            @printf("%-22s %8d %8d   FAILED: %s\n",
                String(model_kind), trunc, nlayers, sprint(showerror, err))
        end
    end
end

main()
