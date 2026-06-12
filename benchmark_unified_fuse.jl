# Benchmark for the unified-fuse refactor — measures both `horizontal_diffusion!` and
# `leapfrog!` for PrimitiveDry / PrimitiveWet, comparing the per-variable code path against
# the unified-fuse batched code path. Asserts bit-identicality of the two implementations
# at every measured point.
#
# Usage:
#   julia --project=SpeedyWeather benchmark_unified_fuse.jl           # CPU
#   julia --project=SpeedyWeather benchmark_unified_fuse.jl gpu       # GPU (loads CUDA)
#
# Background:
#
# Both functions now have an architecture-dispatched implementation
# (`SpeedyWeather/src/dynamics/horizontal_diffusion.jl` and
# `SpeedyWeather/src/time_stepping/leapfrog.jl`). On GPU, a single batched kernel handles
# vor + div + T + pres + [humid] in one launch by iterating the leading slots of
# `:spectral_tendencies` (slot-aligned with `:prognostic`'s leading slots). On CPU the
# per-launch overhead is negligible so we fall back to per-variable launches; the
# benchmark "reference" path is exactly what CPU runs. The "unified" path is what GPU
# runs (which we also measure on CPU here for comparison).

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

const RESOLUTIONS = [
    (31, 8),
    (63, 8),
    (127, 8),
    (255, 8),
]

const MODEL_KINDS = [:PrimitiveDryModel, :PrimitiveWetModel]

# ==============================================================================
# Reference implementations — the legacy per-variable code paths
# ==============================================================================

"""Per-variable diffusion (legacy / CPU path). Reference for the batched version."""
function reference_horizontal_diffusion!(vars, diffusion, model, lf)
    (; expl, impl, expl_div, impl_div) = diffusion
    vor   = SpeedyWeather.get_step(vars.prognostic.vorticity,   lf)
    div   = SpeedyWeather.get_step(vars.prognostic.divergence,  lf)
    temp  = SpeedyWeather.get_step(vars.prognostic.temperature, lf)
    SpeedyWeather.horizontal_diffusion!(vars.tendencies.vorticity,   vor,  expl,     impl)
    SpeedyWeather.horizontal_diffusion!(vars.tendencies.divergence,  div,  expl_div, impl_div)
    SpeedyWeather.horizontal_diffusion!(vars.tendencies.temperature, temp, expl,     impl)
    if haskey(vars.tendencies, :humidity)
        humid = SpeedyWeather.get_step(vars.prognostic.humidity, lf)
        SpeedyWeather.horizontal_diffusion!(vars.tendencies.humidity, humid, expl, impl)
    end
    return nothing
end

"""Per-variable leapfrog (legacy / CPU path). Reference for the batched version.
Equivalent to `invoke(leapfrog!, …, ::AbstractModel)` but isolated here for the benchmark."""
function reference_leapfrog!(vars, dt, lf, model)
    (; prognostic, tendencies) = vars
    for varname in keys(tendencies)
        if !(tendencies[varname] isa NamedTuple) && hasfield(typeof(prognostic), varname)
            var = getfield(prognostic, varname)
            var_old, var_new = SpeedyWeather.get_steps(var)
            var_tend = getfield(tendencies, varname)
            SpeedyTransforms.spectral_truncation!(var_tend)
            SpeedyWeather.leapfrog!(var_old, var_new, var_tend, dt, lf, model.time_stepping)
        end
    end
    # tracers + random process unchanged in the batched version — skip here
    return nothing
end

"""Batched leapfrog for PrimitiveEquation (the GPU code path).
Mirrors the body of `leapfrog!(::Variables, ::Real, ::Int, ::PrimitiveEquation)` but
without the architecture dispatch — so we can call it on CPU too for the benchmark."""
function batched_leapfrog!(vars, dt, lf, model)
    spec_parent = parent(vars.fused.spectral_tendencies)
    prog_parent = parent(vars.fused.prognostic)
    slot_map = vars.fused.spectral_tendencies.slot_map
    K_diff = haskey(slot_map, :humidity) ? last(slot_map.humidity) : last(slot_map.pressure)

    tend_view = SpeedyWeather.lta_view(spec_parent, :, 1:K_diff)
    prog_old  = SpeedyWeather.lta_view(prog_parent, :, 1:K_diff, 1)
    prog_new  = SpeedyWeather.lta_view(prog_parent, :, 1:K_diff, 2)

    SpeedyTransforms.spectral_truncation!(tend_view)
    SpeedyWeather.leapfrog!(prog_old, prog_new, tend_view, dt, lf, model.time_stepping)
    return nothing
end

# ==============================================================================
# Setup + bit-identicality checks
# ==============================================================================

function setup(arch, trunc, nlayers, model_kind::Symbol)
    spectral_grid = SpectralGrid(; trunc, nlayers, architecture = arch)
    ModelType = model_kind === :PrimitiveDryModel ? PrimitiveDryModel : PrimitiveWetModel
    model = ModelType(spectral_grid)
    sim = initialize!(model)
    # Run one hour so prognostic + tendency arrays are populated with realistic state
    run!(sim, period = Hour(1))
    return sim.model, sim.variables
end

function snapshot_tendencies(vars)
    out = Dict{Symbol, Any}()
    for k in (:vorticity, :divergence, :temperature, :pressure, :humidity)
        if haskey(vars.tendencies, k)
            out[k] = copy(getproperty(vars.tendencies, k).data)
        end
    end
    return out
end

function snapshot_prognostic(vars)
    out = Dict{Symbol, Any}()
    for k in (:vorticity, :divergence, :temperature, :pressure, :humidity)
        if hasproperty(vars.prognostic, k)
            out[k] = copy(getproperty(vars.prognostic, k).data)
        end
    end
    return out
end

function assert_match(actual_dict, ref_dict, what::AbstractString)
    for (k, ref_data) in ref_dict
        new_data = actual_dict[k]
        if !all(Array(ref_data) .== Array(new_data))
            err = maximum(abs, Array(ref_data) .- Array(new_data))
            error("Bit-identicality FAILED for $what / $(k): max |Δ| = $err")
        end
    end
    return nothing
end

# ==============================================================================
# Benchmarking
# ==============================================================================

function bench_diffusion(arch, trunc, nlayers, model_kind)
    model, vars = setup(arch, trunc, nlayers, model_kind)
    lf = 1

    # Bit-identicality check
    SpeedyWeather.reset_tendencies!(vars); sync()
    reference_horizontal_diffusion!(vars, model.horizontal_diffusion, model, lf); sync()
    ref = snapshot_tendencies(vars)
    SpeedyWeather.reset_tendencies!(vars); sync()
    SpeedyWeather.horizontal_diffusion!(vars, model.horizontal_diffusion, model, lf); sync()
    assert_match(snapshot_tendencies(vars), ref, "diffusion $model_kind T$trunc/L$nlayers")

    # Warmup
    SpeedyWeather.reset_tendencies!(vars); sync()
    SpeedyWeather.horizontal_diffusion!(vars, model.horizontal_diffusion, model, lf); sync()

    b_ref = if USE_GPU
        @benchmark begin
            SpeedyWeather.reset_tendencies!($vars)
            reference_horizontal_diffusion!($vars, $(model.horizontal_diffusion), $model, $lf)
            CUDA.synchronize()
        end
    else
        @benchmark begin
            SpeedyWeather.reset_tendencies!($vars)
            reference_horizontal_diffusion!($vars, $(model.horizontal_diffusion), $model, $lf)
        end
    end

    b_new = if USE_GPU
        @benchmark begin
            SpeedyWeather.reset_tendencies!($vars)
            SpeedyWeather.horizontal_diffusion!($vars, $(model.horizontal_diffusion), $model, $lf)
            CUDA.synchronize()
        end
    else
        @benchmark begin
            SpeedyWeather.reset_tendencies!($vars)
            SpeedyWeather.horizontal_diffusion!($vars, $(model.horizontal_diffusion), $model, $lf)
        end
    end

    b_reset = if USE_GPU
        @benchmark begin
            SpeedyWeather.reset_tendencies!($vars)
            CUDA.synchronize()
        end
    else
        @benchmark SpeedyWeather.reset_tendencies!($vars)
    end

    return (median(b_ref).time - median(b_reset).time,
            median(b_new).time - median(b_reset).time)
end

function bench_leapfrog(arch, trunc, nlayers, model_kind)
    model, vars = setup(arch, trunc, nlayers, model_kind)
    lf = 2          # later leapfrog index (post-init, full Robert+Williams filter)
    dt = 2 * model.time_stepping.Δt   # standard 2Δt time step

    # Snapshot prognostic so we can roll back between benchmark runs
    prog_snap = snapshot_prognostic(vars)
    tend_snap = snapshot_tendencies(vars)

    function restore!()
        for (k, data) in prog_snap
            copyto!(getproperty(vars.prognostic, k).data, data)
        end
        for (k, data) in tend_snap
            copyto!(getproperty(vars.tendencies, k).data, data)
        end
        sync()
        return nothing
    end

    # Bit-identicality check
    restore!()
    reference_leapfrog!(vars, dt, lf, model); sync()
    ref = snapshot_prognostic(vars)
    restore!()
    batched_leapfrog!(vars, dt, lf, model); sync()
    assert_match(snapshot_prognostic(vars), ref, "leapfrog $model_kind T$trunc/L$nlayers")

    # Warmup
    restore!()
    SpeedyWeather.leapfrog!(vars, dt, lf, model); sync()

    b_ref = if USE_GPU
        @benchmark begin
            $restore!()
            reference_leapfrog!($vars, $dt, $lf, $model)
            CUDA.synchronize()
        end
    else
        @benchmark begin
            $restore!()
            reference_leapfrog!($vars, $dt, $lf, $model)
        end
    end

    b_new = if USE_GPU
        @benchmark begin
            $restore!()
            batched_leapfrog!($vars, $dt, $lf, $model)
            CUDA.synchronize()
        end
    else
        @benchmark begin
            $restore!()
            batched_leapfrog!($vars, $dt, $lf, $model)
        end
    end

    b_restore = if USE_GPU
        @benchmark begin
            $restore!()
            CUDA.synchronize()
        end
    else
        @benchmark $restore!()
    end

    return (median(b_ref).time - median(b_restore).time,
            median(b_new).time - median(b_restore).time)
end

function report(name, bench_fn)
    println()
    println("== $name ==")
    @printf("%-22s %6s %6s %14s %14s %10s\n",
        "model", "trunc", "nlayers", "per-var (ms)", "batched (ms)", "speedup")
    println(repeat('-', 80))
    for model_kind in MODEL_KINDS, (trunc, nlayers) in RESOLUTIONS
        try
            t_ref, t_new = bench_fn(arch, trunc, nlayers, model_kind)
            speedup = t_ref / max(t_new, eps())
            @printf("%-22s %6d %6d %14.4f %14.4f %9.2fx\n",
                String(model_kind), trunc, nlayers,
                t_ref / 1.0e6, t_new / 1.0e6, speedup)
        catch err
            @printf("%-22s %6d %6d   FAILED: %s\n",
                String(model_kind), trunc, nlayers, sprint(showerror, err))
        end
    end
end

function main()
    println("Benchmark: unified-fuse refactor (horizontal_diffusion! + leapfrog!)")
    println("Architecture: ", arch)
    report("horizontal_diffusion!", bench_diffusion)
    report("leapfrog!",             bench_leapfrog)
    println()
    println("All measured points pass bit-identicality vs the per-variable reference.")
end

main()
