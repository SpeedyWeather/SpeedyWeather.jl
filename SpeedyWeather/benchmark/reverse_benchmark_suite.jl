# Reverse-mode (AD) benchmarks with Reactant: wallclock time of n time steps forward and of the
# reverse-mode gradient over the same n time steps. Opt-in with `--reverse`, the driver includes this
# file only then, after `using Reactant, Enzyme`, as the `@compile` macro must be resolvable when this
# file is parsed. Plan: docs/dev/2026-10/2026-10-21-reactant-reverse-benchmarks.md

# The compile options of every compilation in SpeedyWeather plus `raise_first = true` for the
# differentiation, as in test/reactant/differentation.jl. The forward run is compiled with them too, so
# that Gradient/Forward measures the differentiation, not different compile options.
const REACTANT_AD_COMPILE_OPTIONS = SpeedyWeather.reactant_compile_options(; raise_first = true)

"""Scalar loss to differentiate: squared norm of the prognostic vorticity after `n` time steps, as in
test/reactant/differentation.jl. `options` are the keyword arguments of `time_steps!` (`checkpointing`),
a positional argument as `Enzyme.autodiff` takes no keywords."""
function loss_n_steps(vars, model, n, options = (;))
    SpeedyWeather.time_steps!(Simulation(vars, model), n; options...)
    return sum(abs2, vars.prognostic.vorticity.data)
end

@kwdef mutable struct BenchmarkSuiteReverse <: AbstractBenchmarkSuite
    title::String
    nruns::Int = 1
    model::Vector = fill(PrimitiveWetModel, nruns)
    truncation::Vector{Int} = fill(SpeedyWeather.DEFAULT_TRUNCATION, nruns)
    nlayers::Vector{Int} = default_nlayers(model)
    nlat::Vector{Int} = fill(0, nruns)
    # number of time steps of the forward run and the gradient
    nsteps::Vector{Int} = fill(50, nruns)
    # checkpointing of the traced time loop for the gradient, `checkpointing` keyword of `time_steps!`.
    # Revolve with a fixed number of checkpoints: storing every step needs memory ∝ nsteps,
    # Reactant.Periodic doesn't compile for the PrimitiveWetModel (test/reactant/differentation.jl)
    checkpointing::Vector = fill(Reactant.Binomial(32), nruns)
    Δt::Vector{Float64} = fill(0.0, nruns)
    # wallclock times [s] of the forward run and the gradient (minimum over the samples), and their compilation
    forward_time::Vector{Float64} = fill(NaN, nruns)
    gradient_time::Vector{Float64} = fill(NaN, nruns)
    forward_compile_time::Vector{Float64} = fill(NaN, nruns)
    gradient_compile_time::Vector{Float64} = fill(NaN, nruns)
    # "ok", "non-finite or zero gradient" or "failed: <error type>"
    status::Vector{String} = fill("", nruns)
    # timed samples after one warm-up, no further samples after one that took longer than `max_sample_time` [s]
    nsamples::Int = 3
    max_sample_time::Float64 = 60
    architecture::Any = SpeedyWeather.CPU()
end

"""Minimum over the samples of `timed_sample()`, which returns the wallclock time of one sample (and
does its untimed setup), after one warm-up call; see `nsamples` and `max_sample_time` of `suite`."""
function minimum_sample_time(timed_sample, suite::BenchmarkSuiteReverse)
    timed_sample()      # warm-up
    times = Float64[]
    while length(times) < suite.nsamples && (isempty(times) || times[end] <= suite.max_sample_time)
        push!(times, timed_sample())
    end
    return minimum(times)
end

# Wait for the device to have written `array`, an output of a compiled program. Views (the fused
# variables) wait for their parent array.
synchronize_array(array::SubArray) = synchronize_array(parent(array))
synchronize_array(array) = Reactant.synchronize(array)

function run_benchmark!(suite::BenchmarkSuiteReverse, i::Integer)
    Model = suite.model[i]
    truncation = suite.truncation[i]
    nlayers = suite.nlayers[i]
    n = suite.nsteps[i]
    options = (; checkpointing = suite.checkpointing[i])
    architecture = suite.architecture

    spectral_grid = SpectralGrid(; truncation, nlayers, architecture)
    suite.nlat[i] = spectral_grid.nlat
    model = build_model(Model, spectral_grid, architecture)     # with the SpectralTransform
    suite.Δt[i] = model.time_stepping.Δt

    # `initialize!(simulation; steps)` scales the prognostic variables, so initialize once and only
    # run copies of the initialized variables
    simulation = initialize!(model)
    initialize!(simulation; steps = n)
    vars0 = simulation.variables

    # FORWARD: n time steps and the loss, one compiled program
    vars = deepcopy(vars0)
    suite.forward_compile_time[i] = @elapsed forward! = @compile compile_options = REACTANT_AD_COMPILE_OPTIONS loss_n_steps(vars, model, n)
    @info @sprintf("     compiled %d time steps forward in %.1f s", n, suite.forward_compile_time[i])
    suite.forward_time[i] = minimum_sample_time(suite) do
        GC.gc()                                         # free the device memory of earlier samples
        vars = deepcopy(vars0)
        return @elapsed Float64(forward!(vars, model, n))   # reading the loss waits for the device
    end

    # GRADIENT of the loss with respect to all variables before the n time steps (incl. the forward pass)
    vars, dvars = deepcopy(vars0), Enzyme.make_zero(vars0)
    suite.gradient_compile_time[i] = @elapsed gradient! = @compile compile_options = REACTANT_AD_COMPILE_OPTIONS Enzyme.autodiff(
        Reverse, loss_n_steps, Active, Duplicated(vars, dvars), Const(model), Const(n), Const(options),
    )
    @info @sprintf("     compiled the gradient of %d time steps in %.1f s", n, suite.gradient_compile_time[i])
    suite.gradient_time[i] = minimum_sample_time(suite) do
        GC.gc()
        vars, dvars = deepcopy(vars0), Enzyme.make_zero(vars0)     # a new shadow, autodiff accumulates into it
        return @elapsed begin
            gradient!(Reverse, loss_n_steps, Active, Duplicated(vars, dvars), Const(model), Const(n), Const(options))
            synchronize_array(dvars.prognostic.vorticity.data)
        end
    end

    gradient = on_architecture(SpeedyWeather.CPU(), dvars.prognostic.vorticity.data)
    suite.status[i] = all(isfinite, gradient) && any(!iszero, gradient) ? "ok" : "non-finite or zero gradient"
    @info @sprintf("     %d steps: forward %.3g s, gradient %.3g s", n, suite.forward_time[i], suite.gradient_time[i])
    return
end

checkpointing_label(checkpointing::Reactant.Binomial) = "Binomial($(checkpointing.budget))"
checkpointing_label(checkpointing::Reactant.Periodic) = "Periodic($(checkpointing.n))"
checkpointing_label(checkpointing) = string(checkpointing)

# seconds with 3 significant digits, "—" if not available (e.g. a failed run)
format_time(t) = isfinite(t) ? @sprintf("%.3g s", t) : "—"

function write_results(md, suite::BenchmarkSuiteReverse)
    write(md, "\n### $(suite.title)\n\n")
    write(md, "Wallclock time of the time steps forward and of the reverse-mode gradient with respect to the ")
    write(md, "initial state over the same time steps, minimum over up to $(suite.nsamples) runs, ")
    write(md, "compilation excluded. The loss is the squared norm of the vorticity after the time steps.\n\n")

    print_status = any(!=("ok"), suite.status)
    column_header = "| truncation | L | Rings | Δt | Steps | Checkpointing | Forward | Gradient | Gradient/Forward | Compile forward | Compile gradient "
    column_header *= print_status ? "| Status |" : "|"
    ncolumns = length(findall('|', column_header)) - 1
    write(md, "$column_header\n")
    write(md, repeat("| --- ", ncolumns) * "|\n")

    for i in 1:suite.nruns
        ratio = suite.gradient_time[i] / suite.forward_time[i]
        row = "| $(suite.truncation[i]) | $(suite.nlayers[i]) | $(suite.nlat[i]) | $(round(Int, suite.Δt[i])) "
        row *= "| $(suite.nsteps[i]) | $(checkpointing_label(suite.checkpointing[i])) "
        row *= "| $(format_time(suite.forward_time[i])) | $(format_time(suite.gradient_time[i])) "
        row *= "| $(isfinite(ratio) ? @sprintf("%.1f", ratio) : "—") "
        row *= "| $(format_seconds(suite.forward_compile_time[i])) | $(format_seconds(suite.gradient_compile_time[i])) "
        row *= print_status ? "| $(suite.status[i]) |" : "|"
        write(md, "$row\n")
    end
    return
end

"""The numbers of the reverse-mode suite, the `reverse` field of the arch record in the results store."""
function reverse_data(suite::BenchmarkSuiteReverse)
    return Dict(
        "truncation" => collect(suite.truncation),
        "nlayers" => collect(suite.nlayers),
        "nlat" => collect(suite.nlat),
        "nsteps" => collect(suite.nsteps),
        "checkpointing" => checkpointing_label.(suite.checkpointing),
        "dt" => map(json_safe, suite.Δt),
        "forward_time" => map(json_safe, suite.forward_time),
        "gradient_time" => map(json_safe, suite.gradient_time),
        "forward_compile_time" => map(json_safe, suite.forward_compile_time),
        "gradient_compile_time" => map(json_safe, suite.gradient_compile_time),
        "status" => collect(suite.status),
    )
end
