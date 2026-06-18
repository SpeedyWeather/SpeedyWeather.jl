using BenchmarkTools
using SpeedyWeather: synchronize, _jit, ReactantDevice, MatrixSpectralTransform

# Build the model with the right spectral transform for `arch`.
# Reactant requires MatrixSpectralTransform; otherwise honour `transform_kind`
# (`:default` for FFT+Legendre, `:matrix` for MatrixSpectralTransform).
function build_model(Model, spectral_grid, arch; transform_kind::Symbol = :default, kwargs...)
    if arch isa ReactantDevice
        M = MatrixSpectralTransform(spectral_grid)
        return Model(spectral_grid; spectral_transform = M, feedback = nothing, output = nothing, kwargs...)
    elseif transform_kind == :matrix
        M = MatrixSpectralTransform(spectral_grid)
        return Model(spectral_grid; spectral_transform = M, kwargs...)
    elseif transform_kind == :default
        return Model(spectral_grid; kwargs...)
    else
        error("Unknown transform_kind = $transform_kind. Use :default or :matrix.")
    end
end

# Format SYPD: one decimal digit if < 10, integer otherwise. Used by both the
# per-arch tables and the cross-arch overview so the column reads consistently.
function format_sypd(s)
    (s isa Number && isfinite(s)) || return "0"
    return s < 10 ? string(round(s; digits = 1)) : string(round(Int, s))
end

# Resolve a NamedTuple of per-run model kwargs against `spectral_grid`.
# Convention: a `Type` value is instantiated with `spectral_grid` (every
# component has a `T(spectral_grid)` constructor); a
# `Function` value is called with `spectral_grid` (for parameterised cases
# like `sg -> WhichZenith(sg, my_planet)`); everything else is passed through.
_resolve_kwarg_value(v::Type, sg) = v(sg)
_resolve_kwarg_value(v::Function, sg) = v(sg)
_resolve_kwarg_value(v, _) = v
function resolve_model_kwargs(nt::NamedTuple, spectral_grid)
    return NamedTuple{keys(nt)}(map(v -> _resolve_kwarg_value(v, spectral_grid), values(nt)))
end

abstract type AbstractBenchmarkSuite end
@kwdef mutable struct BenchmarkSuite <: AbstractBenchmarkSuite
    title::String
    nruns::Int = 1
    model::Vector = fill(PrimitiveWetModel, nruns)
    NF::Vector = fill(SpeedyWeather.DEFAULT_NF, nruns)
    trunc::Vector{Int} = fill(SpeedyWeather.DEFAULT_TRUNC, nruns)
    nlayers::Vector{Int} = default_nlayers(model)
    Grid::Vector = fill(SpeedyWeather.DEFAULT_GRID, nruns)
    nlat::Vector{Int} = fill(0, nruns)
    dynamics::Vector{Bool} = fill(true, nruns)
    physics::Vector{Bool} = fill(true, nruns)
    spectral_transform::Vector{Symbol} = fill(:default, nruns)
    # Per-run extra kwargs to splat into the model constructor. Values that
    # are `Type` or `Function` are resolved against `spectral_grid` first;
    # Use this to vary arbitrary model components, e.g.
    # `model_kwargs = [(convection = NoConvection,), (convection = BettsMillerConvection,)]`.
    model_kwargs::Vector = fill(NamedTuple(), nruns)
    SYPD::Vector{Float64} = fill(0.0, nruns)
    Δt::Vector{Float64} = fill(0.0, nruns)
    memory::Vector{Int} = fill(0, nruns)
    architecture::Any = SpeedyWeather.CPU()
end

default_nlayers(::Type{<:Barotropic}) = 1
default_nlayers(::Type{<:ShallowWater}) = 1
default_nlayers(::Type{<:PrimitiveEquation}) = 8
default_nlayers(models) = [default_nlayers(model) for model in models]

# this should return number of timesteps so that every simulation
# only takes seconds
n_timesteps(trunc, nlayers) = max(10, round(Int, 4.0e8 / trunc^3 / nlayers^2))

function run_benchmark_suite!(suite::BenchmarkSuite)
    for i in 1:suite.nruns

        # unpack
        Model = suite.model[i]
        NF = suite.NF[i]
        trunc = suite.trunc[i]
        nlayers = suite.nlayers[i]
        Grid = suite.Grid[i]
        dynamics = suite.dynamics[i]
        physics = suite.physics[i]
        transform_kind = suite.spectral_transform[i]
        architecture = suite.architecture

        spectral_grid = SpectralGrid(; NF, trunc, Grid, nlayers, architecture)
        suite.nlat[i] = spectral_grid.nlat

        extra_components = resolve_model_kwargs(suite.model_kwargs[i], spectral_grid)
        model = build_model(Model, spectral_grid, architecture; transform_kind, extra_components..., feedback = Feedback(verbose = true))
        if Model <: PrimitiveEquation
            model.dynamics = dynamics
            model.dynamics_only = !physics
        else
            suite.dynamics[i] = true
            suite.physics[i] = false
        end

        simulation = initialize!(model)
        suite.memory[i] = Base.summarysize(simulation)

        nsteps = n_timesteps(trunc, nlayers)
        period = Second(round(Int, model.time_stepping.Δt_sec * (nsteps + 1)))
        run!(simulation; period)
        synchronize(architecture)

        time_elapsed = model.feedback.progress_meter.tlast - model.feedback.progress_meter.tinit
        sypd = model.time_stepping.Δt_sec * nsteps / (time_elapsed * 365.25)

        suite.Δt[i] = model.time_stepping.Δt_sec
        suite.SYPD[i] = sypd
    end
    return
end

function write_results(md, suite::BenchmarkSuite)

    write(md, "\n### $(suite.title)\n\n")

    print_NF = any(suite.NF .!= suite.NF[1])
    print_Grid = any(suite.Grid .!= suite.Grid[1])
    print_nlat = any(suite.nlat .!= suite.nlat[1])
    print_dynamics = any(suite.dynamics .!= suite.dynamics[1])
    print_physics = any(suite.physics .!= suite.physics[1])
    print_transform = any(suite.spectral_transform .!= suite.spectral_transform[1])

    column_header = "| Model "
    column_header *= print_NF ? "| NF " : ""
    column_header *= "| T "
    column_header *= "| L "
    column_header *= print_Grid ? "| Grid " : ""
    column_header *= print_nlat ? "| Rings " : ""
    column_header *= print_dynamics ? "| Dynamics " : ""
    column_header *= print_physics ? "| Physics " : ""
    column_header *= print_transform ? "| Transform " : ""
    column_header *= "| Δt | SYPD | Memory|"

    ncolumns = length(findall('|', column_header)) - 1
    second_row = repeat("| --- ", ncolumns) * "|"

    write(md, "$column_header\n")
    write(md, "$second_row\n")

    for i in 1:suite.nruns

        row = "| $(suite.model[i]) "
        row *= print_NF ? "| $(suite.NF[i]) " : ""
        row *= "| $(suite.trunc[i]) "
        row *= "| $(suite.nlayers[i]) "
        row *= print_Grid ? "| $(suite.Grid[i]) " : ""
        row *= print_nlat ? "| $(suite.nlat[i]) " : ""
        row *= print_dynamics ? "| $(suite.dynamics[i]) " : ""
        row *= print_physics ? "| $(suite.physics[i]) " : ""
        row *= print_transform ? "| $(suite.spectral_transform[i]) " : ""

        Δt = round(Int, suite.Δt[i])
        SYPD = format_sypd(suite.SYPD[i])
        memory = prettymemory(suite.memory[i])
        row *= "| $Δt | $SYPD | $memory |"

        write(md, "$row\n")
    end
    return
end

abstract type AbstractBenchmarkSuiteTimed <: AbstractBenchmarkSuite end

@kwdef mutable struct BenchmarkSuiteDynamics <: AbstractBenchmarkSuiteTimed
    title::String
    nruns::Int = 1
    model::Vector = fill(PrimitiveWetModel, nruns)
    NF::Vector = fill(SpeedyWeather.DEFAULT_NF, nruns)
    trunc::Vector{Int} = fill(SpeedyWeather.DEFAULT_TRUNC, nruns)
    nlayers::Vector{Int} = default_nlayers(model)
    Grid::Vector = fill(SpeedyWeather.DEFAULT_GRID, nruns)
    nlat::Vector{Int} = fill(0, nruns)
    function_names::Vector{String} = default_function_names()
    time::Vector{Vector{Float64}} = [fill(0.0, length(function_names)) for i in 1:nruns]
    memory::Vector{Vector{Int}} = [fill(0, length(function_names)) for i in 1:nruns]
    allocs::Vector{Vector{Int}} = [fill(0, length(function_names)) for i in 1:nruns]
    architecture::Any = SpeedyWeather.CPU()
end

default_function_names() = ["pressure_gradient_flux!", "linear_virtual_temperature!", "geopotential!", "vertical_integration!", "surface_pressure_tendency!", "vertical_velocity!", "linear_pressure_gradient!", "vertical_advection!", "vordiv_tendencies!", "temperature_tendency!", "humidity_tendency!", "bernoulli_potential!"]

function add_results!(suite::AbstractBenchmarkSuiteTimed, trial::BenchmarkTools.Trial, i_run::Integer, i_func::Integer)

    t = minimum(trial)

    suite.time[i_run][i_func] = t.time
    suite.memory[i_run][i_func] = t.memory
    suite.allocs[i_run][i_func] = t.allocs

    return suite
end

# Run one @benchmark and store the result. On failure, record NaN/0 and warn
# (used for individual function benchmarks that may not be GPU-compatible).
function safe_benchmark!(f, suite::AbstractBenchmarkSuiteTimed, i_run::Integer, i_func::Integer)
    name = suite.function_names[i_func]
    try
        trial = f()
        add_results!(suite, trial, i_run, i_func)
    catch err
        @warn "Benchmark for $name failed on $(suite.architecture); recording N/A" exception = (err, catch_backtrace())
        suite.time[i_run][i_func] = NaN
        suite.memory[i_run][i_func] = 0
        suite.allocs[i_run][i_func] = 0
    end
    return suite
end

function run_benchmark_suite!(suite::BenchmarkSuiteDynamics)
    arch = suite.architecture

    for i in 1:suite.nruns

        Model = suite.model[i]
        NF = suite.NF[i]
        trunc = suite.trunc[i]
        nlayers = suite.nlayers[i]
        Grid = suite.Grid[i]

        spectral_grid = SpectralGrid(; NF, trunc, Grid, nlayers, architecture = arch)
        suite.nlat[i] = spectral_grid.nlat

        model = build_model(Model, spectral_grid, arch)

        simulation = initialize!(model)

        vars, model = SpeedyWeather.unpack(simulation)
        (; orography, geometry, spectral_transform, geopotential, atmosphere, implicit, time_stepping) = model

        # Each benchmark sample also calls `synchronize(arch)` to wait for the device.
        safe_benchmark!(suite, i, 1) do
            @benchmark (_jit($arch, SpeedyWeather.pressure_gradient_flux!, $vars, $spectral_transform, $time_stepping); synchronize($arch))
        end
        safe_benchmark!(suite, i, 2) do
            @benchmark (_jit($arch, SpeedyWeather.linear_virtual_temperature!, $vars, $model); synchronize($arch))
        end
        safe_benchmark!(suite, i, 3) do
            @benchmark (_jit($arch, SpeedyWeather.geopotential!, $vars, $geopotential, $orography); synchronize($arch))
        end
        safe_benchmark!(suite, i, 4) do
            @benchmark (_jit($arch, SpeedyWeather.vertical_integration!, $vars, $geometry, $time_stepping); synchronize($arch))
        end
        safe_benchmark!(suite, i, 5) do
            @benchmark (_jit($arch, SpeedyWeather.surface_pressure_tendency!, $vars, $spectral_transform, $time_stepping); synchronize($arch))
        end
        safe_benchmark!(suite, i, 6) do
            @benchmark (_jit($arch, SpeedyWeather.vertical_velocity!, $vars, $geometry, $time_stepping); synchronize($arch))
        end
        safe_benchmark!(suite, i, 7) do
            @benchmark (_jit($arch, SpeedyWeather.linear_pressure_gradient!, $vars, $atmosphere, $implicit, $time_stepping); synchronize($arch))
        end
        safe_benchmark!(suite, i, 8) do
            @benchmark (_jit($arch, SpeedyWeather.vertical_advection!, $vars, $model); synchronize($arch))
        end
        safe_benchmark!(suite, i, 9) do
            @benchmark (_jit($arch, SpeedyWeather.vordiv_tendencies!, $vars, $model); synchronize($arch))
        end
        safe_benchmark!(suite, i, 10) do
            @benchmark (_jit($arch, SpeedyWeather.temperature_tendency!, $vars, $model); synchronize($arch))
        end
        safe_benchmark!(suite, i, 11) do
            @benchmark (_jit($arch, SpeedyWeather.humidity_tendency!, $vars, $model); synchronize($arch))
        end
        safe_benchmark!(suite, i, 12) do
            @benchmark (_jit($arch, SpeedyWeather.bernoulli_potential!, $vars, $spectral_transform, $time_stepping); synchronize($arch))
        end
    end

    return suite
end

function write_results(md, suite::AbstractBenchmarkSuiteTimed)

    write(md, "\n### $(suite.title)\n\n")

    for i_run in 1:suite.nruns

        title_row = "$(suite.model[i_run]) "
        title_row *= "| $(suite.NF[i_run]) "
        title_row *= "| T$(suite.trunc[i_run]) "
        title_row *= "L$(suite.nlayers[i_run]) "
        title_row *= "| $(suite.Grid[i_run]) "
        title_row *= "| $(suite.nlat[i_run]) Rings"

        write(md, "\n#### $title_row\n\n")

        column_header = "| Function "
        column_header *= "| Time "
        column_header *= "| Memory "
        column_header *= "| Allocations |"

        ncolumns = length(findall('|', column_header)) - 1
        second_row = repeat("| --- ", ncolumns) * "|"

        write(md, "$column_header\n")
        write(md, "$second_row\n")

        for i_func in 1:length(suite.function_names)

            t = suite.time[i_run][i_func]
            if isnan(t)
                time = "N/A"
                memory = "N/A"
                allocs = "N/A"
            else
                time = BenchmarkTools.prettytime(t)
                memory = BenchmarkTools.prettymemory(suite.memory[i_run][i_func])
                allocs = string(suite.allocs[i_run][i_func])
            end

            row = "| $(suite.function_names[i_func]) "
            row *= "| $time"
            row *= "| $memory"
            row *= "| $allocs |"

            write(md, "$row\n")
        end
    end
    return
end
