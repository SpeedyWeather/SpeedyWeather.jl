using BenchmarkTools 

abstract type AbstractBenchmarkSuite end 
@kwdef mutable struct BenchmarkSuite <: AbstractBenchmarkSuite
    title::String
    nruns::Int = 1
    model::Vector = fill(PrimitiveWetModel, nruns)
    NF::Vector = fill(SpeedyWeather.DEFAULT_NF, nruns)
    trunc::Vector{Int} = fill(SpeedyWeather.DEFAULT_TRUNC, nruns)
    nlev::Vector{Int} = default_nlev(model)
    Grid::Vector = fill(SpeedyWeather.DEFAULT_GRID, nruns)
    nlat::Vector{Int} = fill(0, nruns)
    dynamics::Vector{Bool} = fill(true, nruns)
    physics::Vector{Bool} = fill(true, nruns)
    SYPD::Vector{Float64} = fill(0.0, nruns)
    Δt::Vector{Float64} = fill(0.0, nruns)
    memory::Vector{Int} = fill(0, nruns)
end

default_nlev(::Type{<:Barotropic}) = 1
default_nlev(::Type{<:ShallowWater}) = 1
default_nlev(::Type{<:PrimitiveEquation}) = 8
default_nlev(models) = [default_nlev(model) for model in models]

# this should return number of timesteps so that every simulation
# only takes seconds
n_timesteps(trunc, nlev) = max(10, round(Int, 4e8/trunc^3/nlev^2))

function run_benchmark_suite!(suite::BenchmarkSuite)
    for i in 1:suite.nruns

        # unpack 
        Model = suite.model[i]
        NF = suite.NF[i]
        trunc = suite.trunc[i]
        nlayers = suite.nlev[i]
        Grid = suite.Grid[i]
        dynamics = suite.dynamics[i]
        physics = suite.physics[i]

        spectral_grid = SpectralGrid(;NF, trunc, Grid, nlayers)
        suite.nlat[i] = spectral_grid.nlat

        model = Model(;spectral_grid)
        if Model <: PrimitiveEquation
            model.physics = physics
            model.dynamics = dynamics
        else
            suite.dynamics[i] = true
            suite.physics[i] = false
        end

        simulation = initialize!(model)
        suite.memory[i] = Base.summarysize(simulation)

        nsteps = n_timesteps(trunc, nlayers)
        period = Second(round(Int,model.time_stepping.Δt_sec * (nsteps+1)))
        run!(simulation; period)

        time_elapsed = model.feedback.progress_meter.tlast - model.feedback.progress_meter.tinit
        sypd = model.time_stepping.Δt_sec*nsteps / (time_elapsed * 365.25)

        suite.Δt[i] = model.time_stepping.Δt_sec
        suite.SYPD[i] = sypd
    end
end

function write_results(md, suite::BenchmarkSuite)

    write(md, "\n## $(suite.title)\n\n")

    print_NF = any(suite.NF .!= suite.NF[1])
    print_Grid = any(suite.Grid .!= suite.Grid[1])
    print_nlat = any(suite.nlat .!= suite.nlat[1]) 
    print_dynamics = any(suite.dynamics .!= suite.dynamics[1])
    print_physics = any(suite.physics .!= suite.physics[1])

    column_header = "| Model "
    column_header *= print_NF ? "| NF " : ""
    column_header *= "| T "
    column_header *= "| L "
    column_header *= print_Grid ? "| Grid " : ""
    column_header *= print_nlat ? "| Rings " : ""
    column_header *= print_dynamics ? "| Dynamics " : ""
    column_header *= print_physics ? "| Physics " : ""
    column_header *= "| Δt | SYPD | Memory|"

    ncolumns = length(findall('|',column_header)) - 1
    second_row = repeat("| - ", ncolumns) * "|"

    write(md,"$column_header\n")
    write(md,"$second_row\n")

    for i in 1:suite.nruns

        row = "| $(suite.model[i]) "
        row *= print_NF ? "| $(suite.NF[i]) " : ""
        row *= "| $(suite.trunc[i]) "
        row *= "| $(suite.nlev[i]) "
        row *= print_Grid ? "| $(suite.Grid[i]) " : ""
        row *= print_nlat ? "| $(suite.nlat[i]) " : ""
        row *= print_dynamics ? "| $(suite.dynamics[i]) " : ""
        row *= print_physics  ? "| $(suite.physics[i]) " : ""

        Δt = round(Int,suite.Δt[i])
        sypd = suite.SYPD[i]
        SYPD = isfinite(sypd) ? round(Int, sypd) : 0
        memory = prettymemory(suite.memory[i])
        row *= "| $Δt | $SYPD | $memory |"

        write(md,"$row\n")
    end
end 

@kwdef mutable struct BenchmarkSuiteDynamics <: AbstractBenchmarkSuite
    title::String
    nruns::Int = 1
    model::Vector = fill(PrimitiveWetModel, nruns)
    NF::Vector = fill(SpeedyWeather.DEFAULT_NF, nruns)
    trunc::Vector{Int} = fill(SpeedyWeather.DEFAULT_TRUNC, nruns)
    nlev::Vector{Int} = default_nlev(model)
    Grid::Vector = fill(SpeedyWeather.DEFAULT_GRID, nruns)
    nlat::Vector{Int} = fill(0, nruns)
    dynamics::Vector{Bool} = fill(true, nruns)
    physics::Vector{Bool} = fill(true, nruns)
    function_names::Vector{String} = default_function_names()
    time_median::Vector{Vector{Float64}} = [fill(0.0, length(function_names)) for i=1:nruns]
    memory::Vector{Vector{Int}} = [fill(0, length(function_names)) for i=1:nruns]
    allocs::Vector{Vector{Int}} = [fill(0, length(function_names)) for i=1:nruns]
end 

default_function_names() = ["pressure_gradient_flux!", "linear_virtual_temperature!", "temperature_anomaly!", "geopotential!", "vertical_integration!", "surface_pressure_tendency!", "vertical_velocity!", "linear_pressure_gradient!", "vertical_advection!","vordiv_tendencies!", "temperature_tendency!", "humidity_tendency!", "bernoulli_potential!"]

function add_results!(suite::BenchmarkSuiteDynamics, trial::BenchmarkTools.Trial, i_run::Integer, i_func::Integer)

    t = median(trial)

    suite.time_median[i_run][i_func] = t.time 
    suite.memory[i_run][i_func] = t.memory 
    suite.allocs[i_run][i_func] = t.allocs

    return suite
end 

function run_benchmark_suite!(suite::BenchmarkSuiteDynamics)
    for i in 1:suite.nruns

        Model = suite.model[i]
        NF = suite.NF[i]
        trunc = suite.trunc[i]
        nlayers = suite.nlev[i]
        Grid = suite.Grid[i]
        dynamics = suite.dynamics[i]
        physics = suite.physics[i]

        spectral_grid = SpectralGrid(;NF, trunc, Grid, nlayers)
        suite.nlat[i] = spectral_grid.nlat

        model = Model(;spectral_grid)
        if Model <: PrimitiveEquation
            model.physics = physics
            model.dynamics = dynamics
        else
            suite.dynamics[i] = true
            suite.physics[i] = false
        end

        simulation = initialize!(model)

        diagn = simulation.diagnostic_variables
        progn = simulation.prognostic_variables 
        lf = 2 
        (; orography, geometry, spectral_transform, geopotential, atmosphere, implicit) = model
        lf_implicit = implicit.α == 0 ? lf : 1

        b = @benchmark SpeedyWeather.pressure_gradient_flux!($diagn, $progn, $lf, $spectral_transform)
        add_results!(suite, b, i, 1)

        b = @benchmark SpeedyWeather.linear_virtual_temperature!($diagn, $progn, $lf_implicit, $model)
        add_results!(suite, b, i, 2)

        b = @benchmark SpeedyWeather.temperature_anomaly!($diagn, $implicit)
        add_results!(suite, b, i, 3)

        b = @benchmark SpeedyWeather.geopotential!($diagn, $geopotential, $orography)
        add_results!(suite, b, i, 4)

        b = @benchmark SpeedyWeather.vertical_integration!($diagn, $progn, $lf_implicit, $geometry)
        add_results!(suite, b, i, 5)

        b = @benchmark SpeedyWeather.surface_pressure_tendency!($diagn, $spectral_transform)
        add_results!(suite, b, i, 6)

        b = @benchmark SpeedyWeather.vertical_velocity!($diagn, $geometry)
        add_results!(suite, b, i, 7)

        b = @benchmark SpeedyWeather.linear_pressure_gradient!($diagn, $progn, $lf_implicit, $atmosphere, $implicit)
        add_results!(suite, b, i, 8)

        b = @benchmark SpeedyWeather.vertical_advection!($diagn, $model)
        add_results!(suite, b, i, 9)

        b = @benchmark SpeedyWeather.vordiv_tendencies!($diagn, $model)
        add_results!(suite, b, i, 10)

        b = @benchmark SpeedyWeather.temperature_tendency!($diagn, $model)
        add_results!(suite, b, i, 11)

        b = @benchmark SpeedyWeather.humidity_tendency!($diagn, $model)
        add_results!(suite, b, i, 12)

        b = @benchmark SpeedyWeather.bernoulli_potential!($diagn, $spectral_transform)
        add_results!(suite, b, i, 13)
    end

    return suite
end 

function write_results(md, suite::BenchmarkSuiteDynamics)

    write(md, "\n## $(suite.title)\n\n")

    for i_run in 1:suite.nruns

        title_row = "$(suite.model[i_run]) "
        title_row *= "| $(suite.NF[i_run]) " 
        title_row *= "| T$(suite.trunc[i_run]) "
        title_row *= "L$(suite.nlev[i_run]) "
        title_row *= "| $(suite.Grid[i_run]) " 
        title_row *= "| $(suite.nlat[i_run]) Rings" 

        write(md, "\n### $title_row\n\n")

        column_header = "| Function "
        column_header *= "| Time "
        column_header *= "| Memory "
        column_header *= "| Allocations |" 

        ncolumns = length(findall('|',column_header)) - 1
        second_row = repeat("| - ", ncolumns) * "|"

        write(md,"$column_header\n")
        write(md,"$second_row\n")

        for i_func in 1:length(suite.function_names)

            time = BenchmarkTools.prettytime(suite.time_median[i_run][i_func])
            memory = BenchmarkTools.prettymemory(suite.memory[i_run][i_func])

            row = "| $(suite.function_names[i_func]) "
            row *= "| $time"
            row *= "| $memory"
            row *= "| $(suite.allocs[i_run][i_func]) |"

            write(md,"$row\n")
        end
    end
end 