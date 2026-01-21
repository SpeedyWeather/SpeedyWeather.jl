include("../benchmark_suite.jl")

@kwdef mutable struct BenchmarkSuiteTransform <: AbstractBenchmarkSuiteTimed
    title::String
    nruns::Int = 1
    model::Vector = fill(:CPU, nruns)
    NF::Vector = fill(SpeedyWeather.DEFAULT_NF, nruns)
    trunc::Vector{Int} = fill(SpeedyWeather.DEFAULT_TRUNC, nruns)
    nlayers::Vector{Int} = fill(SpeedyWeather.DEFAULT_NLAYERS, nruns)
    Grid::Vector = fill(SpeedyWeather.DEFAULT_GRID, nruns)
    nlat::Vector{Int} = fill(0, nruns)
    function_names::Vector{String} = default_function_names()
    time::Vector{Vector{Float64}} = [fill(0.0, length(function_names)) for i in 1:nruns]
    memory::Vector{Vector{Int}} = [fill(0, length(function_names)) for i in 1:nruns]
    allocs::Vector{Vector{Int}} = [fill(0, length(function_names)) for i in 1:nruns]
end

default_function_names() = ["forward_legendre", "inverse_legendre", "forward_fourier", "inverse_fourier"]


function generate_random_inputs(spectral_grid::SpectralGrid)
    (; spectrum, grid, nlayers, NF, architecture) = spectral_grid

    grids = on_architecture(architecture, rand(NF, grid, nlayers))
    specs = on_architecture(architecture, rand(Complex{NF}, spectrum, nlayers))

    return specs, grids
end

function run_benchmark_suite!(suite::BenchmarkSuiteTransform)

    for i in 1:suite.nruns
        NF = suite.NF[i]
        trunc = suite.trunc[i]
        nlayers = suite.nlayers[i]
        Grid = suite.Grid[i]
        architecture = suite.model[i]

        spectral_grid = SpectralGrid(; NF, trunc, Grid, nlayers, architecture)
        suite.nlat[i] = spectral_grid.nlat
        S = SpectralTransform(spectral_grid)

        specs, grids = generate_random_inputs(spectral_grid)

        println("Running benchmark for architecture=$architecture, trunc=$trunc, nlayers=$nlayers, NF=$NF \n")

        # Forward _legendre
        b = @benchmark CUDA.@sync SpeedyTransforms._legendre!($specs, $S.scratch_memory.north, $S.scratch_memory.south, $S.scratch_memory.column, $S)
        add_results!(suite, b, i, 1)

        # Inverse _legendre
        b = @benchmark CUDA.@sync SpeedyTransforms._legendre!($S.scratch_memory.north, $S.scratch_memory.south, $specs, $S.scratch_memory.column, $S)
        add_results!(suite, b, i, 2)

        # Forward _fourier!(..., grids, S)
        b = @benchmark CUDA.@sync SpeedyTransforms._fourier!($S.scratch_memory.north, $S.scratch_memory.south, $grids, $S)
        add_results!(suite, b, i, 3)

        # Reverse _fourier!(grids, ..., S)
        b = @benchmark CUDA.@sync SpeedyTransforms._fourier!($grids, $S.scratch_memory.north, $S.scratch_memory.south, $S)
        add_results!(suite, b, i, 4)
    end

    return suite
end

@kwdef mutable struct BenchmarkSuiteModel <: AbstractBenchmarkSuiteTimed
    title::String
    nruns::Int = 1
    model::Vector = fill(SpeedyWeather.CPU, nruns)
    NF::Vector = fill(SpeedyWeather.DEFAULT_NF, nruns)
    trunc::Vector{Int} = fill(SpeedyWeather.DEFAULT_TRUNC, nruns)
    nlayers::Vector{Int} = fill(SpeedyWeather.DEFAULT_NLAYERS, nruns)
    Grid::Vector = fill(SpeedyWeather.DEFAULT_GRID, nruns)
    nlat::Vector{Int} = fill(0, nruns)
    function_names::Vector{String} = ["parameterization_tendencies", "dynamics_tendencies", "implicit_correction"]
    time::Vector{Vector{Float64}} = [fill(0.0, length(function_names)) for i in 1:nruns]
    memory::Vector{Vector{Int}} = [fill(0, length(function_names)) for i in 1:nruns]
    allocs::Vector{Vector{Int}} = [fill(0, length(function_names)) for i in 1:nruns]
end

function run_benchmark_suite!(suite::BenchmarkSuiteModel)

    for i in 1:suite.nruns
        NF = suite.NF[i]
        trunc = suite.trunc[i]
        nlayers = suite.nlayers[i]
        Grid = suite.Grid[i]
        architecture = suite.model[i]
        lf = 2

        spectral_grid = SpectralGrid(; NF, trunc, Grid, nlayers, architecture)
        suite.nlat[i] = spectral_grid.nlat
        model = PrimitiveWetModel(spectral_grid)
        simulation = initialize!(model)
        # spin up
        run!(simulation; steps = 10)
        initialize!(simulation)
        progn, diagn, model = SpeedyWeather.unpack(simulation)
        fill!(diagn.tendencies, 0, typeof(model))

        println("Running benchmark for architecture=$architecture, trunc=$trunc, nlayers=$nlayers, NF=$NF \n")

        b = @benchmark CUDA.@sync SpeedyWeather.parameterization_tendencies!($diagn, $progn, $model)
        add_results!(suite, b, i, 1)

        b = @benchmark CUDA.@sync SpeedyWeather.dynamics_tendencies!($diagn, $progn, $lf, $model)
        add_results!(suite, b, i, 2)

        b = @benchmark CUDA.@sync SpeedyWeather.implicit_correction!($diagn, $progn, $model.implicit, $model)
        add_results!(suite, b, i, 3)
    end

    return suite
end

@kwdef mutable struct BenchmarkSuiteDynamics <: AbstractBenchmarkSuiteTimed
    title::String
    nruns::Int = 1
    model::Vector = fill(SpeedyWeather.CPU, nruns)
    NF::Vector = fill(SpeedyWeather.DEFAULT_NF, nruns)
    trunc::Vector{Int} = fill(SpeedyWeather.DEFAULT_TRUNC, nruns)
    nlayers::Vector{Int} = fill(SpeedyWeather.DEFAULT_NLAYERS, nruns)
    Grid::Vector = fill(SpeedyWeather.DEFAULT_GRID, nruns)
    nlat::Vector{Int} = fill(0, nruns)
    function_names::Vector{String} = ["forcing!", "drag!", "pressure_gradient_flux!", "linear_virtual_temperature!", "geopotential!", "vertical_integration!", "surface_pressure_tendency!", "vertical_velocity!", "linear_pressure_gradient!", "vertical_advection!", "vordiv_tendencies!", "temperature_tendency!", "humidity_tendency!", "bernoulli_potential!", "tracer_advection!", "transform! (forward)", "transform! (inverse)"]
    time::Vector{Vector{Float64}} = [fill(0.0, length(function_names)) for i in 1:nruns]
    memory::Vector{Vector{Int}} = [fill(0, length(function_names)) for i in 1:nruns]
    allocs::Vector{Vector{Int}} = [fill(0, length(function_names)) for i in 1:nruns]
end

function run_benchmark_suite!(suite::BenchmarkSuiteDynamics)

    for i in 1:suite.nruns
        NF = suite.NF[i]
        trunc = suite.trunc[i]
        nlayers = suite.nlayers[i]
        Grid = suite.Grid[i]
        architecture = suite.model[i]
        lf = 2

        spectral_grid = SpectralGrid(; NF, trunc, Grid, nlayers, architecture)
        suite.nlat[i] = spectral_grid.nlat
        model = PrimitiveWetModel(spectral_grid)
        simulation = initialize!(model)
        # spin up
        run!(simulation; steps = 10)
        initialize!(simulation)
        progn, diagn, model = SpeedyWeather.unpack(simulation)
        fill!(diagn.tendencies, 0, typeof(model))
        (; orography, geometry, spectral_transform, geopotential, atmosphere, implicit) = model

        println("Running benchmark for architecture=$architecture, trunc=$trunc, nlayers=$nlayers, NF=$NF \n")

        b = @benchmark CUDA.@sync SpeedyWeather.forcing!($diagn, $progn, $lf, $model)
        add_results!(suite, b, i, 1)

        b = @benchmark CUDA.@sync SpeedyWeather.drag!($diagn, $progn, $lf, $model)
        add_results!(suite, b, i, 2)

        lf_implicit = implicit.α == 0 ? lf : 1

        b = @benchmark CUDA.@sync SpeedyWeather.pressure_gradient_flux!($diagn, $progn, $lf, $spectral_transform)
        add_results!(suite, b, i, 3)

        b = @benchmark CUDA.@sync SpeedyWeather.linear_virtual_temperature!($diagn, $progn, $lf_implicit, $model)
        add_results!(suite, b, i, 4)

        diagn.grid.temp_grid.data .-= implicit.temp_profile'

        b = @benchmark CUDA.@sync SpeedyWeather.geopotential!($diagn, $geopotential, $orography)
        add_results!(suite, b, i, 5)

        b = @benchmark CUDA.@sync SpeedyWeather.vertical_integration!($diagn, $progn, $lf_implicit, $geometry)
        add_results!(suite, b, i, 6)

        b = @benchmark CUDA.@sync SpeedyWeather.surface_pressure_tendency!($diagn, $spectral_transform)
        add_results!(suite, b, i, 7)

        b = @benchmark CUDA.@sync SpeedyWeather.vertical_velocity!($diagn, $geometry)
        add_results!(suite, b, i, 8)

        b = @benchmark CUDA.@sync SpeedyWeather.linear_pressure_gradient!($diagn, $progn, $lf_implicit, $atmosphere, $implicit)
        add_results!(suite, b, i, 9)

        b = @benchmark CUDA.@sync SpeedyWeather.vertical_advection!($diagn, $model)
        add_results!(suite, b, i, 10)

        b = @benchmark CUDA.@sync SpeedyWeather.vordiv_tendencies!($diagn, $model)
        add_results!(suite, b, i, 11)

        b = @benchmark CUDA.@sync SpeedyWeather.temperature_tendency!($diagn, $model)
        add_results!(suite, b, i, 12)

        b = @benchmark CUDA.@sync SpeedyWeather.humidity_tendency!($diagn, $model)
        add_results!(suite, b, i, 13)

        b = @benchmark CUDA.@sync SpeedyWeather.bernoulli_potential!($diagn, $spectral_transform)
        add_results!(suite, b, i, 14)

        b = @benchmark CUDA.@sync SpeedyWeather.tracer_advection!($diagn, $model)
        add_results!(suite, b, i, 15)

        diagn.grid.temp_grid.data .+= implicit.temp_profile'

        # now just transforms as comparisions
        b = @benchmark CUDA.@sync SpeedyWeather.transform!($diagn.grid.vor_grid, $progn.vor, $spectral_transform)
        add_results!(suite, b, i, 16)

        # and inverse
        b = @benchmark CUDA.@sync SpeedyWeather.transform!($progn.vor, $diagn.grid.vor_grid, $spectral_transform)
        add_results!(suite, b, i, 17)

    end

    return suite
end
