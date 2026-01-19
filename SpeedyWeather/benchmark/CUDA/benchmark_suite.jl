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
        device = suite.model[i]

        spectral_grid = SpectralGrid(; NF, trunc, Grid, nlayers, device)
        suite.nlat[i] = spectral_grid.nlat
        S = SpectralTransform(spectral_grid)

        specs, grids = generate_random_inputs(spectral_grid)

        println("Running benchmark for device=$device, trunc=$trunc, nlayers=$nlayers, NF=$NF \n")

        # Forward _legendre
        b = @benchmark CUDA.@sync SpeedyTransforms._legendre!($specs, $S.scratch_memory_north, $S.scratch_memory_south, $S)
        add_results!(suite, b, i, 1)

        # Inverse _legendre
        b = @benchmark CUDA.@sync SpeedyTransforms._legendre!($S.scratch_memory_north, $S.scratch_memory_south, $specs, $S)
        add_results!(suite, b, i, 2)

        # Forward _fourier!(..., grids, S)
        b = @benchmark CUDA.@sync SpeedyTransforms._fourier!($S.scratch_memory_north, $S.scratch_memory_south, $grids, $S)
        add_results!(suite, b, i, 3)

        # Reverse _fourier!(grids, ..., S)
        b = @benchmark CUDA.@sync SpeedyTransforms._fourier!($grids, $S.scratch_memory_north, $S.scratch_memory_south, $S)
        add_results!(suite, b, i, 4)
    end

    return suite
end 

@kwdef mutable struct BenchmarkSuiteModel <: AbstractBenchmarkSuiteTimed
    title::String 
    nruns::Int = 1
    architecture::Vector = fill(SpeedyWeather.CPU, nruns)
    NF::Vector = fill(SpeedyWeather.DEFAULT_NF, nruns)
    trunc::Vector{Int} = fill(SpeedyWeather.DEFAULT_TRUNC, nruns)
    nlayers::Vector{Int} = fill(SpeedyWeather.DEFAULT_NLAYERS, nruns)
    Grid::Vector = fill(SpeedyWeather.DEFAULT_GRID, nruns)
    nlat::Vector{Int} = fill(0, nruns)
    function_names::Vector{String} = ["parameterization_tendencies", "dynamics_tendencies", "implicit_correction"]
    time::Vector{Vector{Float64}} = [fill(0.0, length(function_names)) for i=1:nruns]
    memory::Vector{Vector{Int}} = [fill(0, length(function_names)) for i=1:nruns]
    allocs::Vector{Vector{Int}} = [fill(0, length(function_names)) for i=1:nruns]
end 

function run_benchmark_suite!(suite::BenchmarkSuiteModel)

    for i in 1:suite.nruns
        NF = suite.NF[i]
        trunc = suite.trunc[i]
        nlayers = suite.nlayers[i]
        Grid = suite.Grid[i]
        architecture = suite.architecture[i]

        spectral_grid = SpectralGrid(;NF, trunc, Grid, nlayers, architecture)
        suite.nlat[i] = spectral_grid.nlat
        model = PrimitiveWetModel(spectral_grid)
        simulation = initialize!(model)
        # spin up 
        run!(simulation; steps = 10)
        initialize!(simulation)
        progn, diagn, model = SpeedyWeather.unpack(simulation)
        fill!(diagn.tendencies, 0, typeof(model))

        println("Running benchmark for device=$architecture, trunc=$trunc, nlayers=$nlayers, NF=$NF \n")

        b = @benchmark CUDA.@sync SpeedyWeather.parameterization_tendencies!($diagn, $progn, $model)
        add_results!(suite, b, i, 1)

        b = @benchmark CUDA.@sync SpeedyWeather.dynamics_tendencies!($diagn, $progn, $model)
        add_results!(suite, b, i, 2)

        b = @benchmark CUDA.@sync SpeedyWeather.implicit_correction!($diagn, $progn, $model)
        add_results!(suite, b, i, 3)
    end

    return suite
end 

