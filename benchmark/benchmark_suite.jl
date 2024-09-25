Base.@kwdef mutable struct BenchmarkSuite
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