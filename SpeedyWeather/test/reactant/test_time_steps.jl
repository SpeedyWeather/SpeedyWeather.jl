# Several time steps per compiled program: `time_steps!` runs n time steps in one traced loop and
# `time_stepping!` runs one compiled program for the time steps between two outputs (without
# output one program for the whole run)

using NCDatasets

# the same initial conditions in every simulation: the Barotropic model's are random, seed them
function seeded_reactant_model(ModelType::Type)
    ModelType == BarotropicModel || return create_reactant_model(ModelType)
    spectral_grid = SpectralGrid(architecture = ReactantDevice(), nlayers = 1, truncation = TRUNCATION)
    return create_reactant_model(ModelType; initial_conditions = (; vordiv = RandomVelocity(spectral_grid; seed = 1)))
end

@testset "time_steps! vs single compiled time steps: $ModelType" for ModelType in (BarotropicModel, PrimitiveWetModel)
    nsteps = 4
    simulation_single = initialize!(seeded_reactant_model(ModelType))
    simulation_loop = initialize!(seeded_reactant_model(ModelType))
    initialize!(simulation_single, steps = nsteps)
    initialize!(simulation_loop, steps = nsteps)

    vorticity(simulation) = Array(parent(simulation.variables.prognostic.vorticity))
    @test vorticity(simulation_loop) == vorticity(simulation_single)

    r_time_step! = @compile compile_options = SpeedyWeather.reactant_compile_options() SpeedyWeather.time_step!(simulation_single)
    for _ in 1:nsteps
        r_time_step!(simulation_single)
    end

    r_time_steps! = @compile compile_options = SpeedyWeather.reactant_compile_options() SpeedyWeather.time_steps!(simulation_loop, nsteps)
    r_time_steps!(simulation_loop, nsteps)

    # a differently fused program, so not bitwise identical
    for name in filter(in(keys(simulation_loop.variables.prognostic)), (:vorticity, :divergence, :temperature))
        single = Array(parent(getfield(simulation_single.variables.prognostic, name)))
        loop = Array(parent(getfield(simulation_loop.variables.prognostic, name)))
        @test any(!iszero, loop)
        @test loop ≈ single rtol = 1.0e-5
    end

    clock_single, clock_loop = simulation_single.variables.prognostic.clock, simulation_loop.variables.prognostic.clock
    @test DateTime(clock_loop.time) == DateTime(clock_single.time)
    @test Int(clock_loop.step_counter) == Int(clock_single.step_counter) == nsteps
    @test Int(clock_loop.time_step_counter) == Int(clock_single.time_step_counter)    # spin-up steps don't count
end

@testset "run! as one compiled program without output, with feedback" begin
    model = create_reactant_model(BarotropicModel; feedback = Feedback(verbose = false))
    simulation = initialize!(model)
    run!(simulation, steps = 6)
    (; clock) = simulation.variables.prognostic
    @test Int(clock.step_counter) == 6
    @test model.feedback.progress_meter.counter == 6
    @test !model.feedback.nans_detected
    @test all(isfinite, Array(parent(simulation.variables.prognostic.vorticity)))
end

@testset "run! with output: spin-up steps, time steps between outputs, remainder" begin
    # leapfrog with 1 spin-up step, output every 6 hours, 26 hours: 4 outputs and 2 hours remainder
    tmp = mktempdir(pwd(), prefix = "tmp_reactant_time_steps_")
    spectral_grid = SpectralGrid(architecture = ReactantDevice(), truncation = TRUNCATION, nlayers = 8)
    output = NetCDFOutput(spectral_grid, PrimitiveDryModel; path = tmp, interval = Hour(6))
    model = PrimitiveDryModel(
        spectral_grid; spectral_transform = MatrixSpectralTransform(spectral_grid), feedback = nothing,
        output, model_kwargs(PrimitiveDryModel, spectral_grid)...,
    )
    simulation = initialize!(model)
    run!(simulation; period = Hour(26), output = true)

    (; clock) = simulation.variables.prognostic
    n_time_steps = Int(clock.n_time_steps)
    @test n_time_steps % output.output_every_n_steps != 0     # a remainder after the last output
    @test Int(clock.time_step_counter) == n_time_steps
    @test DateTime(clock.time) - DateTime(clock.start) == n_time_steps * model.time_stepping.Δt_millisec
    @test output.output_counter == 5        # initial conditions + 4 outputs

    ds = NCDataset(joinpath(output.run_path, output.filename))
    time = ds["time"][:]
    @test time .- time[1] == Hour.(0:6:24)
    vor = ds["vor"][:, :, :, :]
    @test all(isfinite, vor)
    @test vor[:, :, :, 1] != vor[:, :, :, end]
    close(ds)
end
