@testset "Elapsed time of Leapfrog steps" begin
    spectral_grid = SpectralGrid(truncation = 15, nlayers = 4)
    time_stepping = Leapfrog(spectral_grid)
    clock = Clock()
    (; Δt) = time_stepping

    # spin-up: Euler Δt/2, leapfrog Δt but rewinding the clock by Δt/2, then Δt
    elapsed = map(0:3) do step_counter
        clock.step_counter = step_counter
        SpeedyWeather.elapsed_time(time_stepping, clock)
    end
    @test elapsed == [Δt / 2, Δt / 2, Δt, Δt]
end

@testset "Accumulate diagnostic" begin
    spectral_grid = SpectralGrid(truncation = 15, nlayers = 4)
    model = PrimitiveWetModel(spectral_grid)
    model.feedback.verbose = false
    @test model.diagnostics.rain_large_scale isa Accumulate

    simulation = initialize!(model)
    run!(simulation, steps = 10)
    vars = simulation.variables

    # parameterizations only write rates, the diagnostics accumulate them
    @test !haskey(vars.parameterizations, :rain_large_scale)
    @test haskey(vars.parameterizations, :rain_rate_large_scale)
    @test all(vars.output.rain_large_scale .>= 0)
    @test any(vars.output.rain_large_scale .> 0)
    @test haskey(vars.output.land, :river_runoff)

    # one diagnostic step adds Δt * rate
    vars.parameterizations.rain_rate_large_scale .= 2
    vars.output.rain_large_scale .= 1
    SpeedyWeather.diagnostic!(vars, model.diagnostics.rain_large_scale, model)
    Δt = SpeedyWeather.elapsed_time(model.time_stepping, vars.prognostic.clock)
    @test all(vars.output.rain_large_scale .≈ 1 + 2 * Δt)

    # reset to zero
    SpeedyWeather.reset_diagnostic!(vars, model.diagnostics.rain_large_scale, model)
    @test all(vars.output.rain_large_scale .== 0)

    # the rate is not declared by the diagnostic: warn and skip without convection
    model = PrimitiveWetModel(spectral_grid, convection = nothing)
    model.feedback.verbose = false
    simulation = @test_logs (:warn, r"rain_rate_convection not defined") match_mode = :any initialize!(model)
    @test !haskey(simulation.variables.parameterizations, :rain_rate_convection)
    run!(simulation, steps = 2)
    @test all(simulation.variables.output.rain_convection .== 0)
end

@testset "Output writer resets diagnostics" begin
    spectral_grid = SpectralGrid(truncation = 15, nlayers = 4)
    diagnostics = (; PrecipitationDiagnostics()..., temperature_max = Maximum(), temperature_mean = Mean())

    # with output every 6 hours: precipitation since last output, reset after every output
    output = ArrayOutput(interval = Hour(6), groups = (:output,))
    model = PrimitiveWetModel(spectral_grid; output, diagnostics)
    model.feedback.verbose = false
    simulation = initialize!(model)
    run!(simulation, period = Day(1), output = true)

    snapshots = model.output.output
    @test length(snapshots) == 5                        # initial conditions + 4 outputs
    vars = simulation.variables.output
    @test all(vars.rain_large_scale .== 0)              # reset after the last output
    @test all(==(-Inf), vars.surface_air_temperature_maximum)
    @test all(isnan, vars.surface_air_temperature_mean)
    @test vars.surface_air_temperature_mean_count[] == 0

    last_output = snapshots[end].output
    @test all(isfinite, last_output.surface_air_temperature_maximum)
    @test all(last_output.surface_air_temperature_maximum .>= last_output.surface_air_temperature_mean)
    @test last_output.surface_air_temperature_mean_count[] == 6 * 60 * 60 ÷ Int(model.time_stepping.Δt)

    # summing over all output intervals yields the precipitation of a run without output (no reset)
    rain_since_last_output = sum(snapshot.output.rain_large_scale for snapshot in snapshots[2:end])
    model_no_output = PrimitiveWetModel(spectral_grid; output = ArrayOutput(interval = Hour(6)), diagnostics)
    model_no_output.feedback.verbose = false
    simulation_no_output = initialize!(model_no_output)
    run!(simulation_no_output, period = Day(1), output = false)
    rain_total = simulation_no_output.variables.output.rain_large_scale
    @test rain_since_last_output ≈ rain_total rtol = 1.0e-4
end

@testset "Maximum and Mean diagnostics" begin
    spectral_grid = SpectralGrid(truncation = 15, nlayers = 4)
    diagnostics = (; temperature_max = Maximum(), temperature_mean = Mean())
    model = PrimitiveWetModel(spectral_grid; diagnostics)
    model.feedback.verbose = false
    simulation = initialize!(model)
    vars = simulation.variables
    (; temperature_max, temperature_mean) = model.diagnostics

    @test temperature_max.name == :surface_air_temperature_maximum
    @test temperature_mean.name == :surface_air_temperature_mean
    @test all(==(-Inf), vars.output.surface_air_temperature_maximum)
    @test all(isnan, vars.output.surface_air_temperature_mean)

    surface_air_temperature = vars.parameterizations.surface_air_temperature
    temperature1 = 280 .+ 10 .* rand(spectral_grid.npoints)
    temperature2 = 280 .+ 10 .* rand(spectral_grid.npoints)
    for temperature in (temperature1, temperature2)
        surface_air_temperature.data .= temperature
        SpeedyWeather.diagnostics!(vars, model)
    end
    @test vars.output.surface_air_temperature_maximum.data ≈ max.(temperature1, temperature2)
    @test vars.output.surface_air_temperature_mean.data ≈ (temperature1 .+ temperature2) ./ 2
    @test vars.output.surface_air_temperature_mean_count[] == 2

    # after a reset the first value is taken
    SpeedyWeather.reset_diagnostics!(vars, model)
    @test all(==(-Inf), vars.output.surface_air_temperature_maximum)
    surface_air_temperature.data .= temperature1
    SpeedyWeather.diagnostics!(vars, model)
    @test vars.output.surface_air_temperature_maximum.data ≈ temperature1
    @test vars.output.surface_air_temperature_mean.data ≈ temperature1

    # NaNs propagate, e.g. to see where the model blew up
    surface_air_temperature[1] = NaN
    SpeedyWeather.diagnostics!(vars, model)
    @test isnan(vars.output.surface_air_temperature_maximum[1])
end

@testset "Output variables in restart files" begin
    spectral_grid = SpectralGrid(truncation = 15, nlayers = 4)
    model = PrimitiveWetModel(spectral_grid)
    model.feedback.verbose = false
    simulation = initialize!(model)
    run!(simulation, steps = 4)
    rain = simulation.variables.output.rain_large_scale
    rain .= rand(spectral_grid.npoints)

    path = mktempdir()
    restart = WriteVariablesRestartFile(path = joinpath(path, "run_0001"), write_only_with_output = false)
    SpeedyWeather.finalize!(restart, simulation.variables, model)

    model2 = PrimitiveWetModel(spectral_grid, initial_conditions = StartFromFile(; path))
    simulation2 = initialize!(model2)
    @test simulation2.variables.output.rain_large_scale == rain
end
