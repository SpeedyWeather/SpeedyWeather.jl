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
    @test all(vars.output.rain_large_scale .≈ 1 + 2 * model.time_stepping.Δt)

    # the rate is declared by the diagnostic too, so it also works without convection
    model = PrimitiveWetModel(spectral_grid, convection = nothing)
    model.feedback.verbose = false
    simulation = initialize!(model)
    run!(simulation, steps = 2)
    @test all(simulation.variables.output.rain_convection .== 0)
end

@testset "PeriodicMaximum diagnostic" begin
    spectral_grid = SpectralGrid(truncation = 15, nlayers = 4)
    diagnostics = (; PrecipitationDiagnostics()..., temperature_max = PeriodicMaximum())
    model = PrimitiveWetModel(spectral_grid; diagnostics)
    model.feedback.verbose = false
    simulation = initialize!(model)
    vars = simulation.variables
    maximum_diagnostic = model.diagnostics.temperature_max

    @test maximum_diagnostic.name == :surface_air_temperature_maximum
    @test all(isnan, vars.output.surface_air_temperature_maximum)
    @test all(==(-Inf), vars.output.surface_air_temperature_maximum_running)

    # within the period: running maximum of all values so far
    clock = vars.prognostic.clock
    clock.time_step_counter = 1
    surface_air_temperature = vars.parameterizations.surface_air_temperature
    temperature1 = 280 .+ 10 .* rand(spectral_grid.npoints)
    temperature2 = 280 .+ 10 .* rand(spectral_grid.npoints)
    surface_air_temperature.data .= temperature1
    SpeedyWeather.diagnostic!(vars, maximum_diagnostic, model)
    surface_air_temperature.data .= temperature2
    SpeedyWeather.diagnostic!(vars, maximum_diagnostic, model)
    @test vars.output.surface_air_temperature_maximum_running.data ≈ max.(temperature1, temperature2)
    @test all(isnan, vars.output.surface_air_temperature_maximum)

    # end of the period: pass on the maximum including this step's value, reset the running maximum
    clock.time_step_counter = maximum_diagnostic.schedule.every_n_steps
    surface_air_temperature.data .= 300
    SpeedyWeather.diagnostic!(vars, maximum_diagnostic, model)
    @test all(vars.output.surface_air_temperature_maximum .== 300)
    @test all(==(-Inf), vars.output.surface_air_temperature_maximum_running)

    # other periods via the schedule
    six_hourly = PeriodicMaximum(name = :temperature_max_6h, schedule = PeriodicSchedule(every = Hour(6)))
    @test six_hourly.schedule.every == Second(Hour(6))

    # and within a simulation
    simulation = initialize!(model)
    run!(simulation, period = Hour(30))
    @test all(isfinite, simulation.variables.output.surface_air_temperature_maximum)
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
