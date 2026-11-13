@testset "OutputVariables module" begin
    @test OutputVariables isa Module
    @test OutputVariables.Temperature() isa SpeedyWeather.AbstractOutputVariable
    @test OutputVariables.Temperature(keepbits = 9).keepbits == 9
    @test OutputVariables.Tracer(:co2).name == "co2"
    @test all(v -> v isa SpeedyWeather.AbstractOutputVariable, OutputVariables.All())
    @test OutputVariables.Dynamics() isa Tuple

    # old names forward to the new ones with a deprecation warning
    for (old, new) in SpeedyWeather.DEPRECATED_OUTPUT_VARIABLES
        @test isdefined(OutputVariables, new)
    end
    variable = @test_logs (:warn, r"deprecated") match_mode = :any SpeedyWeather.TemperatureOutput(keepbits = 3)
    @test variable isa OutputVariables.Temperature
    @test variable.keepbits == 3
    dynamics = @test_logs (:warn, r"deprecated") match_mode = :any SpeedyWeather.DynamicsOutput()
    @test typeof.(dynamics) == typeof.(OutputVariables.Dynamics())
end

@testset "Precipitable water" begin
    spectral_grid = SpectralGrid(truncation = 21, nlayers = 8)
    model = PrimitiveWetModel(spectral_grid)
    simulation = initialize!(model)
    vars = simulation.variables
    @test any(v -> v isa OutputVariables.PrecipitableWater, OutputVariables.Dynamics())

    # constant specific humidity: precipitable water is q * pₛ / g as the layers' Δp sum to pₛ
    q = 0.01
    vars.grid.humidity .= q
    vars.parameterizations.surface_pressure .= 1.0e5
    vars.parameterizations.surface_pressure[1] = 0.8e5
    pw = vars.scratch.grid.a_2D
    g = model.planet.gravity
    humidity = SpeedyWeather.get_prognostic_step(vars.grid.humidity, model.time_stepping, model.output)
    SpeedyWeather.launch!(
        SpeedyWeather.architecture(pw), SpeedyWeather.LinearWorkOrder, size(pw),
        OutputVariables.precipitable_water_kernel!, pw, humidity,
        vars.parameterizations.surface_pressure, model.geometry.vertical_coordinates, inv(g),
    )
    @test pw[1] ≈ q * 0.8e5 / g rtol = 1.0e-5
    @test all(pw[2:end] .≈ q * 1.0e5 / g)
end
