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
