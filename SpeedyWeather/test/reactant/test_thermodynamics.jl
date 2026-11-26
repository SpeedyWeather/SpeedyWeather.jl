# Saturation vapor pressure formulations with traced temperatures: the formulation determines the number
# format via `convert_number_format`, as `convert(Float32, ::TracedRNumber)` has no method

@testset "Saturation vapor pressure with traced temperatures (Reactant)" begin
    for SVP in (ClausiusClapeyronConstantL, ClausiusClapeyronLinearL, IFSMixedPhase, MagnusBolton)
        atmosphere = EarthAtmosphere(Float32, saturation = SVP(Float32))
        temperatures = Float32[230, 255, 273.15, 300]
        e = SpeedyWeather.saturation_vapor_pressure.(temperatures, atmosphere)
        de = SpeedyWeather.grad_saturation_vapor_pressure.(temperatures, atmosphere)

        # broadcast over a traced array, each element a TracedRNumber{Float32}
        traced_temperatures = Reactant.to_rarray(temperatures)
        e_traced = @jit SpeedyWeather.saturation_vapor_pressure.(traced_temperatures, atmosphere)
        de_traced = @jit SpeedyWeather.grad_saturation_vapor_pressure.(traced_temperatures, atmosphere)
        @test Array(e_traced) ≈ e rtol = 1.0e-5
        @test Array(de_traced) ≈ de rtol = 1.0e-5

        # a traced Float64 temperature is converted to the formulation's Float32
        e64 = @jit SpeedyWeather.saturation_vapor_pressure(Reactant.ConcreteRNumber(255.0), atmosphere)
        @test e64 isa Reactant.ConcreteRNumber{Float32}
        @test Float32(e64) ≈ e[2] rtol = 1.0e-5
    end
end
