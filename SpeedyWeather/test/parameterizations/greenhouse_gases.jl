@testset "Greenhouse gases: allocation" begin
    spectral_grid = SpectralGrid(trunc = 31, nlayers = 8)

    @testset for CO2type in (CO2, ExponentialCO2)
        co2 = CO2type(spectral_grid)
        model = PrimitiveWetModel(spectral_grid; greenhouse_gases = (; co2))
        simulation = initialize!(model)
        @test haskey(simulation.variables.prognostic.greenhouse_gases, :co2)
        @test simulation.variables.prognostic.greenhouse_gases.co2 isa Ref
    end

    @testset "No greenhouse gases" begin
        model = PrimitiveWetModel(spectral_grid; greenhouse_gases = nothing)
        simulation = initialize!(model)
        @test !haskey(simulation.variables.prognostic, :greenhouse_gases)
    end
end

@testset "Greenhouse gases: time stepping" begin
    spectral_grid = SpectralGrid(trunc = 31, nlayers = 8)

    # ExponentialCO2 increases over time, so concentration after a few days > initial
    co2 = ExponentialCO2(spectral_grid; start = DateTime(1850))
    model = PrimitiveWetModel(spectral_grid; greenhouse_gases = (; co2))
    simulation = initialize!(model)

    co2_initial = simulation.variables.prognostic.greenhouse_gases.co2[]
    run!(simulation, period = Day(5))
    co2_after = simulation.variables.prognostic.greenhouse_gases.co2[]

    @test co2_after > co2_initial
end

@testset "Greenhouse gases: CO2 scenarios" begin
    t_before = DateTime(1999)
    t_after  = DateTime(2001)

    @test TwoTimesCO2()(t_before) == SpeedyWeather.DEFAULT_CO2
    @test TwoTimesCO2()(t_after)  == 2 * SpeedyWeather.DEFAULT_CO2

    @test FourTimesCO2()(t_before) == SpeedyWeather.DEFAULT_CO2
    @test FourTimesCO2()(t_after)  == 4 * SpeedyWeather.DEFAULT_CO2

    # constant CO2
    @test CO2(420)(t_before) == 420
    @test CO2(420)(t_after)  == 420
end
