@testset "Simulation(model) equals initialize!(model)" begin
    spectral_grid = SpectralGrid(truncation = 6, nlayers = 1)

    # Simulation(model) is an alias for initialize!(model)
    model = BarotropicModel(spectral_grid)
    simulation = Simulation(model)
    @test simulation isa SpeedyWeather.Simulation
    @test simulation.model === model

    # keyword arguments are forwarded, e.g. the initial time
    time = DateTime(2020, 5, 1)
    simulation = Simulation(BarotropicModel(spectral_grid), time = time)
    @test simulation.variables.prognostic.clock.time == time

    # both entry points give the same result
    model1 = ShallowWaterModel(spectral_grid)
    model2 = ShallowWaterModel(spectral_grid)
    simulation1 = initialize!(model1, time = time)
    simulation2 = Simulation(model2, time = time)
    @test simulation1.variables.prognostic.vorticity == simulation2.variables.prognostic.vorticity
    @test simulation1.variables.prognostic.clock.time == simulation2.variables.prognostic.clock.time
end
