@testset "run PrimitiveWetModel defaults for a year, no blow up" begin
    spectral_grid = SpectralGrid()
    model = PrimitiveWetModel(spectral_grid)
    model.feedback.verbose = false
    simulation = initialize!(model)
    run!(simulation, period = Year(1))
    @test !SpeedyWeather.nans_detected(simulation.model.feedback)
end
