@testset "run PrimitiveWetModel with MatrixSpectralTransform, no blow up" begin
    spectral_grid = SpectralGrid()
    spectral_transform = MatrixSpectralTransform(spectral_grid)
    model = PrimitiveWetModel(spectral_grid; spectral_transform)
    simulation = initialize!(model)
    run!(simulation, period = Day(20))
    @test !SpeedyWeather.nans_detected(simulation.model.feedback)
end
