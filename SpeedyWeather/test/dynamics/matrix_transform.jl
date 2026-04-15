@testset "run PrimitiveWetModel with MatrixSpectralTransform, no blow up" begin
    spectral_grid = SpectralGrid()
    spectral_transform = MatrixSpectralTransform(spectral_grid)
    model = PrimitiveWetModel(spectral_grid; spectral_transform)
    simulation = initialize!(model)
    run!(simulation, period = Day(20))
    @test simulation.model.feedback.nans_detected == false
end