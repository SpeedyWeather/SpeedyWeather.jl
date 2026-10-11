@testset "GPU ShallowWater" begin
    spectral_grid = SpectralGrid(truncation = 33, nlayers = 1, architecture = SpeedyWeather.GPU())
    spectral_transform = MatrixSpectralTransform(spectral_grid)
    model = ShallowWaterModel(spectral_grid; spectral_transform)
    simulation = initialize!(model)
    run!(simulation, steps = 4)

    @test simulation.model.feedback.nans_detected == false
end
