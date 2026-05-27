@testset "GPU ShallowWater" begin
    spectral_grid = SpectralGrid(trunc = 32, nlayers = 1, architecture = SpeedyWeather.GPU())
    model = ShallowWaterModel(spectral_grid = spectral_grid)
    simulation = initialize!(model)
    run!(simulation, steps = 4)

    @test simulation.feedback.nans_detected == false
end