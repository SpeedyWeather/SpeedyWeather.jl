@testset "GPU PrimitiveDryModel" begin
    arch = SpeedyWeather.GPU()
    spectral_grid = SpectralGrid(trunc=32, nlayers=8, architecture=arch)
    model = PrimitiveWetModel(spectral_grid=spectral_grid)
    simulation = initialize!(model)
    run!(simulation, steps=3)

    @test simulation.model.feedback.nans_detected == false
end
