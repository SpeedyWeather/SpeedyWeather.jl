@testset "GPU PrimitiveDryModel" begin
    spectral_grid = SpectralGrid(trunc=32, nlayers=8, architecture=SpeedyWeather.GPU())
    model = PrimitiveDryModel(spectral_grid=spectral_grid, physics=false)
    simulation = CUDA.@allowscalar initialize!(model)
    run!(simulation, steps=4)

    @test simulation.model.feedback.nans_detected == false
end
