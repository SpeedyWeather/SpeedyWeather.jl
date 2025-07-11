@testset "GPU Barotropic" begin
    spectral_grid = SpectralGrid(trunc=32, nlayers=1, architecture=SpeedyWeather.GPU())
    model = CUDA.@allowscalar BarotropicModel(spectral_grid=spectral_grid)
    simulation = CUDA.@allowscalar initialize!(model)
    run!(simulation, steps=4)

    @test simulation.model.feedback.nars_detected == false
end
