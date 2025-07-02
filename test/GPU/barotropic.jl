@testset "GPU Barotropic" begin
    spectral_grid = SpectralGrid(trunc=32, nlayers=1, architecture=SpeedyWeather.GPU())
    model = BarotropicModel(spectral_grid=spectral_grid)
    CUDA.@allowscalar simulation = initialize!(model)
    run!(simulation, steps=4)

    @test simulation.model.feedback.nars_detected == false
end
