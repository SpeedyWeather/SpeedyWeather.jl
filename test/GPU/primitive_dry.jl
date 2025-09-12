@testset "GPU PrimitiveDryModel" begin
    arch = SpeedyWeather.GPU()
    spectral_grid = SpectralGrid(trunc=32, nlayers=8, architecture=arch)
    model = PrimitiveDryModel(spectral_grid=spectral_grid, physics=false)
    simulation = CUDA.@allowscalar initialize!(model)
    run!(simulation, steps=4)

    @test simulation.model.feedback.nans_detected == false
end
