@testset "GPU PrimitiveWetModel" begin
    arch = SpeedyWeather.GPU()

    # includes particles to test GPU particle advection
    spectral_grid = SpectralGrid(trunc = 32, nlayers = 8, architecture = arch)
    particle_advection = ParticleAdvection2D(spectral_grid, nparticles = 10, layer = 1)
    model = PrimitiveWetModel(spectral_grid; particle_advection)
    simulation = initialize!(model)
    run!(simulation, steps = 3)

    @test simulation.model.feedback.nans_detected == false
end
