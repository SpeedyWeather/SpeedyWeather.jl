@testset "GPU PrimitiveWetModel (with SpectralTransform)" begin
    arch = SpeedyWeather.GPU()
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_gpu_netcdf_")

    # includes particles to test GPU particle advection and output on GPU
    spectral_grid = SpectralGrid(truncation = 33, nlayers = 8, architecture = arch)
    spectral_transform = SpectralTransform(spectral_grid)
    particle_advection = ParticleAdvection2D(spectral_grid, nparticles = 10, layer = 1)
    random_process = SpectralAR1Process(spectral_grid)
    sppt = StochasticallyPerturbedParameterizationTendencies(spectral_grid)
    output = NetCDFOutput(spectral_grid, PrimitiveWet, path = tmp_output_path, id = "gpu-netcdf")
    model = PrimitiveWetModel(spectral_grid; spectral_transform, output, particle_advection, random_process, stochastic_physics = sppt)
    simulation = initialize!(model)
    run!(simulation, steps = 3, output = true)

    @test simulation.model.feedback.nans_detected == false
    @test isfile(joinpath(output.run_path, output.filename))
end
