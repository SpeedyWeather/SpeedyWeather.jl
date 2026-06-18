@testset "GPU PrimitiveWetModel" begin
    arch = SpeedyWeather.GPU()
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_gpu_netcdf_")

    # includes particles to test GPU particle advection and output on GPU
    # LEAPFROG 
    spectral_grid = SpectralGrid(trunc = 32, nlayers = 8, architecture = arch)
    particle_advection = ParticleAdvection2D(spectral_grid, nparticles = 10, layer = 1)
    time_stepping = Leapfrog(spectral_grid) 
    output = NetCDFOutput(spectral_grid, PrimitiveWet, path = tmp_output_path, id = "gpu-netcdf")
    model = PrimitiveWetModel(spectral_grid; output, particle_advection, time_stepping)
    simulation = initialize!(model)
    run!(simulation, steps = 3, output = true)

    @test simulation.model.feedback.nans_detected == false
    @test isfile(joinpath(output.run_path, output.filename))

    # NCYCLE LORENZ 
    # includes particles to test GPU particle advection and output on GPU
    spectral_grid = SpectralGrid(trunc = 32, nlayers = 8, architecture = arch)
    particle_advection = ParticleAdvection2D(spectral_grid, nparticles = 10, layer = 1)
    time_stepping = NCycleLorenz(spectral_grid; steps=4, variant=NCycleLorenzAB())
    output = NetCDFOutput(spectral_grid, PrimitiveWet, path = tmp_output_path, id = "gpu-netcdf")
    model = PrimitiveWetModel(spectral_grid; output, particle_advection, time_stepping)
    simulation = initialize!(model)
    run!(simulation, steps = 3, output = true)

    @test simulation.model.feedback.nans_detected == false
    @test isfile(joinpath(output.run_path, output.filename))
end
