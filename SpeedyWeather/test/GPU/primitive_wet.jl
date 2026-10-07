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

@testset "GPU PrimitiveWetModel (default construction, WhichTransform)" begin
    # No component overrides: at the default truncation this lets `WhichTransform`
    # pick `MatrixSpectralTransform` (as it does for any truncation <= 64 on GPU),
    # unlike the testset above which forces `SpectralTransform` explicitly. This is
    # the path a plain `PrimitiveWetModel(spectral_grid)` on GPU actually takes.
    arch = SpeedyWeather.GPU()
    spectral_grid = SpectralGrid(architecture = arch)
    model = PrimitiveWetModel(spectral_grid)
    @test model.spectral_transform isa SpeedyWeather.SpeedyTransforms.MatrixSpectralTransform
    simulation = Simulation(model)
    run!(simulation, steps = 3)

    @test simulation.model.feedback.nans_detected == false
end

@testset "GPU PrimitiveWetModel with prognostic clouds" begin
    # the fused cloud condensate widens the batched transforms (5L+1 prognostic, 12L+1 tendency
    # batch), the default matrix transform at this resolution has to fit them; the column scheme
    # and cloud radiation run in the fused parameterization kernel
    arch = SpeedyWeather.GPU()
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8, architecture = arch)
    large_scale_condensation = PrognosticCloudCondensation(spectral_grid)
    shortwave = OneBandShortwave(spectral_grid; clouds = PrognosticClouds(spectral_grid))
    longwave = OneBandLongwave(spectral_grid; transmissivity = CloudyLongwaveTransmissivity(spectral_grid))
    radiation = Radiation(spectral_grid; shortwave, longwave)
    stochastic_physics = StochasticallyPerturbedParameterizationTendencies(spectral_grid)
    random_process = SpectralAR1Process(spectral_grid)
    model = PrimitiveWetModel(spectral_grid; large_scale_condensation, radiation, stochastic_physics, random_process)
    simulation = initialize!(model)
    run!(simulation, period = Day(1))

    @test simulation.model.feedback.nans_detected == false
    condensate = Array(simulation.variables.grid.cloud_condensate.data)
    cloud_fraction = Array(simulation.variables.parameterizations.cloud_fraction.data)
    @test any(condensate .> 0)
    @test all(0 .<= cloud_fraction .<= 1)
end
