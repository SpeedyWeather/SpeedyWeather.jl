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

@testset "GPU surface longwave over sea ice matches CPU" begin
    # the ocean surface emission over sea ice dispatches on the sea ice model, which therefore is a core
    # component available inside the GPU kernels: GPU and CPU give the same emission (Float32 round-off)
    function set_points!(field, values)     # values per grid point, the same on every layer and step
        data = Array(field.data)
        data .= reshape(values, :, ntuple(_ -> 1, ndims(data) - 1)...)
        return copyto!(field.data, data)
    end

    function ocean_longwave(arch, Longwave, SeaIce)
        spectral_grid = SpectralGrid(truncation = 32, nlayers = 8, architecture = arch)
        radiation = Radiation(spectral_grid; longwave = Longwave(spectral_grid))
        sea_ice = SeaIce === ThermodynamicSeaIce ? ThermodynamicSeaIce(spectral_grid) : SeaIce()
        model = PrimitiveWetModel(spectral_grid; radiation, sea_ice, parameterizations = (:radiation,))
        initialize!(model.radiation, model)
        vars = Variables(model)

        npoints = spectral_grid.npoints
        x = (0:(npoints - 1)) ./ npoints
        temperature = Array(vars.grid.temperature.data)
        for k in 1:spectral_grid.nlayers
            temperature[:, k, :] .= 210 + 10 * (k - 1)
        end
        copyto!(vars.grid.temperature.data, temperature)
        vars.grid.humidity .= 0.003
        vars.parameterizations.surface_pressure .= 1.0e5
        vars.tendencies.grid.temperature .= 0
        vars.prognostic.land.soil_temperature .= 280
        set_points!(vars.prognostic.ocean.sea_surface_temperature, Float32.(265 .+ 30 .* x))
        set_points!(vars.parameterizations.surface_air_temperature, Float32.(230 .+ 50 .* reverse(x)))
        ice = Float32.(ifelse.(isodd.(1:npoints), x, 0))        # sea ice on every other point
        set_points!(vars.prognostic.ocean.sea_ice_concentration, ice)

        SpeedyWeather.column_parameterizations!(vars, model)
        return Array(vars.parameterizations.ocean.surface_longwave_up), ice
    end

    @testset for Longwave in (OneBandLongwave, JeevanjeeRadiation), SeaIce in (PrescribedSeaIce, ThermodynamicSeaIce)
        cpu, ice = ocean_longwave(SpeedyWeather.CPU(), Longwave, SeaIce)
        gpu, _ = ocean_longwave(SpeedyWeather.GPU(), Longwave, SeaIce)
        @test all(isfinite, gpu)
        @test gpu ≈ cpu rtol = 1.0e-5
    end
end
