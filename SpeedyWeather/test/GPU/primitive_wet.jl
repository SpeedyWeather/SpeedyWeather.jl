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
    # component available inside the GPU kernels: GPU and CPU give the same emission (Float32 round-off),
    # with ThermodynamicSeaIce also the same ice surface temperature
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
        vars.grid.humidity .= 0.003f0             # Float32 literals: no Float64 on Metal
        vars.parameterizations.surface_pressure .= 1.0f5
        vars.tendencies.grid.temperature .= 0.0f0
        vars.prognostic.land.soil_temperature .= 280.0f0
        set_points!(vars.prognostic.ocean.sea_surface_temperature, Float32.(265 .+ 30 .* x))
        set_points!(vars.parameterizations.surface_air_temperature, Float32.(230 .+ 50 .* reverse(x)))
        ice = Float32.(ifelse.(isodd.(1:npoints), x, 0))        # sea ice on every other point
        set_points!(vars.prognostic.ocean.sea_ice_concentration, ice)
        if sea_ice isa ThermodynamicSeaIce      # its prognostic variable is the volume, ice surface balance needs albedo, shortwave
            set_points!(vars.prognostic.ocean.sea_ice_volume, SpeedyWeather.sea_ice_volume.(ice, Ref(sea_ice)))
            vars.parameterizations.ocean.sea_ice_albedo .= 0.6f0
            vars.parameterizations.surface_shortwave_down .= 100.0f0
        end

        SpeedyWeather.column_parameterizations!(vars, model)
        ice_surface_temperature = haskey(vars.parameterizations.ocean, :sea_ice_surface_temperature) ?
            Array(vars.parameterizations.ocean.sea_ice_surface_temperature) : nothing
        return Array(vars.parameterizations.ocean.surface_longwave_up), ice, ice_surface_temperature
    end

    @testset for Longwave in (OneBandLongwave, JeevanjeeRadiation), SeaIce in (PrescribedSeaIce, ThermodynamicSeaIce)
        cpu, ice, T_i_cpu = ocean_longwave(SpeedyWeather.CPU(), Longwave, SeaIce)
        gpu, _, T_i_gpu = ocean_longwave(SpeedyWeather.GPU(), Longwave, SeaIce)
        @test all(isfinite, gpu)
        @test gpu ≈ cpu rtol = 1.0e-5
        if SeaIce === ThermodynamicSeaIce       # ice surface temperature from the ice surface energy balance (Newton)
            @test T_i_gpu[ice .> 0] ≈ T_i_cpu[ice .> 0] rtol = 1.0e-5
        end
    end
end
