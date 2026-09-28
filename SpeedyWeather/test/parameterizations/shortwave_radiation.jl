function init_shortwave_state!(vars, model)
    vars.grid.temperature .= 280
    vars.grid.humidity .= 1.0e-3
    vars.grid.pressure .= 1.0e5
    vars.parameterizations.cloud_top .= model.spectral_grid.nlayers + 1
    vars.parameterizations.rain_rate .= 0
    vars.parameterizations.ocean.albedo .= 0.5
    vars.parameterizations.land.albedo .= 0.3
    vars.parameterizations.surface_pressure .= 1.0e5
    return nothing
end

@testset "Shortwave radiation" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
    @testset for SW in (Nothing, TransparentShortwave, OneBandShortwave, OneBandGreyShortwave)
        sw = SW(spectral_grid)
        model = PrimitiveWetModel(spectral_grid; shortwave_radiation = sw)

        initialize!(model.shortwave_radiation, model)

        vars = Variables(model)
        init_shortwave_state!(vars, model)

        SpeedyWeather.parameterization!(vars, model.solar_zenith, model)

        for ij in 1:model.spectral_grid.npoints
            SpeedyWeather.parameterization!(ij, vars, model.shortwave_radiation, model)
        end

        # top of atmosphere radiation down
        trd = model.planet.solar_constant * vars.parameterizations.cos_zenith

        if !(sw isa Nothing)
            osr = vars.parameterizations.outgoing_shortwave
            ssrd = vars.parameterizations.surface_shortwave_down
            @test all(0 .<= osr .<= ssrd .<= trd)
            @test all(isfinite.(osr))
            @test all(isfinite.(ssrd))
        end

        if !(sw isa Nothing)
            @test any(vars.parameterizations.outgoing_shortwave .> 0)
            @test any(vars.parameterizations.surface_shortwave_down .> 0)
        else
            @test !haskey(vars.parameterizations, :outgoing_shortwave)
            @test all(vars.parameterizations.surface_shortwave_down .== 0)
        end
    end
end

@testset "Shortwave radiation transmissivity" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
    @testset for T in (TransparentShortwaveTransmissivity, BackgroundShortwaveTransmissivity)
        transmissivity = T(spectral_grid)
        sw = OneBandShortwave(spectral_grid; transmissivity = transmissivity)
        model = PrimitiveWetModel(spectral_grid; shortwave_radiation = sw)

        initialize!(model.shortwave_radiation, model)

        vars = Variables(model)
        init_shortwave_state!(vars, model)
        vars.parameterizations.cos_zenith .= 1

        clouds = SpeedyWeather.clouds!(1, vars, model.shortwave_radiation.clouds, model)
        t = SpeedyWeather.transmissivity!(1, vars, clouds, model.shortwave_radiation.transmissivity, model)
        for ij in 1:model.spectral_grid.npoints
            clouds = SpeedyWeather.clouds!(ij, vars, model.shortwave_radiation.clouds, model)
            t = SpeedyWeather.transmissivity!(ij, vars, clouds, model.shortwave_radiation.transmissivity, model)
        end
        @test all(0 .< t .<= 1)
    end
end

@testset "Shortwave radiation clouds" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
    @testset for C in (DiagnosticClouds, NoClouds)
        clouds = C(spectral_grid)
        sw = OneBandShortwave(spectral_grid; clouds = clouds)
        model = PrimitiveWetModel(spectral_grid; shortwave_radiation = sw)

        initialize!(model.shortwave_radiation, model)

        vars = Variables(model)
        init_shortwave_state!(vars, model)
        vars.parameterizations.cos_zenith .= 1

        ij = rand(1:model.spectral_grid.npoints)
        clouds_state = SpeedyWeather.clouds!(ij, vars, model.shortwave_radiation.clouds, model)
        @test 0 <= clouds_state.cloud_cover <= 1
        @test 1 <= clouds_state.cloud_top <= model.spectral_grid.nlayers + 1
        SpeedyWeather.parameterization!(ij, vars, model.shortwave_radiation, model)
    end
end

@testset "Diagnostic cloud top at level of maximum relative humidity" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
    model = PrimitiveWetModel(spectral_grid; shortwave_radiation = OneBandShortwave(spectral_grid))
    initialize!(model.shortwave_radiation, model)
    vars = Variables(model)
    init_shortwave_state!(vars, model)

    # relative humidity per layer, all above rh_min but maximum in layer 5
    ij = 1
    for k in 1:8    # nonzero geopotential for cloud top height, increasing upwards
        vars.dynamics.geopotential[ij, k] = 1000 * (9 - k)
    end
    relative_humidity = [0.9, 0.5, 0.6, 0.7, 0.95, 0.8, 0.5, 0.99]
    for k in 1:8
        p = model.geometry.σ_levels_full[k] * 1.0e5
        qsat = SpeedyWeather.saturation_humidity(280.0f0, p, model.atmosphere)
        vars.grid.humidity[ij, k, :] .= relative_humidity[k] * qsat
    end

    # layer 1 (top) and 8 (surface) excluded, so maximum in layer 5 not the highest layer above rh_min
    clouds_state = SpeedyWeather.clouds!(ij, vars, model.shortwave_radiation.clouds, model)
    @test clouds_state.cloud_top == 5
    @test clouds_state.cloud_cover ≈ ((0.95 - 0.3) / (1 - 0.3))^2 rtol = 1.0e-4

    # cloud cover and cloud top height [m] from geopotential stored for output
    @test vars.parameterizations.cloud_cover[ij] == clouds_state.cloud_cover
    g = model.planet.gravity
    @test vars.parameterizations.cloud_top_height[ij] == vars.dynamics.geopotential[ij, 5] / g

    # precipitation cloud top higher than humidity cloud top wins
    vars.parameterizations.cloud_top[ij] = 3
    @test SpeedyWeather.clouds!(ij, vars, model.shortwave_radiation.clouds, model).cloud_top == 3
    @test vars.parameterizations.cloud_top_height[ij] == vars.dynamics.geopotential[ij, 3] / g

    # no cloud: height 0
    vars.grid.humidity[ij, :, :] .= 0
    vars.parameterizations.cloud_top[ij] = 9
    @test SpeedyWeather.clouds!(ij, vars, model.shortwave_radiation.clouds, model).cloud_top == 9
    @test vars.parameterizations.cloud_top_height[ij] == 0
end
