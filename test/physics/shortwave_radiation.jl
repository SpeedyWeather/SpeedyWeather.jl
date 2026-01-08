function init_shortwave_state!(diagn, model)
    diagn.grid.temp_grid_prev .= 280
    diagn.grid.humid_grid_prev .= 1e-3
    diagn.grid.pres_grid_prev .= 1e5
    diagn.physics.cloud_top .= model.spectral_grid.nlayers + 1
    diagn.physics.rain_rate .= 0
    diagn.physics.ocean.albedo .= 0.5
    diagn.physics.land.albedo .= 0.3
    return nothing
end

@testset "Shortwave radiation" begin
    spectral_grid = SpectralGrid(trunc = 31, nlayers = 8)
    @testset for SW in (Nothing, TransparentShortwave, OneBandShortwave, OneBandGreyShortwave)
        sw = SW(spectral_grid)
        model = PrimitiveWetModel(spectral_grid; shortwave_radiation = sw)

        initialize!(model.shortwave_radiation, model)

        progn = PrognosticVariables(model)
        diagn = DiagnosticVariables(model)
        init_shortwave_state!(diagn, model)

        (; time) = progn.clock
        SpeedyWeather.cos_zenith!(diagn, time, model)

        for ij in 1:model.spectral_grid.npoints
            SpeedyWeather.parameterization!(ij, diagn, progn, model.shortwave_radiation, model)

            # top of atmosphere radiation down
            trd = model.planet.solar_constant * diagn.physics.cos_zenith[ij]

            if !(sw isa Nothing)
                osr = diagn.physics.outgoing_shortwave[ij]
                ssrd = diagn.physics.surface_shortwave_down[ij]
                @test 0 <= osr <= ssrd <= trd
                @test isfinite(osr)
                @test isfinite(ssrd)
            end
        end

        if !(sw isa Nothing)
            @test any(diagn.physics.outgoing_shortwave .> 0)
            @test any(diagn.physics.surface_shortwave_down .> 0)
        else
            @test !haskey(diagn.physics, :outgoing_shortwave)
            @test all(diagn.physics.surface_shortwave_down .== 0)
        end
    end
end

@testset "Shortwave radiation transmissivity" begin
    spectral_grid = SpectralGrid(trunc = 31, nlayers = 8)
    @testset for T in (TransparentShortwaveTransmissivity, BackgroundShortwaveTransmissivity)
        transmissivity = T(spectral_grid)
        sw = OneBandShortwave(spectral_grid; transmissivity = transmissivity)
        model = PrimitiveWetModel(spectral_grid; shortwave_radiation = sw)

        initialize!(model.shortwave_radiation, model)

        progn = PrognosticVariables(model)
        diagn = DiagnosticVariables(model)
        init_shortwave_state!(diagn, model)
        diagn.physics.cos_zenith .= 1

        for ij in 1:model.spectral_grid.npoints
            clouds = SpeedyWeather.clouds!(ij, diagn, progn, model.shortwave_radiation.clouds, model)
            t = SpeedyWeather.transmissivity!(ij, diagn, progn, clouds, model.shortwave_radiation.transmissivity, model)
            @test all(0 .<= t[ij, :] .<= 1)
        end
    end
end

@testset "Shortwave radiation clouds" begin
    spectral_grid = SpectralGrid(trunc = 31, nlayers = 8)
    @testset for C in (DiagnosticClouds, NoClouds)
        clouds = C(spectral_grid)
        sw = OneBandShortwave(spectral_grid; clouds = clouds)
        model = PrimitiveWetModel(spectral_grid; shortwave_radiation = sw)

        initialize!(model.shortwave_radiation, model)

        progn = PrognosticVariables(model)
        diagn = DiagnosticVariables(model)
        init_shortwave_state!(diagn, model)
        diagn.physics.cos_zenith .= 1

        ij = rand(1:model.spectral_grid.npoints)
        clouds_state = SpeedyWeather.clouds!(ij, diagn, progn, model.shortwave_radiation.clouds, model)
        @test 0 <= clouds_state.cloud_cover <= 1
        @test 1 <= clouds_state.cloud_top <= model.spectral_grid.nlayers + 1
        SpeedyWeather.parameterization!(ij, diagn, progn, model.shortwave_radiation, model)
    end
end
