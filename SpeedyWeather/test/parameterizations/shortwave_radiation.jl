function init_shortwave_state!(vars, model)
    vars.grid.temp_prev .= 280
    vars.grid.humid_prev .= 1.0e-3
    vars.grid.pres_prev .= 1.0e5
    vars.parameterizations.cloud_top .= model.spectral_grid.nlayers + 1
    vars.parameterizations.rain_rate .= 0
    vars.parameterizations.ocean.albedo .= 0.5
    vars.parameterizations.land.albedo .= 0.3
    return nothing
end

@testset "Shortwave radiation" begin
    spectral_grid = SpectralGrid(trunc = 31, nlayers = 8)
    @testset for SW in (Nothing, TransparentShortwave, OneBandShortwave, OneBandGreyShortwave)
        sw = SW(spectral_grid)
        model = PrimitiveWetModel(spectral_grid; shortwave_radiation = sw)

        initialize!(model.shortwave_radiation, model)

        vars = Variables(model)
        init_shortwave_state!(vars, model)

        SpeedyWeather.parameterization!(vars, model.solar_zenith, model)

        for ij in 1:model.spectral_grid.npoints
            SpeedyWeather.parameterization!(ij, vars, model.shortwave_radiation, model)

            # top of atmosphere radiation down
            trd = model.planet.solar_constant * vars.parameterizations.cos_zenith[ij]

            if !(sw isa Nothing)
                osr = vars.parameterizations.outgoing_shortwave[ij]
                ssrd = vars.parameterizations.surface_shortwave_down[ij]
                @test 0 <= osr <= ssrd <= trd
                @test isfinite(osr)
                @test isfinite(ssrd)
            end
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
    spectral_grid = SpectralGrid(trunc = 31, nlayers = 8)
    @testset for T in (TransparentShortwaveTransmissivity, BackgroundShortwaveTransmissivity)
        transmissivity = T(spectral_grid)
        sw = OneBandShortwave(spectral_grid; transmissivity = transmissivity)
        model = PrimitiveWetModel(spectral_grid; shortwave_radiation = sw)

        initialize!(model.shortwave_radiation, model)

        vars = Variables(model)
        init_shortwave_state!(vars, model)
        vars.parameterizations.cos_zenith .= 1

        for ij in 1:model.spectral_grid.npoints
            clouds = SpeedyWeather.clouds!(ij, vars, model.shortwave_radiation.clouds, model)
            t = SpeedyWeather.transmissivity!(ij, vars, clouds, model.shortwave_radiation.transmissivity, model)
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
