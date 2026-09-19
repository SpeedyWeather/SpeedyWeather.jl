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
    @testset for SW in (Nothing, TransparentShortwave, OneBandShortwave, OneBandGreyShortwave, TwoBandShortwave)
        sw = SW(spectral_grid)
        model = PrimitiveWetModel(spectral_grid; shortwave_radiation = sw)

        initialize!(model.shortwave_radiation, model)

        vars = Variables(model)
        init_shortwave_state!(vars, model)

        SpeedyWeather.parameterization!(vars, model.solar_zenith, model)
        SpeedyWeather.parameterization!(vars, model.shortwave_radiation, model)

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

@testset "Two-band shortwave transmissivity" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
    model = PrimitiveWetModel(spectral_grid; shortwave_radiation = TwoBandShortwave(spectral_grid))
    initialize!(model.shortwave_radiation, model)

    vars = Variables(model)
    init_shortwave_state!(vars, model)
    vars.parameterizations.cos_zenith .= 1

    for ij in 1:model.spectral_grid.npoints
        clouds = SpeedyWeather.clouds!(ij, vars, model.shortwave_radiation.clouds, model)
        t = SpeedyWeather.transmissivity!(ij, vars, clouds, model.shortwave_radiation.transmissivity, model)
        @test t.zenith_factor == 1
    end
    @test all(0 .< vars.scratch.grid.a .<= 1)   # visible
    @test all(0 .< vars.scratch.grid.b .<= 1)   # near infrared
    # near-infrared band is more strongly absorbed by water vapor
    @test all(vars.scratch.grid.b .< vars.scratch.grid.a)
end

@testset "Ozone" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
    @testset for ozone in (SeasonalOzone(spectral_grid), NoOzone(spectral_grid))
        sw = TwoBandShortwave(spectral_grid; ozone)
        model = PrimitiveWetModel(spectral_grid; shortwave_radiation = sw)
        simulation = initialize!(model)
        vars = simulation.variables
        SpeedyWeather.parameterization!(vars, model.shortwave_radiation, model)

        nlayers = spectral_grid.nlayers
        σ_half = model.geometry.σ_levels_half
        absorption = [SpeedyWeather.ozone_absorption(ij, k, vars, ozone, model) for ij in 1:spectral_grid.npoints, k in 1:nlayers]
        @test all(0 .<= absorption .< 1)

        if ozone isa NoOzone
            @test all(absorption .== 0)
        else
            # only absorbed in the stratosphere, total is upper + lower stratosphere
            @test all(all(absorption[:, k] .== 0) for k in 1:nlayers if σ_half[k] >= ozone.σ_lower)
            lower = vars.parameterizations.ozone_absorption_lower.data
            @test vec(sum(absorption, dims = 2)) ≈ ozone.absorption * ozone.upper_fraction .+ lower

            # more ozone absorption towards the poles than at the equator
            @test maximum(lower) > minimum(lower) >= 0
        end
    end
end
