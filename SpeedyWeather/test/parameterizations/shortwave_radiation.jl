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

@testset "Ozone seasonal cycle synchronized with solar declination" begin
    spectral_grid = SpectralGrid(truncation = 31, nlayers = 8)

    # ratio of lower stratospheric ozone absorption in the northernmost over southernmost ring
    function north_south_ratio(time; kwargs...)
        planet = Earth(spectral_grid; kwargs...)
        model = PrimitiveWetModel(spectral_grid; planet, shortwave_radiation = TwoBandShortwave(spectral_grid))
        simulation = initialize!(model; time)
        vars = simulation.variables
        SpeedyWeather.parameterization!(vars, model.shortwave_radiation, model)
        lower = vars.parameterizations.ozone_absorption_lower
        return lower[1] / lower[end]
    end

    @test north_south_ratio(DateTime(2000, 1, 1)) > 1                   # northern winter
    @test north_south_ratio(DateTime(2000, 7, 1)) ≈ 1                   # northern summer, symmetric
    @test north_south_ratio(DateTime(2000, 1, 1), axial_tilt = 0) ≈ 1   # no seasons without tilt
    @test north_south_ratio(DateTime(2000, 1, 1), axial_tilt = -23.4) ≈ 1   # flipped seasons

    # seasons follow the planet's equinox, shifted by half a year here
    @test north_south_ratio(DateTime(2000, 7, 1), equinox = DateTime(2000, 9, 20)) > 1
end

@testset "Diagnostic clouds in cold, dry air and from snow" begin
    spectral_grid = SpectralGrid(truncation = 31, nlayers = 8)
    model = PrimitiveWetModel(spectral_grid)
    clouds = model.shortwave_radiation.clouds
    vars = Variables(model)
    init_shortwave_state!(vars, model)
    nlayers = spectral_grid.nlayers
    ij = 1

    # cold air with humidity below the q threshold but RH > threshold: clouds form above the surface layer
    vars.grid.temperature .= 240
    vars.grid.humidity .= clouds.specific_humidity_threshold_min / 2
    state = SpeedyWeather.clouds!(ij, vars, clouds, model)
    @test state.cloud_cover > 0
    @test state.cloud_top == nlayers - 1

    # completely dry air, snow alone creates clouds
    vars.grid.humidity .= 0
    vars.parameterizations.cloud_top .= nlayers + 1
    @test SpeedyWeather.clouds!(ij, vars, clouds, model).cloud_cover == 0
    vars.parameterizations.snow_rate .= 1 / (86400 * 1000)     # 1 mm/day
    @test SpeedyWeather.clouds!(ij, vars, clouds, model).cloud_cover ≈ clouds.precipitation_weight
end

@testset "Free troposphere layers for any vertical resolution" begin
    for nlayers in (1, 2, 3, 4, 8, 16, 32, 64)
        spectral_grid = SpectralGrid(truncation = 21, nlayers = nlayers)
        model = PrimitiveWetModel(spectral_grid)
        clouds = model.shortwave_radiation.clouds
        layer_top, cloud_base = SpeedyWeather.free_troposphere_layers(clouds, model)
        σ = model.geometry.σ_levels_full

        @test 1 <= layer_top <= cloud_base <= nlayers
        nlayers > 1 && @test cloud_base < nlayers       # clouds never in the surface layer

        # σ boundaries are respected unless that would leave no layer for clouds
        if nlayers >= 4
            @test σ[layer_top] >= clouds.σ_tropopause
            @test σ[cloud_base] <= clouds.σ_boundary_layer
            @test σ[cloud_base + 1] > clouds.σ_boundary_layer || cloud_base == nlayers - 1
        end
        nlayers == 8 && @test (layer_top, cloud_base) == (2, 7)
    end
end
