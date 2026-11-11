using Dates

@testset "Ocean and sea ice models" begin

    spectral_grid = SpectralGrid(truncation = 32, nlayers = 5)

    # just test that these parameters can be set
    SlabOcean(spectral_grid)
    ThermodynamicSeaIce(spectral_grid, freezing_temperature = 273.15)

    @testset for OceanModel in (
            SeasonalOceanClimatology,
            ConstantOceanClimatology,
            AquaPlanet,
            SlabOcean,
            PrescribedOcean,
        )

        @testset for SeaIceModel in (
                ThermodynamicSeaIce,
                PrescribedSeaIce,
                Nothing,
            )

            ocean = OceanModel(spectral_grid)
            sea_ice = SeaIceModel(spectral_grid)
            albedo = OceanLandAlbedo(ocean = OceanSeaIceAlbedo(spectral_grid), land = AlbedoClimatology(spectral_grid))

            model = PrimitiveDryModel(spectral_grid; ocean, sea_ice, albedo)
            model.feedback.verbose = false
            simulation = initialize!(model, time = DateTime(2000, 5, 1))

            # PrescribedOcean does not initialize SST set to something reasonable if it is all zeros
            if all(simulation.variables.prognostic.ocean.sea_surface_temperature .== 0)
                simulation.variables.prognostic.ocean.sea_surface_temperature .= 285
            end

            run!(simulation, period = Day(1))

            @test simulation.model.feedback.nans_detected == false
            @test haskey(simulation.variables.prognostic.ocean, :sea_surface_temperature)

            @test ArrayDimensions.hastime(simulation.variables.prognostic.ocean.sea_surface_temperature) == (ocean isa SpeedyWeather.AbstractDynamicOcean)
            sst = get_step(simulation.variables.prognostic.ocean.sea_surface_temperature, 1)
            @test ndims(sst) == 1
            @test all(0 .<= sst .<= 330)

            if sea_ice isa Nothing
                @test !haskey(simulation.variables.prognostic.ocean, :sea_ice_concentration)
            else
                @test ArrayDimensions.hastime(simulation.variables.prognostic.ocean.sea_ice_concentration) == (sea_ice isa SpeedyWeather.AbstractDynamicSeaIce)
                sic = get_step(simulation.variables.prognostic.ocean.sea_ice_concentration, 1)
                @test ndims(sic) == 1
                @test all(0 .<= sic .<= 1)
            end
        end
    end
end

@testset "Monthly interpolation of climatologies" begin
    # months and weight as computed with Dates' month arithmetic, 2000 is a leap year
    times = DateTime(2000, 1, 1):Hour(13):DateTime(2002, 1, 1)
    @test all(times) do time
        this_month, next_month, weight = SpeedyWeather.monthly_interpolation(time)
        this_month == Dates.month(time) && next_month == mod1(this_month + 1, 12) &&
            weight == Dates.days(time - Dates.firstdayofmonth(time)) / Dates.daysinmonth(time)
    end
    @test SpeedyWeather.monthly_interpolation(DateTime(2000, 2, 29, 12)) == (2, 3, 28 / 29)
    @test SpeedyWeather.monthly_interpolation(DateTime(2001, 12, 31, 23)) == (12, 1, 30 / 31)
    @test SpeedyWeather.monthly_interpolation(DateTime(2100, 2, 15)) == (2, 3, 14 / 28)  # not a leap year

    # weight in the precision of the interpolated variable, as it's a kernel argument (no Float64 on Metal)
    @test SpeedyWeather.monthly_interpolation(DateTime(2000, 2, 29, 12), Float32) === (2, 3, 28.0f0 / 29.0f0)
end

@testset "Coastal ocean initial conditions without land fallback" begin
    # the SST climatology is NaN over land, coastal ocean points have to interpolate from the
    # surrounding ocean, not fall back to the land value. Excluding the Antarctic ice shelves
    # (Ross, Weddell Sea), ocean in the land-sea mask but without any SST data nearby
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
    model = PrimitiveWetModel(spectral_grid)
    simulation = initialize!(model)
    (; sea_surface_temperature) = simulation.variables.prognostic.ocean
    _, latds = RingGrids.get_londlatds(spectral_grid.grid)
    mostly_ocean = (model.land_sea_mask.land_fraction.data .<= 0.5) .& (latds .> -75)
    @test !any(sea_surface_temperature.data[mostly_ocean, :] .== model.ocean.land_temperature)
end

@testset "Snow melting into the slab ocean" begin
    spectral_grid = SpectralGrid(truncation = 21, nlayers = 5)
    model = PrimitiveWetModel(spectral_grid; ocean = SlabOcean(spectral_grid))
    model.feedback.verbose = false
    simulation = initialize!(model, time = DateTime(2000, 5, 1))
    vars = simulation.variables
    ocean = model.ocean
    ocean_stepper = SpeedyWeather.time_stepper(model.time_stepping, :ocean)
    land_stepper = SpeedyWeather.time_stepper(model.time_stepping, :land)
    land_fraction = model.land_sea_mask.land_fraction.data
    all_ocean = land_fraction .== 0

    dsst = SpeedyWeather.get_tendency_step(vars.tendencies.ocean.sea_surface_temperature, ocean_stepper, ocean)
    sst = SpeedyWeather.get_prognostic_step(vars.prognostic.ocean.sea_surface_temperature, ocean_stepper, ocean)
    snow_depth = SpeedyWeather.get_prognostic_step(vars.prognostic.land.snow_depth, land_stepper, ocean)

    # reference tendency without snow
    vars.parameterizations.snow_rate .= 0
    snow_depth.data[all_ocean] .= 0
    SpeedyWeather.timestep!(vars, ocean, model)
    dsst_reference = copy(dsst.data)

    # snowfall cools by melting and warming the meltwater
    (; latent_heat_fusion) = model.atmosphere
    water_density = ocean.density
    s = 1.0e-7                                          # snow rate [m/s]
    vars.parameterizations.snow_rate .= s
    SpeedyWeather.timestep!(vars, ocean, model)
    Q = water_density * s .* (latent_heat_fusion .+ ocean.specific_heat_capacity .* (sst.data .- ocean.melting_temperature))
    expected = dsst_reference .- Q ./ ocean.heat_capacity_mixed_layer
    ocean_points = land_fraction .< 1
    @test dsst.data[ocean_points] ≈ expected[ocean_points]

    # snow on the ground over 100% ocean also melts into the ocean, then removed by filter!
    vars.parameterizations.snow_rate .= 0
    snow_depth.data[all_ocean] .= 0.2
    SpeedyWeather.timestep!(vars, ocean, model)
    @test all(dsst.data[all_ocean] .< dsst_reference[all_ocean])
    @test dsst.data[.!all_ocean] ≈ dsst_reference[.!all_ocean]

    filter!(vars, model.land.snow, model)
    @test all(vars.prognostic.land.snow_depth.data[all_ocean, :] .== 0)

    # and a short run from global snow stays finite without snow over ocean
    vars.prognostic.land.snow_depth .= 0.2
    run!(simulation, period = Day(1))
    @test !model.feedback.nans_detected
    @test all(vars.prognostic.land.snow_depth.data[all_ocean, :] .== 0)
end
