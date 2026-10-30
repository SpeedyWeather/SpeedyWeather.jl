import Statistics: mean

@testset "Single Albedos" begin
    @testset for Model in (PrimitiveWetModel, PrimitiveDryModel)
        spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)

        @testset for AlbedoType in (GlobalConstantAlbedo, ManualAlbedo, AlbedoClimatology)

            albedo = AlbedoType(spectral_grid)
            AlbedoType == ManualAlbedo && set!(albedo, 0.06)

            model = Model(spectral_grid; albedo)

            initialize!(albedo, model)
            initialize!(model.radiation, model)  # precomputes the free troposphere layers for clouds
            vars = Variables(model)
            SpeedyWeather.parameterization_tendencies!(vars, model)

            a, b = extrema(vars.parameterizations.albedo)
            m = mean(vars.parameterizations.albedo)
            @test a >= 0
            @test b <= 1
            @test 0 < m < 1
        end
    end
end

@testset "Albedo composites" begin
    for Model in (PrimitiveWetModel, PrimitiveDryModel)
        spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)

        @testset for OceanAlbedo in (GlobalConstantAlbedo, ManualAlbedo, AlbedoClimatology, OceanSeaIceAlbedo)
            @testset for LandAlbedo in (GlobalConstantAlbedo, ManualAlbedo, AlbedoClimatology, LandSnowAlbedo)

                albedo = OceanLandAlbedo(
                    ocean = OceanAlbedo(spectral_grid),
                    land = LandAlbedo(spectral_grid)
                )

                OceanAlbedo == ManualAlbedo && set!(albedo.ocean, 0.06)
                LandAlbedo == ManualAlbedo && set!(albedo.land, 0.06)

                model = Model(spectral_grid; albedo)

                initialize!(albedo, model)
                initialize!(model.radiation, model)  # precomputes the free troposphere layers for clouds
                vars = Variables(model)
                SpeedyWeather.parameterization_tendencies!(vars, model)

                a, b = extrema(vars.parameterizations.land.albedo)
                m = mean(vars.parameterizations.land.albedo)
                @test a >= 0
                @test b <= 1
                @test 0 < m < 1

                a, b = extrema(vars.parameterizations.ocean.albedo)
                m = mean(vars.parameterizations.ocean.albedo)
                @test a >= 0
                @test b <= 1
                @test 0 < m < 1
            end
        end
    end
end

@testset "LearnedLandAlbedo" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
    land_albedo = LearnedLandAlbedo(spectral_grid)
    @test isbitstype(typeof(land_albedo))       # GPU-compatible

    # coefficients of all learned expressions are exposed as parameters
    params = stripparams(parameters(land_albedo))
    components = (
        :visible_isotropic, :visible_volumetric, :visible_geometric,
        :near_infrared_isotropic, :near_infrared_volumetric, :near_infrared_geometric,
    )
    @test all(haskey(params, component) for component in components)
    @test params.visible_isotropic.c1 == land_albedo.visible_isotropic.c1

    # dry model has no snow, soil moisture or precipitation, so also covers the predictor fallbacks
    @testset for Model in (PrimitiveWetModel, PrimitiveDryModel)
        model = Model(spectral_grid; albedo = OceanLandAlbedo(spectral_grid; land = land_albedo))
        simulation = initialize!(model)
        vars = simulation.variables
        has_snow = haskey(vars.prognostic.land, :snow_depth)

        # put snow on land to exercise the snow terms, run long enough to catch instabilities
        has_snow && (vars.prognostic.land.snow_depth.data .= 0.05)
        run!(simulation, period = Day(2))

        land_points = model.land_sea_mask.land_fraction.data .> 0
        (; white_sky_albedo, black_sky_albedo) = vars.parameterizations.land
        for land_surface_albedo in (vars.parameterizations.land.albedo, white_sky_albedo, black_sky_albedo)
            @test all(isfinite, land_surface_albedo)

            # check that albedo falls within the physical bounds
            @test all(0 .<= land_surface_albedo.data[land_points] .<= 1)
        end

        # snow brightens the land surface on average
        if has_snow
            snowy_albedo = mean(white_sky_albedo.data[land_points])
            vars.prognostic.land.snow_depth.data .= 0
            SpeedyWeather.parameterization_tendencies!(vars, model)
            @test mean(white_sky_albedo.data[land_points]) < snowy_albedo
        end
    end
end
