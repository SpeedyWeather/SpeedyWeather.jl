import Statistics: mean

@testset "Single Albedos" begin
    @testset for Model in (PrimitiveWetModel, PrimitiveDryModel)
        spectral_grid = SpectralGrid(trunc = 31, nlayers = 8)

        @testset for AlbedoType in (GlobalConstantAlbedo, ManualAlbedo, AlbedoClimatology)

            albedo = AlbedoType(spectral_grid)
            AlbedoType == ManualAlbedo && set!(albedo, 0.06)

            model = Model(spectral_grid; albedo)

            initialize!(albedo, model)
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
        spectral_grid = SpectralGrid(trunc = 31, nlayers = 8)

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
