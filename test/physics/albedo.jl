import Statistics: mean

@testset "Single Albedos" begin
    @testset for Model in (PrimitiveWetModel, PrimitiveDryModel)
        spectral_grid = SpectralGrid(trunc=31, nlayers=8)

        @testset for AlbedoType in (GlobalConstantAlbedo, ManualAlbedo, AlbedoClimatology)

            albedo = AlbedoType(spectral_grid)
            AlbedoType == ManualAlbedo && set!(albedo, 0.06)

            model = Model(spectral_grid; albedo)

            initialize!(albedo, model)
            progn = PrognosticVariables(model)
            diagn = DiagnosticVariables(model)
            SpeedyWeather.parameterization_tendencies!(diagn, progn, model)

            a, b = extrema(diagn.physics.albedo)
            m = mean(diagn.physics.albedo)
            @test a >= 0
            @test b <= 1
            @test 0 < m < 1
        end
    end
end

@testset "Albedo composites" begin
    for Model in (PrimitiveWetModel, PrimitiveDryModel)
        spectral_grid = SpectralGrid(trunc=31, nlayers=8)

        @testset for OceanAlbedo in (GlobalConstantAlbedo, ManualAlbedo, AlbedoClimatology, OceanSeaIceAlbedo)
            @testset for LandAlbedo in (GlobalConstantAlbedo, ManualAlbedo, AlbedoClimatology, LandSnowAlbedo)

                albedo = OceanLandAlbedo(
                    ocean=OceanAlbedo(spectral_grid),
                    land=LandAlbedo(spectral_grid))

                OceanAlbedo == ManualAlbedo && set!(albedo.ocean, 0.06)
                LandAlbedo == ManualAlbedo && set!(albedo.land, 0.06)

                model = Model(spectral_grid; albedo)

                initialize!(albedo, model)
                progn = PrognosticVariables(model)
                diagn = DiagnosticVariables(model)
                SpeedyWeather.parameterization_tendencies!(diagn, progn, model)

                a, b = extrema(diagn.physics.land.albedo)
                m = mean(diagn.physics.land.albedo)
                @test a >= 0
                @test b <= 1
                @test 0 < m < 1

                a, b = extrema(diagn.physics.ocean.albedo)
                m = mean(diagn.physics.ocean.albedo)
                @test a >= 0
                @test b <= 1
                @test 0 < m < 1                
            end
        end
    end
end