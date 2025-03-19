import Statistics: mean

@testset "Single Albedos" begin
    for Model in (PrimitiveWetModel, PrimitiveDryModel)
        spectral_grid = SpectralGrid(trunc=31, nlayers=8)

        for AlbedoType in (GlobalConstantAlbedo, ManualAlbedo, AlbedoClimatology)

            albedo = AlbedoType(spectral_grid)
            AlbedoType == ManualAlbedo && set!(albedo, 0.06)

            model = Model(spectral_grid; albedo)
            simulation = initialize!(model)
            run!(simulation, steps=1)

            a, b = extrema(simulation.diagnostic_variables.physics.albedo)
            m = mean(simulation.diagnostic_variables.physics.albedo)
            @test a >= 0
            @test b <= 1
            @test 0 < m < 1
        end
    end
end

@testset "Albedo composites" begin
    for Model in (PrimitiveWetModel, PrimitiveDryModel)
        spectral_grid = SpectralGrid(trunc=31, nlayers=8)

        for AlbedoType1 in (GlobalConstantAlbedo, ManualAlbedo, AlbedoClimatology)
            for AlbedoType2 in (GlobalConstantAlbedo, ManualAlbedo, AlbedoClimatology)

                albedo = Albedo(
                    ocean=AlbedoType1(spectral_grid),
                    land=AlbedoType2(spectral_grid))

                AlbedoType1 == ManualAlbedo && set!(albedo.ocean, 0.06)
                AlbedoType2 == ManualAlbedo && set!(albedo.land, 0.06)

                model = Model(spectral_grid; albedo)
                simulation = initialize!(model)
                run!(simulation, steps=1)

                a, b = extrema(simulation.diagnostic_variables.physics.albedo)
                m = mean(simulation.diagnostic_variables.physics.albedo)
                @test a >= 0
                @test b <= 1
                @test 0 < m < 1
            end
        end
    end
end