import Statistics: mean

@testset "Albedo" begin
    for Model in (PrimitiveWetModel, PrimitiveDryModel)
        spectral_grid = SpectralGrid(trunc=31, nlayers=8)

        # albedo climatology for both ocean and land
        albedo1 = AlbedoClimatology(spectral_grid)

        # constant albedo for ocean and land
        albedo2 = Albedo(
            ocean=GlobalConstantAlbedo(spectral_grid, albedo=0.1),
            land=GlobalConstantAlbedo(spectral_grid, albedo=0.3))

        # constant albedo for ocean, climatology for land
        albedo3 = Albedo(
            ocean=GlobalConstantAlbedo(spectral_grid, albedo=0.1),
            land=AlbedoClimatology(spectral_grid))

        for albedo in (albedo1, albedo2, albedo3)

            model = Model(spectral_grid; albedo)
            simulation = initialize!(model)
            run!(simulation, period=Day(1))

            a, b = extrema(simulation.diagnostic_variables.physics.albedo)
            m = mean(simulation.diagnostic_variables.physics.albedo)
            @test a >= 0
            @test b <= 1
            @test 0 < m < 1
        end
    end
end