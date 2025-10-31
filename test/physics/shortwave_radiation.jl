@testset "Shortwave radiation" begin
    spectral_grid = SpectralGrid(trunc=31, nlayers=8)
    
    @testset for SW in (TransparentShortwave, OneBandShortwave)

            sw = SW(spectral_grid)
            model = PrimitiveWetModel(spectral_grid; shortwave_radiation=sw)
            simulation = initialize!(model)
            run!(simulation, period=Day(5))
            ssrd = simulation.diagnostic_variables.physics.surface_shortwave_down
            TRD = model.planet.solar_constant * simulation.diagnostic_variables.physics.cos_zenith
            @test all(0 .<= ssrd .<= TRD)
    end
end