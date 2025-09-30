@testset "Longwave radiation" begin
    spectral_grid = SpectralGrid(trunc=31, nlayers=8)
    @testset for LW in (Nothing, UniformCooling, JeevanjeeRadiation, )#NBandRadiation)
        longwave_radiation = LW(spectral_grid)
        optical_depth = FriersonOpticalDepth(spectral_grid)

        # with no longwave radiation or uniform cooling there are no 
        # longwave surface radiative fluxes cooling the soil
        # which will heat up like crazy, use prescribed instead
        if longwave_radiation isa Union{Nothing, UniformCooling}
            soil_temperature = SeasonalLandTemperature(spectral_grid)
        else
            soil_temperature = LandBucketTemperature(spectral_grid)
        end

        land = LandModel(spectral_grid, temperature=soil_temperature)

        model = PrimitiveWetModel(spectral_grid; land, optical_depth, longwave_radiation)
        simulation = initialize!(model)
        run!(simulation, period=Day(5))
        
        temp = simulation.diagnostic_variables.grid.temp_grid[:, end]
        
        # just test that the surface temperature isn't completely off
        @test all(200 .< temp .< 330)
    end
end