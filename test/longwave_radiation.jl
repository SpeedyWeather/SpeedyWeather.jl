@testset "Longwave radiation" begin
    spectral_grid = SpectralGrid(trunc=31, nlayers=8)
    @testset for LW in (UniformCooling, JeevanjeeRadiation, NBandRadiation)
        longwave_radiation = LW(spectral_grid)
        optical_depth = FriersonOpticalDepth(spectral_grid)
        model = PrimitiveWetModel(spectral_grid; optical_depth, longwave_radiation)
        simulation = initialize!(model)
        run!(simulation, period=Day(5))
        
        temp = simulation.diagnostic_variables.grid.temp_grid[:, end]

        for ij in eachindex(temp)
            # just check that surface temperature isn't completely off
            @test 200 < temp[ij] < 330
        end
    end
end