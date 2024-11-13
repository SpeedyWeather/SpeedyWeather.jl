@testset "Ocean models" begin
    
    spectral_grid = SpectralGrid(trunc=31, nlayers=5)

    for OceanModel in ( SeasonalOceanClimatology, 
                        ConstantOceanClimatology,
                        AquaPlanet)

        ocean = OceanModel(spectral_grid)
        model = PrimitiveWetModel(spectral_grid; ocean)
        model.feedback.verbose = false
        simulation = initialize!(model, time=DateTime(2000,5,1))
        run!(simulation, period=Day(5))     # run for 5 days as ocean time step is 3 days by default
    end
end