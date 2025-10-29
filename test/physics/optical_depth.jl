@testset "Transmittance" begin
    spectral_grid = SpectralGrid(trunc=31, nlayers=8)

    @testset for TR in (TransparentTransmittance, FriersonTransmittance)
        transmittance = TR(spectral_grid)
        longwave_radiation = NBandRadiation(spectral_grid)

        model = PrimitiveWetModel(spectral_grid; transmittance, longwave_radiation)
        simulation = initialize!(model)
        run!(simulation, period=Day(1))
    
        band = longwave_radiation.nbands
        t = simulation.diagnostic_variables.column.transmittance_longwave[:, band]

        (; nlayers) = spectral_grid

        for k in 2:nlayers
            # transmittance has to be between 0 and 1
            # because optical depth non-negative
            @test 1 >= t[k] >= 0
        end
    end
end