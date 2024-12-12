@testset "Optical depth" begin
    spectral_grid = SpectralGrid(trunc=31, nlayers=8)

    for OD in (ZeroOpticalDepth, FriersonOpticalDepth)
        optical_depth = OD(spectral_grid)
        longwave_radiation = NBandRadiation(spectral_grid)
        
        model = PrimitiveWetModel(spectral_grid; optical_depth, longwave_radiation)
        simulation = initialize!(model)
        run!(simulation, period=Day(1))
    
        band = longwave_radiation.nbands
        dτ = simulation.diagnostic_variables.column.optical_depth_longwave[:, band]

        (; nlayers) = spectral_grid

        for k in 2:nlayers
            # optical depth cannot be negative
            @test dτ[k] >= 0
        end
    end
end