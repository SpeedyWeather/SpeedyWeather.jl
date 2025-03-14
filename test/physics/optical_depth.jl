@testset "Optical depth" begin
    spectral_grid = SpectralGrid(trunc=31, nlayers=8)

    for OD in (ZeroOpticalDepth, FriersonOpticalDepth)
        optical_depth = OD(spectral_grid)
        longwave_radiation = NBandRadiation(spectral_grid)
        
        model = PrimitiveWetModel(spectral_grid; optical_depth, longwave_radiation)
        simulation = initialize!(model)
        run!(simulation, steps=1)
    
        band = longwave_radiation.nbands
        t = simulation.diagnostic_variables.column.transmittance_longwave[:, band]
        @test all(0 .<= t .<= 1)
    end
end