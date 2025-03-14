@testset "Optical depth" begin
    spectral_grid = SpectralGrid(trunc=31, nlayers=8)

    for OD in (ZeroOpticalDepth, FriersonOpticalDepth)
        optical_depth = OD(spectral_grid)
        longwave_radiation = NBandRadiation(spectral_grid)
        
        model = PrimitiveWetModel(spectral_grid; optical_depth, longwave_radiation)
        simulation = initialize!(model)
        run!(simulation, steps=1)
    
        t = simulation.diagnostic_variables.column.transmittance_longwave
        @test all(0 .<= t .<= 1)
        @test any(t .> 0)
        
        t = simulation.diagnostic_variables.column.transmittance_shortwave
        @test all(0 .<= t .<= 1)
    end
end