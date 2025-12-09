@testset "Transmissivity" begin
    spectral_grid = SpectralGrid(trunc=31, nlayers=8)

    @testset "Longwave Transmissivity (via Longwave Radiation)" begin
        @testset for LW in (JeevanjeeRadiation,)
            longwave_radiation = LW(spectral_grid)

            model = PrimitiveWetModel(spectral_grid; longwave_radiation)
            simulation = initialize!(model)
            run!(simulation, period=Day(1))
        
            # Test that transmissivity values are physically reasonable
            t = simulation.diagnostic_variables.column.transmissivity_longwave

            (; nlayers) = spectral_grid
            nbands = size(t, 2)

            for band in 1:nbands
                for k in 1:nlayers
                    # transmissivity has to be between 0 and 1
                    # because optical depth non-negative
                    @test 1 >= t[k, band] >= 0
                end
            end
        end
    end

    @testset "Shortwave Transmissivity (via Shortwave Radiation)" begin
        @testset for SW in (TransparentShortwave, OneBandGreyShortwave)
            shortwave_radiation = SW(spectral_grid)

            model = PrimitiveWetModel(spectral_grid; shortwave_radiation)
            simulation = initialize!(model)
            run!(simulation, period=Day(1))
        
            # Test shortwave transmissivity
            t = simulation.diagnostic_variables.column.transmissivity_shortwave

            (; nlayers) = spectral_grid
            nbands = size(t, 2)

            for band in 1:nbands
                for k in 1:nlayers
                    # transmissivity has to be between 0 and 1
                    # because optical depth non-negative
                    @test 1 >= t[k, band] >= 0
                end
            end
        end
    end
end