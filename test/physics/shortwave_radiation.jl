@testset "Shortwave radiation" begin
    spectral_grid = SpectralGrid(trunc=31, nlayers=8)
    
    @testset for SW in (TransparentShortwave, OneBandShortwave, OneBandGreyShortwave)
            sw = SW(spectral_grid)
            model = PrimitiveWetModel(spectral_grid; shortwave_radiation=sw)
            simulation = initialize!(model)
            run!(simulation, period=Day(5))
            ssrd = simulation.diagnostic_variables.physics.surface_shortwave_down
            TRD = model.planet.solar_constant * simulation.diagnostic_variables.physics.cos_zenith
            @test all(0 .<= ssrd .<= TRD)
    end
    
    @testset "OneBandShortwave component testing" begin
        """Test different cloud and transmittance combinations using convenience constructors."""
        
        # Test the full wet model version (with clouds)
        @testset "OneBandShortwave (wet model with clouds)" begin
            sw = OneBandShortwave(spectral_grid)
            model = PrimitiveWetModel(spectral_grid; shortwave_radiation=sw)
            sim = initialize!(model)
            run!(sim, period=Day(3))
            
            osr = sim.diagnostic_variables.physics.outgoing_shortwave_radiation
            ssrd = sim.diagnostic_variables.physics.surface_shortwave_down
            @test all(osr .>= 0)
            @test all(ssrd .>= 0)
            @test !any(isnan, osr)
            @test !any(isnan, ssrd)
        end
        
        # Test the dry model version (no clouds)
        @testset "OneBandGreyShortwave (dry model, no clouds)" begin
            sw = OneBandGreyShortwave(spectral_grid)
            model = PrimitiveWetModel(spectral_grid; shortwave_radiation=sw)
            sim = initialize!(model)
            run!(sim, period=Day(3))
            
            osr = sim.diagnostic_variables.physics.outgoing_shortwave_radiation
            ssrd = sim.diagnostic_variables.physics.surface_shortwave_down
            @test all(osr .>= 0)
            @test all(ssrd .>= 0)
            @test !any(isnan, osr)
            @test !any(isnan, ssrd)
        end
        
        # Test transparent shortwave (reference case)
        @testset "TransparentShortwave (transparent atmosphere)" begin
            sw = TransparentShortwave(spectral_grid)
            model = PrimitiveWetModel(spectral_grid; shortwave_radiation=sw)
            sim = initialize!(model)
            run!(sim, period=Day(3))
            
            osr = sim.diagnostic_variables.physics.outgoing_shortwave_radiation
            ssrd = sim.diagnostic_variables.physics.surface_shortwave_down
            @test all(osr .>= 0)
            @test all(ssrd .>= 0)
            @test !any(isnan, osr)
            @test !any(isnan, ssrd)
            
            # For transparent atmosphere, surface shortwave down should equal outgoing shortwave
            # (since it's just surface albedo reflection)
            @test all(ssrd .>= osr)  # Surface down >= reflected (outgoing)
        end
        
        # Test cloud parameterization with different settings
        @testset "DiagnosticClouds parameter variations" begin
            # Test with different cloud parameters
            clouds = DiagnosticClouds(
                spectral_grid,
                cloud_albedo = 0.6,  # Higher cloud albedo
                use_stratocumulus = false  # Disable stratocumulus
            )
            sw = OneBandShortwave(
                clouds,
                TransparentShortwaveTransmittance(spectral_grid),
                OneBandShortwave(spectral_grid).radiative_transfer
            )
            model = PrimitiveWetModel(spectral_grid; shortwave_radiation=sw)
            sim = initialize!(model)
            run!(sim, period=Day(3))
            
            osr = sim.diagnostic_variables.physics.outgoing_shortwave_radiation
            ssrd = sim.diagnostic_variables.physics.surface_shortwave_down
            @test all(osr .>= 0)
            @test all(ssrd .>= 0)
            @test !any(isnan, osr)
            @test !any(isnan, ssrd)
        end
    end
end