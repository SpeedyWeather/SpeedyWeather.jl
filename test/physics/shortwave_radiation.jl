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
    
    @testset "OneBandShortwave branches coverage" begin
        """Minimal tests to ensure all conditional branches in OneBandShortwave are covered."""
        
        # Branch: cloud_top_reflection = true (default, applies cloud reflection)
        @testset "cloud_top_reflection=true" begin
            sw = OneBandShortwave(spectral_grid; cloud_top_reflection=true)
            model = PrimitiveWetModel(spectral_grid; shortwave_radiation=sw)
            sim = initialize!(model)
            run!(sim, period=Day(3))
            
            osr = sim.diagnostic_variables.physics.outgoing_shortwave_radiation
            @test all(osr .>= 0)
            @test !any(isnan, osr)
        end
        
        # Branch: cloud_top_reflection = false (no cloud reflection)
        @testset "cloud_top_reflection=false" begin
            sw = OneBandShortwave(spectral_grid; cloud_top_reflection=false)
            model = PrimitiveWetModel(spectral_grid; shortwave_radiation=sw)
            sim = initialize!(model)
            run!(sim, period=Day(3))
            
            osr = sim.diagnostic_variables.physics.outgoing_shortwave_radiation
            @test all(osr .>= 0)
            @test !any(isnan, osr)
        end
        
        # Branch: use_stratocumulus = true (default, stratocumulus parameterization active)
        @testset "use_stratocumulus=true" begin
            sw = OneBandShortwave(spectral_grid; use_stratocumulus=true)
            model = PrimitiveWetModel(spectral_grid; shortwave_radiation=sw)
            sim = initialize!(model)
            run!(sim, period=Day(3))
            
            ssrd = sim.diagnostic_variables.physics.surface_shortwave_down
            @test all(ssrd .>= 0)
            @test !any(isnan, ssrd)
        end
        
        # Branch: use_stratocumulus = false (no stratocumulus)
        @testset "use_stratocumulus=false" begin
            sw = OneBandShortwave(spectral_grid; use_stratocumulus=false)
            model = PrimitiveWetModel(spectral_grid; shortwave_radiation=sw)
            sim = initialize!(model)
            run!(sim, period=Day(3))
            
            ssrd = sim.diagnostic_variables.physics.surface_shortwave_down
            @test all(ssrd .>= 0)
            @test !any(isnan, ssrd)
        end
        
        # Branch: use_speedy_transmittance = true (default, SPEEDY-style transmittance)
        @testset "use_speedy_transmittance=true" begin
            sw = OneBandShortwave(spectral_grid; use_speedy_transmittance=true)
            model = PrimitiveWetModel(spectral_grid; shortwave_radiation=sw)
            sim = initialize!(model)
            run!(sim, period=Day(3))
            
            ssrd = sim.diagnostic_variables.physics.surface_shortwave_down
            osr = sim.diagnostic_variables.physics.outgoing_shortwave_radiation
            @test all(ssrd .>= 0)
            @test all(osr .>= 0)
        end
        
        # Branch: use_speedy_transmittance = false (simpler transmittance)
        @testset "use_speedy_transmittance=false" begin
            sw = OneBandShortwave(spectral_grid; use_speedy_transmittance=false)
            model = PrimitiveWetModel(spectral_grid; shortwave_radiation=sw)
            sim = initialize!(model)
            run!(sim, period=Day(3))
            
            ssrd = sim.diagnostic_variables.physics.surface_shortwave_down
            osr = sim.diagnostic_variables.physics.outgoing_shortwave_radiation
            @test all(ssrd .>= 0)
            @test all(osr .>= 0)
        end
    end
end