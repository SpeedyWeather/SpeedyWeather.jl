@testset "Large-scale condensation" begin
    spectral_grid = SpectralGrid(trunc=31, nlayers=8)
    @testset for r in (0.8, 0.9, 1.0)
        large_scale_condensation = ImplicitCondensation(spectral_grid, relative_humidity_threshold=r)
        model = PrimitiveWetModel(spectral_grid; large_scale_condensation)
        simulation = initialize!(model)
        run!(simulation, period=Day(5))
        
        precip = simulation.diagnostic_variables.physics.precip_large_scale

        for ij in eachindex(precip)
            @test precip[ij] >= 0   # precipitation should always be non-negative
        end
    end
end