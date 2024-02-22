@testset "Land sea masks" begin
    @testset for Mask in (LandSeaMask, AquaPlanetMask)
        spectral_grid = SpectralGrid(trunc=31,nlev=8)
        mask = Mask(spectral_grid)
        model = PrimitiveWetModel(;spectral_grid,land_sea_mask=mask)
        simulation = initialize!(model)
        run!(simulation,period=Day(5))
        @test simulation.model.feedback.nars_detected == false
    end
end