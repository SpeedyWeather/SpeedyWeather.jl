@testset "Land sea masks" begin
    @testset for Mask in (LandSeaMask, AquaPlanetMask, RockyPlanetMask)
        spectral_grid = SpectralGrid(trunc = 31, nlayers = 8)
        mask = Mask(spectral_grid)
        model = PrimitiveWetModel(spectral_grid; land_sea_mask = mask)
        simulation = initialize!(model; feedback = Feedback(verbose = false))
        run!(simulation, period = Day(5))
        @test simulation.feedback.nans_detected == false
    end
end
