@testset "Orographies" begin
    @testset for Orography in (EarthOrography, ZonalRidge)
        spectral_grid = SpectralGrid(trunc=31, nlayers=8)
        orography = Orography(spectral_grid)
        model = PrimitiveWetModel(spectral_grid; orography)
        simulation = initialize!(model)
        run!(simulation, period=Day(5))
        @test simulation.model.feedback.nans_detected == false
    end
end