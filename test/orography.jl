@testset "Orographies" begin
    @testset for Orography in (EarthOrography, ZonalRidge)
        spectral_grid = SpectralGrid(trunc=31, nlev=8)
        orography = Orography(spectral_grid)
        model = PrimitiveWetModel(;spectral_grid, orography)
        simulation = initialize!(model)
        run!(simulation, period=Day(5))
        @test simulation.model.feedback.nars_detected == false
    end
end