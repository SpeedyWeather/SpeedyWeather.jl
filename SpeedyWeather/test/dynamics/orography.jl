@testset "Orographies" begin
    @testset for Orography in (EarthOrography, ZonalRidge)
        spectral_grid = SpectralGrid(trunc = 31, nlayers = 8)
        orography = Orography(spectral_grid)
        model = PrimitiveWetModel(spectral_grid; orography)
        initialize!(model.orography, model)

        @test any(model.orography.orography .!= 0)
        @test any(model.orography.surface_geopotential .!= 0)
    end
end
