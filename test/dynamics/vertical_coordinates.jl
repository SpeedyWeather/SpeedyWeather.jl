@testset "Initialize sigma layers manually" begin
    # automatic levels
    spectral_grid = SpectralGrid(nlayers=4)
    G = Geometry(spectral_grid)
    @test length(G.σ_levels_half) == 5
    @test length(G.σ_levels_full) == 4

    # manual levels
    σ = SigmaCoordinates([0, 0.4, 0.6, 1])
    spectral_grid = SpectralGrid(nlayers=σ.nlayers, vertical_coordinates=σ)
    G = Geometry(spectral_grid)
    @test spectral_grid.nlayers == 3
    @test length(G.σ_levels_half) == 4
    @test length(G.σ_levels_full) == 3

    # specify both
    σ = SigmaCoordinates([0, 0.4, 0.6, 1])
    spectral_grid = SpectralGrid(nlayers=3, vertical_coordinates=σ)
    G = Geometry(spectral_grid)
    @test spectral_grid.nlayers == 3
    @test length(G.σ_levels_half) == 4
    @test length(G.σ_levels_full) == 3
end