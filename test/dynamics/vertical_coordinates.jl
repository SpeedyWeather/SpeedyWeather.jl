@testset "Initialize sigma layers manually" begin

    # constructors
    SigmaCoordinates()
    SigmaCoordinates(3)
    SigmaCoordinates([0, 0.4, 0.6, 1])
    SigmaCoordinates(0:0.25:1)
    SigmaCoordinates(3, [0, 0.4, 0.6, 1])
    SigmaCoordinates{Float32, Vector{Float32}}(4, 0:0.25:1)

    # automatic levels
    spectral_grid = SpectralGrid(nlayers=4)
    G = Geometry(spectral_grid)
    @test length(G.σ_levels_half) == 5
    @test length(G.σ_levels_full) == 4

    # manual levels
    σ = SigmaCoordinates([0, 0.4, 0.6, 1])
    spectral_grid = SpectralGrid(nlayers=σ.nlayers)
    G = Geometry(spectral_grid, vertical_coordinates=σ)
    @test spectral_grid.nlayers == 3
    @test length(G.σ_levels_half) == 4
    @test length(G.σ_levels_full) == 3

    # specify both
    σ = SigmaCoordinates([0, 0.4, 0.6, 1])
    spectral_grid = SpectralGrid(nlayers=3)
    G = Geometry(spectral_grid, vertical_coordinates=σ)
    @test spectral_grid.nlayers == 3
    @test length(G.σ_levels_half) == 4
    @test length(G.σ_levels_full) == 3
end