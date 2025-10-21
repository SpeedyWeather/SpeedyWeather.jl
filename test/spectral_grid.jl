@testset "SpectralGrid" begin
    @testset for trunc in (31, 63, 127)
        spectral_grid = SpectralGrid(trunc=trunc, dealiasing=2, Grid=OctahedralGaussianGrid)
        
        spectral_grid_2 = SpectralGrid(spectral_grid.grid, spectral_grid.NF, 2)

        @test spectral_grid_2.trunc == spectral_grid.trunc
        @test spectral_grid_2.nlat_half == spectral_grid.grid.nlat_half
    end
end