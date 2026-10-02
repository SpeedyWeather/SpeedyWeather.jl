@testset "SpectralGrid" begin
    @testset for truncation in (32, 64, 128)
        spectral_grid = SpectralGrid(truncation = truncation, dealiasing = 2, Grid = OctahedralGaussianGrid)

        spectral_grid_2 = SpectralGrid(spectral_grid.grid; NF = spectral_grid.NF, dealiasing = 2)

        @test spectral_grid_2.truncation == spectral_grid.truncation
        @test spectral_grid_2.nlat_half == spectral_grid.grid.nlat_half
    end

    @testset "trunc kwarg deprecated in favour of truncation" begin
        spectral_grid_trunc = SpectralGrid(trunc = 63, dealiasing = 2, Grid = OctahedralGaussianGrid)
        spectral_grid_truncation = SpectralGrid(truncation = 64, dealiasing = 2, Grid = OctahedralGaussianGrid)

        @test spectral_grid_trunc.trunc == spectral_grid_truncation.trunc
        @test spectral_grid_trunc.truncation == spectral_grid_truncation.truncation
    end
end
