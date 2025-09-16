using Metal

@testset "Test Metal extension" begin

    spectral_grid = SpectralGrid(architecture=GPU())
    @test spectral_grid.architecture isa MetalGPU


end
