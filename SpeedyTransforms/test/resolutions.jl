@testset "Create transform at uncommon resolutions" begin
    @testset for trunc in rand(15:200, 10)
        spectrum = Spectrum(trunc)
        grid = FullClenshawGrid(SpeedyTransforms.get_nlat_half(trunc))
        S = SpectralTransform(spectrum, grid)
    end
end
