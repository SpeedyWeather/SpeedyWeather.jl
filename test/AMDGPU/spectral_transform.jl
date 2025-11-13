@testset "AMDGPU spectral transform rountrip " begin
    spectral_grid = SpectralGrid(trunc=41, nlayers=8, architecture=GPU())
    S = SpectralTransform(spectral_grid)
    
    # first roundtrip
    L = randn(ComplexF32, spectral_grid.spectrum, 8)
    field = zeros(Float32, spectral_grid.grid, 8)
    transform!(field, L, S)
    transform!(L, field, S)

    # 2nd roundtrip
    L2 = deepcopy(L)
    transform!(field, L, S)
    transform!(L, field, S)

    @test L ≈ L2
end
