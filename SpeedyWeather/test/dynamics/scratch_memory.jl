@testset "Scratch memory allocation" begin
    SG = SpectralGrid(trunc = 21, nlayers = 3)

    # created via transform or directly
    S = SpectralTransform(SG)
    SM = SpeedyTransforms.ScratchMemory(SG.NF, SG.architecture, SG.grid, 3)

    # created via model initialization
    model = PrimitiveDryModel(SG, spectral_transform = S)
    variables = Variables(model)

    # change the memory of the transform
    S.scratch_memory.north[1] = 1.23

    # same as in variables?
    @test variables.scratch.transform_memory == S.scratch_memory
end
