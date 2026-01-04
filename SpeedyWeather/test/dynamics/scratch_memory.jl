@testset "Scratch memory allocation" begin
    SG = SpectralGrid(trunc = 21, nlayers = 3)
    S = SpectralTransform(SG)
    SpeedyTransforms.ScratchMemory(SG.NF, SG.architecture, SG.grid, 3)
    DynamicsVariables(SG)
    DynamicsVariables(SG, spectral_transform = S)
end