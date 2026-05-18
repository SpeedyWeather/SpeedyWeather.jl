using Metal

@testset "Test Metal extension" begin

    spectral_grid = SpectralGrid(architecture = GPU())
    @test spectral_grid.architecture.device isa MetalBackend
    @test spectral_grid.VectorType <: MtlArray
    @test RingGrids.array_type(spectral_grid.GridVariable2D) <: MtlArray
    @test RingGrids.array_type(spectral_grid.SpectralVariable2D) <: MtlArray

    # allocate variables with Metal arrays
    model = PrimitiveWetModel(spectral_grid)
    vars = Variables(model)
    @test vars.prognostic.vorticity.data isa MtlArray
    @test vars.prognostic.ocean.sea_surface_temperature.data isa MtlArray
    @test vars.grid.vorticity.data isa MtlArray
end
