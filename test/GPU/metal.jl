using Metal

@testset "Test Metal extension" begin

    spectral_grid = SpectralGrid(architecture=GPU())
    @test spectral_grid.architecture.device isa MetalBackend
    @test spectral_grid.VectorType <: MtlArray
    @test RingGrids.array_type(spectral_grid.GridVariable2D) <: MtlArray
    @test RingGrids.array_type(spectral_grid.SpectralVariable2D) <: MtlArray

    # allocate variables with Metal arrays
    progn = PrognosticVariables(spectral_grid)
    @test progn.vor.data isa MtlArray
    @test progn.ocean.sea_surface_temperature.data isa MtlArray
    
    diagn = DiagnosticVariables(spectral_grid)
    @test diagn.grid.vor_grid.data isa MtlArray
end
