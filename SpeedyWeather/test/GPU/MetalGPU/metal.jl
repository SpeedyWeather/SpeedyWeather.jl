using Metal

@testset "Test Metal extension" begin

    spectral_grid = SpectralGrid(architecture = GPU())
    @test spectral_grid.architecture.device isa MetalBackend
    @test spectral_grid.VectorType <: MtlArray
    @test RingGrids.array_type(spectral_grid.GridVariable2D) <: MtlArray
    @test RingGrids.array_type(spectral_grid.SpectralVariable2D) <: MtlArray

    # allocate variables with Metal arrays
    # `.data` may be a `SubArray` into a fused parent buffer (see Variables system), so
    # check architecture rather than array type directly, which correctly unwraps views
    model = PrimitiveWetModel(spectral_grid)
    vars = Variables(model)
    @test architecture(vars.prognostic.vorticity.data) isa SpeedyWeather.GPU
    @test architecture(vars.prognostic.ocean.sea_surface_temperature.data) isa SpeedyWeather.GPU
    @test architecture(vars.grid.vorticity.data) isa SpeedyWeather.GPU
end
