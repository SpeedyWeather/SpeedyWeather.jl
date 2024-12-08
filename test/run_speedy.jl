@testset "run_speedy no errors, no blowup" begin
    # Barotropic
    spectral_grid = SpectralGrid(nlayers=1)
    model = BarotropicModel(spectral_grid)
    simulation = initialize!(model)
    run!(simulation, period=Day(10))
    @test simulation.model.feedback.nars_detected == false
    
    # ShallowWater
    spectral_grid = SpectralGrid(nlayers=1)
    model = ShallowWaterModel(spectral_grid)
    simulation = initialize!(model)
    run!(simulation, period=Day(10))
    @test simulation.model.feedback.nars_detected == false

    # PrimitiveDry
    spectral_grid = SpectralGrid()
    model = PrimitiveDryModel(spectral_grid)
    simulation = initialize!(model)
    run!(simulation, period=Day(10))
    @test simulation.model.feedback.nars_detected == false

    # PrimitiveWet
    spectral_grid = SpectralGrid()
    model = PrimitiveWetModel(spectral_grid)
    simulation = initialize!(model)
    run!(simulation, period=Day(10))
    @test simulation.model.feedback.nars_detected == false
end