@testset "run_speedy no errors, no blowup" begin
    # Barotropic
    spectral_grid = SpectralGrid(Barotropic)
    model = Model(;spectral_grid)
    simulation = initialize!(model)
    run!(simulation,n_days=10)
    @test simulation.model.feedback.nars_detected == false
    
    # ShallowWater
    spectral_grid = SpectralGrid(ShallowWater)
    model = Model(;spectral_grid)
    simulation = initialize!(model)
    run!(simulation,n_days=10)
    @test simulation.model.feedback.nars_detected == false

    # PrimitiveDry
    spectral_grid = SpectralGrid(PrimitiveDry)
    model = Model(;spectral_grid)
    simulation = initialize!(model)
    run!(simulation,n_days=10)
    @test simulation.model.feedback.nars_detected == false

    # PrimitiveWet
    spectral_grid = SpectralGrid(PrimitiveWet)
    model = Model(;spectral_grid)
    simulation = initialize!(model)
    run!(simulation,n_days=10)
    @test simulation.model.feedback.nars_detected == false
end