@testset "run PrimitiveWetModel defaults for a year, no blow up" begin
    spectral_grid = SpectralGrid()
    model = PrimitiveWetModel(spectral_grid)
    simulation = initialize!(model)
    run!(simulation, period = Year(1))
    @test simulation.model.feedback.nans_detected == false
end