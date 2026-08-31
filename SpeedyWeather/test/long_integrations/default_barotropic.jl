@testset "run BarotropicModel defaults for a year, no blow up" begin
    spectral_grid = SpectralGrid(nlayers = 1)
    model = BarotropicModel(spectral_grid)
    model.feedback.verbose = false
    simulation = initialize!(model)
    run!(simulation, period = Year(1))
    @test simulation.model.feedback.nans_detected == false
end
