@testset "Single time step at uncommon resolutions" begin
    @testset for trunc in rand(15:200, 10)
        spectral_grid = SpectralGrid(; trunc, Grid=FullClenshawGrid)
        model = PrimitiveWetModel(spectral_grid)
        simulation = initialize!(model)
        run!(simulation, period=Day(0))
        @test simulation.model.feedback.nars_detected == false
    end
end