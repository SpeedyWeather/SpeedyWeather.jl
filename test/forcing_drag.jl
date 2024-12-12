@testset "Forcing and drag" begin
    spectral_grid = SpectralGrid(trunc=31, nlayers=1)
    drag = JetDrag(spectral_grid, time_scale=Day(6))
    forcing = StochasticStirring(spectral_grid)
    initial_conditions = StartFromRest()

    for Model in (BarotropicModel, ShallowWaterModel, PrimitiveDryModel, PrimitiveWetModel)
        model = Model(spectral_grid; initial_conditions, forcing, drag)
        simulation = initialize!(model)

        run!(simulation, period=Day(5))
        @test simulation.model.feedback.nars_detected == false
    end
end