@testset "Forcing and drag: 2D" begin
    # 2D models
    spectral_grid = SpectralGrid(trunc=31, nlayers=1)
    drag = JetDrag(spectral_grid, time_scale=Day(6))
    forcing = StochasticStirring(spectral_grid)
    random_process = SpectralAR1Process(spectral_grid)
    initial_conditions = StartFromRest()

    @testset for Model in (BarotropicModel, ShallowWaterModel)
        model = Model(spectral_grid; initial_conditions, forcing, drag, random_process)
        simulation = initialize!(model)

        run!(simulation, period=Day(15), output=true)
        @test simulation.model.feedback.nars_detected == false
    end
end

@testset "Forcing and drag: 3D" begin

    # 3D models
    spectral_grid = SpectralGrid(trunc=31, nlayers=8)
    drag = JetDrag(spectral_grid, time_scale=Day(6))
    forcing = StochasticStirring(spectral_grid)
    random_process = SpectralAR1Process(spectral_grid)
    initial_conditions = StartFromRest()

    @testset for Model in (PrimitiveDryModel, PrimitiveWetModel)
        model = Model(spectral_grid; initial_conditions, forcing, drag, random_process)
        simulation = initialize!(model)

        run!(simulation, period=Day(15), output=true)
        @test simulation.model.feedback.nars_detected == false
    end
end