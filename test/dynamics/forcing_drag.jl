@testset "Forcing and drag: 2D" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_testruns_")  # Cleaned up when the process exits

    # 2D models
    spectral_grid = SpectralGrid(trunc=31, nlayers=1)
    output = NetCDFOutput(spectral_grid, path=tmp_output_path)
    drag = JetDrag(spectral_grid, time_scale=Day(6))
    forcing = StochasticStirring(spectral_grid)
    random_process = SpectralAR1Process(spectral_grid)
    initial_conditions = StartFromRest()

    @testset for Model in (BarotropicModel, ShallowWaterModel)
        model = Model(spectral_grid; initial_conditions, forcing, drag, random_process, output)
        simulation = initialize!(model)

        run!(simulation, period=Day(15), output=true)
        @test simulation.model.feedback.nars_detected == false
    end
end

@testset "Forcing and drag: 3D" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_testruns_")  # Cleaned up when the process exits

    # 3D models
    spectral_grid = SpectralGrid(trunc=31, nlayers=8)
    output = NetCDFOutput(spectral_grid, path=tmp_output_path)
    drag = JetDrag(spectral_grid, time_scale=Day(6))
    forcing = StochasticStirring(spectral_grid)
    random_process = SpectralAR1Process(spectral_grid)
    initial_conditions = StartFromRest()

    @testset for Model in (PrimitiveDryModel, PrimitiveWetModel)
        model = Model(spectral_grid; initial_conditions, forcing, drag, random_process, output)
        simulation = initialize!(model)

        run!(simulation, period=Day(15), output=true)
        @test simulation.model.feedback.nars_detected == false
    end
end