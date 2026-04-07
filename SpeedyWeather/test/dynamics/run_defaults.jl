@testset "run defaults for a year, no blow up" begin

    @testset "BarotropicModel" begin
        spectral_grid = SpectralGrid(nlayers = 1)
        model = BarotropicModel(spectral_grid)
        simulation = initialize!(model)
        run!(simulation, period = Year(1))
        @test simulation.model.feedback.nans_detected == false
    end

    @testset "ShallowWaterModel" begin
        spectral_grid = SpectralGrid(nlayers = 1)
        model = ShallowWaterModel(spectral_grid)
        simulation = initialize!(model)
        run!(simulation, period = Year(1))
        @test simulation.model.feedback.nans_detected == false
    end

    @testset "PrimitiveDryModel" begin
        spectral_grid = SpectralGrid()
        model = PrimitiveDryModel(spectral_grid)
        simulation = initialize!(model)
        run!(simulation, period = Year(1))
        @test simulation.model.feedback.nans_detected == false
    end

    @testset "PrimitiveWetModel" begin
        spectral_grid = SpectralGrid()
        model = PrimitiveWetModel(spectral_grid)
        simulation = initialize!(model)
        run!(simulation, period = Year(1))
        @test simulation.model.feedback.nans_detected == false
    end
end

@testset "run defaults with MatrixSpectralTransform, no blow up" begin
    # PrimitiveWet, with MatrixSpectralTransform
    spectral_grid = SpectralGrid()
    spectral_transform = MatrixSpectralTransform(spectral_grid)
    model = PrimitiveWetModel(spectral_grid; spectral_transform)
    simulation = initialize!(model)
    run!(simulation, period = Day(20))
    @test simulation.model.feedback.nans_detected == false
end
