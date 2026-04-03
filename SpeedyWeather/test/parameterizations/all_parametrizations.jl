@testset "Test all parametrizations with boundschecks" begin
    model_classes = (PrimitiveWetModel, PrimitiveDryModel)

    spectral_grid = SpectralGrid()
    for model_class in model_classes

        model = model_class(spectral_grid)

        simulation = initialize!(model)
        run!(simulation, steps = 4)
        vars, model = SpeedyWeather.unpack(simulation)

        # call parameterization without @inbounds to check for boundserrors
        SpeedyWeather._column_parameterizations_cpu!(vars, SpeedyWeather.get_parameterizations(model), model)

        for key in keys(vars.parameterizations)
            if key != :land && key != :ocean
                @test !any(isnan.(vars.parameterizations[key]))
            end
        end
    end
end
