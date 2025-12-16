@testset "Test all parametrizations with boundschecks" begin 
    model_classes = (PrimitiveWetModel, PrimitiveDryModel)

    spectral_grid = SpectralGrid()
    for model_class in model_classes
        
        model = model_class(spectral_grid)

        simulation = initialize!(model)
        run!(simulation, steps=4)
        progn, diagn, model = SpeedyWeather.unpack(simulation)
        
        # call parameterization without @inbounds to check for boundserrors 
        SpeedyWeather._call_parameterizations_cpu!(diagn, progn, SpeedyWeather.get_parameterizations(model), model)
        
        for key in keys(diagn.physics)
            if key != :land && key != :ocean
                @test !any(isnan.(diagn.physics[key]))
            end
        end
    end
end