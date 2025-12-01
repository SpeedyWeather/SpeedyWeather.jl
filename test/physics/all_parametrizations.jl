@testset "Test all parametrizations with boundschecks" begin 
    model_classes = (PrimitiveWetModel, PrimitiveDryModel)

    spectral_grid = SpectralGrid()
    for model_class in model_classes
        
        model = model_class(spectral_grid)

        simulation = initialize!(model)
        run!(simulation, steps=4)
        
        # call parameterization without @inbounds to check for boundserrors 
        _call_parameterizations_cpu!(diagn, progn, get_parameterizations(model), model)
        
        for key in keys(diagn.physics)
            @test !any(isnan.(diagn.physics[key]))
        end
    end
end