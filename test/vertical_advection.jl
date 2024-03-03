@testset "Vertical advection runs" begin
    spectral_grid = SpectralGrid()
    model_types = (PrimitiveDryModel, PrimitiveWetModel)
    advection_schems = (SpeedyWeather.WENOVerticalAdvection,
                        SpeedyWeather.CenteredVerticalAdvection,
                        SpeedyWeather.UpwindVerticalAdvection)

    for Model in model_types, VerticalAdvection in advection_schems
        model = Model(; spectral_grid,
                        vertical_advection = VerticalAdvection(spectral_grid),
                        physics=false)
        simulation = initialize!(model)
        run!(simulation,period=Day(10))
        @test simulation.model.feedback.nars_detected == false    
    end
end