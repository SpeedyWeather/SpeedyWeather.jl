using SpeedyWeather: WENOVerticalAdvection, CenteredVerticalAdvection, UpwindVerticalAdvection
using SpeedyWeather: vertical_advection!

@testset "Vertical advection runs" begin
    spectral_grid = SpectralGrid()
    model_types = (PrimitiveDryModel, PrimitiveWetModel)
    advection_schems = (WENOVerticalAdvection, CenteredVerticalAdvection, UpwindVerticalAdvection)

    for Model in model_types, VerticalAdvection in advection_schems
        model = Model(; spectral_grid, vertical_advection = VerticalAdvection(spectral_grid))
        simulation = initialize!(model)
        run!(simulation,n_days=10)
        @test simulation.model.feedback.nars_detected == false    
    end
end