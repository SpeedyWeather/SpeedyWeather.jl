@testset "Vertical advection stencils" begin
    @test (1, 1, 2) == SpeedyWeather.SpeedyWeather.retrieve_stencil(1, 8, SpeedyWeather.CenteredVerticalAdvection{Float32, 1}())
    @test (1, 2, 3) == SpeedyWeather.retrieve_stencil(2, 8, SpeedyWeather.CenteredVerticalAdvection{Float32, 1}())
    @test (2, 3, 4) == SpeedyWeather.retrieve_stencil(3, 8, SpeedyWeather.CenteredVerticalAdvection{Float32, 1}())
    @test (7, 8, 8) == SpeedyWeather.retrieve_stencil(8, 8, SpeedyWeather.CenteredVerticalAdvection{Float32, 1}())

    @test (1, 1, 2) == SpeedyWeather.retrieve_stencil(1, 5, SpeedyWeather.CenteredVerticalAdvection{Float32, 1}())
    @test (1, 2, 3) == SpeedyWeather.retrieve_stencil(2, 5, SpeedyWeather.CenteredVerticalAdvection{Float32, 1}())
    @test (2, 3, 4) == SpeedyWeather.retrieve_stencil(3, 5, SpeedyWeather.CenteredVerticalAdvection{Float32, 1}())
    @test (4, 5, 5) == SpeedyWeather.retrieve_stencil(5, 5, SpeedyWeather.CenteredVerticalAdvection{Float32, 1}())

    @test (1, 1, 1, 2, 3) == SpeedyWeather.retrieve_stencil(1, 8, SpeedyWeather.CenteredVerticalAdvection{Float32, 2}())
    @test (1, 1, 2, 3, 4) == SpeedyWeather.retrieve_stencil(2, 8, SpeedyWeather.CenteredVerticalAdvection{Float32, 2}())
    @test (1, 2, 3, 4, 5) == SpeedyWeather.retrieve_stencil(3, 8, SpeedyWeather.CenteredVerticalAdvection{Float32, 2}())
    @test (6, 7, 8, 8, 8) == SpeedyWeather.retrieve_stencil(8, 8, SpeedyWeather.CenteredVerticalAdvection{Float32, 2}())
end

@testset "Vertical advection runs" begin
    spectral_grid = SpectralGrid(trunc=31, nlayers=8)

    advection_schemes = (SpeedyWeather.WENOVerticalAdvection,
                        SpeedyWeather.CenteredVerticalAdvection,
                        SpeedyWeather.UpwindVerticalAdvection)

    for VerticalAdvection in advection_schemes
        model = PrimitiveWetModel(spectral_grid;
                        vertical_advection = VerticalAdvection(spectral_grid),
                        physics=false)
        simulation = initialize!(model)
        run!(simulation, period=Day(1))
        @test simulation.model.feedback.nars_detected == false    
    end
end