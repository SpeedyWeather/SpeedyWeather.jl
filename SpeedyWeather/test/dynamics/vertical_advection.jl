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

@testset "Vertical advection face consistency" begin
    # the flux-form rewrite of the vertical advection kernel relies on the "+" face of
    # layer k (tail of its stencil) being bit-identical to the "-" face of layer k+1
    # (front of its stencil), for every k including the clamped boundaries, so that each
    # interior face can be reconstructed once and shared between neighbouring layers.
    NF = Float32
    advection_schemes = (
        SpeedyWeather.CenteredVerticalAdvection{NF, 1}(),
        SpeedyWeather.CenteredVerticalAdvection{NF, 2}(),
        SpeedyWeather.UpwindVerticalAdvection{NF, 1}(),
        SpeedyWeather.UpwindVerticalAdvection{NF, 2}(),
        SpeedyWeather.UpwindVerticalAdvection{NF, 3}(),
        SpeedyWeather.WENOVerticalAdvection{NF}(),
    )

    for adv in advection_schemes, nlayers in (2, 3, 5, 8, 24)
        for k in 1:(nlayers - 1)
            stencil_k = SpeedyWeather.retrieve_stencil(k, nlayers, adv)
            stencil_kp1 = SpeedyWeather.retrieve_stencil(k + 1, nlayers, adv)
            @test Base.tail(stencil_k) == Base.front(stencil_kp1)
        end
    end
end

@testset "Vertical advection runs" begin
    spectral_grid = SpectralGrid(trunc = 31, nlayers = 8)

    advection_schemes = (
        SpeedyWeather.WENOVerticalAdvection,
        SpeedyWeather.CenteredVerticalAdvection,
        SpeedyWeather.UpwindVerticalAdvection,
    )

    for VerticalAdvection in advection_schemes
        model = PrimitiveWetModel(
            spectral_grid;
            vertical_advection = VerticalAdvection(spectral_grid),
            dynamics_only = true
        )
        model.feedback.verbose = false
        simulation = initialize!(model)
        run!(simulation, period = Day(1))
        @test simulation.model.feedback.nans_detected == false
    end
end
