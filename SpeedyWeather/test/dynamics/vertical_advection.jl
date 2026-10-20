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
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)

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

@testset "Vertical advection of the reference temperature profile" begin
    # The semi-implicit reference profile Tₖ is subtracted from the grid temperature
    # for the explicit tendencies, but the vertical advection has to see the full
    # temperature T = T' + Tₖ as implicit_correction! expects the vertical advection
    # of Tₖ in the explicit tendency (issue #1285). Changing Tₖ may then only change
    # the grid temperature tendency through the T'D term, by exactly -ΔTₖ⋅D.
    spectral_grid = SpectralGrid(truncation = 31, nlayers = 8)

    advection_schemes = (
        SpeedyWeather.CenteredVerticalAdvection,
        SpeedyWeather.UpwindVerticalAdvection,
        SpeedyWeather.WENOVerticalAdvection,
    )

    for VerticalAdvection in advection_schemes
        model = PrimitiveDryModel(
            spectral_grid;
            vertical_advection = VerticalAdvection(spectral_grid),
            dynamics_only = true
        )
        model.feedback.verbose = false
        simulation = initialize!(model)
        run!(simulation, period = Hour(6))  # develop a non-zero vertical velocity
        (; variables) = simulation
        (; time_stepping) = model

        function grid_temperature_tendency()
            SpeedyWeather.reset_tendencies!(variables, time_stepping)
            SpeedyWeather.dynamics_tendencies!(variables, model)
            temperature_tendency = variables.tendencies.grid.temperature
            return copy(SpeedyWeather.get_tendency_step(temperature_tendency, time_stepping, SpeedyWeather.DynamicalCore()).data)
        end

        tendency = grid_temperature_tendency()

        # change the stratification, not just the mean, of the reference profile
        profile_change = 10 .* model.geometry.σ_levels_full
        model.implicit.temp_profile .+= profile_change
        tendency_changed_profile = grid_temperature_tendency()

        divergence = SpeedyWeather.get_prognostic_step(
            variables.grid.divergence, time_stepping, SpeedyWeather.DynamicalCore()
        ).data
        @test tendency_changed_profile ≈ tendency .- profile_change' .* divergence
    end
end
