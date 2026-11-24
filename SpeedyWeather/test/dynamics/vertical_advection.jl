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

@testset "Vertical advection order of accuracy" begin
    # reconstruct sin(2πx) at the faces from its cell averages on N uniform cells and
    # check that the error converges at the nominal order of each scheme, for both signs of u
    cell_averages(N) = reshape([(cos(2π * (j - 1) / N) - cos(2π * j / N)) / (2π / N) for j in 1:N], 1, N, 1)

    function face_error(scheme, N, u)
        ξ = cell_averages(N)
        return maximum((N ÷ 4):(3N ÷ 4)) do k    # interior faces k+1/2, away from the clamped boundaries
            stencil = Base.tail(SpeedyWeather.retrieve_stencil(k, N, scheme))
            abs(SpeedyWeather.reconstructed_at_face(ξ, 1, 1, stencil, u, scheme) - sin(2π * k / N))
        end
    end

    for (scheme, order) in (
            (SpeedyWeather.UpwindVerticalAdvection{Float64, 1}(), 1),
            (SpeedyWeather.UpwindVerticalAdvection{Float64, 2}(), 3),
            (SpeedyWeather.UpwindVerticalAdvection{Float64, 3}(), 5),
            (SpeedyWeather.CenteredVerticalAdvection{Float64, 1}(), 2),
            (SpeedyWeather.CenteredVerticalAdvection{Float64, 2}(), 4),
            (SpeedyWeather.WENOVerticalAdvection{Float64}(), 5),
        ), u in (1.0, -1.0)
        observed_order = log2(face_error(scheme, 32, u) / face_error(scheme, 64, u))
        @test observed_order ≈ order atol = 0.2
    end
end

@testset "Vertical advection runs" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)

    advection_schemes = (
        SpeedyWeather.WENOVerticalAdvection,
        SpeedyWeather.CenteredVerticalAdvection,
        SpeedyWeather.UpwindVerticalAdvection,
        SpeedyWeather.SubcycledVerticalAdvection,
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
        SpeedyWeather.SubcycledVerticalAdvection,
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

@testset "Subcycled vertical advection" begin
    # a single column: compare with WENO for small Courant numbers, stay bounded for large ones
    nlayers = 40
    Δσ = fill(1 / nlayers, nlayers)
    ξ = zeros(1, nlayers, 2)
    ξ[1, :, 1] .= [exp(-((k - 15) / 4)^2) for k in 1:nlayers]     # Gaussian bump, step 1 is advected
    function tendency(adv, w, interval)
        ξ_tend = zeros(1, nlayers, 2)
        if adv isa SpeedyWeather.SubcycledVerticalAdvection
            adv.interval[] = interval
            S, A, B = adv.S, adv.A, adv.B
            SpeedyWeather.launch!(
                SpeedyWeather.CPU(), SpeedyWeather.LinearWorkOrder, (1,), SpeedyWeather.subcycled_vertical_advection_kernel!,
                ξ_tend, 1, w, ξ, 1, Δσ, nlayers, S, A, B, interval, adv.courant_target, adv.max_substeps,
                SpeedyWeather.WENOVerticalAdvection{Float64}(),
            )
        else
            SpeedyWeather.launch!(
                SpeedyWeather.CPU(), SpeedyWeather.RingGridWorkOrder, (1, nlayers), SpeedyWeather.vertical_advection_kernel!,
                ξ_tend, 1, w, ξ, 1, Δσ, nlayers, adv,
            )
        end
        return ξ_tend[1, :, 1]
    end

    spectral_grid = SpectralGrid(truncation = 5; nlayers, NF = Float64)
    @test SubcycledVerticalAdvection(spectral_grid) isa SpeedyWeather.SubcycledVerticalAdvection{Float64}
    sub = SpeedyWeather.SubcycledVerticalAdvection{Float64, Matrix{Float64}, Base.RefValue{Float64}}(
        1.0, 64, zeros(1, nlayers), zeros(1, nlayers), zeros(1, nlayers), Ref(0.0))
    weno = SpeedyWeather.WENOVerticalAdvection{Float64}()

    w = zeros(1, nlayers)
    w[1, 1:(nlayers - 1)] .= 1.0e-4                     # uniform descent (σ increasing), w = 0 at the surface
    # tiny interval: one substep, SSP-RK3 converges to the instantaneous WENO tendency
    @test tendency(sub, w, 1.0e-3) ≈ tendency(weno, w, 1.0e-3) rtol = 1.0e-3

    # Courant number 10 over the interval (bump moves ~10 layers, away from the boundaries):
    # subcycled stays finite and bounded to within 1 % (WENO is not strictly monotone)
    interval = 10 * Δσ[1] / 1.0e-4
    ξ_end = ξ[1, :, 1] .+ interval .* tendency(sub, w, interval)
    @test all(isfinite, ξ_end)
    @test maximum(ξ_end) <= maximum(ξ[1, :, 1]) + 1.0e-2
    @test minimum(ξ_end) >= -1.0e-2
    @test argmax(ξ_end) - argmax(ξ[1, :, 1]) in 9:11   # the bump moved down ~10 layers (towards larger σ)
end
