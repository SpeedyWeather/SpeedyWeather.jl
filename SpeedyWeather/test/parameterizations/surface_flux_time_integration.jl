@testset "Surface flux time integration: damping factor" begin
    explicit, implicit, exponential, rational = ExplicitFluxes(), ImplicitFluxes(), ExponentialFluxes(), RationalFluxes()
    for NF in (Float32, Float64)
        # z → 0: all schemes converge to the explicit flux
        for scheme in (explicit, implicit, exponential, rational)
            @test SpeedyWeather.flux_damping(scheme, zero(NF)) === one(NF)
            @test SpeedyWeather.flux_damping(scheme, NF(1.0e-6)) ≈ 1 atol = 1.0e-5
        end

        z = NF.(10 .^ range(-3, 3, 20))
        φ_implicit = SpeedyWeather.flux_damping.(implicit, z)
        φ_exponential = SpeedyWeather.flux_damping.(exponential, z)
        @test all(SpeedyWeather.flux_damping.(explicit, z) .== 1)
        @test φ_implicit ≈ 1 ./ (1 .+ z)
        @test φ_exponential ≈ (1 .- exp.(-z)) ./ z

        # monotonically decreasing, in (0, 1], implicit damps more than exponential
        @test all(diff(φ_implicit) .< 0)
        @test all(diff(φ_exponential) .< 0)
        @test all(0 .< φ_implicit .<= φ_exponential .<= 1)

        # rational approximates exponential: 2nd order for small z, 1/z for large z
        φ_rational = SpeedyWeather.flux_damping.(rational, z)
        @test φ_rational ≈ (6 .+ z) ./ (6 .+ 4z .+ z .^ 2)
        @test all(diff(φ_rational) .< 0)
        @test all(abs.(φ_rational .- φ_exponential) .<= 0.11 .* φ_exponential)
        moderate = z .<= 5     # more accurate than implicit for moderate z = Δt/τ
        @test all((abs.(φ_rational .- φ_exponential) .<= abs.(φ_implicit .- φ_exponential))[moderate])
    end
    @test SpeedyWeather.flux_damping(rational, 0.01) ≈ SpeedyWeather.flux_damping(exponential, 0.01) rtol = 1.0e-6
    @test 1000 * SpeedyWeather.flux_damping(rational, 1000.0) ≈ 1 rtol = 0.01

    # single relaxation step x⁺ = x + z*φ(z)*(x_surface - x): implicit and exponential stay
    # between x and x_surface for any z, explicit overshoots for z > 1 and is unstable for z > 2
    x, x_surface = 1.0, 0.0
    for z in (0.5, 1.5, 5.0, 100.0)
        for scheme in (implicit, exponential)
            x_new = x + z * SpeedyWeather.flux_damping(scheme, z) * (x_surface - x)
            @test x_surface <= x_new <= x
        end
        # rational: stable, overshoots x_surface by at most ~10% for z > 3
        x_rational = x + z * SpeedyWeather.flux_damping(rational, z) * (x_surface - x)
        @test x_surface - 0.1 <= x_rational <= x
        x_explicit = x + z * (x_surface - x)
        z > 1 && @test x_explicit < x_surface
        z > 2 && @test abs(x_explicit - x_surface) > abs(x - x_surface)
    end
    # exponential is exact
    z = 3.0
    @test x + z * SpeedyWeather.flux_damping(exponential, z) * (x_surface - x) ≈ exp(-z)
end

@testset "Surface flux time integration: model" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
    model = PrimitiveWetModel(spectral_grid)

    # Defaults for all relaxation surface fluxes
    @test model.surface_momentum_flux.time_integration isa SpeedyWeather.DEFAULT_FLUX_INTEGRATION
    @test model.surface_heat_flux.ocean.time_integration isa SpeedyWeather.DEFAULT_FLUX_INTEGRATION
    @test model.surface_heat_flux.land.time_integration isa SpeedyWeather.DEFAULT_FLUX_INTEGRATION
    @test model.surface_humidity_flux.ocean.time_integration isa SpeedyWeather.DEFAULT_FLUX_INTEGRATION
    @test model.surface_humidity_flux.land.time_integration isa SpeedyWeather.DEFAULT_FLUX_INTEGRATION

    # option chosen per flux
    for scheme in (ExplicitFluxes(), ImplicitFluxes(), ExponentialFluxes())
        @test SurfaceMomentumFlux(spectral_grid; time_integration = scheme).time_integration === scheme
        @test SurfaceOceanHeatFlux(spectral_grid; time_integration = scheme).time_integration === scheme
        @test SurfaceLandHeatFlux(spectral_grid; time_integration = scheme).time_integration === scheme
        @test SurfaceOceanHumidityFlux(spectral_grid; time_integration = scheme).time_integration === scheme
        @test SurfaceLandHumidityFlux(spectral_grid; time_integration = scheme).time_integration === scheme
    end

    # damping uses the lowermost layer's pressure thickness and the time stepper's time step
    pₛ = 1.0f5
    conductance = 0.05f0
    Δp = SpeedyWeather.pressure_thickness(8, pₛ, model.geometry.vertical_coordinates)
    z = SpeedyWeather.default_time_step(model.time_stepping) * model.planet.gravity * conductance / Δp
    @test SpeedyWeather.surface_flux_damping(ImplicitFluxes(), conductance, pₛ, model) ≈ 1 / (1 + z)
    @test SpeedyWeather.surface_flux_damping(ExplicitFluxes(), conductance, pₛ, model) == 1
    flux = 2.0f0
    @test SpeedyWeather.surface_flux_to_tendency(ImplicitFluxes(), flux, conductance, pₛ, model) ≈
        SpeedyWeather.surface_flux_to_tendency(flux, pₛ, model) / (1 + z)
end

@testset "Surface flux time integration: Frierson sigma coordinates" begin
    # thin lowermost layer made explicit surface fluxes unstable at default Δt (#1335)
    # FriersonSigmaCoordinates are the default vertical coordinates
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
    @test Geometry(spectral_grid).vertical_coordinates.σ_half ≈ FriersonSigmaCoordinates(spectral_grid).σ_half
    for Model in (PrimitiveDryModel, PrimitiveWetModel)
        model = Model(spectral_grid)
        simulation = initialize!(model)
        run!(simulation, period = Day(5))
        @test model.feedback.nans_detected == false
    end
end
