@testset "Large-scale condensation" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
    @testset for r in (0.9, 0.95)
        @testset for snow in (true, false)
            @testset for reevaporation in (0, 10, 30)
                large_scale_condensation = ImplicitCondensation(
                    spectral_grid;
                    relative_humidity_threshold = r, snow, reevaporation
                )
                model = PrimitiveWetModel(spectral_grid; large_scale_condensation)
                model.feedback.verbose = false
                simulation = initialize!(model)
                run!(simulation, period = Day(1))

                rain_fall = simulation.variables.parameterizations.rain_large_scale
                snow_fall = simulation.variables.parameterizations.snow_large_scale

                @test model.feedback.nans_detected == false
                @test all(rain_fall .>= 0)      # precipitation should always be non-negative

                if snow
                    @test all(snow_fall .>= 0)  # snow should always be non-negative
                    @test any(snow_fall .> 0)   # should have some snow
                    @test any(snow_fall .== 0)  # should have some areas without snow
                else
                    @test all(snow_fall .== 0)  # no snow should be produced
                end
            end
        end
    end
end

@testset "Large-scale condensation energy budget" begin
    # prescribed column: snow forms in cold supersaturated layers aloft and melts in warm layers
    # below, the column enthalpy change has to match latent heat of net condensation plus
    # latent heat of fusion for snow reaching the surface, see #1276
    spectral_grid = SpectralGrid(truncation = 15, nlayers = 8, NF = Float64)
    model = PrimitiveWetModel(spectral_grid)
    model.feedback.verbose = false
    simulation = initialize!(model)

    (; geometry, planet, atmosphere, time_stepping) = model
    condensation = model.large_scale_condensation
    vars = simulation.variables
    temp = SpeedyWeather.get_prognostic_step(vars.grid.temperature, time_stepping, condensation)
    humid = SpeedyWeather.get_prognostic_step(vars.grid.humidity, time_stepping, condensation)
    temp_tend = SpeedyWeather.get_tendency_step(vars.tendencies.grid.temperature, time_stepping, condensation)
    humid_tend = SpeedyWeather.get_tendency_step(vars.tendencies.grid.humidity, time_stepping, condensation)

    ij = 1
    nlayers = spectral_grid.nlayers
    surface_pressure = 1.0e5
    vars.parameterizations.surface_pressure[ij] = surface_pressure
    coord = geometry.vertical_coordinates
    for k in 1:nlayers
        temp[ij, k] = k <= 4 ? 250 : 285       # snow forms aloft, melts below
        pₖ = SpeedyWeather.pressure(k, surface_pressure, coord)
        sat_humid = SpeedyWeather.saturation_humidity(temp[ij, k], pₖ, atmosphere)
        humid[ij, k] = (k <= 4 ? 1.2 : 0.8) * sat_humid
    end

    temp_tend[ij, :] .= 0
    humid_tend[ij, :] .= 0
    vars.parameterizations.rain_rate[ij] = 0
    vars.parameterizations.snow_rate[ij] = 0
    SpeedyWeather.large_scale_condensation!(ij, vars, condensation, geometry, planet, atmosphere, time_stepping)

    g = planet.gravity
    ρ = atmosphere.water_density
    cₚ = atmosphere.heat_capacity
    Lᵥ = atmosphere.latent_heat_condensation
    Lᵢ = atmosphere.latent_heat_fusion
    Δp = [SpeedyWeather.pressure_thickness(k, surface_pressure, coord) for k in 1:nlayers]

    snow_surface = ρ * vars.parameterizations.snow_rate_large_scale[ij]    # [kg/m²/s]
    @test snow_surface >= 0
    @test vars.parameterizations.rain_rate_large_scale[ij] > 0      # rain only from melted snow

    heating = cₚ / g * sum(temp_tend[ij, k] * Δp[k] for k in 1:nlayers)          # [W/m²]
    net_condensation = -sum(humid_tend[ij, k] * Δp[k] for k in 1:nlayers) / g     # [kg/m²/s]
    @test heating ≈ Lᵥ * net_condensation + Lᵢ * snow_surface rtol = 1.0e-10

    # implicit correction in the top layer (no incoming precipitation), Frierson et al. 2006 eq. (21)
    # with Clausius-Clapeyron dqsat/dT = qsat Lᵥ / (Rᵥ T²) [1/K]
    (; relative_humidity_threshold, time_scale) = condensation
    p₁ = SpeedyWeather.pressure(1, surface_pressure, coord)
    sat_humid = SpeedyWeather.saturation_humidity(temp[ij, 1], p₁, atmosphere)
    dqsat_dT = relative_humidity_threshold * sat_humid * Lᵥ / (atmosphere.R_vapor * temp[ij, 1]^2)
    δq_cond = relative_humidity_threshold * sat_humid - humid[ij, 1]
    @test humid_tend[ij, 1] ≈ δq_cond / ((1 + Lᵥ / cₚ * dqsat_dT) * time_scale * time_stepping.Δt)
end
