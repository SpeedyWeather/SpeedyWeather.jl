@testset "Neutral transfer coefficient" begin
    κ = 0.4f0
    z = 30.0f0      # height of the lowermost layer [m]

    # present surfaces: κ²/(log(z/z₀M) log(z/z₀X)), capped one order of magnitude above the roughness length
    @test SpeedyWeather.neutral_transfer_coefficient(z, 1.0f-3, 1.0f-4, κ) ≈ κ^2 / (log(z / 1.0f-3) * log(z / 1.0f-4))
    @test SpeedyWeather.neutral_transfer_coefficient(z, 10.0f0, 10.0f0, κ) ≈ κ^2 / log(10.0f0)^2

    # absent surface (zero roughness): finite instead of log(z/0) = Inf
    @test isfinite(SpeedyWeather.neutral_transfer_coefficient(z, 0.0f0, 0.0f0, κ))
    @test isfinite(SpeedyWeather.neutral_transfer_coefficient(z, 1.0f-3, 0.0f0, κ))
end

@testset "BulkRichardsonDrag at land- and ocean-only points" begin
    @testset for Model in (PrimitiveWetModel, PrimitiveDryModel)
        spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
        model = Model(spectral_grid)
        simulation = initialize!(model)
        run!(simulation, steps = 3)

        (; land, ocean) = simulation.variables.parameterizations
        (; drag_min) = model.boundary_layer.drag
        for drag in (
                land.boundary_layer_drag_momentum, land.boundary_layer_drag_heat, land.boundary_layer_drag_humidity,
                ocean.boundary_layer_drag_momentum, ocean.boundary_layer_drag_heat, ocean.boundary_layer_drag_humidity,
            )
            @test all(isfinite, drag)
            @test all(>=(drag_min), drag)
        end
    end
end

@testset "BulkRichardsonDrag stability from the surface temperature" begin
    @testset for Model in (PrimitiveWetModel, PrimitiveDryModel)
        spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
        model = Model(spectral_grid)
        simulation = initialize!(model)
        run!(simulation, steps = 3)

        vars = simulation.variables
        (; land, ocean) = vars.parameterizations
        (; drag) = model.boundary_layer
        (; drag_min, von_Karman) = drag

        # lowermost layer temperature as read by the drag, and the neutral heat transfer coefficients
        air_temperature = SpeedyWeather.get_prognostic_step(vars.grid.temperature, model.time_stepping, drag).data[:, end]
        z = vars.dynamics.geopotential.data[:, end] ./ model.planet.gravity .- model.orography.orography.data
        C_land = SpeedyWeather.neutral_transfer_coefficient.(z, land.momentum_roughness.data, land.heat_roughness.data, von_Karman)
        C_ocean = SpeedyWeather.neutral_transfer_coefficient.(z, ocean.momentum_roughness.data, ocean.heat_roughness.data, von_Karman)

        # set soil and sea surface temperature to the air temperature + ΔT and recompute the drag
        function heat_drag(ΔT, surface_wind_speed)
            vars.prognostic.land.soil_temperature.data .= air_temperature .+ ΔT
            vars.prognostic.ocean.sea_surface_temperature.data .= air_temperature .+ ΔT
            vars.parameterizations.surface_wind_speed.data .= surface_wind_speed
            for ij in eachindex(air_temperature)
                SpeedyWeather.boundary_layer_drag!(ij, vars, drag, model)
            end
            return copy(land.boundary_layer_drag_heat.data), copy(ocean.boundary_layer_drag_heat.data)
        end

        # surface warmer than the air: unstable (Ri < 0), neutral transfer coefficient
        land_drag, ocean_drag = heat_drag(20, 5)
        @test land_drag ≈ C_land
        @test ocean_drag ≈ C_ocean

        # surface colder than the air: stable (0 < Ri < Ri_c), reduced transfer coefficient
        land_drag, ocean_drag = heat_drag(-10, 10)
        @test all(drag_min .< land_drag .< C_land)
        @test all(drag_min .< ocean_drag .< C_ocean)

        # surface much colder than the air at weak wind: strongly stable (Ri ≥ Ri_c), minimum drag
        land_drag, ocean_drag = heat_drag(-20, 0.5)
        @test all(==(drag_min), land_drag)
        @test all(==(drag_min), ocean_drag)
    end
end
