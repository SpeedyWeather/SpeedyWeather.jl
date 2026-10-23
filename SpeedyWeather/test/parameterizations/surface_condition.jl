@testset "LearnedWindSlowdown roughness option" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)

    for use_computed_roughness in (false, true)
        wind_slowdown = LearnedWindSlowdown(spectral_grid; use_computed_roughness)
        surface_condition = SurfaceCondition(spectral_grid; wind_slowdown)
        boundary_layer = BoundaryLayer(spectral_grid; surface_condition)
        model = PrimitiveWetModel(spectral_grid; boundary_layer)
        simulation = initialize!(model)
        run!(simulation, steps = 2)

        # recompute the surface wind with different roughness left over from a "previous" time step
        (; variables) = simulation
        (; parameterizations) = variables
        fill!(parameterizations.momentum_roughness, 1.0e-3)
        SpeedyWeather.parameterization_tendencies!(variables, model)
        surface_wind_speed = copy(parameterizations.surface_wind_speed)
        fill!(parameterizations.momentum_roughness, 0.5)
        SpeedyWeather.parameterization_tendencies!(variables, model)

        @test all(isfinite, parameterizations.surface_wind_speed)

        # fixed roughness by default: no dependency on the previous time step
        # computed roughness: surface wind depends on the previous time step's roughness
        independent = surface_wind_speed == parameterizations.surface_wind_speed
        @test independent == !use_computed_roughness
    end
end
