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
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
    model = PrimitiveWetModel(spectral_grid)
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
