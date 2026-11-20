@testset "ClipNegatives stores the fill" begin
    spectral_grid = SpectralGrid(truncation = 16, nlayers = 4)
    humidity = rand(Float32, spectral_grid.grid, spectral_grid.nlayers) .- 0.5f0
    humidity_clipped = copy(humidity)
    humidity_fill = zero(humidity)
    SpeedyWeather.clip_negatives!(humidity_clipped, humidity_fill)

    @test all(humidity_clipped .>= 0)
    @test all(humidity_fill .>= 0)
    @test humidity_clipped .- humidity_fill == humidity     # fill = clipped - unclipped, exactly
    @test humidity_clipped == max.(humidity, 0)

    # clipping a field without storing the fill
    humidity_clipped2 = copy(humidity)
    SpeedyWeather.hole_filling!(humidity_clipped2, ClipNegatives())
    @test humidity_clipped2 == humidity_clipped
end

@testset "Hole filling tendency fills the spectral humidity" begin
    # constant negative humidity: clipped to zero on the grid, so that all humidity tendencies
    # but the hole filling one are zero (no physics, constant humidity is not advected).
    # Without the hole filling tendency the spectral humidity would stay at humidity0
    spectral_grid = SpectralGrid(truncation = 16, nlayers = 4)
    humidity0 = -1.0f-5

    model = PrimitiveWetModel(spectral_grid; dynamics_only = true)
    model.feedback.verbose = false
    simulation = initialize!(model)
    set!(simulation, humidity = humidity0)
    run!(simulation, period = Day(1))

    S = model.spectral_transform
    humidity_spectral = SpeedyWeather.get_step(simulation.variables.prognostic.humidity, 2)
    humidity_unclipped = transform(humidity_spectral, S)
    global_mean = real(humidity_spectral[1, :]) ./ real(transform(ones(Float32, spectral_grid.grid), S)[1])

    # the fill is added to the spectral humidity, so that its grid humidity is no longer negative.
    # Not exactly: the Robert-Asselin-Williams filter sees the unfilled previous step and moves
    # ~2.5% of the fill back, refilled two steps later, which leaves ~1% of the hole as water
    @test all(abs.(global_mean) .< 0.05 * abs(humidity0))
    @test all(humidity_unclipped .> -0.05 * abs(humidity0))
    @test all(simulation.variables.grid.humidity .>= 0)
end

@testset "Hole filling: grid humidity and fill consistent with spectral humidity" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
    model = PrimitiveWetModel(spectral_grid)
    model.feedback.verbose = false
    simulation = initialize!(model)
    run!(simulation, period = Day(1))

    (; prognostic, grid, dynamics) = simulation.variables
    humidity_unclipped = transform(SpeedyWeather.get_step(prognostic.humidity, 2), model.spectral_transform)
    @test SpeedyWeather.get_step(grid.humidity, 2) == max.(humidity_unclipped, 0)
    @test SpeedyWeather.get_step(dynamics.humidity_hole_fill, 2) == max.(.-humidity_unclipped, 0)
    @test all(dynamics.humidity_hole_fill .>= 0)
end
