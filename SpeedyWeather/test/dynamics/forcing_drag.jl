@testset "Forcing and drag: 2D" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_testruns_")  # Cleaned up when the process exits

    # 2D models
    spectral_grid = SpectralGrid(truncation = 43, nlayers = 1)
    output = NetCDFOutput(spectral_grid, path = tmp_output_path)
    add!(output, SpeedyWeather.RandomPatternOutput())
    drag = LinearVorticityDrag(spectral_grid)
    forcing = StochasticStirring(spectral_grid)
    random_process = SpectralAR1Process(spectral_grid, time_scale = Day(2))
    initial_conditions = StartFromRest(spectral_grid)

    @testset for Model in (BarotropicModel, ShallowWaterModel)
        model = Model(spectral_grid; initial_conditions, forcing, drag, random_process, output)
        model.feedback.verbose = false
        simulation = initialize!(model)

        run!(simulation, period = Day(15), output = true)
        @test simulation.model.feedback.nans_detected == false
    end
end

@testset "Forcing and drag: 3D" begin
    tmp_output_path = mktempdir(pwd(), prefix = "tmp_testruns_")  # Cleaned up when the process exits

    # 3D models
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
    output = NetCDFOutput(spectral_grid, path = tmp_output_path)
    add!(output, SpeedyWeather.RandomPatternOutput())
    drag = LinearVorticityDrag(spectral_grid)
    forcing = StochasticStirring(spectral_grid)
    random_process = SpectralAR1Process(spectral_grid, time_scale = Day(2))
    initial_conditions = StartFromRest(spectral_grid)

    @testset for Model in (PrimitiveDryModel, PrimitiveWetModel)
        model = Model(spectral_grid; initial_conditions, forcing, drag, random_process, output)
        model.feedback.verbose = false
        simulation = initialize!(model)

        run!(simulation, period = Day(15), output = true)
        @test simulation.model.feedback.nans_detected == false
    end
end

@testset "SpeedLimitDrag is antiparallel to the flow" begin
    # (Fu, Fv) = -c*max(0, |(u,v)| - speed_limit)^2 * (u,v)/|(u,v)|, i.e. it decelerates the
    # flow without rotating it. The previous (sign(u), sign(v)) form pointed along the quadrant
    # diagonal instead and was √2 too large for a 45° flow.
    spectral_grid = SpectralGrid(trunc = 15, nlayers = 2)
    c, speed_limit = 3.0e-6, 80.0
    drag = SpeedLimitDrag(spectral_grid; drag = c, speed_limit)
    model = PrimitiveDryModel(spectral_grid; drag)
    simulation = initialize!(model)
    (; variables) = simulation
    time_stepping = model.time_stepping

    u = SpeedyWeather.get_prognostic_step(variables.grid.u, time_stepping, drag)
    v = SpeedyWeather.get_prognostic_step(variables.grid.v, time_stepping, drag)
    Fu = SpeedyWeather.get_tendency_step(variables.tendencies.grid.u, time_stepping, drag)
    Fv = SpeedyWeather.get_tendency_step(variables.tendencies.grid.v, time_stepping, drag)

    # a few flow states: below the limit, purely zonal, purely meridional, and 45°
    below, zonal, meridional, diagonal = (60.0, 0.0), (100.0, 0.0), (0.0, -100.0), (70.0, 70.0)
    @testset for (uv, expected_excess) in (
            (below, 0.0), (zonal, 20.0), (meridional, 20.0), (diagonal, sqrt(2) * 70 - 80),
        )
        u .= uv[1]
        v .= uv[2]
        Fu .= 0
        Fv .= 0
        SpeedyWeather.drag!(variables, drag, model)

        speed = hypot(uv...)
        expected = c * expected_excess^2
        fu, fv = Array(parent(Fu))[1], Array(parent(Fv))[1]

        # magnitude is c*excess^2 regardless of flow direction
        @test hypot(fu, fv) ≈ expected rtol = 1.0e-5
        if expected > 0     # and it is exactly antiparallel to (u, v)
            @test fu ≈ -expected * uv[1] / speed rtol = 1.0e-5
            @test fv ≈ -expected * uv[2] / speed rtol = 1.0e-5
        end
    end

    # no drag below the limit, and no NaN for zero wind (the 0/0 guard)
    u .= 0
    v .= 0
    Fu .= 0
    Fv .= 0
    SpeedyWeather.drag!(variables, drag, model)
    @test all(iszero, Array(parent(Fu)))
    @test all(iszero, Array(parent(Fv)))
end
