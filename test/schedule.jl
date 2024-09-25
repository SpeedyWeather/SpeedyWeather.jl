@testset "Periodic schedule" begin
    # especially T170 uses 337500 milliseconds time steps
    # not representable as seconds
    for trunc in (31,42,63,85,127,170,255,341)
        spectral_grid = SpectralGrid(trunc=trunc, nlayers=1)
        time_stepping = Leapfrog(spectral_grid)

        clock = Clock()
        period = Day(1)
        SpeedyWeather.set_period!(clock, period)
        initialize!(clock, time_stepping)

        hour = Hour(1)
        schedule = Schedule(every=hour)
        initialize!(schedule, clock)

        # adapted schedule time step should be within 20%
        @test schedule.every.value ≈ Second(hour).value rtol=2e-1
        @test schedule.steps ≈ period/hour rtol=2e-1
    end
end