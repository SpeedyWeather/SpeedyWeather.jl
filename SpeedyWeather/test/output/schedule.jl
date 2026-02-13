@testset "Periodic schedule" begin
    # especially T170 uses 337500 milliseconds time steps
    # not representable as seconds
    for trunc in (31, 42, 63, 85, 127, 170, 255, 341)
        spectral_grid = SpectralGrid(trunc = trunc, nlayers = 1)
        time_stepping = Leapfrog(spectral_grid)

        clock = Clock()
        period = Day(1)
        initialize!(clock, time_stepping, period)

        hour = Hour(2)
        schedule = Schedule(every = hour)
        initialize!(schedule, clock)

        # adapted schedule time step should be within 20%
        @test schedule.every.value ≈ Second(hour).value rtol = 2.0e-1
        @test schedule.steps ≈ period / hour rtol = 2.0e-1

        for _ in 1:clock.n_timesteps
            SpeedyWeather.timestep!(clock, clock.Δt)
            isscheduled(schedule, clock)
        end

        @test schedule.counter == sum(schedule.schedule) > 1

        # reinitialize clock
        clock = Clock()
        period = Day(1)
        initialize!(clock, time_stepping, period)

        schedule = Schedule(clock.time + Hour(6))
        initialize!(schedule, clock)
        @test sum(schedule.schedule) == 1
        @test schedule.counter == 0     # no execution of isscheduled yet

        for _ in 1:clock.n_timesteps
            SpeedyWeather.timestep!(clock, clock.Δt)
            isscheduled(schedule, clock)
        end
        @test schedule.counter == 1
    end
end
