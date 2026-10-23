@testset "Periodic schedule" begin
    # especially T170 uses 337500 milliseconds time steps
    # not representable as seconds
    for truncation in (32, 43, 64, 86, 128, 171, 256, 342)
        spectral_grid = SpectralGrid(truncation = truncation, nlayers = 1)
        for TS in (Leapfrog, NCycleLorenz)
            time_stepping = TS(spectral_grid)

            clock = Clock()
            period = Day(1)
            initialize!(clock, time_stepping, period)

            hour = Hour(2)
            schedule = Schedule(every = hour)
            initialize!(schedule, clock)

            # adapted schedule time step should be within 20%
            @test schedule.every.value ≈ Second(hour).value rtol = 2.0e-1
            @test schedule.steps ≈ period / hour rtol = 2.0e-1

            for _ in 1:clock.n_time_steps
                SpeedyWeather.time_step!(clock, clock.Δt)
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

            for _ in 1:clock.n_time_steps
                SpeedyWeather.time_step!(clock, clock.Δt)
                isscheduled(schedule, clock)
            end
            @test schedule.counter == 1
        end
    end
end

@testset "PeriodicSchedule and StepSchedule" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 1)
    time_stepping = Leapfrog(spectral_grid)
    clock = Clock()
    initialize!(clock, time_stepping, Day(1))

    periodic = PeriodicSchedule(Hour(6))
    initialize!(periodic, time_stepping)
    n = periodic.every_n_steps
    @test n == round(Int, Millisecond(Hour(6)) / time_stepping.Δt_millisec)

    never = initialize!(PeriodicSchedule(every = Second(0)), time_stepping)
    @test never.every_n_steps == 0

    step = StepSchedule(3)
    @test initialize!(step, time_stepping) === step

    executions = 0
    for _ in 1:clock.n_time_steps
        @test !isscheduled(never, clock)
        @test isscheduled(step, clock) == (clock.time_step_counter > 0 && clock.time_step_counter % 3 == 0)
        executions += isscheduled(periodic, clock)
        SpeedyWeather.time_step!(clock, clock.Δt)
    end
    @test executions == periodic.counter == clock.n_time_steps ÷ n - (clock.n_time_steps % n == 0)
end
