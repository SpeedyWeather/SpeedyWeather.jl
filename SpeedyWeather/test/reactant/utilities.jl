# Test other Reactant related utilities

@testset "copy!(::Clock, ::Clock) with Reactant" begin
    # create a standard CPU clock with non-default values
    clock_cpu = Clock(
        time = DateTime(2020, 6, 15), start = DateTime(2020, 6, 1), period = Day(14),
        step_counter = 5, time_step_counter = 5, n_steps = 10, n_time_steps = 10,
        Δt = Millisecond(3600_000)
    )

    # create a Reactant clock (has ConcreteRNumber fields via track_numbers)
    clock_reactant = SpeedyWeather.Clock(ReactantDevice())

    # copy CPU -> Reactant
    copy!(clock_reactant, clock_cpu)

    @test Int(clock_reactant.time_step_counter) == clock_cpu.time_step_counter
    @test Int(clock_reactant.step_counter) == clock_cpu.step_counter
    @test Int(clock_reactant.n_time_steps) == clock_cpu.n_time_steps
    @test Int(clock_reactant.n_steps) == clock_cpu.n_steps
    @test DateTime(clock_reactant.time) == clock_cpu.time
    @test DateTime(clock_reactant.start) == clock_cpu.start

    # copy Reactant -> CPU (fresh clock)
    clock_cpu2 = Clock()
    copy!(clock_cpu2, clock_reactant)

    @test clock_cpu2.step_counter == clock_cpu.step_counter
    @test clock_cpu2.time_step_counter == clock_cpu.time_step_counter
    @test clock_cpu2.n_time_steps == clock_cpu.n_time_steps
    @test clock_cpu2.n_steps == clock_cpu.n_steps
    @test clock_cpu2.time == clock_cpu.time
    @test clock_cpu2.start == clock_cpu.start
    @test clock_cpu2.period == clock_cpu.period
    @test clock_cpu2.Δt == clock_cpu.Δt
end
