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

@testset "Views into concrete Reactant arrays" begin
    x = Reactant.to_rarray(ones(Float32, 6))
    v = view(x, 2:4)

    # broadcasting into a view writes into the parent
    v .*= 3
    @test Array(x) == Float32[1, 3, 3, 3, 1, 1]

    # deepcopy keeps the view aliased with the copied parent, as for Array
    copied = deepcopy((; parent = x, view = v))
    @test parent(copied.view) === copied.parent
    fill!(copied.parent, 0)
    @test all(iszero, Array(copied.view))
    @test Array(x) == Float32[1, 3, 3, 3, 1, 1]     # original untouched

    # to the CPU: the parent is transferred and the view copied out on the host (instead of compiling
    # a program that materializes the view on the device on every call)
    y = Reactant.to_rarray(reshape(Float32.(1:24), 4, 3, 2))
    w = view(y, :, 2:3, :)
    w_cpu = on_architecture(SpeedyWeather.CPU(), w)
    @test w_cpu isa Array{Float32, 3}
    @test w_cpu == Array(w)
end
