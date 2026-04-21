# Helper functions (GPU-compatible, no Reactant dependency)

"""Synchronize prognostic variables from GPU to CPU simulation
using copy_variables! which handles arrays (via copyto!), Ref values,
and mutable structs (Clock) recursively."""
function sync_variables!(sim_cpu, sim_gpu)
    vars_cpu, _ = SpeedyWeather.unpack(sim_cpu)
    vars_gpu, _ = SpeedyWeather.unpack(sim_gpu)
    SpeedyWeather.copy!(vars_cpu, vars_gpu)
    return
end

"""Compare all array entries in two NamedTuples, skipping non-array entries
and keys not present in both. Returns a Dict mapping variable names to
comparison statistics."""
function compare_arrays(nt_cpu::NamedTuple, nt_gpu::NamedTuple; rtol = RTOL, atol = ATOL)
    results = Dict{Symbol, NamedTuple}()

    for key in keys(nt_cpu)
        getfield(nt_cpu, key) isa AbstractArray || continue
        haskey(nt_gpu, key) || continue

        arr_cpu = Array(getfield(nt_cpu, key))
        arr_gpu = Array(getfield(nt_gpu, key))
        abs_diff = abs.(arr_cpu .- arr_gpu)
        rel_diff = abs_diff ./ max.(abs.(arr_cpu), abs.(arr_gpu), eps(eltype(real(arr_cpu))))
        results[key] = (
            max_abs_diff = maximum(abs_diff),
            mean_abs_diff = mean(abs_diff),
            max_rel_diff = maximum(rel_diff),
            mean_rel_diff = mean(rel_diff),
            matches = isapprox(arr_cpu, arr_gpu, rtol = rtol, atol = atol),
        )
    end

    return results
end

"""Compare that the clock is running in the same way."""
function compare_clock(sim_cpu, sim_gpu)
    clock_cpu = sim_cpu.variables.prognostic.clock
    clock_gpu = sim_gpu.variables.prognostic.clock

    @test clock_cpu.n_timesteps == clock_gpu.n_timesteps
    @test clock_cpu.timestep_counter == clock_gpu.timestep_counter
    return @test DateTime(clock_cpu.time) == DateTime(clock_gpu.time)
end

"""Test tendencies after running both simulations identically."""
function test_tendencies!(sim_cpu, sim_gpu, model_name; nsteps = 1, rtol = RTOL, atol = ATOL)
    println("\n" * "-"^60)
    println("Testing tendencies ($nsteps steps)")
    println("-"^60)

    # Sync GPU state to CPU so both start from identical state
    sync_variables!(sim_cpu, sim_gpu)

    # Initialize both
    initialize!(sim_cpu; steps = nsteps)
    initialize!(sim_gpu; steps = nsteps)

    # Sync again so CPU has the exact same state as GPU post-initialize
    sync_variables!(sim_cpu, sim_gpu)

    # Run one step on each
    SpeedyWeather.time_stepping!(sim_cpu)
    SpeedyWeather.time_stepping!(sim_gpu)

    SpeedyWeather.finalize!(sim_gpu)
    SpeedyWeather.finalize!(sim_cpu)
    println("  ✓ Tendencies computed")

    # Compare tendencies
    vars_cpu, _ = SpeedyWeather.unpack(sim_cpu)
    vars_gpu, _ = SpeedyWeather.unpack(sim_gpu)
    tend_results = compare_arrays(vars_cpu.tendencies, vars_gpu.tendencies; rtol, atol)

    println("\nTendency comparison:")
    for (name, r) in tend_results
        println("  $name:")
        println("    Max absolute difference:  $(r.max_abs_diff)")
        println("    Mean absolute difference: $(r.mean_abs_diff)")
        println("    Max relative difference:  $(r.max_rel_diff)")
        println("    Mean relative difference: $(r.mean_rel_diff)")
    end

    @testset "$model_name Tendency Comparison" begin
        for (name, r) in tend_results
            @test r.matches
        end
    end

    return tend_results
end

"""Test prognostic and grid variables after running for nsteps on already-initialized simulations."""
function test_time_stepping!(sim_cpu, sim_gpu, model_name; nsteps = NSTEPS, rtol = RTOL, atol = ATOL)
    println("\n" * "-"^60)
    println("Testing time stepping ($nsteps steps)")
    println("-"^60)

    # Sync GPU state to CPU so both start from identical state
    sync_variables!(sim_cpu, sim_gpu)

    # Run time stepping
    println("  Running CPU model...")
    SpeedyWeather.run!(sim_cpu; steps = nsteps)

    println("  Running GPU model...")
    initialize!(sim_gpu; steps = nsteps)
    SpeedyWeather.time_stepping!(sim_gpu)
    SpeedyWeather.finalize!(sim_gpu)
    println("  ✓ Time stepping completed")

    # Compare results
    vars_cpu, _ = SpeedyWeather.unpack(sim_cpu)
    vars_gpu, _ = SpeedyWeather.unpack(sim_gpu)
    progn_results = compare_arrays(vars_cpu.prognostic, vars_gpu.prognostic; rtol, atol)
    grid_results = compare_arrays(vars_cpu.grid, vars_gpu.grid; rtol, atol)

    println("\nPrognostic variable comparison after $nsteps steps:")
    for (name, r) in progn_results
        println("  $name:")
        println("    Max absolute difference:  $(r.max_abs_diff)")
        println("    Mean absolute difference: $(r.mean_abs_diff)")
        println("    Max relative difference:  $(r.max_rel_diff)")
        println("    Mean relative difference: $(r.mean_rel_diff)")
    end

    println("\nGrid variable comparison after $nsteps steps:")
    for (name, r) in grid_results
        println("  $name:")
        println("    Max absolute difference:  $(r.max_abs_diff)")
        println("    Mean absolute difference: $(r.mean_abs_diff)")
        println("    Max relative difference:  $(r.max_rel_diff)")
        println("    Mean relative difference: $(r.mean_rel_diff)")
    end

    @testset "$model_name Time Stepping" begin
        @testset "Clock" begin
            compare_clock(sim_cpu, sim_gpu)
        end

        @testset "Prognostic variables" begin
            for (name, r) in progn_results
                @test r.matches
            end
        end

        @testset "Grid variables" begin
            for (name, r) in grid_results
                @test r.matches
            end
        end
    end

    return (progn_results, grid_results)
end

"""Run all correctness tests for a given model type (CPU vs GPU)."""
function test_model(ModelType::Type; trunc = TRUNC, nsteps = NSTEPS, rtol = RTOL, atol = ATOL)
    model_name = string(ModelType)

    println("="^60)
    println("$model_name: CPU vs GPU Correctness Tests")
    println("="^60)

    println("\n[1/3] Setting up CPU model...")
    model_cpu = create_cpu_model(ModelType; trunc)
    simulation_cpu = initialize!(model_cpu)
    println("  ✓ CPU model initialized (T$trunc)")

    println("\n[2/3] Setting up GPU model...")
    model_gpu = create_gpu_model(ModelType; trunc)
    simulation_gpu = initialize!(model_gpu)
    println("  ✓ GPU model initialized")

    println("\n[3/3] Spinning up models...")
    run!(simulation_cpu; period = Day(1))
    run!(simulation_gpu; period = Day(10))
    println("  ✓ Models spun up")

    @testset "$model_name CPU vs GPU" begin
        test_tendencies!(simulation_cpu, simulation_gpu, model_name; rtol, atol)
        test_time_stepping!(simulation_cpu, simulation_gpu, model_name; nsteps, rtol, atol)
    end

    println("\n" * "="^60)
    println("$model_name tests completed!")
    println("="^60)

    return nothing
end

# Run tests
test_model(PrimitiveWetModel)
