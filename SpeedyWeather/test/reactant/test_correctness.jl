# Helper functions

"""Synchronize prognostic variables from Reactant to CPU simulation
using copy_variables! which handles arrays (via copyto!), Ref values,
and mutable structs (Clock) recursively."""
function sync_variables!(sim_cpu, sim_reactant)
    vars_cpu, _ = SpeedyWeather.unpack(sim_cpu)
    vars_reactant, _ = SpeedyWeather.unpack(sim_reactant)

    # Copy all variables
    SpeedyWeather.copy!(vars_cpu, vars_reactant)
    return
end

"""Look up the (rtol, atol) tolerance pair for a variable from a namedtuple of
tolerances. Variables not explicitly listed fall back to `tolerances.default`.
A trailing `_prev` (used for Leapfrog copies) is stripped before lookup so
e.g. `temperature_prev` shares the tolerance of `temperature`."""
function tolerance_for(name::Symbol, tolerances::NamedTuple)
    base = name
    sname = String(name)
    if endswith(sname, "_prev")
        base = Symbol(sname[1:(end - length("_prev"))])
    end
    return get(tolerances, base, tolerances.default)
end

"""Compare all array entries in two NamedTuples, skipping non-array entries
and keys not present in both. Returns a Dict mapping variable names to
comparison statistics. Tolerances are looked up per-variable from
`tolerances`; unknown variables get `tolerances.default`."""
function compare_arrays(nt_cpu::NamedTuple, nt_reactant::NamedTuple; tolerances::NamedTuple = TOLERANCES)
    results = Dict{Symbol, NamedTuple}()

    for key in keys(nt_cpu)
        getfield(nt_cpu, key) isa AbstractArray || continue
        haskey(nt_reactant, key) || continue

        arr_cpu = Array(getfield(nt_cpu, key))
        arr_reactant = Array(getfield(nt_reactant, key))
        abs_diff = abs.(arr_cpu .- arr_reactant)
        rel_diff = abs_diff ./ max.(abs.(arr_cpu), abs.(arr_reactant), eps(eltype(real(arr_cpu))))
        tol = tolerance_for(key, tolerances)
        results[key] = (
            max_abs_diff = maximum(abs_diff),
            mean_abs_diff = mean(abs_diff),
            max_rel_diff = maximum(rel_diff),
            mean_rel_diff = mean(rel_diff),
            rtol = tol.rtol,
            atol = tol.atol,
            matches = isapprox(arr_cpu, arr_reactant, rtol = tol.rtol, atol = tol.atol),
        )
    end

    return results
end

"""Compare that the clock is running in the same way."""
function compare_clock(sim_cpu, sim_reactant)
    clock_cpu = sim_cpu.variables.prognostic.clock
    clock_reactant = sim_reactant.variables.prognostic.clock

    @test clock_cpu.n_timesteps == clock_reactant.n_timesteps
    @test clock_cpu.timestep_counter == clock_reactant.timestep_counter
    # convert to DateTime to compare because Reactant ReactantDatetime might be used
    return @test DateTime(clock_cpu.time) == DateTime(clock_reactant.time)
end

"""Test tendencies after running both simulations identically.

Tendencies are computed from grid-point variables which are produced by
`transform!` (spectral → grid). Because the CPU and Reactant backends
execute transform! slightly differently (identical matrices but different
numerical paths through @jit), the resulting grid variables — and therefore
the tendencies — can differ at the level of the spectral transform
discrepancy. We therefore compare tendencies after running both models in
the same self-consistent way, mirroring `test_time_stepping!`."""
function test_tendencies!(sim_cpu, sim_reactant, model_name, r_first! = nothing, r_later! = nothing, mode=:reactant; nsteps = 1, tolerances::NamedTuple = TOLERANCES)
    println("\n" * "-"^60)
    println("Testing tendencies ($nsteps steps)")
    println("-"^60)

    # Pre-compile Reactant functions if needed (@compile may mutate sim_reactant
    # as a side effect, so compile first, then finalize + sync to restore clean state)
    if isnothing(r_first!) && mode == :reactant
        initialize!(sim_reactant; steps = nsteps)
        r_first! = @compile first_timesteps!(sim_reactant)
        SpeedyWeather.finalize!(sim_reactant)
    end
    if isnothing(r_later!) && mode == :reactant
        initialize!(sim_reactant; steps = nsteps)
        r_later! = @compile later_timestep!(sim_reactant)
        SpeedyWeather.finalize!(sim_reactant)
    end

    # Sync from Reactant to CPU (both now in unscaled post-finalize state)
    sync_variables!(sim_cpu, sim_reactant)

    # Initialize both (scales prognostic vars, transforms to grid)
    initialize!(sim_cpu; steps = nsteps)
    initialize!(sim_reactant; steps = nsteps)

    # Sync again so CPU has the exact same state as Reactant
    sync_variables!(sim_cpu, sim_reactant)

    # Run one step to compute tendencies
    SpeedyWeather.time_stepping!(sim_cpu)

    # Run one step on Reactant with pre-compiled functions
    if mode == :reactant
        SpeedyWeather.time_stepping!(sim_reactant, r_first!, r_later!)
    else 
        SpeedyWeather.time_stepping!(sim_reactant)
    end 

    SpeedyWeather.finalize!(sim_reactant)
    SpeedyWeather.finalize!(sim_cpu)
    println("  ✓ Tendencies computed")

    # Compare tendencies
    vars_cpu, _ = SpeedyWeather.unpack(sim_cpu)
    vars_reactant, _ = SpeedyWeather.unpack(sim_reactant)
    tend_results = compare_arrays(vars_cpu.tendencies, vars_reactant.tendencies; tolerances)

    println("\nTendency comparison:")
    for (name, r) in tend_results
        println("  $name:")
        println("    Max absolute difference:  $(r.max_abs_diff)")
        println("    Mean absolute difference: $(r.mean_abs_diff)")
        println("    Max relative difference:  $(r.max_rel_diff)")
        println("    Mean relative difference: $(r.mean_rel_diff)")
        println("    Tolerance (rtol, atol):   ($(r.rtol), $(r.atol))")
    end

    @testset "$model_name Tendency Comparison" begin
        for (name, r) in tend_results
            @test r.matches
        end
    end

    return tend_results
end

"""Test prognostic and grid variables after running for nsteps on already-initialized simulations."""
function test_time_stepping!(sim_cpu, sim_reactant, model_name, r_first! = nothing, r_later! = nothing, mode=:reactant; nsteps = NSTEPS, tolerances::NamedTuple = TOLERANCES)
    println("\n" * "-"^60)
    println("Testing time stepping ($nsteps steps)")
    println("-"^60)

    # Pre-compile Reactant functions if needed (this may mutate sim_reactant as a side effect)
    if isnothing(r_first!) && mode == :reactant
        initialize!(sim_reactant; steps = nsteps)
        r_first! = @compile first_timesteps!(sim_reactant)
    end
    if isnothing(r_later!) && mode == :reactant
        initialize!(sim_reactant; steps = nsteps)
        r_later! = @compile later_timestep!(sim_reactant)
    end

    # Sync after @compile to undo any side effects on sim_reactant
    sync_variables!(sim_cpu, sim_reactant)

    # Run time stepping
    println("  Running CPU model...")
    SpeedyWeather.run!(sim_cpu; steps = nsteps)

    println("  Running Reactant model...")
    initialize!(sim_reactant; steps = nsteps)
    if mode == :reactant
        SpeedyWeather.time_stepping!(sim_reactant, r_first!, r_later!)
    else 
        SpeedyWeather.time_stepping!(sim_reactant)
    end 
    SpeedyWeather.finalize!(sim_reactant)
    println("  ✓ Time stepping completed")

    # Compare results
    vars_cpu, _ = SpeedyWeather.unpack(sim_cpu)
    vars_reactant, _ = SpeedyWeather.unpack(sim_reactant)
    progn_results = compare_arrays(vars_cpu.prognostic, vars_reactant.prognostic; tolerances)
    grid_results = compare_arrays(vars_cpu.grid, vars_reactant.grid; tolerances)

    println("\nPrognostic variable comparison after $nsteps steps:")
    for (name, r) in progn_results
        println("  $name:")
        println("    Max absolute difference:  $(r.max_abs_diff)")
        println("    Mean absolute difference: $(r.mean_abs_diff)")
        println("    Max relative difference:  $(r.max_rel_diff)")
        println("    Mean relative difference: $(r.mean_rel_diff)")
        println("    Tolerance (rtol, atol):   ($(r.rtol), $(r.atol))")
    end

    println("\nGrid variable comparison after $nsteps steps:")
    for (name, r) in grid_results
        println("  $name:")
        println("    Max absolute difference:  $(r.max_abs_diff)")
        println("    Mean absolute difference: $(r.mean_abs_diff)")
        println("    Max relative difference:  $(r.max_rel_diff)")
        println("    Mean relative difference: $(r.mean_rel_diff)")
        println("    Tolerance (rtol, atol):   ($(r.rtol), $(r.atol))")
    end

    @testset "$model_name Time Stepping" begin
        @testset "Clock" begin
            compare_clock(sim_cpu, sim_reactant)
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

"""Run all correctness tests for a given model type."""
function test_model(ModelType::Type; trunc = TRUNC, nsteps = NSTEPS, tolerances::NamedTuple = TOLERANCES, precompile = false)
    model_name = string(ModelType)

    println("="^60)
    println("$model_name: CPU vs Reactant Correctness Tests")
    println("="^60)

    # Setup CPU model
    println("\n[1/3] Setting up CPU model...")
    model_cpu = create_cpu_model(ModelType; trunc)
    simulation_cpu = initialize!(model_cpu)
    println("  ✓ CPU model initialized (T$trunc)")

    # Setup Reactant model
    println("\n[2/3] Setting up Reactant model...")
    model_reactant = create_reactant_model(ModelType; trunc)
    simulation_reactant = initialize!(model_reactant)
    println("  ✓ Reactant model initialized")

    # spin up models a bit
    println("\n[3/3] Spinning up models...")
    run!(simulation_cpu; period = Day(1)) # dummy steps so that `simulation_cpu` also continues a prior simulation, we could also manually set e.g. `continue_with_leapfrog` arguments instead, but those might change in the future
    run!(simulation_reactant; period = Day(10)) # we copy from Reactant to cpu later so only there we need a spin up, we need a long spin up, because with some ICs we get mostly zonal flow otherwise
    println("  ✓ Models spun up")

    # Run tests
    # Pre-compile Reactant functions once
    if precompile
        println("\n[4/4] Pre-compiling Reactant functions...")
        initialize!(simulation_reactant; steps = nsteps)
        r_first! = @compile first_timesteps!(simulation_reactant)
        r_later! = @compile later_timestep!(simulation_reactant)
        println("  ✓ Reactant functions compiled")
    else
        r_first! = nothing
        r_later! = nothing
    end

    @testset "$model_name CPU vs Reactant" begin
        tend_results = test_tendencies!(simulation_cpu, simulation_reactant, model_name, r_first!, r_later!; tolerances)
        stepping_results = test_time_stepping!(simulation_cpu, simulation_reactant, model_name, r_first!, r_later!; nsteps, tolerances)
    end

    println("\n" * "="^60)
    println("$model_name tests completed!")
    println("="^60)

    return nothing
end

# Run tests
test_model(PrimitiveWetModel)
#test_model(BarotropicModel)
