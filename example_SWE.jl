using Test
using SpeedyWeather
using Dates

function test_shallow_water_lorenz_readiness(spectral_grid, time_stepping, implicit)
    println("Checking shallow water Lorenz N-cycle readiness...")
    all_passed = true

    # Build model and simulation once — reused across all checks
    model      = ShallowWaterModel(spectral_grid; time_stepping, implicit)
    simulation = initialize!(model)

    progn, diagn, model = SpeedyWeather.unpack(simulation)
    println(pointer(diagn.tendencies.div_tend.data))
    println(pointer(SpeedyWeather.get_step(progn.div, 2).data))
    println(diagn.tendencies.div_tend.data === SpeedyWeather.get_step(progn.div, 2).data)

    progn, diagn, _ = SpeedyWeather.unpack(simulation)

    # ── 1. prognostic_variables includes :pres ───────────────────────────────
    pvars = SpeedyWeather.prognostic_variables(model)
    result = :pres in pvars
    @test result
    if result
        println("  OK 1: prognostic_variables = $pvars")
    else
        println("  FAIL 1: :pres missing from prognostic_variables(ShallowWater), got: $pvars")
        all_passed = false
    end

    # ── 2. tendency names exist as fields on diagn.tendencies ────────────────
    tnames = SpeedyWeather.tendency_names(model)
    for tname in tnames
        result = hasfield(typeof(diagn.tendencies), tname)
        @test result
        if result
            println("  OK 2: diagn.tendencies.$tname exists")
        else
            println("  FAIL 2: diagn.tendencies has no field :$tname (got: $tnames)")
            all_passed = false
        end
    end

    # ── 3. get_step works on progn.pres (2D, no layer dimension) ────────────
    try
        x = SpeedyWeather.get_step(progn.pres, 1)
        G = SpeedyWeather.get_step(progn.pres, 2)
        @test x !== nothing
        @test G !== nothing
        println("  OK 3: get_step(progn.pres, 1/2) works, size=$(size(x))")
    catch e
        @test false
        println("  FAIL 3: get_step(progn.pres, ...) threw: $e")
        all_passed = false
    end

    # ── 4. α is consistent between time_stepping and implicit ────────────────
    result = time_stepping.α == implicit.α
    @test result
    if result
        println("  OK 4: α=$(time_stepping.α) consistent between time_stepping and implicit")
    else
        println("  FAIL 4: α mismatch — time_stepping.α=$(time_stepping.α), implicit.α=$(implicit.α)")
        all_passed = false
    end

    # ── 5. initialize! stores ξ = α*Δt (radius-scaled), not α*Δt_sec ────────
    # ImplicitShallowWater.initialize! signature is (implicit, dt, args...) — no diagn/model needed
    SpeedyWeather.initialize!(model.implicit, time_stepping.Δt)
    expected_ξ = implicit.α * time_stepping.Δt
    wrong_ξ    = implicit.α * time_stepping.Δt_sec   # what you'd get if Δt_sec was passed by mistake
    result = implicit.time_step ≈ expected_ξ
    @test result
    if result
        println("  OK 5: implicit.time_step = α*Δt = $expected_ξ (radius-scaled)")
    else
        println("  FAIL 5: implicit.time_step=$(implicit.time_step)")
        println("         expected α*Δt=$expected_ξ (radius-scaled)")
        if implicit.time_step ≈ wrong_ξ
            println("         Got α*Δt_sec instead — Δt_sec was passed rather than Δt")
        end
        all_passed = false
    end

    # ── 6. Single lorenz_ncycle_step! runs without error ────────────────────
    try
        fill!(diagn.tendencies, 0, SpeedyWeather.ShallowWater)
        SpeedyWeather.dynamics_tendencies!(diagn, progn, 1, model)
        SpeedyWeather.horizontal_diffusion!(diagn, progn, model.horizontal_diffusion, model)
        SpeedyWeather.lorenz_ncycle_step!(progn, diagn, diagn.tendencies, time_stepping.Δt, model)
        @test true
        println("  OK 6: single lorenz_ncycle_step! completed without error")
    catch e
        @test false
        println("  FAIL 6: lorenz_ncycle_step! threw: $e")
        println("         ", sprint(showerror, e))
        all_passed = false
    end

    # ── 7. lorenz_implicit_correction! runs standalone without error ─────────
    try
        SpeedyWeather.lorenz_implicit_correction!(diagn, progn, model.implicit, model)
        @test true
        println("  OK 7: lorenz_implicit_correction! completed without error")
    catch e
        @test false
        println("  FAIL 7: lorenz_implicit_correction! threw: $e")
        println("         ", sprint(showerror, e))
        all_passed = false
    end

    # ── 8. Short integration stays finite ────────────────────────────────────
    try
        run!(simulation, period=Day(1))
        # Check that key prognostic fields are still finite after 1 day
        vor_ok  = all(isfinite, SpeedyWeather.get_step(progn.vor,  1).data)
        div_ok  = all(isfinite, SpeedyWeather.get_step(progn.div,  1).data)
        pres_ok = all(isfinite, SpeedyWeather.get_step(progn.pres, 1).data)
        @test vor_ok
        @test div_ok
        @test pres_ok
        if vor_ok && div_ok && pres_ok
            println("  OK 8: 1-day integration stays finite (vor, div, pres all finite)")
        else
            fields = filter(x -> !last(x), [:vor=>vor_ok, :div=>div_ok, :pres=>pres_ok])
            println("  FAIL 8: NaNs/Infs detected in: $(first.(fields))")
            all_passed = false
        end
    catch e
        @test false
        println("  FAIL 8: run! threw: $e")
        println("         ", sprint(showerror, e))
        all_passed = false
    end

    println(all_passed ? "\n✓ All checks passed — safe to run full SW simulation" :
                         "\n✗ Some checks failed — fix before running")
    return all_passed
end

# ── Run it ────────────────────────────────────────────────────────────────────
spectral_grid = SpectralGrid(trunc=42, nlayers=1)
time_stepping = NCycleLorenz(spectral_grid; cycles=4, α=0.5)
implicit      = ImplicitShallowWater(spectral_grid; α=0.5)



test_shallow_water_lorenz_readiness(spectral_grid, time_stepping, implicit)