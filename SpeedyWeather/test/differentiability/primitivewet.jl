import Random
Random.seed!(123)

#
# GROUP 1 — DYNAMICAL CORE (dynamics_only model, fast)
#
@testset "Differentiability: PrimitiveWet dynamics_tendencies!" begin
    spectral_grid = SpectralGrid(trunc = 8, nlayers = 4)
    model = PrimitiveWetModel(; spectral_grid, time_stepping = Leapfrog(spectral_grid), dynamics_only = true)
    simulation = initialize_with_spinup!(model)
    adsim = ADSimulation(simulation)

    vars, dvars = ADseed(adsim, :tendencies)

    @info "Running reverse-mode AD"
    @time autodiff(set_runtime_activity(Reverse), SpeedyWeather.dynamics_tendencies!, Const, Duplicated(vars, dvars), Const(model))

    dvec, _ = to_vec(dvars)
    @test all(isfinite.(dvec))
    @test any(abs.(dvec) .> 0)

    vars2, dvars2 = ADseed(adsim, :tendencies)
    function dynamics_tendencies(vars, model)
        vars_new = deepcopy(vars)
        SpeedyWeather.dynamics_tendencies!(vars_new, deepcopy(model))
        return vars_new
    end
    @info "Running finite differences"
    fd_vjp = @time FiniteDifferences.j′vp(central_fdm(5, 1), x -> dynamics_tendencies(x, model), dvars2, vars2)
    @test_broken all(isapprox.(to_vec(fd_vjp[1])[1], to_vec(dvars)[1], rtol = 1.0e-4, atol = 1.0e-1))
end

@testset "Differentiability: PrimitiveWet implicit_correction!" begin
    spectral_grid = SpectralGrid(trunc = 8, nlayers = 4)
    model = PrimitiveWetModel(; spectral_grid, time_stepping = Leapfrog(spectral_grid), dynamics_only = true)
    simulation = initialize_with_spinup!(model)
    adsim = ADSimulation(simulation)

    vars, dvars = ADseed(adsim, :tendencies)

    @info "Running reverse-mode AD"
    @time autodiff(
        set_runtime_activity(Reverse), SpeedyWeather.implicit_correction!, Const,
        Duplicated(vars, dvars), Const(model.implicit), Const(model.time_stepping), Const(model),
    )

    dvec, _ = to_vec(dvars)
    @test all(isfinite.(dvec))

    vars2, dvars2 = ADseed(adsim, :tendencies)
    function implicit_correction(vars, model)
        vars_new = deepcopy(vars)
        m = deepcopy(model)
        SpeedyWeather.implicit_correction!(vars_new, m.implicit, m.time_stepping, m)
        return vars_new
    end
    @info "Running finite differences"
    fd_vjp = @time FiniteDifferences.j′vp(central_fdm(9, 1), x -> implicit_correction(x, model), dvars2, vars2)
    @test all(isapprox.(to_vec(fd_vjp[1])[1], to_vec(dvars)[1], rtol = 1.0e-4, atol = 1.0e-1))
end

#
# GROUP 2 — PHYSICS, PARAMETER AD
#
@testset "Differentiability: PrimitiveWet parameterization_tendencies! (parameter AD)" begin
    # Regression guard: with Duplicated(model), the CPU column-parameterization loop
    # must not materialize the get_parameterizations NamedTuple inside a non-inlined
    # function — that broke Enzyme's reverse-mode primal on Julia 1.12
    # ("Vararg length is negative"); fixed by the direct model-field access in
    # _column_parameterizations_cpu!(vars, model).
    spectral_grid = SpectralGrid(trunc = 5, nlayers = 1)
    model = PrimitiveWetModel(; spectral_grid)
    simulation = initialize!(model)
    initialize!(simulation)
    run!(simulation, period = Hour(6))          # spin-up for nonzero fields
    adsim = ADSimulation(simulation)
    m = adsim.model
    dm = make_zero(m)

    vars, dvars = ADseed(adsim, :tendencies)

    @info "Running reverse-mode AD (parameter AD)"
    @time autodiff(
        set_runtime_activity(Reverse), SpeedyWeather.parameterization_tendencies!, Const,
        Duplicated(vars, dvars), Duplicated(m, dm),
    )

    dvec, _ = to_vec(dvars)
    @test all(isfinite.(dvec))
end

@testset "Differentiability: PrimitiveWet transform!(::Variables)" begin
    spectral_grid = SpectralGrid(trunc = 8, nlayers = 4)
    model = PrimitiveWetModel(; spectral_grid, time_stepping = Leapfrog(spectral_grid), dynamics_only = true)
    simulation = initialize_with_spinup!(model)
    adsim = ADSimulation(simulation)

    vars, dvars = ADseed(adsim, :grid)

    @info "Running reverse-mode AD"
    @time autodiff(set_runtime_activity(Reverse), SpeedyWeather.transform!, Const, Duplicated(vars, dvars), Const(model))

    dvec, _ = to_vec(dvars)
    @test all(isfinite.(dvec))
    @test any(abs.(dvec) .> 0)

    vars2, dvars2 = ADseed(adsim, :grid)
    function transform_step(vars, model)
        vars_new = deepcopy(vars)
        SpeedyWeather.transform!(vars_new, deepcopy(model))
        return vars_new
    end
    @info "Running finite differences"
    fd_vjp = @time FiniteDifferences.j′vp(central_fdm(5, 1), x -> transform_step(x, model), dvars2, vars2)
    @test all(isapprox.(to_vec(fd_vjp[1].prognostic)[1], to_vec(dvars.prognostic)[1], rtol = 1.0e-3, atol = 1.0e-3))
end
