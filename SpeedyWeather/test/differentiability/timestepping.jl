# High-level tests whether time stepping in all models is differentiable.

@testset "Differentiability: Timestepping ($(nameof(model_type)))" for model_type in
        (ShallowWaterModel, PrimitiveDryModel, PrimitiveWetModel)

    # FiniteDifferences struggles with the NaN when we have a land-sea mask, so we test on AquaPlanets
    spectral_grid = SpectralGrid(trunc = 8, nlayers = 1, NF = Float64)
    model = model_type(; spectral_grid)
    simulation = initialize!(model)
    initialize!(simulation)
    run!(simulation, period = Day(3))
    (; variables, model) = simulation

    vars0 = deepcopy(variables)
    model0 = deepcopy(model)
    fresh_model() = deepcopy(model0)

    # random tangent (input) and cotangent (output), both on the prognostic variables
    function seed!(shadow, which)
        rng = Random.MersenneTwister(which)
        n = 0
        for name in (:vorticity, :divergence, :temperature, :humidity, :pressure)
            haskey(shadow.prognostic, name) || continue
            d = getfield(shadow.prognostic, name).data
            R = real(eltype(d))
            d .= eltype(d) <: Complex ?
                complex.(randn(rng, R, size(d)), randn(rng, R, size(d))) :
                randn(rng, R, size(d))
            n += length(d)
        end
        return n
    end
    v = make_zero(vars0); w = make_zero(vars0)
    @test seed!(v, 1) > 0
    seed!(w, 2)

    function run_ad(mode, seed)
        vars_new = make_zero(vars0)
        dvars_new = mode === Reverse ? deepcopy(seed) : make_zero(vars0)
        vars_old = deepcopy(vars0)
        dvars_old = mode === Reverse ? make_zero(vars0) : deepcopy(seed)
        m = fresh_model()
        autodiff(
            set_runtime_activity(mode), timestep_oop!, Const,
            Duplicated(vars_new, dvars_new), Duplicated(vars_old, dvars_old), Const(m),
        )
        return mode === Reverse ? dvars_old : dvars_new
    end

    jvp = run_ad(Forward, v)        # J v      (tangent of the new state)
    vjp = run_ad(Reverse, w)        # Jᵀ w     (gradient w.r.t. the initial state)

    # differentiating w.r.t. the initial conditions gives a finite, non-trivial gradient
    gvec = to_vec(vjp)[1]
    @test all(isfinite, gvec)
    @test any(!iszero, gvec)

    # forward and reverse must agree exactly
    lhs = sum(to_vec(jvp)[1] .* to_vec(w)[1])
    rhs = sum(to_vec(v)[1] .* gvec)
    @test lhs != 0
    @test isapprox(lhs, rhs; rtol = 1.0e-8)

    # ...and w.r.t. a model parameter, checked against a scalar finite difference in gravity.
    # `timestep_oop!(vars_new, vars_old, model, p)` rebuilds the model from `p` via `reconstruct`.
    let p = vec(parameters(model0)), dp = zero(vec(parameters(model0)))
        vars_new = make_zero(vars0)
        autodiff(
            set_runtime_activity(Reverse), timestep_oop!, Const,
            Duplicated(vars_new, deepcopy(w)), Duplicated(deepcopy(vars0), make_zero(vars0)),
            Const(fresh_model()), Duplicated(p, dp),
        )

        # the same projection, as a function of gravity alone
        function project(g)
            q = copy(p)
            q.planet.gravity = g
            out = timestep_oop(deepcopy(vars0), fresh_model(), q)
            return sum(to_vec(out)[1] .* to_vec(w)[1])
        end
        fd_gravity = central_fdm(5, 1)(project, p.planet.gravity)
        @test isapprox(dp.planet.gravity, fd_gravity; rtol = 1.0e-4)
    end
end
