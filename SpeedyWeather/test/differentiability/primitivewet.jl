import Random
Random.seed!(123)

#    NO `j′vp`. `FiniteDifferences._j′vp` is `transpose(jacobian(fdm, f, x)) * ȳ`, i.e. it
#    MATERIALISES the full Jacobian. `to_vec` of a PrimitiveWet `Variables` here is 112_532 long, so
#    that Jacobian is 112_532² ≈ 50–100 GB and the process is OOM-killed before printing anything.
#    (The barotropic tests get away with it: 4_434² ≈ 157 MB.)
#
#    Instead we use the adjoint identity
#
#        <J v, w>  ==  <v, Jᵀ w>
#
#    which needs one forward-mode and one reverse-mode call and no Jacobian at all. It is exact —
#    no step size to tune — and it pins reverse mode against an independently computed forward mode.
#    As an arbiter that does not involve AD at all, we add a DIRECTIONAL finite difference of the
#    same scalar, which costs a handful of primal evaluations rather than one per input dimension.

_arrays(v, group) = AbstractArray[
    getfield(getfield(v, group), k) for k in keys(getfield(v, group))
        if getfield(getfield(v, group), k) isa AbstractArray
]

# seed every array of one variable group with random values
function _seed_group!(shadow, group, seed; only = nothing)
    rng = Random.MersenneTwister(seed)
    n = 0
    for k in keys(getfield(shadow, group))
        (only === nothing || k in only) || continue
        f = getfield(getfield(shadow, group), k)
        f isa AbstractArray || continue
        R = real(eltype(f))
        f .= eltype(f) <: Complex ?
            complex.(randn(rng, R, size(f)), randn(rng, R, size(f))) :
            randn(rng, R, size(f))
        n += length(f)
    end
    return n
end

# real inner product restricted to one variable group. `v` is zero outside `ingroup` and `w` zero
# outside `outgroup`, so restricting the sums this way is exact, and it keeps unrelated fields out.
_gdot(a, b, group) = sum(
    sum(real.(x) .* real.(y)) + (eltype(x) <: Complex ? sum(imag.(x) .* imag.(y)) : 0)
        for (x, y) in zip(_arrays(a, group), _arrays(b, group))
)

"""
Check `<J v, w> == <v, Jᵀ w>` for the in-place `f!(vars, args...)`, and cross-check the same scalar
against a directional finite difference. `fd_only` restricts the FD direction to a subset of the
input variables (used to keep a non-differentiable path out of the finite difference).
"""
function _adjoint_check(f!, vars0, args...; ingroup, outgroup, rtol_ad, rtol_fd, fd_only = nothing)
    v = make_zero(vars0)
    w = make_zero(vars0)
    @test _seed_group!(v, ingroup, 1) > 0
    @test _seed_group!(w, outgroup, 2) > 0

    jv = deepcopy(v)
    x1 = deepcopy(vars0)
    autodiff(set_runtime_activity(Forward), f!, Const, Duplicated(x1, jv), args...)

    jtw = deepcopy(w)
    x2 = deepcopy(vars0)
    autodiff(set_runtime_activity(Reverse), f!, Const, Duplicated(x2, jtw), args...)

    @test all(isfinite, to_vec(jtw)[1])
    @test any(!iszero, to_vec(jtw)[1])

    lhs = _gdot(jv, w, outgroup)        # <J v, w>   from forward mode
    rhs = _gdot(v, jtw, ingroup)        # <v, Jᵀ w>  from reverse mode
    @test lhs != 0
    @test isapprox(lhs, rhs; rtol = rtol_ad)

    # independent arbiter: directional finite difference of eps -> <f(x + eps*v), w>
    vfd = v
    if fd_only !== nothing
        vfd = make_zero(vars0)
        _seed_group!(vfd, ingroup, 1; only = fd_only)
    end
    function scalar(eps)
        x = deepcopy(vars0)
        for (xa, va) in zip(_arrays(x, ingroup), _arrays(vfd, ingroup))
            xa .+= eps .* va
        end
        f!(x, map(a -> a.val, args)...)
        return _gdot(x, w, outgroup)
    end
    jv_fd = deepcopy(vfd)
    x3 = deepcopy(vars0)
    autodiff(set_runtime_activity(Forward), f!, Const, Duplicated(x3, jv_fd), args...)
    @test isapprox(central_fdm(5, 1)(scalar, 0.0), _gdot(jv_fd, w, outgroup); rtol = rtol_fd)
    return nothing
end

# one day of spin-up: five days is unstable at this resolution (see the note at the top)
_spinup(model) = initialize_with_spinup!(model, Day(1), Hour(6))

#
# GROUP 1 — DYNAMICAL CORE (dynamics_only model, fast)
#
@testset "Differentiability: PrimitiveWet dynamics_tendencies!" begin
    spectral_grid = SpectralGrid(truncation = 9, nlayers = 4)
    model = PrimitiveWetModel(; spectral_grid, time_stepping = Leapfrog(spectral_grid), dynamics_only = true)
    simulation = _spinup(model)
    vars0 = deepcopy(simulation.variables)
    @test all(isfinite, to_vec(vars0)[1])   # the spin-up must not have blown up

    _adjoint_check(
        SpeedyWeather.dynamics_tendencies!, vars0, Const(model);
        ingroup = :grid, outgroup = :tendencies, rtol_ad = 1.0e-5, rtol_fd = 1.0e-2,
    )
end

@testset "Differentiability: PrimitiveWet implicit_correction!" begin
    spectral_grid = SpectralGrid(truncation = 9, nlayers = 4)
    model = PrimitiveWetModel(; spectral_grid, time_stepping = Leapfrog(spectral_grid), dynamics_only = true)
    simulation = _spinup(model)
    vars0 = deepcopy(simulation.variables)

    _adjoint_check(
        SpeedyWeather.implicit_correction!, vars0,
        Const(model.implicit), Const(model.time_stepping), Const(model);
        ingroup = :prognostic, outgroup = :tendencies, rtol_ad = 1.0e-5, rtol_fd = 1.0e-2,
    )
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
    spectral_grid = SpectralGrid(truncation = 9, nlayers = 4)
    model = PrimitiveWetModel(; spectral_grid, time_stepping = Leapfrog(spectral_grid), dynamics_only = true)
    simulation = _spinup(model)
    vars0 = deepcopy(simulation.variables)

    # The finite-difference direction leaves `humidity` out. `transform!` calls `hole_filling!` to
    # remove negative humidity, which is not differentiable: perturbing humidity across that kink
    # makes FD disagree with AD by ~56%, and the disagreement does NOT shrink with the step size
    # (5.6e-1 at the default range, 5.7e-1 at max_range=1e-3) — the signature of a genuine kink
    # rather than truncation error. The adjoint identity below still covers humidity in full, and
    # forward vs reverse agree there to ~1e-7.
    _adjoint_check(
        SpeedyWeather.transform!, vars0, Const(model);
        ingroup = :prognostic, outgroup = :grid, rtol_ad = 1.0e-5, rtol_fd = 3.0e-2,
        fd_only = (:vorticity, :divergence, :temperature, :pressure),
    )
end
