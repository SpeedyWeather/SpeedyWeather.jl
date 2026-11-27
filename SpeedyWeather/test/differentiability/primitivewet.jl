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
# Float64, not the default Float32: in Float32 the directional finite difference is limited by
# round-off, not by the gradient. FiniteDifferences sizes its step assuming one ULP of error in the
# scalar, but a whole Float32 tendency evaluation is ~10x noisier, so the FD estimate scatters by ~1%
# (dynamics_tendencies!: median 7e-3, max 4e-2 over 20 random directions; implicit_correction!'s
# adjoint identity up to 2e-5) and the outcome flipped between CI runners. In Float64 the same
# checks agree to ≤ 2e-10 (FD) and ≤ 4e-14 (adjoint identity), so the tolerances below are tight
# with a margin of > 1000x. Float32 AD is still covered by the parameter-AD testset below and by
# the Barotropic tests. See docs/dev/2026-11/2026-11-18-extended-diff-fd-noise.md.
_group1_spectral_grid() = SpectralGrid(; NF = Float64, truncation = 9, nlayers = 4)

@testset "Differentiability: PrimitiveWet dynamics_tendencies!" begin
    spectral_grid = _group1_spectral_grid()
    model = PrimitiveWetModel(; spectral_grid, time_stepping = Leapfrog(spectral_grid), dynamics_only = true)
    simulation = _spinup(model)
    vars0 = deepcopy(simulation.variables)
    @test all(isfinite, to_vec(vars0)[1])   # the spin-up must not have blown up

    _adjoint_check(
        SpeedyWeather.dynamics_tendencies!, vars0, Const(model);
        ingroup = :grid, outgroup = :tendencies, rtol_ad = 1.0e-10, rtol_fd = 1.0e-6,
    )
end

@testset "Differentiability: PrimitiveWet implicit_correction!" begin
    spectral_grid = _group1_spectral_grid()
    model = PrimitiveWetModel(; spectral_grid, time_stepping = Leapfrog(spectral_grid), dynamics_only = true)
    simulation = _spinup(model)
    vars0 = deepcopy(simulation.variables)

    _adjoint_check(
        SpeedyWeather.implicit_correction!, vars0,
        Const(model.implicit), Const(model.time_stepping), Const(model);
        ingroup = :prognostic, outgroup = :tendencies, rtol_ad = 1.0e-10, rtol_fd = 1.0e-6,
    )
end

#
# GROUP 2 — PHYSICS, PARAMETER AD
#
@testset "Differentiability: PrimitiveWet parameterization_tendencies! (parameter AD)" begin
    # Regression guard for reverse mode with Duplicated(model) through the CPU column loops:
    # - the get_parameterizations NamedTuple (or its elements) must not be materialized as arguments
    #   of a non-inlined function (Julia 1.12: "Vararg length is negative"; EnzymeNoTypeError)
    # - llvm/llvm-project#219885: this testset aborted with "Broken module found" in Enzyme's post-AD
    #   GVN, worked around by SPLIT_COLUMN_LOOPS in parameterizations/tendencies.jl. Whether the
    #   workaround is still needed: set SPLIT_COLUMN_LOOPS = false and run this testset.
    spectral_grid = SpectralGrid(truncation = 6, nlayers = 1)
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
    spectral_grid = _group1_spectral_grid()
    model = PrimitiveWetModel(; spectral_grid, time_stepping = Leapfrog(spectral_grid), dynamics_only = true)
    simulation = _spinup(model)
    vars0 = deepcopy(simulation.variables)

    # The finite-difference direction leaves `humidity` out. `transform!` calls `hole_filling!` to
    # remove negative humidity, which is not differentiable: perturbing humidity across that kink
    # makes FD disagree with AD by 9e-2 at the default step and still by 4e-3 at max_range = 1e-6,
    # where the same FD without humidity agrees to 3e-7. The error only shrinks as the stencil stops
    # crossing near-zero humidity values, the signature of a kink rather than truncation error.
    # The adjoint identity below still covers humidity in full.
    _adjoint_check(
        SpeedyWeather.transform!, vars0, Const(model);
        ingroup = :prognostic, outgroup = :grid, rtol_ad = 1.0e-10, rtol_fd = 1.0e-6,
        fd_only = (:vorticity, :divergence, :temperature, :pressure),
    )
end
