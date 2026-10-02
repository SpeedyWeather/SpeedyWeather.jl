# Forward mode over a whole time step, checked against reverse mode.
#
# For the same time step, forward and reverse must satisfy the dot-product identity
#
#     <J v, w>  ==  <v, Jᵀ w>
#
# for any input tangent `v` and output cotangent `w`. This is exact: there is no finite-difference
# step size to tune and no requirement that the model be smooth, which makes it a far sharper check
# than comparing either mode against finite differences. It is also what settled the grid `u`/`v`
# question — a fixed-step central difference disagreed with the JVP by ~5e-3 there, but the two AD
# modes agreed to ~1e-15 and a 9th-order FD rule converged onto the AD value.

function _seed_prognostic!(shadow, vars0, seed) # seed only prognostic atmospheric variables with rms
    rng = Random.MersenneTwister(seed)
    n = 0
    for name in (:vorticity, :divergence, :temperature, :humidity, :pressure)
        haskey(shadow.prognostic, name) || continue
        d = getfield(shadow.prognostic, name).data
        x = getfield(vars0.prognostic, name).data
        rms = sqrt(sum(abs2, x) / max(length(x), 1))
        rms = iszero(rms) ? one(rms) : rms
        R = real(eltype(d))
        d .= eltype(d) <: Complex ?
            rms .* complex.(randn(rng, R, size(d)), randn(rng, R, size(d))) :
            rms .* randn(rng, R, size(d))
        n += length(d)
    end
    return n
end

@testset "Differentiability: forward vs reverse time_step! ($(nameof(MT)))" for (MT, nlayers) in
    ((BarotropicModel, 1), (ShallowWaterModel, 1), (PrimitiveWetModel, 2))

    spectral_grid = SpectralGrid(; truncation = 6, nlayers, NF = Float64)
    simulation = initialize!(MT(; spectral_grid))
    initialize!(simulation)
    run!(simulation, period = Hour(6))
    vars0 = deepcopy(simulation.variables)
    model0 = deepcopy(simulation.model)

    v = make_zero(vars0)
    w = make_zero(vars0)
    @test _seed_prognostic!(v, vars0, 1) > 0
    _seed_prognostic!(w, vars0, 2)

    # `make_zero` must produce a genuinely zero shadow: in forward mode the shadow IS the tangent,
    # so a leftover value anywhere (scratch, a scalar Ref) is a bogus seed that corrupts the JVP.
    @test iszero(sum(abs, to_vec(make_zero(vars0))[1]))

    function run_ad(mode, seed)
        vars = deepcopy(vars0)
        dvars = deepcopy(seed)
        model = deepcopy(model0)
        dmodel = make_zero(model)
        autodiff(
            set_runtime_activity(mode), SpeedyWeather.time_step!, Const,
            Duplicated(vars, dvars),
            Duplicated(model.time_stepping, dmodel.time_stepping),
            Duplicated(model, dmodel),
        )
        return dvars
    end

    jvp = run_ad(Forward, v)        # J v
    vjp = run_ad(Reverse, w)        # Jᵀ w

    jv = to_vec(jvp)[1]
    @test all(isfinite, jv)
    @test any(!iszero, jv)

    lhs = sum(jv .* to_vec(w)[1])                   # <J v, w>
    rhs = sum(to_vec(v)[1] .* to_vec(vjp)[1])       # <v, Jᵀ w>
    @test lhs != 0
    @test isapprox(lhs, rhs; rtol = 1.0e-10)

    # A GRID-SPACE cotangent as well, not just prognostic. Seeding only prognostic variables leaves
    # the grid-space block of the Jacobian untested
    let w_grid = make_zero(vars0)
        rng = Random.MersenneTwister(3)
        w_grid.grid.u.data .= randn(rng, eltype(w_grid.grid.u.data), size(w_grid.grid.u.data))

        vjp_grid = run_ad(Reverse, w_grid)
        jtw = to_vec(vjp_grid)[1]
        @test any(!iszero, jtw)                      # the adjoint path must not vanish

        lhs_g = sum(jv .* to_vec(w_grid)[1])         # <J v, w_grid>
        rhs_g = sum(to_vec(v)[1] .* jtw)             # <v, Jᵀ w_grid>
        @test lhs_g != 0
        @test isapprox(lhs_g, rhs_g; rtol = 1.0e-10)
    end
end
