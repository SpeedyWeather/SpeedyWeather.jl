# `@parameterized` expands to code that refers to `ParameterEditing` unqualified, so the name has
# to be bound here — `using SpeedyWeather` alone does not bring it into scope.
import SpeedyWeather: ParameterEditing

@testset "Differentiability of reconstruct" begin
    SpeedyWeather.@parameterized @kwdef struct TestStruct{T}
        @param a::T = 1.0
        @param b::T = 2.0
        c::T = 3.0
    end

    function testfunc!(y, x, p, obj::TestStruct)
        new_obj = SpeedyWeather.reconstruct(obj, p)
        y[1] = new_obj.a * x[1]^2
        y[2] = new_obj.b * x[2]^3 + new_obj.c
        return nothing
    end

    y, dy = zeros(2), ones(2)
    x, dx = 2 * ones(2), zeros(2)
    obj = TestStruct()
    p, dp = vec(parameters(obj)), zero(vec(parameters(obj)))
    @test length(p) == length(dp) == 2
    autodiff(Reverse, testfunc!, Const, Duplicated(y, dy), Duplicated(x, dx), Duplicated(p, dp), Const(obj))
    # check derivatives
    @test all(dp .== [x[1]^2, x[2]^3])
    @test all(dx .== [2 * p.a * x[1], 3 * p.b * x[2]^2])

    # FORWARD mode on a small struct: no rule involved (the rule is restricted to `AbstractModel`),
    # Enzyme differentiates `reconstruct` natively here. y1 = a*x1², y2 = b*x2³ + c, so seeding both
    # parameters with 1 gives dy = (da*x1², db*x2³) = (x1², x2³).
    @testset "forward mode, small struct (no rule)" begin
        let y = zeros(2), dy = zeros(2), dp = zero(p)
            dp .= 1
            autodiff(
                Forward, testfunc!, Const,
                Duplicated(y, dy), Const(x), Duplicated(p, dp), Const(obj),
            )
            @test dy ≈ [x[1]^2, x[2]^3]
        end

        # width > 1: batch members stay independent
        let y = zeros(2), d1 = zeros(2), d2 = zeros(2)
            p1 = zero(p); p1[1] = 1     # only `a`
            p2 = zero(p); p2[2] = 1     # only `b`
            autodiff(
                Forward, testfunc!, Const,
                BatchDuplicated(y, (d1, d2)), Const(x), BatchDuplicated(p, (p1, p2)), Const(obj),
            )
            @test d1 ≈ [x[1]^2, 0]
            @test d2 ≈ [0, x[2]^3]
        end
    end
end

# The case the `EnzymeRules.forward` rule in SpeedyWeatherEnzymeExt exists for: `reconstruct` on a
# full model type, which without it overruns Enzyme's type-analysis size budget
# (`EnzymeNoTypeError`). `LinearVorticityDrag` gives an exact oracle — seeding every parameter with
# 1 makes the tangent of any single reconstructed parameter exactly 1, and a directional seed
# reproduces itself.
#
# `PrimitiveWetModel` is included deliberately: its type is far larger than `BarotropicModel`
# (91 parameters vs 15), so it is what shows the rule removes the dependence on type size rather
# than just raising a threshold.
@testset "Forward-mode reconstruct on $(nameof(MT))" for (MT, nlayers) in
        ((BarotropicModel, 1), (PrimitiveWetModel, 2))

    spectral_grid = SpectralGrid(; truncation = 9, nlayers, NF = Float64)
    model = MT(; spectral_grid, drag = LinearVorticityDrag(spectral_grid))
    p = vec(parameters(model))
    @test length(p) > 0

    readout(m, q) = SpeedyWeather.reconstruct(m, q).drag.c
    @test readout(model, p) == model.drag.c

    # all-ones seed: the tangent of any single reconstructed parameter is 1
    dp = zero(p); dp .= 1
    @test autodiff(set_runtime_activity(Forward), readout, Duplicated, Const(model), Duplicated(p, dp))[1] ≈ 1

    # directional seed on drag.c alone reproduces itself, and does not pick up other parameters
    dp2 = zero(p); dp2.drag.c = 2.5
    @test autodiff(set_runtime_activity(Forward), readout, Duplicated, Const(model), Duplicated(p, dp2))[1] ≈ 2.5

    # a seed that misses drag.c must give a zero tangent
    dp3 = zero(p); dp3.planet .= 1
    @test autodiff(set_runtime_activity(Forward), readout, Duplicated, Const(model), Duplicated(p, dp3))[1] ≈ 0

    # and with an active model shadow rather than Const
    dmodel = make_zero(model)
    @test autodiff(set_runtime_activity(Forward), readout, Duplicated, Duplicated(model, dmodel), Duplicated(p, dp2))[1] ≈ 2.5
end

@testset "Reverse-mode reconstruct on $(nameof(MT))" for (MT, nlayers) in
        ((BarotropicModel, 1), (PrimitiveWetModel, 2))

    spectral_grid = SpectralGrid(; truncation = 9, nlayers, NF = Float64)
    model = MT(; spectral_grid, drag = LinearVorticityDrag(spectral_grid))
    p = vec(parameters(model))

    # reading one parameter back out has gradient e_{drag.c}: exactly 1 there, exactly 0 elsewhere
    readout(m, q) = SpeedyWeather.reconstruct(m, q).drag.c
    let dp = zero(p)
        autodiff(set_runtime_activity(Reverse), readout, Active, Const(model), Duplicated(p, dp))
        @test dp.drag.c ≈ 1
        # no leakage into any other parameter
        @test sum(abs2, vec(dp)) - dp.drag.c^2 ≈ 0 atol = 1.0e-20
    end

    # an active model shadow exercises the aliasing path in the rule (the output shadow's
    # non-parameter fields alias `obj`'s shadow) and must give the same gradient
    let dp = zero(p), dmodel = make_zero(model)
        autodiff(set_runtime_activity(Reverse), readout, Active, Duplicated(model, dmodel), Duplicated(p, dp))
        @test dp.drag.c ≈ 1
    end

    # a readout mixing several parameters recovers each coefficient independently
    mixed(m, q) = (mm = SpeedyWeather.reconstruct(m, q); 3 * mm.drag.c + 2 * mm.planet.gravity)
    let dp = zero(p)
        autodiff(set_runtime_activity(Reverse), mixed, Active, Const(model), Duplicated(p, dp))
        @test dp.drag.c ≈ 3
        @test dp.planet.gravity ≈ 2
    end

    # the rule must ACCUMULATE into a pre-seeded cotangent, not overwrite it
    let dp = zero(p)
        dp.drag.c = 10
        autodiff(set_runtime_activity(Reverse), readout, Active, Const(model), Duplicated(p, dp))
        @test dp.drag.c ≈ 11
    end
end

# The tracer-registry `iterate` rules (SpeedyWeatherEnzymeExt). The forward rule unblocks forward
# mode (Enzyme fails with an internal error compiling `iterate(::Dict)`); the reverse rule exists
# because a function carrying any custom rule must have one for whichever mode is running, which
# NESTED AD (reverse over forward) needs. Iterating the registry yields only a `Symbol` name and a
# `Tracer` config, so neither rule may contribute a derivative — the differentiable tracer data is
# reached separately via `vars.*.tracers[name]`.
@testset "Tracer registry iterate rules" begin
    spectral_grid = SpectralGrid(truncation = 9, nlayers = 1, NF = Float64)
    model = BarotropicModel(; spectral_grid, drag = LinearVorticityDrag(spectral_grid))

    # sum over the registry, scaled by an active parameter: the loop must not add any derivative of
    # its own, so d/dx is just (number of tracers) * 1 — 0 for an empty registry.
    function count_tracers(m, x)
        n = 0.0
        for (name, tracer) in m.tracers
            n += tracer.active ? 1.0 : 0.0
        end
        return n * x[1]
    end

    @test count_tracers(model, [2.0]) == 0          # default registry is empty
    let dx = zeros(1)
        autodiff(set_runtime_activity(Reverse), count_tracers, Active, Const(model), Duplicated([2.0], dx))
        @test dx[1] ≈ 0                              # empty registry -> no contribution
    end
    let dx = zeros(1)
        d = autodiff(set_runtime_activity(Forward), count_tracers, Duplicated, Const(model), Duplicated([2.0], [1.0]))[1]
        @test d ≈ 0
    end

    # with tracers registered the loop count enters linearly
    add!(model, Tracer(:tracer_a), Tracer(:tracer_b))
    @test count_tracers(model, [1.0]) == 2
    let dx = zeros(1)
        autodiff(set_runtime_activity(Reverse), count_tracers, Active, Const(model), Duplicated([2.0], dx))
        @test dx[1] ≈ 2                              # d/dx of 2x
    end

    # FORWARD mode over a POPULATED tracer registry is a known gap: it throws `EnzymeNoShadowError`.
    # this is just a temporary workaround
    @test_broken autodiff(
        set_runtime_activity(Forward), count_tracers, Duplicated,
        Const(model), Duplicated([2.0], [1.0]),
    )[1] ≈ 2
end
