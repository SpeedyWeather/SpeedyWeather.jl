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

    spectral_grid = SpectralGrid(; trunc = 8, nlayers, NF = Float64)
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
