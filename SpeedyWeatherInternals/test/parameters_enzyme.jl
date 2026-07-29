# Tests for the Enzyme rules in SpeedyWeatherInternalsEnzymeExt.
#
# `reconstruct` cannot be differentiated in Enzyme FORWARD mode once `obj` is a large struct: it is
# `@generated` and expands to a `setproperties` over the whole type, which overruns Enzyme's
# type-analysis size budget (`EnzymeNoTypeError`). The extension provides an `EnzymeRules.forward`
# rule instead, which is exact because `reconstruct` does no arithmetic — it is a pure structural
# scatter (parameter fields from `values`, everything else from `obj`), hence linear, so its tangent
# is the same scatter over the tangents. Reverse mode already works and gets no rule.
using Enzyme
using EnzymeTestUtils
using FiniteDifferences
using SpeedyWeatherInternals.ParameterEditing
using SpeedyWeatherInternals.ParameterEditing: reconstruct

ParameterEditing.@parameterized @kwdef struct EnzymeTestStruct{T}
    @param a::T = 1.0
    @param b::T = 2.0
    c::T = 3.0                      # not a parameter: must carry no tangent from `values`
end

# y1 = a*x1², y2 = b*x2³ + c
function enzyme_testfunc!(y, x, p, obj::EnzymeTestStruct)
    new_obj = reconstruct(obj, p)
    y[1] = new_obj.a * x[1]^2
    y[2] = new_obj.b * x[2]^3 + new_obj.c
    return nothing
end

@testset "reconstruct Enzyme rules" begin
    obj = EnzymeTestStruct()
    x = 2 * ones(2)
    p = vec(parameters(obj))
    @test length(p) == 2

    @testset "reverse mode" begin
        y, dy = zeros(2), ones(2)
        dx, dp = zeros(2), zero(p)
        autodiff(
            Reverse, enzyme_testfunc!, Const,
            Duplicated(y, dy), Duplicated(x, dx), Duplicated(p, dp), Const(obj),
        )
        @test all(dp .== [x[1]^2, x[2]^3])
        @test all(dx .== [2 * p.a * x[1], 3 * p.b * x[2]^2])
    end

    @testset "forward mode" begin
        # analytic oracle: seeding both parameters with 1 gives dy = (da*x1², db*x2³) = (x1², x2³)
        let y = zeros(2), dy = zeros(2), dp = zero(p)
            dp .= 1
            autodiff(
                Forward, enzyme_testfunc!, Const,
                Duplicated(y, dy), Const(x), Duplicated(p, dp), Const(obj),
            )
            @test dy ≈ [x[1]^2, x[2]^3]
        end

        # a seed on one parameter must not leak into the other, and the non-parameter field `c`
        # contributes no tangent at all
        let y = zeros(2), dy = zeros(2), dp = zero(p)
            dp.a = 2.5
            autodiff(
                Forward, enzyme_testfunc!, Const,
                Duplicated(y, dy), Const(x), Duplicated(p, dp), Const(obj),
            )
            @test dy ≈ [2.5 * x[1]^2, 0]
        end

        # width > 1: batch members stay independent
        let y = zeros(2), d1 = zeros(2), d2 = zeros(2)
            p1 = zero(p); p1[1] = 1     # only `a`
            p2 = zero(p); p2[2] = 1     # only `b`
            autodiff(
                Forward, enzyme_testfunc!, Const,
                BatchDuplicated(y, (d1, d2)), Const(x), BatchDuplicated(p, (p1, p2)), Const(obj),
            )
            @test d1 ≈ [x[1]^2, 0]
            @test d2 ≈ [0, x[2]^3]
        end

        # against finite differences, for both activities of `obj`
        @testset "vs finite differences, obj::$Tobj" for Tobj in (Const, Duplicated)
            test_forward(
                reconstruct, Duplicated, (obj, Tobj), (p, Duplicated);
                fdm = FiniteDifferences.central_fdm(5, 1), rtol = 1.0e-6, atol = 1.0e-6,
                runtime_activity = true,
            )
        end

        # a nested struct, i.e. the recursive `reconstruct` path the generated method takes
        ParameterEditing.@parameterized @kwdef struct EnzymeNested{T}
            @param inner::T = EnzymeTestStruct() (group = :g,)
        end
        let nested = EnzymeNested(), q = vec(parameters(EnzymeNested()))
            readout(o, v) = reconstruct(o, v).inner.b
            dq = zero(q); dq.inner.b = 3.0
            @test autodiff(Forward, readout, Duplicated, Const(nested), Duplicated(q, dq))[1] ≈ 3.0
        end
    end
end
