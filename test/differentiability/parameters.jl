@testset "Differentiability of reconstruct" begin
    SpeedyWeather.@parameterized @kwdef struct TestStruct{T}
        @param a::T = 1.0
        @param b::T = 2.0
        c::T = 3.0
    end

    function testfunc!(y, x, p, obj::TestStruct)
        new_obj = SpeedyWeather.reconstruct(obj, p)
        y[1] = new_obj.a*x[1]^2
        y[2] = new_obj.b*x[2]^3 + new_obj.c
        return nothing
    end

    y, dy = zeros(2), ones(2)
    x, dx = 2*ones(2), zeros(2)
    obj = TestStruct()
    p, dp = vec(parameters(obj)), zero(vec(parameters(obj)))
    @test length(p) == length(dp) == 2
    autodiff(Reverse, testfunc!, Const, Duplicated(y, dy), Duplicated(x, dx), Duplicated(p, dp), Const(obj))
    # check derivatives
    @test all(dp .== [x[1]^2, x[2]^3])
    @test all(dx .== [2*p.a*x[1], 3*p.b*x[2]^2])
end
