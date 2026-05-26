using ComponentArrays: ComponentVector
using Dates
using DomainSets: Domain, RealLine, UnitInterval
using SpeedyWeatherInternals.ParameterEditing
using SpeedyWeatherInternals.ParameterEditing: NumberParam

import ModelParameters: ModelParameters, Model, Param, params, update

@testset "NumberParam" begin
    # Test basic construction and getter functions
    p = NumberParam(42.0, bounds = RealLine(), desc = "Test parameter", unit = "m", category = "generic")
    @test ParameterEditing.value(p) == p.val == 42.0
    @test bounds(p) == p.bounds == RealLine()
    @test description(p) == p.desc == "Test parameter"
    @test attributes(p) == (unit = "m", category = "generic")

    # Test ModelParameters interface
    @test params(p) == (p,)
    @test Model(p)[:val] == (42.0,)
    @test ParameterEditing.value(update(p, (1.0,))) == 1.0
    @test stripparams(p) == 42.0
end

@testset "@parameterized" begin
    # test single parameter, no kwdef
    ParameterEditing.@parameterized struct TestType1{T}
        "non parameter"
        x::T
        "parameter"
        @param y::T
    end
    ps = parameters(TestType1(0.0, 1.0))
    @test length(ps) == 1
    @test vec(ps).y == 1.0
    @test ps[:desc] == ("parameter",)

    # test two parameters, no kwdef, also with docstring
    """Docstring for `TestType2`"""
    ParameterEditing.@parameterized struct TestType2{TX, TY, TZ}
        "parameter"
        @param x::TX
        "non-parameter"
        y::TY
        @param z::TZ
    end
    ps = parameters(TestType2(0.0, 1.0, 2.0))
    @test length(ps) == 2
    @test vec(ps) == ComponentVector(x = 0.0, z = 2.0)
    @test ps[:desc] == ("parameter", "")

    # test one parameter, with kwdef and docstrin
    """Docstring for `TestType3`"""
    ParameterEditing.@parameterized @kwdef struct TestType3{T}
        "parameter"
        @param x::T = 1.0
        "non-parameter"
        y::T = 2.0
    end
    ps = parameters(TestType3())
    @test length(ps) == 1
    @test vec(ps) == ComponentVector(x = 1.0)
    @test ps[:desc] == ("parameter",)

    # test multiple parameters, with kwdef and docstring
    """Docstring for `TestType4`"""
    ParameterEditing.@parameterized @kwdef struct TestType4{T1, T2}
        "parameter 1"
        @param x::T1 = 1.0
        "parameter 2"
        @param y::T1 = 2.0
        @param z::T2 = 3.0
        date::DateTime = DateTime(0)
    end
    ps = parameters(TestType4())
    @test length(ps) == 3
    @test vec(ps) == ComponentVector(x = 1.0, y = 2.0, z = 3.0)
    @test ps[:desc] == ("parameter 1", "parameter 2", "")

    # test multiple parameters with heterogeneous attributes
    """Docstring for `TestType5`"""
    ParameterEditing.@parameterized @kwdef struct TestType5{T1, T2}
        @param x::T1 = 1.0 (a = 1,)
        @param y::T1 = 2.0 (b = 2,)
        @param z::T2 = 3.0 (b = 3,)
        date::DateTime = DateTime(0)
    end
    ps = parameters(TestType5());
    @test haskey(ps, :a)
    @test haskey(ps, :b)
    @test !haskey(ps, :c)
    @test length(ps) == 3
    @test vec(ps) == ComponentVector(x = 1.0, y = 2.0, z = 3.0)

    # test parameters for nested type
    ParameterEditing.@parameterized @kwdef struct MyModel{T}
        @param component::T = TestType4() (group = :group1,)
    end
    ps = parameters(MyModel())
    @test length(ps) == 3
    @test haskey(parent(ps), :component)
    @test vec(ps) == ComponentVector(component = (x = 1.0, y = 2.0, z = 3.0))
    @test all(map(==(:group1), ps[:group]))
end
