using ComponentArrays: ComponentVector
using DomainSets: Domain, RealLine, UnitInterval
using SpeedyWeather: value, bounds, description, attributes, parameters, reconstruct, stripparams

import ModelParameters: ModelParameters, Model, Param, params, update

@testset "SpeedyParam" begin
    # Test basic construction and getter functions
    p = SpeedyParam(42.0, bounds=RealLine(), desc="Test parameter", unit="m", category="generic")
    @test value(p) == p.val == 42.0
    @test bounds(p) == p.bounds == RealLine()
    @test description(p) == p.desc == "Test parameter"
    @test attributes(p) == (unit="m", category="generic")

    # Test ModelParameters interface
    @test params(p) == (p,)
    @test Model(p)[:val] == (42.0,)
    @test value(update(p, (1.0,))) == 1.0
    @test stripparams(p) == 42.0
end

@testset "parameters" begin
    # check basic function of parameters method
    @test parameters(1.0) == SpeedyParam(1.0)
    @test parameters(Param, 1.0) == Param(1.0)
    @test parameters(SpeedyParam, 1.0, bounds=UnitInterval()) == SpeedyParam(1.0, bounds=UnitInterval())
    # test for non-nested type
    earth = Earth(Float32)
    ps = parameters(earth, category="planet")
    pvals = (earth.rotation, earth.gravity, earth.axial_tilt, earth.solar_constant)
    @test isa(ps, SpeedyParams) && SpeedyParams <: ModelParameters.AbstractModel
    @test values(stripparams(ps)) == pvals
    @test ps[:val] == pvals
    # this fails with an indexing error... looks like a bug in ModelParameters.jl;
    # @test all(map(p -> p.category == "planet", ps))
    @test all(map(p -> p.category == "planet", params(ps)))
    # test for nested model type
    spectral_grid = SpectralGrid(trunc=31, nlayers=1)   # define resolution
    model = BarotropicModel(spectral_grid)
    model_ps = parameters(model)
    new_model = reconstruct(model, 2*vec(model_ps))
    new_model_ps = parameters(new_model)
    @test all(vec(new_model_ps) .== 2*vec(model_ps))
    # test parameter subsets
    planet_ps = model_ps["planet"]
    @test isa(planet_ps, SpeedyParams)
    @test length(planet_ps) == length(parameters(earth))
    @test vec(planet_ps).planet == vec(parameters(earth))
    planet_hc_ps = model_ps[["planet","atmosphere.heat_capacity"]]
    @test isa(planet_hc_ps, SpeedyParams)
    @test length(planet_hc_ps) == length(parameters(earth))+1
    @test haskey(vec(planet_hc_ps), :planet) && haskey(vec(planet_hc_ps), :atmosphere)
    # test reconstruction from subset
    new_model_ps2 = 2*vec(planet_hc_ps)
    new_model2 = reconstruct(model, new_model_ps2)
    @test vec(parameters(new_model2.planet)) == new_model_ps2.planet
    @test vec(parameters(new_model2.atmosphere)).heat_capacity == new_model_ps2.atmosphere.heat_capacity
    ## check if functional constraints are preserved
    @test new_model2.atmosphere.κ ≈ model.atmosphere.κ / 2
end

@testset "@parameterized" begin
    # test single parameter, no kwdef
    SpeedyWeather.@parameterized struct TestType1{T}
        "non parameter"
        x::T
        "parameter"
        @param y::T
    end
    ps = parameters(TestType1(0.0,1.0))
    @test length(ps) == 1
    @test vec(ps).y == 1.0
    @test ps[:desc] == ("parameter",)

    # test two parameters, no kwdef
    SpeedyWeather.@parameterized struct TestType2{TX,TY,TZ}
        "parameter"
        @param x::TX
        "non-parameter"
        y::TY
        @param z::TZ
    end
    ps = parameters(TestType2(0.0,1.0,2.0))
    @test length(ps) == 2
    @test vec(ps) == ComponentVector(x=0.0, z=2.0)
    @test ps[:desc] == ("parameter","")

    # test one parameter, with kwdef
    SpeedyWeather.@parameterized @kwdef struct TestType3{T}
        "parameter"
        @param x::T = 1.0
        "non-parameter"
        y::T = 2.0
    end
    ps = parameters(TestType3())
    @test length(ps) == 1
    @test vec(ps) == ComponentVector(x=1.0)
    @test ps[:desc] == ("parameter",)

    # test multiple parameters, with kwdef
    SpeedyWeather.@parameterized @kwdef struct TestType4{T1,T2}
        "parameter 1"
        @param x::T1 = 1.0
        "parameter 2"
        @param y::T1 = 2.0
        @param z::T2 = 3.0
        date::DateTime = DateTime(0)
    end
    ps = parameters(TestType4())
    @test length(ps) == 3
    @test vec(ps) == ComponentVector(x=1.0,y=2.0,z=3.0)
    @test ps[:desc] == ("parameter 1","parameter 2","")

    # test parameters for nested type
    SpeedyWeather.@parameterized @kwdef struct MyModel{T}
        @param component::T = TestType4() (group=:group1,)
    end
    ps = parameters(MyModel())
    @test length(ps) == 3
    @test haskey(parent(ps), :component)
    @test vec(ps) == ComponentVector(component=(x=1.0,y=2.0,z=3.0))
    @test all(map(==(:group1), ps[:group]))
end
