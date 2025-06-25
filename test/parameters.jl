using ModelParameters
using SpeedyWeather: value, bounds, description, attributes, parameters, reconstruct
using DomainSets: Domain, RealLine

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
    earth = Earth(Float64)
    ps = parameters(earth, category="planet")
    @test values(stripparams(ps)) == (earth.rotation, earth.gravity, earth.axial_tilt, earth.solar_constant)
    @test all(map(p -> p.category == "planet", ps))
end

spectral_grid = SpectralGrid(trunc=31, nlayers=1)   # define resolution
model = Barotropic(spectral_grid)