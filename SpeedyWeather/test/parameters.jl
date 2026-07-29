using ComponentArrays: ComponentVector
using DomainSets: Domain, RealLine, UnitInterval
using SpeedyWeather.SpeedyWeatherInternals.ParameterEditing
using SpeedyWeather.SpeedyWeatherInternals.ParameterEditing: NumberParam

import ModelParameters: ModelParameters, Model, Param, params, update

@testset "parameters" begin
    # check basic function of parameters method
    @test parameters(1.0) == NumberParam(1.0)
    @test parameters(Param, 1.0) == Param(1.0)
    @test parameters(NumberParam, 1.0, bounds = UnitInterval()) == NumberParam(1.0, bounds = UnitInterval())
    # test for non-nested type
    earth = Earth(Float32)
    ps = parameters(earth, category = "planet")
    pvals = (earth.radius, earth.rotation, earth.gravity, earth.axial_tilt, earth.solar_constant)
    @test isa(ps, ParameterTable) && ParameterTable <: ModelParameters.AbstractModel
    @test values(stripparams(ps)) == pvals
    @test ps[:val] == pvals
    # this fails with an indexing error... looks like a bug in ModelParameters.jl;
    # @test all(map(p -> p.category == "planet", ps))
    @test all(map(p -> p.category == "planet", params(ps)))
    # test for nested model type
    spectral_grid = SpectralGrid(trunc = 31, nlayers = 1)   # define resolution
    model = BarotropicModel(spectral_grid)
    model_ps = parameters(model)
    new_ps = 2 * vec(model_ps)
    new_model = @inferred reconstruct(model, new_ps)
    new_model_ps = parameters(new_model)
    @test all(vec(new_model_ps) .== 2 * vec(model_ps))
    # test parameter subsets
    planet_ps = model_ps["planet"]
    @test isa(planet_ps, ParameterTable)
    @test length(planet_ps) == length(parameters(earth))
    @test vec(planet_ps).planet == vec(parameters(earth))
    planet_hc_ps = model_ps[["planet", "atmosphere.heat_capacity"]]
    @test isa(planet_hc_ps, ParameterTable)
    @test length(planet_hc_ps) == length(parameters(earth)) + 1
    @test haskey(vec(planet_hc_ps), :planet) && haskey(vec(planet_hc_ps), :atmosphere)
    # test reconstruction from subset
    new_model_ps2 = 2 * vec(planet_hc_ps)
    new_model2 = @inferred reconstruct(model, new_model_ps2)
    @test vec(parameters(new_model2.planet)) == new_model_ps2.planet
    @test vec(parameters(new_model2.atmosphere)).heat_capacity == new_model_ps2.atmosphere.heat_capacity
    ## κ is a parameter, so it gets updated like any other parameter in reconstruction
    ## Partial reconstruction: κ not in patch, so it retains original value
    @test new_model2.atmosphere.κ ≈ model.atmosphere.κ  # κ unchanged in partial reconstruction
    # Full reconstruction: κ is in the parameter vector, so it gets doubled like other parameters
    new_model_ps_full = 2 * vec(parameters(model))
    new_model_full = @inferred reconstruct(model, new_model_ps_full)
    @test new_model_full.atmosphere.κ ≈ 2 * model.atmosphere.κ  # κ doubled in full reconstruction
    # regression test for 20aa240611:
    # we set all parameters to identical values and check that parameter subsets still work
    model_ps_zero = parameters(reconstruct(model, zero(vec(model_ps))))
    model_ps_zero_subset = model_ps_zero[["planet"]]
    @test all(==(:planet), model_ps_zero_subset[:component])
    @test all(iszero, vec(model_ps_zero_subset))
    @test length(model_ps_zero_subset) == length(vec(model_ps).planet)
end
