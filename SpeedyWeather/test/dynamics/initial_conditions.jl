@testset "RandomVorticity" begin
    spectral_grid = SpectralGrid(nlayers = 1)

    # Test with default parameters
    ic = RandomVorticity(spectral_grid)
    model = BarotropicModel(spectral_grid)
    vars = SpeedyWeather.Variables(model)
    initialize!(vars, ic, model)
    vor = vars.prognostic.vor

    # Check vorticity is non-zero and has no NaNs
    @test !all(vor .== 0)
    @test !any(isnan.(vor))

    # Test different parameters produce different results``
    ic_custom = RandomVorticity(spectral_grid; power = -2, amplitude = 1.0e-5, max_wavenumber = 15)
    model_custom = BarotropicModel(spectral_grid)
    vars_custom = SpeedyWeather.Variables(model_custom)
    initialize!(vars_custom, ic_custom, model_custom)
    vor_custom = vars_custom.prognostic.vor

    @test !any(isnan.(vor_custom))
    @test vor_custom != vor  # Different parameters should give different results
end

@testset "RandomVelocity" begin
    spectral_grid = SpectralGrid(nlayers = 1)

    # Test with default parameters
    ic = RandomVelocity(spectral_grid)
    model = BarotropicModel(spectral_grid)
    vars = Variables(model)
    initialize!(vars, ic, model)
    vor = vars.prognostic.vor

    # Check vorticity is non-zero and has no NaNs
    @test !all(vor .== 0)
    @test !any(isnan.(vor))

    # Test different parameters produce different results
    ic_custom = RandomVelocity(spectral_grid; max_speed = 40, truncation = 10)
    model_custom = BarotropicModel(spectral_grid)
    vars_custom = Variables(model_custom)
    initialize!(vars_custom, ic_custom, model_custom)
    vor_custom = vars_custom.prognostic.vor

    @test !any(isnan.(vor_custom))
    @test vor_custom != vor  # Different parameters should give different results
end

@testset "ZonalJet" begin
    spectral_grid = SpectralGrid(nlayers = 1)

    # Test with default parameters (Galewsky jet)
    ic = ZonalJet(spectral_grid)
    model = ShallowWaterModel(spectral_grid)
    vars = Variables(model)
    initialize!(vars, ic, model)
    vor = vars.prognostic.vor
    η = vars.prognostic.η

    # Check vorticity and pressure are initialized
    @test !all(vor .== 0)
    @test !any(isnan.(vor))
    @test !all(η .== 0)
    @test !any(isnan.(η))

    # Test different parameters produce different results
    ic_custom = ZonalJet(spectral_grid; latitude = 30, width = 15, umax = 60)
    model_custom = ShallowWaterModel(spectral_grid)
    vars_custom = Variables(model_custom)
    initialize!(vars_custom, ic_custom, model_custom)
    vor_custom = vars_custom.prognostic.vor

    @test !any(isnan.(vor_custom))
    @test vor_custom != vor  # Different parameters should give different results
end

@testset "ZonalWind" begin
    spectral_grid = SpectralGrid(nlayers = 8)

    # Test with default parameters (Jablonowski)
    ic = ZonalWind(spectral_grid)
    model = PrimitiveDryModel(spectral_grid)
    vars = Variables(model)
    initialize!(vars, ic, model)
    vor = vars.prognostic.vor
    div = vars.prognostic.div

    # Check vorticity and divergence are initialized
    @test !all(vor .== 0)
    @test !any(isnan.(vor))
    @test !all(div .== 0)
    @test !any(isnan.(div))

    # Test different parameters produce different results
    ic_custom = ZonalWind(spectral_grid; u₀ = 40, perturb_uₚ = 2)
    model_custom = PrimitiveDryModel(spectral_grid)
    vars_custom = Variables(model_custom)
    initialize!(vars_custom, ic_custom, model_custom)
    vor_custom = vars_custom.prognostic.vor

    @test !any(isnan.(vor_custom))
    @test vor_custom != vor  # Different parameters should give different results
end

@testset "RossbyHaurwitzWave" begin
    spectral_grid = SpectralGrid(nlayers = 1)

    # Test with ShallowWater model
    ic = RossbyHaurwitzWave(spectral_grid)
    model = ShallowWaterModel(spectral_grid)
    vars = Variables(model)
    initialize!(vars, ic, model)
    vor = vars.prognostic.vor
    η = vars.prognostic.η

    # Check vorticity and pressure are initialized
    @test !all(vor .== 0)
    @test !any(isnan.(vor))
    @test !all(η .== 0)
    @test !any(isnan.(η))

    # Test different parameters produce different results
    ic_custom = RossbyHaurwitzWave(spectral_grid; m = 6)
    model_custom = ShallowWaterModel(spectral_grid)
    vars_custom = Variables(model_custom)
    initialize!(vars_custom, ic_custom, model_custom)
    vor_custom = vars_custom.prognostic.vor

    @test !any(isnan.(vor_custom))
    @test vor_custom != vor  # Different wavenumber should give different results
end

@testset "JablonowskiTemperature" begin
    spectral_grid = SpectralGrid(nlayers = 8)

    # Test with default parameters
    ic = JablonowskiTemperature(spectral_grid)
    model = PrimitiveWetModel(spectral_grid)
    vars = Variables(model)
    initialize!(vars, ic, model)
    temp = vars.prognostic.temp

    # Check temperature is initialized
    @test !all(temp .== 0)
    @test !any(isnan.(temp))

    # Test different parameters produce different results
    ic_custom = JablonowskiTemperature(spectral_grid; ΔT = 4.8e5)
    model_custom = PrimitiveWetModel(spectral_grid)
    vars_custom = Variables(model_custom)
    initialize!(vars_custom, ic_custom, model_custom)
    temp_custom = vars_custom.prognostic.temp

    @test !any(isnan.(temp_custom))
    @test temp_custom != temp  # Different parameters should give different results
end

@testset "PressureOnOrography" begin
    spectral_grid = SpectralGrid(nlayers = 8)

    # Test with PrimitiveDry model
    ic = PressureOnOrography(spectral_grid)
    model = PrimitiveDryModel(spectral_grid)
    vars = Variables(model)
    initialize!(vars, ic, model)
    pres = vars.prognostic.pres

    # Check pressure is initialized
    @test !all(pres .== 0)
    @test !any(isnan.(pres))

    # Test with PrimitiveWet model
    model_wet = PrimitiveWetModel(spectral_grid)
    vars_wet = Variables(model_wet)
    initialize!(vars_wet, ic, model_wet)
    pres_wet = vars_wet.prognostic.pres

    @test !all(pres_wet .== 0)
    @test !any(isnan.(pres_wet))

    # Different models should give same results for pressure
    @test pres_wet == pres
end

@testset "ConstantPressure" begin
    spectral_grid = SpectralGrid(nlayers = 8)

    # Test with PrimitiveDry model
    ic = ConstantPressure(spectral_grid)
    model = PrimitiveDryModel(spectral_grid)
    vars = Variables(model)
    initialize!(vars, ic, model)
    pres = vars.prognostic.pres

    # Check pressure is initialized to constant value
    @test !any(isnan.(pres))

    # Test with ShallowWater model (should do nothing)
    spectral_grid_sw = SpectralGrid(nlayers = 1)
    ic_sw = ConstantPressure(spectral_grid_sw)
    model_sw = ShallowWaterModel(spectral_grid_sw)
    vars_sw = Variables(model_sw)
    initialize!(vars_sw, ic_sw, model_sw)
    η = vars_sw.prognostic.η

    @test all(η .≈ 0)
    @test !any(isnan.(η))
end

@testset "ConstantRelativeHumidity" begin
    spectral_grid = SpectralGrid(nlayers = 8)

    # Test with default relative humidity
    # we need a nonzero temperature to set the humidity, so we use the Jablonowski temperature
    # important to have temp therefore first in the named tuple!
    ic = (; temp = JablonowskiTemperature(spectral_grid), humid = ConstantRelativeHumidity(spectral_grid))
    model = PrimitiveWetModel(spectral_grid)
    vars = Variables(model)
    initialize!(vars, ic, model)
    humid = vars.prognostic.humid

    # Check humidity is initialized
    @test !all(humid .== 0)
    @test !any(isnan.(humid))

    # Test different parameters produce different results
    ic_custom = (; temp = JablonowskiTemperature(spectral_grid), humid = ConstantRelativeHumidity(spectral_grid; relhumid_ref = 0.5))
    model_custom = PrimitiveWetModel(spectral_grid)
    vars_custom = Variables(model_custom)
    initialize!(vars_custom, ic_custom, model_custom)
    humid_custom = vars_custom.prognostic.humid

    @test !any(isnan.(humid_custom))
    @test humid_custom != humid  # Different relative humidity should give different results
end

@testset "RandomWaves" begin
    spectral_grid = SpectralGrid(nlayers = 1)

    # Test with default parameters
    ic = RandomWaves(spectral_grid)
    model = ShallowWaterModel(spectral_grid)
    vars = Variables(model)
    initialize!(vars, ic, model)
    η = vars.prognostic.η

    # Check pressure is initialized
    @test !all(η .== 0)
    @test !any(isnan.(η))

    # Test different parameters produce different results
    ic_custom = RandomWaves(spectral_grid; amplitude = 1000, lmin = 5, lmax = 20)
    model_custom = ShallowWaterModel(spectral_grid)
    vars_custom = Variables(model_custom)
    initialize!(vars_custom, ic_custom, model_custom)
    η_custom = vars_custom.prognostic.η

    @test !any(isnan.(η_custom))
    @test η_custom != η   # Different parameters should give different results
end

@testset "StartFromRest" begin
    spectral_grid = SpectralGrid(nlayers = 8)

    # Test with PrimitiveDry model
    ic = StartFromRest(spectral_grid)
    model = PrimitiveDryModel(spectral_grid)
    vars = Variables(model)
    initialize!(vars, ic, model)
    pres = vars.prognostic.pres
    temp = vars.prognostic.temp

    # Check pressure and temperature are initialized
    @test !all(iszero.(pres))
    @test !all(iszero.(temp))

    # Test with PrimitiveWet model
    model_wet = PrimitiveWetModel(spectral_grid)
    vars_wet = Variables(model_wet)
    initialize!(vars_wet, ic, model_wet)
    humid = vars_wet.prognostic.humid

    @test !all(iszero.(humid))
end
