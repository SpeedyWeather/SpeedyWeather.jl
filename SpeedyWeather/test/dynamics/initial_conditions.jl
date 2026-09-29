@testset "RandomVorticity" begin
    spectral_grid = SpectralGrid(nlayers = 1)

    # Test with default parameters
    ic = RandomVorticity(spectral_grid)
    model = BarotropicModel(spectral_grid)
    vars = SpeedyWeather.Variables(model)
    initialize!(vars, ic, model)
    vor = vars.prognostic.vorticity

    # Check vorticity is non-zero and has no NaNs
    @test !all(vor .== 0)
    @test !any(isnan.(vor))

    # Test different parameters produce different results``
    ic_custom = RandomVorticity(spectral_grid; power = -2, amplitude = 1.0e-5, max_wavenumber = 15)
    model_custom = BarotropicModel(spectral_grid)
    vars_custom = SpeedyWeather.Variables(model_custom)
    initialize!(vars_custom, ic_custom, model_custom)
    vor_custom = vars_custom.prognostic.vorticity

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
    vor = vars.prognostic.vorticity

    # Check vorticity is non-zero and has no NaNs
    @test !all(vor .== 0)
    @test !any(isnan.(vor))

    # Test different parameters produce different results
    ic_custom = RandomVelocity(spectral_grid; max_speed = 40, truncation = 10)
    model_custom = BarotropicModel(spectral_grid)
    vars_custom = Variables(model_custom)
    initialize!(vars_custom, ic_custom, model_custom)
    vor_custom = vars_custom.prognostic.vorticity

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
    vor = vars.prognostic.vorticity
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
    vor_custom = vars_custom.prognostic.vorticity

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
    vor = vars.prognostic.vorticity
    div = vars.prognostic.divergence

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
    vor_custom = vars_custom.prognostic.vorticity

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
    vor = vars.prognostic.vorticity
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
    vor_custom = vars_custom.prognostic.vorticity

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
    temp = vars.prognostic.temperature

    # Check temperature is initialized
    @test !all(temp .== 0)
    @test !any(isnan.(temp))

    # Test different parameters produce different results
    ic_custom = JablonowskiTemperature(spectral_grid; ΔT = 4.8e5)
    model_custom = PrimitiveWetModel(spectral_grid)
    vars_custom = Variables(model_custom)
    initialize!(vars_custom, ic_custom, model_custom)
    temp_custom = vars_custom.prognostic.temperature

    @test !any(isnan.(temp_custom))
    @test temp_custom != temp  # Different parameters should give different results
end

@testset "BalancedZonalState" begin
    spectral_grid = SpectralGrid(truncation = 31, nlayers = 8)
    model = PrimitiveWetModel(spectral_grid)
    @test model.initial_conditions.zonal isa BalancedZonalState    # default for PrimitiveWet

    ic = BalancedZonalState(spectral_grid)
    vars = Variables(model)
    initialize!(vars, ic, model)
    S = model.spectral_transform
    temp = transform(SpeedyWeather.get_step(vars.prognostic.temperature, 1), S)
    vor = transform(SpeedyWeather.get_step(vars.prognostic.vorticity, 1), S)
    @test !any(isnan.(temp))
    @test !any(isnan.(vor))

    (; R_dry) = model.atmosphere
    (; gravity, radius) = model.planet
    σ = model.geometry.σ_levels_full
    φ = model.geometry.latds

    # area-weighted global mean via zonal means and the quadrature weights of the rings
    weights = RingGrids.get_quadrature_weights(spectral_grid.grid)
    temp_zonal_mean = zonal_mean(temp)

    for k in 1:spectral_grid.nlayers
        # global mean temperature is the prescribed profile T̄(σ)
        T̄ = max(ic.T₀ * σ[k]^(R_dry * ic.lapse_rate / gravity), ic.Tmin)
        @test sum(temp_zonal_mean[:, k] .* weights) / sum(weights) ≈ T̄ rtol = 1.0e-3

        # vorticity of u = U(σ)sin²(2φ) in the southern hemisphere away from the perturbation
        U = ic.u₀ * tanh(-log(σ[k]) / ic.H)
        ij = findfirst(φ .< -30)
        ζ = -4U / radius * sind(φ[ij]) * cosd(φ[ij]) * (2 - 5sind(φ[ij])^2)
        @test vor[ij, k] ≈ ζ rtol = 0.02
    end

    # warmer tropics than poles near the surface
    tropics, poles = abs.(φ) .< 10, abs.(φ) .> 80
    @test sum(temp[tropics, end]) / count(tropics) > sum(temp[poles, end]) / count(poles) + 30
end

@testset "PressureOnOrography" begin
    spectral_grid = SpectralGrid(nlayers = 8)

    # Test with PrimitiveDry model
    ic = PressureOnOrography(spectral_grid)
    model = PrimitiveDryModel(spectral_grid)
    vars = Variables(model)
    initialize!(vars, ic, model)
    pres = vars.prognostic.pressure

    # Check pressure is initialized
    @test !all(pres .== 0)
    @test !any(isnan.(pres))

    # Test with PrimitiveWet model
    model_wet = PrimitiveWetModel(spectral_grid)
    vars_wet = Variables(model_wet)
    initialize!(vars_wet, ic, model_wet)
    pres_wet = vars_wet.prognostic.pressure

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
    pres = vars.prognostic.pressure

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
    humid = vars.prognostic.humidity

    # Check humidity is initialized
    @test !all(humid .== 0)
    @test !any(isnan.(humid))

    # Test different parameters produce different results
    ic_custom = (; temp = JablonowskiTemperature(spectral_grid), humid = ConstantRelativeHumidity(spectral_grid; relhumid_ref = 0.5))
    model_custom = PrimitiveWetModel(spectral_grid)
    vars_custom = Variables(model_custom)
    initialize!(vars_custom, ic_custom, model_custom)
    humid_custom = vars_custom.prognostic.humidity

    @test !any(isnan.(humid_custom))
    @test humid_custom != humid  # Different relative humidity should give different results
end

@testset "ConstantRelativeHumidity vertical profile" begin
    spectral_grid = SpectralGrid(truncation = 31, nlayers = 8)
    model = PrimitiveWetModel(spectral_grid)
    humid_ic = ConstantRelativeHumidity(spectral_grid; relhumid_ref = 0.7, σ_moist = 0.3, σ_dry = 0.1)
    ic = (; zonal = BalancedZonalState(spectral_grid), pres = ConstantPressure(spectral_grid), humid = humid_ic)
    vars = Variables(model)
    initialize!(vars, ic, model)

    S = model.spectral_transform
    temp = transform(SpeedyWeather.get_step(vars.prognostic.temperature, 1), S)
    humid = transform(SpeedyWeather.get_step(vars.prognostic.humidity, 1), S)
    pₛ = model.atmosphere.reference_pressure
    σ = model.geometry.σ_levels_full

    for k in 1:spectral_grid.nlayers
        RH = [humid[ij, k] / SpeedyWeather.saturation_humidity(temp[ij, k], σ[k] * pₛ, model.atmosphere) for ij in eachindex(temp[:, k])]
        RH_mean = sum(RH) / length(RH)
        if σ[k] >= 0.3
            @test RH_mean ≈ 0.7 rtol = 0.02
        elseif σ[k] <= 0.1
            @test maximum(abs, humid[:, k]) < 1.0e-6
        else    # linearly decreasing
            @test RH_mean ≈ 0.7 * (σ[k] - 0.1) / (0.3 - 0.1) rtol = 0.05
        end
    end
end

@testset "No initial precipitation burst" begin
    spectral_grid = SpectralGrid(truncation = 31, nlayers = 8)
    model = PrimitiveWetModel(spectral_grid)
    model.feedback.verbose = false
    simulation = initialize!(model)
    run!(simulation, period = Day(1))

    P = simulation.variables.parameterizations
    rain = P.rain_convection .+ P.rain_large_scale .+ P.snow_large_scale     # accumulated [m]
    weights = RingGrids.get_quadrature_weights(spectral_grid.grid)
    global_mean_rain = sum(zonal_mean(rain) .* weights) / sum(weights) * 1000  # [mm] in 1 day
    @test global_mean_rain < 3
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
    pres = vars.prognostic.pressure
    temp = vars.prognostic.temperature

    # Check pressure and temperature are initialized
    @test !all(iszero.(pres))
    @test !all(iszero.(temp))

    # Test with PrimitiveWet model
    model_wet = PrimitiveWetModel(spectral_grid)
    vars_wet = Variables(model_wet)
    initialize!(vars_wet, ic, model_wet)
    humid = vars_wet.prognostic.humidity

    @test !all(iszero.(humid))
end
