# Both schemes construct, models build with them, one column update and a few time steps
# give finite, physically signed output.

spectral_grid = SpectralGrid(truncation = 16, nlayers = 8)
NF = spectral_grid.NF
nlayers, npoints = spectral_grid.nlayers, spectral_grid.npoints

# idealised column state: temperature increasing towards the surface (k = 1 is the top),
# moist, warm ocean and land surfaces, sun at 60° zenith angle
function set_test_state!(variables, model)
    variables.parameterizations.surface_pressure .= 100000        # [Pa]
    for k in 1:model.spectral_grid.nlayers
        variables.grid.temperature[:, k, :] .= 220 + 9 * (k - 1)   # all time steps
        variables.grid.humidity[:, k, :] .= 0.005
    end
    variables.prognostic.ocean.sea_surface_temperature .= 295
    variables.prognostic.land.soil_temperature .= 285
    variables.parameterizations.ocean.albedo .= 0.06
    variables.parameterizations.land.albedo .= 0.3
    variables.parameterizations.cos_zenith .= 0.5
    variables.tendencies.grid.temperature .= 0
    return nothing
end

@testset "Extension is active" begin
    @test Base.get_extension(SpeedyWeather, :SpeedyWeatherNumericalRadiationExt) !== nothing
end

@testset "AnalyticBandLongwave" begin
    longwave = AnalyticBandLongwave(spectral_grid)
    @test longwave isa NumericalRadiation.AnalyticBandLongwave{NF}
    @test AnalyticBandLongwave(spectral_grid; diffusivity = 1.5).diffusivity == NF(1.5)

    # longwave only, no other parameterizations
    radiation = Radiation(spectral_grid; shortwave = nothing, longwave)
    model = PrimitiveWetModel(spectral_grid; radiation, parameterizations = (:radiation,))
    initialize!(model.radiation, model)
    variables = Variables(model)
    set_test_state!(variables, model)
    SpeedyWeather.column_parameterizations!(variables, model)

    @test all(isfinite, variables.tendencies.grid.temperature)
    @test any(!=(zero(NF)), variables.tendencies.grid.temperature)
    @test all(>(zero(NF)), variables.parameterizations.outgoing_longwave)

    # without a CO₂ component the scheme uses 280 ppm: same result as prescribing it
    olr = copy(variables.parameterizations.outgoing_longwave)
    model_co2 = PrimitiveWetModel(spectral_grid; radiation, parameterizations = (:radiation,),
                                  greenhouse_gases = (; co2 = CO2(spectral_grid, 280)))
    initialize!(model_co2.radiation, model_co2)
    variables_co2 = Variables(model_co2)
    set_test_state!(variables_co2, model_co2)
    variables_co2.prognostic.greenhouse_gases.co2[] = 280
    SpeedyWeather.column_parameterizations!(variables_co2, model_co2)
    @test variables_co2.parameterizations.outgoing_longwave == olr
end

@testset "ClearSkyEcCKDRadiation" begin
    radiation = ClearSkyEcCKDRadiation(spectral_grid)     # the reference 32x32 tables (ecrad_data artifact)
    @test radiation isa NumericalRadiation.ClearSkyEcCKDRadiation{NF}
    @test eltype(radiation.gas_optics) === NF
    @test radiation.surface_emissivity == NF(0.98)
    @test ClearSkyEcCKDRadiation(spectral_grid; mole_fractions = (; co2 = 400e-6)).mole_fractions.co2 == NF(400e-6)

    model = PrimitiveWetModel(spectral_grid; radiation, parameterizations = (:radiation,))
    initialize!(model.radiation, model)
    variables = Variables(model)
    W = variables.parameterizations.ecckd
    @test size(W.longwave_optical_depth) == (npoints, nlayers, 32)
    @test size(W.shortwave_optical_depth) == (npoints, nlayers, 32)
    @test size(W.longwave_up) == (npoints, nlayers + 1)

    set_test_state!(variables, model)
    SpeedyWeather.column_parameterizations!(variables, model)

    P = variables.parameterizations
    @test all(isfinite, variables.tendencies.grid.temperature)
    @test all(>(zero(NF)), P.outgoing_longwave)
    @test all(>(zero(NF)), P.outgoing_shortwave)                         # Rayleigh + surface reflection
    @test all(P.surface_shortwave_down .< model.planet.solar_constant * 0.5)
    @test all(P.outgoing_shortwave .< P.surface_shortwave_down)
end

@testset "Full model time steps" begin
    model = PrimitiveWetModel(spectral_grid; radiation = ClearSkyEcCKDRadiation(spectral_grid))
    simulation = initialize!(model)
    run!(simulation, steps = 2)
    @test all(isfinite, simulation.variables.parameterizations.outgoing_longwave)
    @test all(isfinite, simulation.variables.prognostic.temperature)
end
