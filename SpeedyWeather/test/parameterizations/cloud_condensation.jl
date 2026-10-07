using SpeedyWeather: get_prognostic_step, get_tendency_step, pressure, pressure_thickness,
    saturation_humidity, ice_fraction, mixed_phase_saturation_humidity, xu_randall_cloud_fraction,
    liquid_autoconversion, default_time_step, cloud_condensation!, large_scale_condensation!

# a component that only declares the fused cloud condensate (no tendencies)
struct CondensateOnly <: SpeedyWeather.AbstractParameterization end
SpeedyWeather.variables(::CondensateOnly, model::SpeedyWeather.AbstractModel) =
    SpeedyWeather.cloud_condensate_variables(SpeedyWeather.get_nsteps(model.time_stepping, model))

# a model with the prognostic cloud condensation and its variables, Float64 for tight budgets
function cloud_test_model(; NF = Float64, nlayers = 8, kwargs...)
    spectral_grid = SpectralGrid(; truncation = 15, nlayers, NF)
    large_scale_condensation = PrognosticCloudCondensation(spectral_grid; kwargs...)
    model = PrimitiveWetModel(spectral_grid; large_scale_condensation)
    model.feedback.verbose = false
    simulation = initialize!(model)
    return model, simulation.variables
end

# prescribe column ij (lagged step that the physics reads): temperature [K], relative humidity
# over liquid [1] and condensate [kg/kg], all per layer; zero the tendencies and precipitation
function set_cloud_column!(vars, model, ij; temperature, relative_humidity, condensate, surface_pressure = 1.0e5)
    scheme = model.large_scale_condensation
    (; time_stepping, geometry, atmosphere) = model
    T = get_prognostic_step(vars.grid.temperature, time_stepping, scheme)
    q = get_prognostic_step(vars.grid.humidity, time_stepping, scheme)
    qc = get_prognostic_step(vars.grid.cloud_condensate, time_stepping, scheme)
    vars.parameterizations.surface_pressure[ij] = surface_pressure
    for k in eachindex(temperature)
        p = pressure(k, surface_pressure, geometry.vertical_coordinates)
        T[ij, k] = temperature[k]
        q[ij, k] = relative_humidity[k] * saturation_humidity(temperature[k], p, atmosphere)
        qc[ij, k] = condensate[k]
    end
    for tend in (vars.tendencies.grid.temperature, vars.tendencies.grid.humidity, vars.tendencies.grid.cloud_condensate)
        get_tendency_step(tend, time_stepping, scheme)[ij, :] .= 0
    end
    vars.parameterizations.rain_rate[ij] = 0
    vars.parameterizations.snow_rate[ij] = 0
    vars.parameterizations.cloud_top[ij] = length(temperature) + 1
    return nothing
end

function run_cloud_column!(vars, model, ij)
    (; geometry, planet, atmosphere, land_sea_mask, time_stepping) = model
    cloud_condensation!(ij, vars, model.large_scale_condensation, geometry, planet, atmosphere, land_sea_mask, time_stepping)
    return nothing
end

# column tendencies of column ij as vectors
function cloud_column_tendencies(vars, model, ij)
    scheme = model.large_scale_condensation
    ts = model.time_stepping
    dT = get_tendency_step(vars.tendencies.grid.temperature, ts, scheme)[ij, :]
    dq = get_tendency_step(vars.tendencies.grid.humidity, ts, scheme)[ij, :]
    dqc = get_tendency_step(vars.tendencies.grid.cloud_condensate, ts, scheme)[ij, :]
    return dT, dq, dqc
end

@testset "Prognostic cloud condensation: variables fused with the atmosphere" begin
    model, vars = cloud_test_model(NF = Float32)
    nlayers = model.spectral_grid.nlayers

    # condensate is fused into the existing parents, aligned between spectral and grid
    @test haskey(vars.fused.prognostic.slot_map, :cloud_condensate)
    @test vars.fused.prognostic.slot_map.cloud_condensate == vars.fused.grid.slot_map.cloud_condensate
    @test vars.fused.spectral_tendencies.slot_map.cloud_condensate == vars.fused.grid_tendencies.slot_map.cloud_condensate
    @test vars.fused.spectral_tendencies.slot_map.uqc == vars.fused.grid_tendencies.slot_map.uqc
    @test length(vars.fused.prognostic.slot_map.cloud_condensate) == nlayers

    # grid copy with layers and both leapfrog steps, one tendency step
    @test size(vars.grid.cloud_condensate) == (model.spectral_grid.npoints, nlayers, 2)
    @test size(vars.tendencies.grid.cloud_condensate) == (model.spectral_grid.npoints, nlayers, 1)
    @test :cloud_condensate in SpeedyWeather.tendency_names(vars)

    # the transform can do the widest batch in one call
    @test size(parent(vars.fused.grid_tendencies).data, 2) <= SpeedyWeather.max_transform_batch(nlayers)

    # without the scheme there is no condensate
    model_default = PrimitiveWetModel(model.spectral_grid)
    @test !haskey(Variables(model_default).prognostic, :cloud_condensate)
end

@testset "Prognostic cloud condensation: thermodynamic helpers" begin
    model, vars = cloud_test_model()
    (; atmosphere) = model
    T₀ = atmosphere.temperature_freezing

    @test ice_fraction(280.0, T₀, 253.15) == 0
    @test ice_fraction(250.0, T₀, 253.15) == 1
    @test ice_fraction((T₀ + 253.15) / 2, T₀, 253.15) ≈ 0.5

    # saturation over ice equals liquid at freezing, is lower below, derivative by finite differences
    p = 50000.0
    q_liquid = saturation_humidity(T₀, p, atmosphere)
    @test mixed_phase_saturation_humidity(T₀, p, 1.0, atmosphere, true)[1] ≈ q_liquid
    @test mixed_phase_saturation_humidity(250.0, p, 1.0, atmosphere, true)[1] < saturation_humidity(250.0, p, atmosphere)
    @test mixed_phase_saturation_humidity(250.0, p, 1.0, atmosphere, false)[1] ≈ saturation_humidity(250.0, p, atmosphere)
    for f in (0.0, 0.3, 1.0)
        q_sat, dq_sat_dT = mixed_phase_saturation_humidity(250.0, p, f, atmosphere, true)
        δ = 1.0e-4
        fd = (
            mixed_phase_saturation_humidity(250.0 + δ, p, f, atmosphere, true)[1] -
                mixed_phase_saturation_humidity(250.0 - δ, p, f, atmosphere, true)[1]
        ) / 2δ
        @test dq_sat_dT ≈ fd rtol = 1.0e-6
    end

    # Xu-Randall: no cloud without condensate, more cloud with more condensate, never above 1
    q_sat = 0.005
    @test xu_randall_cloud_fraction(0.0, 0.9, q_sat, p, 2000.0, 0.001) == 0
    covers = [xu_randall_cloud_fraction(qc, 0.9, q_sat, p, 2000.0, 0.001) for qc in (1.0e-5, 1.0e-4, 1.0e-3)]
    @test issorted(covers)
    @test all(0 .< covers .<= 1)
    @test xu_randall_cloud_fraction(1.0e-3, 1.0, q_sat, p, 2000.0, 0.001) ≈ 1
end

@testset "Prognostic cloud condensation: processes on a column" begin
    model, vars = cloud_test_model()
    scheme = model.large_scale_condensation
    nlayers = model.spectral_grid.nlayers
    Δt_prognostic = default_time_step(model.time_stepping)
    @test Δt_prognostic == 2 * model.time_stepping.Δt
    ij = 1

    # supersaturated layers condense into the condensate, not into precipitation directly
    set_cloud_column!(
        vars, model, ij;
        temperature = fill(285.0, nlayers), relative_humidity = fill(1.1, nlayers), condensate = zeros(nlayers)
    )
    run_cloud_column!(vars, model, ij)
    dT, dq, dqc = cloud_column_tendencies(vars, model, ij)
    @test all(dqc .> 0)
    @test all(dq .< 0)
    @test all(dT .> 0)
    @test all(vars.parameterizations.cloud_fraction[ij, :] .> 0)
    @test vars.parameterizations.cloud_top[ij] == 1

    # condensate in subsaturated air evaporates, not more than available, cools
    qc = fill(1.0e-4, nlayers)
    set_cloud_column!(
        vars, model, ij;
        temperature = fill(285.0, nlayers), relative_humidity = fill(0.5, nlayers), condensate = qc
    )
    run_cloud_column!(vars, model, ij)
    dT, dq, dqc = cloud_column_tendencies(vars, model, ij)
    @test all(dqc .< 0)
    @test all(qc .+ Δt_prognostic .* dqc .>= -eps())
    @test all(dq .> 0)
    @test all(dT .< 0)

    # autoconversion at the threshold humidity (no condensation, no evaporation) in a warm top
    # layer with nothing above: the condensate decays by the Sundqvist rate, exponential in time
    set_cloud_column!(
        vars, model, ij;
        temperature = fill(285.0, nlayers), relative_humidity = fill(scheme.relative_humidity_threshold, nlayers),
        condensate = fill(1.0e-3, nlayers)
    )
    run_cloud_column!(vars, model, ij)
    dT, dq, dqc = cloud_column_tendencies(vars, model, ij)
    p₁ = pressure(1, 1.0e5, model.geometry.vertical_coordinates)
    cover = vars.parameterizations.cloud_fraction[ij, 1]
    expected = liquid_autoconversion(1.0e-3, cover, 0.0, 285.0, scheme.minimum_condensate * p₁ / 1.0e5, scheme, Δt_prognostic)
    @test dqc[1] ≈ -expected / Δt_prognostic rtol = 1.0e-6
    @test 0 < expected < 1.0e-3
    @test dq[1] ≈ 0 atol = 1.0e-15                         # no phase change of vapour
    @test vars.parameterizations.rain_rate_large_scale[ij] > 0

    # negative condensate (spectral transport) is filled from vapour with latent heating
    set_cloud_column!(
        vars, model, ij;
        temperature = fill(285.0, nlayers), relative_humidity = fill(0.9, nlayers), condensate = fill(-1.0e-5, nlayers)
    )
    run_cloud_column!(vars, model, ij)
    dT, dq, dqc = cloud_column_tendencies(vars, model, ij)
    @test all(abs.(-1.0e-5 .+ Δt_prognostic .* dqc) .< 1.0e-18)
    @test all(dq .< 0)
    @test all(dT .> 0)
end

@testset "Prognostic cloud condensation: non-negative condensate at extreme rates" begin
    model, vars = cloud_test_model(
        autoconversion_rate = 1.0e3, ice_autoconversion_rate = 1.0e3, evaporation_time_scale = 0.01
    )
    nlayers = model.spectral_grid.nlayers
    Δt_prognostic = default_time_step(model.time_stepping)
    ij = 2
    qc = [0.0, 1.0e-3, 2.0e-3, 5.0e-4, 1.0e-3, 3.0e-3, 1.0e-4, 1.0e-3]
    set_cloud_column!(
        vars, model, ij;
        temperature = [215.0, 230.0, 245.0, 258.0, 268.0, 280.0, 288.0, 295.0],
        relative_humidity = [0.2, 1.3, 0.4, 1.2, 0.3, 1.1, 0.5, 0.95], condensate = qc
    )
    run_cloud_column!(vars, model, ij)
    dT, dq, dqc = cloud_column_tendencies(vars, model, ij)
    @test all(qc .+ Δt_prognostic .* dqc .>= -1.0e-15)
    @test all(isfinite, dT)
end

@testset "Prognostic cloud condensation: water and enthalpy budgets" begin
    model, vars = cloud_test_model()
    scheme = model.large_scale_condensation
    (; geometry, planet, atmosphere, time_stepping) = model
    nlayers = model.spectral_grid.nlayers
    g = planet.gravity
    ρ = atmosphere.water_density
    cₚ = atmosphere.heat_capacity
    Lᵥ = atmosphere.latent_heat_condensation
    Lᵢ = atmosphere.latent_heat_fusion
    T₀ = atmosphere.temperature_freezing

    columns = (
        # snow forms in cold supersaturated layers aloft and melts below, negative condensate
        # is filled, rain reevaporates in a dry layer
        (
            temperature = [220.0, 235.0, 250.0, 262.0, 272.0, 282.0, 290.0, 295.0],
            relative_humidity = [0.5, 1.1, 1.2, 1.1, 0.9, 1.05, 0.6, 0.8],
            condensate = [0.0, 2.0e-5, 3.0e-4, -1.0e-5, 5.0e-4, 2.0e-4, 1.0e-4, 0.0],
        ),
        # cold column, snow reaches the surface
        (
            temperature = [215.0, 225.0, 235.0, 242.0, 248.0, 252.0, 255.0, 258.0],
            relative_humidity = [0.4, 1.2, 1.3, 1.1, 0.8, 1.05, 0.9, 0.7],
            condensate = [0.0, 1.0e-4, 2.0e-4, 1.0e-4, 0.0, 3.0e-4, 1.0e-4, 0.0],
        ),
    )

    @testset for (ij, column) in enumerate(columns)
        set_cloud_column!(vars, model, ij; column...)
        run_cloud_column!(vars, model, ij)
        dT, dq, dqc = cloud_column_tendencies(vars, model, ij)

        Δp = [pressure_thickness(k, 1.0e5, geometry.vertical_coordinates) for k in 1:nlayers]
        f_ice = ice_fraction.(column.temperature, T₀, scheme.ice_temperature)
        rain = ρ * vars.parameterizations.rain_rate_large_scale[ij]    # [kg/m²/s]
        snow = ρ * vars.parameterizations.snow_rate_large_scale[ij]
        @test rain + snow > 0

        # water: vapour and condensate lost by the column fall out as precipitation
        @test sum((dq .+ dqc) .* Δp) / g + rain + snow ≈ 0 atol = 1.0e-12

        # enthalpy: heating = latent heat of vapour lost + latent heat of fusion of the frozen water
        # gained (ice in the condensate at the layer's ice fraction, snow at the surface)
        heating = cₚ / g * sum(dT .* Δp)
        latent = Lᵥ * (-sum(dq .* Δp) / g) + Lᵢ * (sum(f_ice .* dqc .* Δp) / g + snow)
        @test heating ≈ latent rtol = 1.0e-10
    end
    @test ρ * vars.parameterizations.snow_rate_large_scale[2] > 0  # cold column: snow at the surface
end

@testset "Prognostic cloud condensation: limit of ImplicitCondensation" begin
    # with instantaneous autoconversion, no reevaporation and a warm column the condensate is
    # rained out in the same step: identical to ImplicitCondensation (without snow)
    model, vars = cloud_test_model(
        autoconversion_rate = 1.0e10, autoconversion_water = 1.0e-12, minimum_condensate = 0, reevaporation = 0
    )
    nlayers = model.spectral_grid.nlayers
    ij = 3
    column = (
        temperature = [280.0, 282.0, 285.0, 288.0, 291.0, 293.0, 296.0, 299.0],
        relative_humidity = [0.5, 1.1, 1.2, 0.8, 1.05, 0.9, 1.3, 0.7],
        condensate = zeros(nlayers),
    )
    set_cloud_column!(vars, model, ij; column...)
    run_cloud_column!(vars, model, ij)
    dT, dq, dqc = cloud_column_tendencies(vars, model, ij)
    rain = vars.parameterizations.rain_rate_large_scale[ij]

    implicit = ImplicitCondensation(model.spectral_grid; reevaporation = 0, snow = false)
    set_cloud_column!(vars, model, ij; column...)
    (; geometry, planet, atmosphere, time_stepping) = model
    large_scale_condensation!(ij, vars, implicit, geometry, planet, atmosphere, time_stepping)
    dT_implicit, dq_implicit, _ = cloud_column_tendencies(vars, model, ij)

    @test dqc ≈ zeros(nlayers) atol = 1.0e-15
    @test dq ≈ dq_implicit rtol = 1.0e-10
    @test dT ≈ dT_implicit rtol = 1.0e-10
    @test rain ≈ vars.parameterizations.rain_rate_large_scale[ij] rtol = 1.0e-10
end

@testset "Cloud radiation: PrognosticClouds and CloudyLongwaveTransmissivity" begin
    spectral_grid = SpectralGrid(truncation = 15, nlayers = 8)
    shortwave = OneBandShortwave(spectral_grid; clouds = PrognosticClouds(spectral_grid))
    longwave = OneBandLongwave(spectral_grid; transmissivity = CloudyLongwaveTransmissivity(spectral_grid))
    radiation = Radiation(spectral_grid; shortwave, longwave)
    model = PrimitiveWetModel(spectral_grid; radiation)     # no prognostic cloud scheme
    model.feedback.verbose = false
    simulation = initialize!(model)
    vars = simulation.variables
    P = vars.parameterizations
    ij = 1
    P.surface_pressure .= 1.0e5                 # only set from the spectral state when time stepping

    # without a cloud scheme the cloud state is zero: no clouds, clear-sky transmissivity
    clouds = SpeedyWeather.clouds!(ij, vars, model.radiation.shortwave.clouds, model)
    @test clouds.cloud_cover == 0
    @test clouds.cloud_top == spectral_grid.nlayers + 1
    t_clear = copy(SpeedyWeather.transmissivity!(ij, vars, model.radiation.longwave.transmissivity.clear_sky, model)[ij, :])
    t = SpeedyWeather.transmissivity!(ij, vars, model.radiation.longwave.transmissivity, model)[ij, :]
    @test t == t_clear

    # overlap: maximum for adjacent cloudy layers, random for separated ones
    P.cloud_liquid_water[ij, :] .= 0
    P.cloud_ice_water[ij, :] .= 0
    P.cloud_liquid_effective_radius[ij, :] .= 10.0e-6
    P.cloud_ice_effective_radius[ij, :] .= 50.0e-6
    P.cloud_fraction[ij, :] .= [0, 0, 0.5, 0.5, 0, 0, 0, 0]
    @test SpeedyWeather.clouds!(ij, vars, model.radiation.shortwave.clouds, model).cloud_cover ≈ 0.5
    P.cloud_fraction[ij, :] .= [0, 0, 0.5, 0, 0.5, 0, 0, 0]
    clouds = SpeedyWeather.clouds!(ij, vars, model.radiation.shortwave.clouds, model)
    @test clouds.cloud_cover ≈ 0.75
    @test clouds.cloud_top == 3

    # cloud albedo increases with liquid water, longwave transmissivity decreases in cloudy layers
    P.cloud_liquid_water[ij, 3] = 1.0e-5
    albedo_thin = SpeedyWeather.clouds!(ij, vars, model.radiation.shortwave.clouds, model).cloud_albedo
    P.cloud_liquid_water[ij, 3] = 1.0e-4
    albedo_thick = SpeedyWeather.clouds!(ij, vars, model.radiation.shortwave.clouds, model).cloud_albedo
    @test 0 < albedo_thin < albedo_thick < 1

    P.cloud_liquid_water[ij, 3] = 1.0e-3     # optically thick
    t = SpeedyWeather.transmissivity!(ij, vars, model.radiation.longwave.transmissivity, model)[ij, :]
    @test t[3] < t_clear[3]
    @test t[3] ≈ t_clear[3] * (1 - 0.5) rtol = 1.0e-3          # optically thick: emissivity ≈ 1
    @test t[[1, 2, 4, 6, 7, 8]] == t_clear[[1, 2, 4, 6, 7, 8]]
end

@testset "SPPT perturbs the cloud condensate tendency" begin
    spectral_grid = SpectralGrid(truncation = 15, nlayers = 8)
    random_process = SpectralAR1Process(spectral_grid, seed = 1)
    sppt = StochasticallyPerturbedParameterizationTendencies(spectral_grid)
    large_scale_condensation = PrognosticCloudCondensation(spectral_grid)
    model = PrimitiveWetModel(spectral_grid; random_process, stochastic_physics = sppt, large_scale_condensation)
    initialize!(model.stochastic_physics, model)
    vars = Variables(model)

    ij = 1
    vars.grid.random_pattern[ij] = 0.5      # perturbation factor 1 + 0.5 (no vertical tapering)
    vars.tendencies.grid.cloud_condensate[ij, :, 1] .= 1.0f-6
    vars.tendencies.grid.humidity[ij, :, 1] .= -1.0f-6
    SpeedyWeather.parameterization!(ij, vars, model.stochastic_physics, model)

    # same perturbation for humidity and condensate, so the water budget is kept
    @test all(vars.tendencies.grid.cloud_condensate[ij, :, 1] .≈ 1.5f-6)
    @test all(vars.tendencies.grid.humidity[ij, :, 1] .≈ -1.5f-6)
end

@testset "Prognostic clouds in the full model" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
    large_scale_condensation = PrognosticCloudCondensation(spectral_grid)
    shortwave = OneBandShortwave(spectral_grid; clouds = PrognosticClouds(spectral_grid))
    longwave = OneBandLongwave(spectral_grid; transmissivity = CloudyLongwaveTransmissivity(spectral_grid))
    radiation = Radiation(spectral_grid; shortwave, longwave)
    model = PrimitiveWetModel(spectral_grid; large_scale_condensation, radiation)
    model.feedback.verbose = false
    simulation = initialize!(model)
    run!(simulation, period = Day(2))

    vars = simulation.variables
    @test model.feedback.nans_detected == false
    @test any(vars.grid.cloud_condensate .> 0)
    @test all(0 .<= vars.parameterizations.cloud_fraction .<= 1)
    @test any(vars.parameterizations.cloud_fraction .> 0)
    @test all(0 .<= vars.parameterizations.cloud_cover .<= 1)
    @test all(vars.parameterizations.liquid_water_path .>= 0)
    @test all(vars.parameterizations.rain_large_scale .>= 0)
end

@testset "Fused condensate leaves the rest of the model unchanged" begin
    # with a component that only declares the fused condensate (zero tendencies) all other prognostic
    # variables agree with the model without it, to rounding (the batched transforms are wider)
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
    function run_steps(custom_parameterization)
        model = PrimitiveWetModel(spectral_grid; custom_parameterization)
        model.feedback.verbose = false
        simulation = initialize!(model)
        run!(simulation, steps = 10)
        return simulation.variables
    end
    vars_default = run_steps(nothing)
    vars_condensate = run_steps(CondensateOnly())

    @test all(iszero, vars_condensate.prognostic.cloud_condensate)
    for name in (:vorticity, :divergence, :temperature, :humidity, :pressure)
        a = getfield(vars_default.prognostic, name)
        b = getfield(vars_condensate.prognostic, name)
        @test a.data ≈ b.data rtol = 1.0e-5
    end
end
