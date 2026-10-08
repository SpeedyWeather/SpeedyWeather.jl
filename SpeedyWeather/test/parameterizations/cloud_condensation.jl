using SpeedyWeather: get_prognostic_step, get_tendency_step, pressure, pressure_thickness,
    saturation_humidity, ice_fraction, mixed_phase_saturation_humidity, xu_randall_cloud_fraction,
    liquid_autoconversion, default_time_step, cloud_condensation!, large_scale_condensation!

# a component that only declares the fused cloud condensate (no tendencies)
struct CondensateOnly <: SpeedyWeather.AbstractParameterization end
SpeedyWeather.variables(::CondensateOnly, model::SpeedyWeather.AbstractModel) =
    SpeedyWeather.cloud_condensate_variables(SpeedyWeather.get_nsteps(model.time_stepping, model))

# a model with the prognostic cloud condensation and its variables, Float64 for tight budgets
function cloud_test_model(; NF = Float64, nlayers = 8, closure = RelaxationClosure, closure_kwargs = (;), kwargs...)
    spectral_grid = SpectralGrid(; truncation = 15, nlayers, NF)
    large_scale_condensation = PrognosticCloudCondensation(spectral_grid; closure = closure(spectral_grid; closure_kwargs...), kwargs...)
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
        temperature = fill(285.0, nlayers), relative_humidity = fill(scheme.closure.relative_humidity_threshold, nlayers),
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
        autoconversion_rate = 1.0e3, ice_autoconversion_rate = 1.0e3, closure_kwargs = (; evaporation_time_scale = 0.01)
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

# column water and enthalpy budgets of column ij: (water imbalance [kg/m²/s], heating [W/m²], latent heat [W/m²])
function cloud_column_budgets(vars, model, ij, temperature)
    scheme = model.large_scale_condensation
    (; geometry, planet, atmosphere) = model
    nlayers = model.spectral_grid.nlayers
    g = planet.gravity
    ρ = atmosphere.water_density
    dT, dq, dqc = cloud_column_tendencies(vars, model, ij)
    Δp = [pressure_thickness(k, vars.parameterizations.surface_pressure[ij], geometry.vertical_coordinates) for k in 1:nlayers]
    f_ice = ice_fraction.(temperature, atmosphere.temperature_freezing, scheme.ice_temperature)
    rain = ρ * vars.parameterizations.rain_rate_large_scale[ij]
    snow = ρ * vars.parameterizations.snow_rate_large_scale[ij]
    water = sum((dq .+ dqc) .* Δp) / g + rain + snow
    heating = atmosphere.heat_capacity / g * sum(dT .* Δp)
    latent = atmosphere.latent_heat_condensation * (-sum(dq .* Δp) / g) +
        atmosphere.latent_heat_fusion * (sum(f_ice .* dqc .* Δp) / g + snow)
    return water, heating, latent
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

@testset "Sundqvist closure" begin
    using SpeedyWeather: critical_relative_humidity, sundqvist_cloud_fraction
    model, vars = cloud_test_model(closure = SundqvistClosure, autoconversion_rate = 0, ice_autoconversion_rate = 0)
    scheme = model.large_scale_condensation
    closure = scheme.closure
    (; atmosphere, geometry) = model
    nlayers = model.spectral_grid.nlayers
    Δt_prognostic = default_time_step(model.time_stepping)
    references = vars.prognostic.clouds
    ij = 4

    # critical relative humidity and cloud fraction
    @test critical_relative_humidity(closure, 1.0e5, 1.0e5) ≈ closure.critical_relative_humidity_surface
    closure_profile = SundqvistClosure(model.spectral_grid; critical_relative_humidity_surface = 0.95, critical_relative_humidity_top = 0.7)
    @test critical_relative_humidity(closure_profile, 1.0e5, 1.0e5) ≈ 0.95
    @test critical_relative_humidity(closure_profile, 1.0e3, 1.0e5) ≈ 0.7 atol = 1.0e-6
    @test sundqvist_cloud_fraction(0.8, 0.9) == 0
    @test sundqvist_cloud_fraction(0.95, 0.9) ≈ 1 - sqrt(0.5)
    @test sundqvist_cloud_fraction(1.2, 0.9) ≈ 1 atol = 1.0e-2

    # no reference state yet: no supply, no condensation; condensate in a clear cell evaporates
    @test all(iszero, references.surface_pressure_reference)
    qc = fill(1.0e-4, nlayers)
    set_cloud_column!(vars, model, ij; temperature = fill(285.0, nlayers), relative_humidity = fill(0.95, nlayers), condensate = zeros(nlayers))
    run_cloud_column!(vars, model, ij)
    dT, dq, dqc = cloud_column_tendencies(vars, model, ij)
    @test all(dqc .== 0)
    set_cloud_column!(vars, model, ij; temperature = fill(285.0, nlayers), relative_humidity = fill(0.5, nlayers), condensate = qc)
    run_cloud_column!(vars, model, ij)
    dT, dq, dqc = cloud_column_tendencies(vars, model, ij)
    @test all(dqc .< 0)
    @test all(qc .+ Δt_prognostic .* dqc .>= -eps())

    # references: the last call's state moves to "two calls before", this call's is stored
    @test all(references.surface_pressure_reference[ij, :] .== 1.0e5)      # two calls, both stored
    T = get_prognostic_step(vars.grid.temperature, model.time_stepping, scheme)
    q = get_prognostic_step(vars.grid.humidity, model.time_stepping, scheme)
    @test references.temperature_reference[ij, :, 2] ≈ T[ij, :] .+ Δt_prognostic .* atmosphere.latent_heat_condensation ./ atmosphere.heat_capacity .* (-dq)
    @test references.humidity_reference[ij, :, 2] ≈ q[ij, :] .+ Δt_prognostic .* dq

    # a moisture supply in a partly cloudy cell (RH 0.95 > u = 0.9) condenses β M / (1 + L/cₚ f ∂q*/∂T),
    # β = b without condensate
    set_cloud_column!(vars, model, ij; temperature = fill(285.0, nlayers), relative_humidity = fill(0.95, nlayers), condensate = zeros(nlayers))
    supply = 1.0e-9                                             # [kg/kg/s]
    references.temperature_reference[ij, :, 1] .= T[ij, :]
    references.humidity_reference[ij, :, 1] .= q[ij, :] .- supply * Δt_prognostic
    references.surface_pressure_reference[ij, 1] = 1.0e5
    run_cloud_column!(vars, model, ij)
    dT, dq, dqc = cloud_column_tendencies(vars, model, ij)
    for k in 1:nlayers
        p = pressure(k, 1.0e5, geometry.vertical_coordinates)
        q_sat, dq_sat_dT = mixed_phase_saturation_humidity(285.0, p, 0.0, atmosphere, true)
        f = q[ij, k] / q_sat
        b = sundqvist_cloud_fraction(f, critical_relative_humidity(closure, p, 1.0e5))
        expected = b * supply / (1 + atmosphere.latent_heat_condensation / atmosphere.heat_capacity * f * dq_sat_dT)
        @test dqc[k] ≈ expected rtol = 1.0e-6
    end

    # water and enthalpy budgets close with the Sundqvist closure too
    model, vars = cloud_test_model(closure = SundqvistClosure)
    references = vars.prognostic.clouds
    column = (
        temperature = [220.0, 235.0, 250.0, 262.0, 272.0, 282.0, 290.0, 295.0],
        relative_humidity = [0.5, 0.97, 1.05, 0.95, 0.85, 0.99, 0.6, 0.8],
        condensate = [0.0, 2.0e-5, 3.0e-4, -1.0e-5, 5.0e-4, 2.0e-4, 1.0e-4, 0.0],
    )
    set_cloud_column!(vars, model, ij; column...)
    T = get_prognostic_step(vars.grid.temperature, model.time_stepping, model.large_scale_condensation)
    q = get_prognostic_step(vars.grid.humidity, model.time_stepping, model.large_scale_condensation)
    references.temperature_reference[ij, :, 1] .= T[ij, :] .- 1.0e-5 * Δt_prognostic     # cooling and
    references.humidity_reference[ij, :, 1] .= q[ij, :] .- 1.0e-9 * Δt_prognostic        # moistening supply
    references.surface_pressure_reference[ij, 1] = 1.0e5
    run_cloud_column!(vars, model, ij)
    water, heating, latent = cloud_column_budgets(vars, model, ij, column.temperature)
    @test any(cloud_column_tendencies(vars, model, ij)[3] .> 0)        # something condensed
    @test water ≈ 0 atol = 1.0e-12
    @test heating ≈ latent rtol = 1.0e-10
end

@testset "Convective detrainment of condensate" begin
    spectral_grid = SpectralGrid(truncation = 15, nlayers = 8, NF = Float64)
    model = PrimitiveWetModel(spectral_grid; large_scale_condensation = PrognosticCloudCondensation(spectral_grid))
    model.feedback.verbose = false
    simulation = initialize!(model)
    run!(simulation, steps = 20)    # spin up so some columns convect
    vars = simulation.variables
    P = vars.parameterizations
    (; time_stepping, geometry, planet, atmosphere) = model
    dq = get_tendency_step(vars.tendencies.grid.humidity, time_stepping, model.convection)
    dqc = get_tendency_step(vars.tendencies.grid.cloud_condensate, time_stepping, model.convection)

    # the same state with and without detrainment
    function convect!(detrainment)
        dq .= 0
        dqc .= 0
        P.rain_convection .= 0
        P.snow_convection .= 0
        convection = BettsMillerConvection(spectral_grid; detrainment)
        SpeedyWeather._column_parameterizations_cpu!(vars, (; convection), model)
        return copy(dq), copy(dqc), P.rain_convection .+ P.snow_convection
    end
    dq₀, dqc₀, precip₀ = convect!(0.0)
    dq₁, dqc₁, precip₁ = convect!(0.3)

    @test any(precip₀ .> 0)
    @test all(dqc₀ .== 0)
    @test dq₁ == dq₀                                            # humidity tendency unchanged
    @test precip₁ ≈ 0.7 .* precip₀                              # 30 % of the precipitation detrained

    # detrained condensate equals the precipitation removed, all in one layer per column
    nlayers = spectral_grid.nlayers
    Δt = time_stepping.Δt
    for ij in 1:spectral_grid.npoints
        Δp = [pressure_thickness(k, P.surface_pressure[ij], geometry.vertical_coordinates) for k in 1:nlayers]
        detrained = sum(dqc₁[ij, :] .* Δp) / planet.gravity     # [kg/m²/s]
        @test detrained ≈ 0.3 * precip₀[ij] * atmosphere.water_density / Δt atol = 1.0e-12
        @test count(!iszero, dqc₁[ij, :]) <= 1
    end

    # without a prognostic condensate detrainment does nothing
    model_default = PrimitiveWetModel(spectral_grid; convection = BettsMillerConvection(spectral_grid; detrainment = 0.3))
    @test SpeedyWeather.convective_detrainment(model_default.convection, Variables(model_default)) == 0
end

@testset "Two-stream cloud layer" begin
    two_stream = SpeedyWeather.two_stream_diffuse_layer
    for NF in (Float32, Float64)
        g = NF(0.85)
        γ₁ = sqrt(NF(3)) / 2 * (1 - g)
        for τ in NF.((0, 0.01, 1, 10, 100))
            # non-absorbing: R + T = 1 and R = γ₁τ/(1 + γ₁τ)
            R, T = two_stream(τ, one(NF), g)
            @test R + T ≈ 1 rtol = 1.0e-4
            @test R ≈ γ₁ * τ / (1 + γ₁ * τ) atol = 1.0e-4
            # absorbing: R + T < 1, both within [0, 1]
            R, T = two_stream(τ, NF(0.999), g)
            @test 0 <= R <= 1 && 0 <= T <= 1
            @test R + T <= 1 + eps(NF)
        end
    end
    R, T = two_stream(20.0, 0.999, 0.85)
    @test 0.01 < 1 - R - T < 0.1                            # a few percent absorption by a thick cloud
end

@testset "Per-layer shortwave clouds" begin
    spectral_grid = SpectralGrid(truncation = 15, nlayers = 8, NF = Float64)
    radiation = Radiation(spectral_grid; shortwave = OneBandCloudyShortwave(spectral_grid), longwave = nothing)
    model = PrimitiveWetModel(spectral_grid; radiation, parameterizations = (:radiation,))
    model.feedback.verbose = false
    simulation = initialize!(model)
    vars = simulation.variables
    P = vars.parameterizations
    shortwave = model.radiation.shortwave
    (; geometry, planet, atmosphere) = model
    nlayers = spectral_grid.nlayers
    ij = 1
    P.surface_pressure .= 1.0e5
    P.cos_zenith .= 0.7
    P.cloud_liquid_effective_radius .= 10.0e-6
    P.cloud_ice_effective_radius .= 50.0e-6
    Δp = [pressure_thickness(k, 1.0e5, geometry.vertical_coordinates) for k in 1:nlayers]
    dTdt = vars.tendencies.grid.temperature

    function shortwave_column!(; cloud_fraction, liquid, ice = zeros(nlayers), albedo = 0.2, radiative_transfer = shortwave.radiative_transfer)
        P.cloud_fraction[ij, :] .= cloud_fraction
        P.cloud_liquid_water[ij, :] .= liquid
        P.cloud_ice_water[ij, :] .= ice
        P.ocean.albedo[ij] = albedo
        P.land.albedo[ij] = albedo
        dTdt[ij, :, 1] .= 0
        clouds = SpeedyWeather.clouds!(ij, vars, shortwave.clouds, model)
        t = SpeedyWeather.transmissivity!(ij, vars, clouds, shortwave.transmissivity, model)
        SpeedyWeather.shortwave_radiative_transfer!(ij, vars, t, clouds, radiative_transfer, model)
        absorbed = atmosphere.heat_capacity / planet.gravity * sum(dTdt[ij, :, 1] .* Δp)   # [W/m²]
        return (; absorbed, surface_down = P.surface_shortwave_down[ij], outgoing = P.outgoing_shortwave[ij])
    end
    D_toa = planet.solar_constant * 0.7

    # energy conservation: absorbed in the atmosphere + at the surface + reflected to space = incoming
    cloudy = (
        cloud_fraction = [0, 0.3, 0.6, 0, 0.5, 0.8, 0.2, 0], liquid = [0, 1.0e-5, 1.0e-4, 0, 2.0e-4, 3.0e-4, 1.0e-5, 0],
        ice = [0, 5.0e-5, 2.0e-5, 0, 0, 0, 0, 0],
    )
    for albedo in (0.0, 0.2, 0.8)
        fluxes = shortwave_column!(; cloudy..., albedo)
        @test fluxes.absorbed + (1 - albedo) * fluxes.surface_down + fluxes.outgoing ≈ D_toa rtol = 1.0e-10
    end

    # clouds reflect: more outgoing, less at the surface than clear sky
    clear = shortwave_column!(; cloud_fraction = zeros(nlayers), liquid = zeros(nlayers))
    cloudy_fluxes = shortwave_column!(; cloudy...)
    @test cloudy_fluxes.outgoing > clear.outgoing
    @test cloudy_fluxes.surface_down < clear.surface_down

    # an overcast, optically thick layer reflects most of the sunlight
    overcast = shortwave_column!(; cloud_fraction = [0, 0, 0, 0, 1, 0, 0, 0], liquid = [0, 0, 0, 0, 2.0e-3, 0, 0, 0])
    @test overcast.outgoing > 0.5 * D_toa
    @test overcast.surface_down < 0.3 * D_toa

    # clear sky, no ozone, black surface: the same as the one-band transfer without clouds
    no_ozone = SpeedyWeather.CloudyShortwaveRadiativeTransfer(spectral_grid; ozone_absorption = 0)
    one_band = OneBandShortwaveRadiativeTransfer(spectral_grid; ozone_absorption = 0)
    a = shortwave_column!(; cloud_fraction = zeros(nlayers), liquid = zeros(nlayers), albedo = 0.0, radiative_transfer = no_ozone)
    heating_cloudy = copy(dTdt[ij, :, 1])
    b = shortwave_column!(; cloud_fraction = zeros(nlayers), liquid = zeros(nlayers), albedo = 0.0, radiative_transfer = one_band)
    @test a.surface_down ≈ b.surface_down rtol = 1.0e-10
    @test heating_cloudy ≈ dTdt[ij, :, 1] rtol = 1.0e-10
    @test a.outgoing ≈ 0 atol = 1.0e-12
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

@testset "Prognostic clouds with the Sundqvist closure in the full model" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
    large_scale_condensation = PrognosticCloudCondensation(spectral_grid; closure = SundqvistClosure(spectral_grid))
    radiation = Radiation(spectral_grid; shortwave = OneBandCloudyShortwave(spectral_grid), longwave = OneBandCloudyLongwave(spectral_grid))
    model = PrimitiveWetModel(spectral_grid; large_scale_condensation, radiation)
    model.feedback.verbose = false
    simulation = initialize!(model)
    run!(simulation, period = Day(2))

    vars = simulation.variables
    @test model.feedback.nans_detected == false
    @test all(vars.prognostic.clouds.surface_pressure_reference .> 0)
    @test any(vars.grid.cloud_condensate .> 0)
    @test all(0 .<= vars.parameterizations.cloud_fraction .<= 1)
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
