# Differentiability of the prognostic cloud scheme and the cloudy one-band radiation, the
# prerequisite for tuning them with SpeedyCalibration.jl (single-step gradients with respect to
# parameters). Reverse-mode gradients of a scalar with respect to all parameters of a scheme, via
# `parameters` and `reconstruct`, are compared with central finite differences, first per column
# (cloud scheme with both closures, shortwave, longwave), then through all parameterizations of a
# time step with the model `Duplicated`, which is how SpeedyCalibration.jl differentiates.
using SpeedyWeather: get_prognostic_step, get_tendency_step, pressure, saturation_humidity,
    cloud_condensation!, parameterization!, reconstruct

# central finite-difference gradient of f with respect to the vector p, relative steps
function central_difference_gradient(f, p; relative_step = 1.0e-6)
    g = zero(p)
    for i in eachindex(p)
        h = relative_step * max(abs(p[i]), 1.0e-3)
        p₊ = copy(p); p₊[i] += h
        p₋ = copy(p); p₋[i] -= h
        g[i] = (f(p₊) - f(p₋)) / 2h
    end
    return g
end

# a column with supersaturated cold layers aloft, cloud in mixed-phase and warm layers, negative
# condensate, rain falling into dry layers below; the lagged step that the physics reads
function set_cloud_ad_column!(vars, model, ij)
    scheme = model.large_scale_condensation
    (; time_stepping, geometry, atmosphere) = model
    temperature = [220.0, 235.0, 250.0, 262.0, 272.0, 282.0, 290.0, 295.0]
    relative_humidity = [0.5, 1.05, 1.1, 0.97, 0.93, 1.02, 0.6, 0.8]
    condensate = [0.0, 2.0e-5, 3.0e-4, -1.0e-5, 5.0e-4, 2.0e-4, 1.0e-4, 0.0]
    T = get_prognostic_step(vars.grid.temperature, time_stepping, scheme)
    q = get_prognostic_step(vars.grid.humidity, time_stepping, scheme)
    qc = get_prognostic_step(vars.grid.cloud_condensate, time_stepping, scheme)
    vars.parameterizations.surface_pressure[ij] = 1.0e5
    for k in eachindex(temperature)
        p = pressure(k, 1.0e5, geometry.vertical_coordinates)
        T[ij, k] = temperature[k]
        q[ij, k] = relative_humidity[k] * saturation_humidity(temperature[k], p, atmosphere)
        qc[ij, k] = condensate[k]
    end
    return nothing
end

# scalar of column ij after the cloud scheme: weighted tendencies, cloud state and precipitation.
# Resets everything the scheme accumulates into or shifts (the Sundqvist references) so that
# every evaluation starts from the same state
function cloud_column_loss(scheme, vars, model, ij, references)
    (; time_stepping, geometry, planet, atmosphere, land_sea_mask) = model
    temp_tend = get_tendency_step(vars.tendencies.grid.temperature, time_stepping, scheme)
    humid_tend = get_tendency_step(vars.tendencies.grid.humidity, time_stepping, scheme)
    condensate_tend = get_tendency_step(vars.tendencies.grid.cloud_condensate, time_stepping, scheme)
    nlayers = size(temp_tend, 2)
    for k in 1:nlayers
        temp_tend[ij, k] = 0
        humid_tend[ij, k] = 0
        condensate_tend[ij, k] = 0
    end
    P = vars.parameterizations
    P.rain_rate[ij] = 0
    P.snow_rate[ij] = 0
    P.cloud_top[ij] = nlayers + 1
    if haskey(vars.prognostic, :clouds)     # Sundqvist closure: a prescribed supply
        clouds = vars.prognostic.clouds
        for k in 1:nlayers
            clouds.temperature_reference[ij, k, 1] = references.temperature[k]
            clouds.humidity_reference[ij, k, 1] = references.humidity[k]
        end
        clouds.surface_pressure_reference[ij, 1] = 1.0e5
    end

    cloud_condensation!(ij, vars, scheme, geometry, planet, atmosphere, land_sea_mask, time_stepping)

    loss = zero(eltype(temp_tend))
    for k in 1:nlayers
        loss += 1.0e3 * temp_tend[ij, k] + 1.0e7 * (humid_tend[ij, k] + 2 * condensate_tend[ij, k]) +
            P.cloud_fraction[ij, k] + 1.0e3 * (P.cloud_liquid_water[ij, k] + P.cloud_ice_water[ij, k])
    end
    return loss + 1.0e5 * (P.rain_rate_large_scale[ij] + P.snow_rate_large_scale[ij])
end

# Outgoing radiation and the heating profile of column `ij` after one radiation stream. Defined at
# top level: as a closure inside the testset, assigning `P` would rebind the testset's `P` (a boxed
# capture), and Enzyme would read the outgoing fluxes through the constant closure, without gradient
function radiation_loss(scheme, vars, model, ij)
    dTdt = vars.tendencies.grid.temperature
    nlayers = size(dTdt, 2)
    for k in 1:nlayers
        dTdt[ij, k, 1] = 0
    end
    parameterization!(ij, vars, scheme, model)
    P = vars.parameterizations
    heating = zero(eltype(P.outgoing_shortwave))
    for k in 1:nlayers
        heating += k * dTdt[ij, k, 1]
    end
    return P.outgoing_shortwave[ij] + P.outgoing_longwave[ij] + 1.0e5 * heating
end

@testset "Differentiability: prognostic cloud column (parameter AD)" begin
    spectral_grid = SpectralGrid(truncation = 15, nlayers = 8, NF = Float64)
    @testset for closure in (RelaxationClosure(spectral_grid), SundqvistClosure(spectral_grid))
        scheme = PrognosticCloudCondensation(spectral_grid; closure)
        model = PrimitiveWetModel(spectral_grid; large_scale_condensation = scheme)
        model.feedback.verbose = false
        vars = initialize!(model).variables
        ij = 1
        set_cloud_ad_column!(vars, model, ij)

        # a moistening and cooling supply for the Sundqvist closure
        T = get_prognostic_step(vars.grid.temperature, model.time_stepping, scheme)
        q = get_prognostic_step(vars.grid.humidity, model.time_stepping, scheme)
        Δt_prognostic = SpeedyWeather.default_time_step(model.time_stepping)
        references = (temperature = T[ij, :] .- 1.0e-5 * Δt_prognostic, humidity = q[ij, :] .- 2.0e-9 * Δt_prognostic)

        p = vec(parameters(scheme))
        @test length(p) > 15
        loss(p, vars) = cloud_column_loss(reconstruct(scheme, p), vars, model, ij, references)

        dp = zero(p)
        autodiff(
            set_runtime_activity(Reverse), Const(loss), Active,
            Duplicated(p, dp), Duplicated(vars, make_zero(vars)),
        )
        dp_fd = central_difference_gradient(q -> loss(q, vars), p)

        @test all(isfinite, dp)
        @test count(!iszero, dp) >= 8                      # most parameters matter in this column
        @test dp ≈ dp_fd rtol = 1.0e-4
    end
end

@testset "Differentiability: cloudy one-band radiation column (parameter and cloud AD)" begin
    spectral_grid = SpectralGrid(truncation = 15, nlayers = 8, NF = Float64)
    # without the bulk cloud absorption as in OneBandCloudyShortwave, but off the kink of
    # min(absorptivity_cloud_base * q, absorptivity_cloud_limit) at 0 = 0 for finite differences
    transmissivity = BackgroundShortwaveTransmissivity(spectral_grid; absorptivity_cloud_base = 0, absorptivity_cloud_limit = 1)
    shortwave = OneBandCloudyShortwave(spectral_grid; transmissivity)
    radiation = Radiation(spectral_grid; shortwave, longwave = OneBandCloudyLongwave(spectral_grid))
    model = PrimitiveWetModel(spectral_grid; radiation, parameterizations = (:radiation,))
    model.feedback.verbose = false
    vars = initialize!(model).variables
    P = vars.parameterizations
    nlayers = spectral_grid.nlayers
    ij = 1
    P.surface_pressure .= 1.0e5
    P.cos_zenith .= 0.6
    P.ocean.albedo .= 0.1
    P.land.albedo .= 0.3
    P.cloud_fraction[ij, :] .= [0, 0.3, 0.6, 0, 0.5, 0.8, 0.2, 0]
    P.cloud_liquid_water[ij, :] .= [0, 1.0e-5, 1.0e-4, 0, 2.0e-4, 3.0e-4, 1.0e-5, 0]
    P.cloud_ice_water[ij, :] .= [0, 5.0e-5, 2.0e-5, 0, 0, 0, 0, 0]
    P.cloud_liquid_effective_radius .= 10.0e-6
    P.cloud_ice_effective_radius .= 50.0e-6
    for k in 1:nlayers
        vars.grid.temperature[ij, k, :] .= 220 + 10 * (k - 1)
        vars.grid.humidity[ij, k, :] .= 1.0e-3 * k
    end

    @testset for stream in (:shortwave, :longwave)
        scheme = getproperty(model.radiation, stream)
        p = vec(parameters(scheme))
        loss(p, vars) = radiation_loss(reconstruct(scheme, p), vars, model, ij)
        dp = zero(p)
        autodiff(set_runtime_activity(Reverse), Const(loss), Active, Duplicated(p, dp), Duplicated(vars, make_zero(vars)))
        dp_fd = central_difference_gradient(q -> loss(q, vars), p)
        @test all(isfinite, dp)
        @test dp ≈ dp_fd rtol = 1.0e-4
    end

    # with respect to the cloud state, as the cloud scheme hands it to radiation
    dvars = make_zero(vars)
    autodiff(
        set_runtime_activity(Reverse), Const((vars) -> radiation_loss(model.radiation.shortwave, vars, model, ij)), Active,
        Duplicated(vars, dvars),
    )
    for k in (3, 6)
        h = 1.0e-9
        P.cloud_liquid_water[ij, k] += h
        loss₊ = radiation_loss(model.radiation.shortwave, vars, model, ij)
        P.cloud_liquid_water[ij, k] -= 2h
        loss₋ = radiation_loss(model.radiation.shortwave, vars, model, ij)
        P.cloud_liquid_water[ij, k] += h
        @test dvars.parameterizations.cloud_liquid_water[ij, k] ≈ (loss₊ - loss₋) / 2h rtol = 1.0e-4
    end
end

@testset "Differentiability: all parameterizations of a step with prognostic clouds (model AD)" begin
    # the single-step gradient SpeedyCalibration.jl uses: the outgoing radiation after all
    # parameterizations with respect to parameters of the model, with the model `Duplicated`.
    # Radiation reads the cloud state the condensation writes in the same step, so cloud
    # microphysics parameters reach the outgoing shortwave within one step
    spectral_grid = SpectralGrid(truncation = 15, nlayers = 8, NF = Float64)
    radiation = Radiation(spectral_grid; shortwave = OneBandCloudyShortwave(spectral_grid), longwave = OneBandCloudyLongwave(spectral_grid))
    model = PrimitiveWetModel(spectral_grid; radiation, large_scale_condensation = PrognosticCloudCondensation(spectral_grid))
    model.feedback.verbose = false
    simulation = initialize!(model)
    run!(simulation, period = Day(2))           # spin up for clouds
    vars = simulation.variables
    @test any(vars.parameterizations.cloud_liquid_water .> 0)

    function outgoing_radiation(vars, model)
        SpeedyWeather.parameterization_tendencies!(vars, model)
        P = vars.parameterizations
        return (sum(P.outgoing_shortwave) + sum(P.outgoing_longwave)) / length(P.outgoing_shortwave)
    end

    dmodel = make_zero(model)
    autodiff(
        set_runtime_activity(Reverse), Const(outgoing_radiation), Active,
        Duplicated(vars, make_zero(vars)), Duplicated(model, dmodel),
    )

    # finite differences, reconstructing the component with the perturbed parameter
    function directional(path)
        component, field_path = first(path), Base.tail(path)
        original = getproperty(model, component)
        value = foldl(getproperty, field_path; init = original)
        perturb!(x) = setproperty!(model, component, reconstruct(original, foldr((name, inner) -> NamedTuple{(name,)}((inner,)), field_path; init = x)))
        h = 1.0e-5 * value
        perturb!(value + h)
        loss₊ = outgoing_radiation(vars, model)
        perturb!(value - h)
        loss₋ = outgoing_radiation(vars, model)
        setproperty!(model, component, original)
        return (loss₊ - loss₋) / 2h
    end

    for (gradient, path) in (
            (dmodel.large_scale_condensation.autoconversion_rate, (:large_scale_condensation, :autoconversion_rate)),
            (dmodel.large_scale_condensation.cloud_fraction_coefficient, (:large_scale_condensation, :cloud_fraction_coefficient)),
            (dmodel.radiation.shortwave.radiative_transfer.asymmetry_factor, (:radiation, :shortwave, :radiative_transfer, :asymmetry_factor)),
        )
        @test isfinite(gradient)
        @test gradient != 0
        @test gradient ≈ directional(path) rtol = 1.0e-3
    end
end
