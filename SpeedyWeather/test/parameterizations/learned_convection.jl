@testset "Learned convection trigger" begin
    @testset "Projection across vertical resolutions" begin
        for NF in (Float32, Float64), nlayers in (7, 8, 16, 32)
            spectral_grid = SpectralGrid(; truncation = 7, nlayers, NF)
            trigger = LearnedConvectionTrigger(spectral_grid)
            convection = BettsMillerConvection(spectral_grid; trigger)
            # Non-default spacing checks that initialization uses the model's geometry.
            coordinate = SigmaCoordinates(spectral_grid, σ -> σ^2)
            model = PrimitiveWetModel(spectral_grid; convection, geometry = Geometry(spectral_grid; vertical_coordinates = coordinate))
            SpeedyWeather.initialize!(convection, model)
            σ = model.geometry.σ_levels_full
            V = SpeedyWeather.chebvander(2 .* Float64.(σ) .- 1, 6)
            # A nonzero c₆ catches accidentally replacing the degree-six fit with degree five.
            coefficients = [250, 40, -5, -2, 1, -0.2, 3]
            profile = NF.(permutedims(hcat(V * coefficients, V * (2 .* coefficients))))
            atol = NF == Float32 ? 2.0e-3 : 1.0e-10
            for ij in 1:2
                fitted = SpeedyWeather.chebyshev_projection(ij, profile, trigger.projection)
                @test collect(fitted) ≈ ij .* coefficients[1:6] atol = atol
                @test (@inferred SpeedyWeather.chebyshev_projection(ij, profile, trigger.projection)) isa NTuple{6, NF}
            end
            adapted = SpeedyWeather.Adapt.adapt(Array, convection)
            @test adapted.trigger.projection == trigger.projection
        end
        spectral_grid = SpectralGrid(truncation = 7, nlayers = 6)
        convection = BettsMillerConvection(spectral_grid; trigger = LearnedConvectionTrigger(spectral_grid))
        model = PrimitiveWetModel(spectral_grid; convection)
        @test_throws ArgumentError SpeedyWeather.initialize!(convection, model)
    end

    @testset "Exported standardisation and logit" begin
        # Metadata in the exact shuffled order provided by training, independent of the
        # implementation's t/q coefficient order. Expected logits evaluated separately.
        features = (
            :t_c2, :t_c1, :t_c0, :t_c3, :q_c2, :q_c0,
            :q_c1, :t_c4, :q_c5, :t_c5, :q_c3, :q_c4,
        )
        means = (
            -5.468058282704101, 41.82535069647731, 245.06016575811805,
            -2.583673702119638, 0.0020090165867202524, 0.0027977503305968415,
            0.0034685521790958898, 1.8256633228157173, -0.00023787588637710037,
            -0.1829393592716407, -4.118595780749046e-5, 0.0001369833605204781,
        )
        stds = (
            13.50237977263115, 20.193051958292635, 10.742484108475601,
            8.761061887086523, 0.002113051901874061, 0.002218712083825191,
            0.0031575491893514125, 4.79222373752904, 0.0007737742755982443,
            2.3348637844543667, 0.0008579501524223445, 0.0008291460913819926,
        )
        for NF in (Float32, Float64)
            for (standardised, expected) in (
                    (zeros(12), 0.2661257804716962),
                    ([0.4, -0.7, 0.3, 1.2, -0.8, -0.1, 0.6, -0.2, 0.9, 0.5, 0.2, -0.3], 12.454697284602949),
                )
                raw = NamedTuple{features}(Tuple(NF.(means .+ stds .* standardised)))
                t = ntuple(i -> getproperty(raw, Symbol("t_c", i - 1)), 6)
                q = ntuple(i -> getproperty(raw, Symbol("q_c", i - 1)), 6)
                logit = @inferred SpeedyWeather.learned_convection_logit(t, q)
                @test logit isa NF
                @test logit ≈ expected rtol = 1.0e-4
                # At the training means the logit is between 0 and 0.5: thresholding
                # the logit at 0.5 instead of zero would incorrectly reject this column.
                # Use identity projection so the supplied profiles are coefficient fixtures.
                trigger = LearnedConvectionTrigger(Matrix{NF}(SpeedyWeather.LinearAlgebra.I, 6, 6))
                @test SpeedyWeather.convecting(1, reshape(collect(t), 1, :), reshape(collect(q), 1, :), trigger)
            end
        end
    end

    @testset "ML replaces the physical on/off criteria" begin
        spectral_grid = SpectralGrid(truncation = 7, nlayers = 24, NF = Float64)
        geometry = Geometry(spectral_grid)
        atmosphere = EarthAtmosphere(spectral_grid)
        model = (; geometry, atmosphere, planet = Earth(spectral_grid), time_stepping = EulerForward(spectral_grid))
        nlayers = spectral_grid.nlayers
        σ = geometry.σ_levels_full
        Δσ = geometry.σ_levels_thick
        # The wet fixture has PT == 0 and Pq > 0. The second fixture has PT < 0:
        # moisture makes the lifted parcel virtually buoyant despite its lower T.
        for wet in (true, false)
            T = fill(300.0, nlayers)
            q = zeros(nlayers)
            if wet
                T[1:end-1] .= 301
                q .= 0.03
                q[end] = 0.04
            else
                q[end] = 0.02
            end
            vars = (
                grid = (; temperature = reshape(copy(T), 1, nlayers, 1), humidity = reshape(copy(q), 1, nlayers, 1)),
                tendencies = (; grid = (; temperature = zeros(1, nlayers, 1), humidity = zeros(1, nlayers, 1))),
                dynamics = (; geopotential = reshape(-atmosphere.R_dry .* T .* log.(σ), 1, nlayers)),
                scratch = (; grid = (; a = zeros(1, nlayers), b = zeros(1, nlayers))),
                parameterizations = (
                    surface_pressure = [1.0e5], rain_convection = [0.0], snow_convection = [0.0],
                    rain_rate_convection = [0.0], snow_rate_convection = [0.0],
                    rain_rate = [0.0], snow_rate = [0.0], cloud_top = [nlayers + 1],
                ),
            )
            # Deliberately map these columns to a known accepted feature vector.
            # This tests the coupling independently of the classifier's skill.
            t_mean = [245.06016575811805, 41.82535069647731, -5.468058282704101,
                -2.583673702119638, 1.8256633228157173, -0.1829393592716407]
            q_mean = [0.0027977503305968415, 0.0034685521790958898, 0.0020090165867202524,
                -4.118595780749046e-5, 0.0001369833605204781, -0.00023787588637710037]
            projection = zeros(6, nlayers)
            projection[:, end-1:end] .= hcat(t_mean, q_mean) / hcat(T[end-1:end], q[end-1:end])
            trigger = LearnedConvectionTrigger(projection)
            learned = BettsMillerConvection(spectral_grid; trigger, snow = false)
            original = BettsMillerConvection(spectral_grid; snow = false)
            temp = SpeedyWeather.get_prognostic_step(vars.grid.temperature, model.time_stepping, learned)
            humid = SpeedyWeather.get_prognostic_step(vars.grid.humidity, model.time_stepping, learned)
            @test SpeedyWeather.convecting(1, temp, humid, trigger)
            SpeedyWeather.pseudo_adiabat!(1, vars.scratch.grid.a, vars.scratch.grid.b, temp, humid,
                vars.dynamics.geopotential, 1.0e5, σ, atmosphere, learned.relative_humidity, learned.entrainment)
            Pq = sum((q .- vec(vars.scratch.grid.b)) .* Δσ)
            PT = -sum((T .- vec(vars.scratch.grid.a)) .* Δσ)
            @test PT <= 0
            @test (Pq > 0) == wet
            # Default Betts-Miller must reject these columns without adding tendencies.
            SpeedyWeather.convection!(1, vars, original, model)
            @test all(iszero, vars.tendencies.grid.temperature)
            @test all(iszero, vars.tendencies.grid.humidity)
            @test iszero(vars.parameterizations.rain_convection[1])
            # The ML decision must activate an adjustment in both cases, with rain
            # selected by the moisture budget rather than the old PT > 0 criterion.
            SpeedyWeather.convection!(1, vars, learned, model)
            dT = vec(vars.tendencies.grid.temperature)
            dq = vec(vars.tendencies.grid.humidity)
            rain = vars.parameterizations.rain_convection[1]
            @test all(isfinite, dT) && all(isfinite, dq) && isfinite(rain)
            @test any(!iszero, dT) || any(!iszero, dq)
            @test wet ? rain > 0 : iszero(rain)
            water = sum(dq .* Δσ) + rain / model.time_stepping.Δt * model.planet.gravity * atmosphere.water_density / 1.0e5
            energy = atmosphere.heat_capacity * sum(dT .* Δσ) + atmosphere.latent_heat_condensation * sum(dq .* Δσ)
            @test water ≈ 0 atol = 1e-15
            @test energy ≈ 0 atol = 1e-10
        end
    end

    @testset "Gate suppresses tendencies and precipitation" begin
        spectral_grid = SpectralGrid(truncation = 64, nlayers = 24)
        model = PrimitiveWetModel(spectral_grid)
        model.feedback.verbose = false
        simulation = initialize!(model)
        run!(simulation, steps = 20)
        vars = simulation.variables
        baseline_rain = copy(vars.parameterizations.rain_convection)
        @test any(baseline_rain .> 0)

        convection = BettsMillerConvection(spectral_grid; trigger = LearnedConvectionTrigger(spectral_grid))
        SpeedyWeather.initialize!(convection, model)
        # A zero projection supplies zero raw coefficients, giving a negative logit.
        # This deliberately rejects every column to isolate gate integration from the
        # classifier's meteorological skill on the spun-up state.
        convection.trigger.projection .= 0
        @test SpeedyWeather.learned_convection_logit(ntuple(_ -> 0.0f0, 6), ntuple(_ -> 0.0f0, 6)) < 0
        temp_tend = SpeedyWeather.get_tendency_step(vars.tendencies.grid.temperature, model.time_stepping, convection)
        humid_tend = SpeedyWeather.get_tendency_step(vars.tendencies.grid.humidity, model.time_stepping, convection)
        temp_tend .= 0
        humid_tend .= 0
        vars.parameterizations.rain_convection .= 0
        vars.parameterizations.snow_convection .= 0
        SpeedyWeather._column_parameterizations_cpu!(vars, (; convection), model)
        @test all(iszero, temp_tend)
        @test all(iszero, humid_tend)
        @test all(iszero, vars.parameterizations.rain_convection)
        @test all(iszero, vars.parameterizations.snow_convection)

        learned_model = PrimitiveWetModel(
            spectral_grid; convection = BettsMillerConvection(
                spectral_grid; trigger = LearnedConvectionTrigger(spectral_grid)
            )
        )
        learned_model.feedback.verbose = false
        learned_simulation = initialize!(learned_model)
        run!(learned_simulation, steps = 3)
        @test all(isfinite, learned_simulation.variables.prognostic.temperature)
    end
end
