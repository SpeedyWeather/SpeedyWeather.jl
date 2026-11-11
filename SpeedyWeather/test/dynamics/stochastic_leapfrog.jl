using Statistics, Random

# default initial conditions of some models are random (global RNG): seed it for comparisons
seeded_initialize!(model) = (Random.seed!(1234); initialize!(model))

# all prognostic variables (incl ocean, land namespaces) and the clock equal?
function same_state(a, b)
    pa, pb = a.variables.prognostic, b.variables.prognostic
    ok = true
    for name in (:vorticity, :divergence, :temperature, :humidity, :pressure, :η)
        haskey(pa, name) && (ok &= pa[name] == pb[name])
    end
    for ns in (:ocean, :land), name in keys(get(pa, ns, (;)))
        ok &= pa[ns][name] == pb[ns][name]
    end
    return ok & (pa.clock.time == pb.clock.time)
end

const STOCHASTIC_TEST_MODELS = (BarotropicModel, ShallowWaterModel, PrimitiveDryModel, PrimitiveWetModel)

@testset "StochasticLeapfrog with zero amplitudes is Leapfrog" begin
    for Model in STOCHASTIC_TEST_MODELS
        nlayers = Model <: Union{BarotropicModel, ShallowWaterModel} ? 1 : 4
        spectral_grid = SpectralGrid(; truncation = 21, nlayers)
        zero_noise = AdditiveWienerNoise(
            spectral_grid; seed = 1, vorticity = 0, divergence = 0, temperature = 0,
            humidity = 0, pressure = 0, sea_surface_temperature = 0,
        )
        sim_lf = seeded_initialize!(Model(spectral_grid; time_stepping = Leapfrog(spectral_grid)))
        sim_sl = seeded_initialize!(Model(spectral_grid; time_stepping = StochasticLeapfrog(spectral_grid; noise = zero_noise)))
        run!(sim_lf, steps = 10)
        run!(sim_sl, steps = 10)
        @test same_state(sim_lf, sim_sl)
    end
end

@testset "StochasticLeapfrog: reproducible, seed-dependent, stable" begin
    for Model in STOCHASTIC_TEST_MODELS
        nlayers = Model <: Union{BarotropicModel, ShallowWaterModel} ? 1 : 4
        spectral_grid = SpectralGrid(; truncation = 21, nlayers)
        # 2D models: a stronger vorticity noise so that 10 steps clearly differ from Leapfrog
        kw = Model <: Union{BarotropicModel, ShallowWaterModel} ? (; vorticity = 1.0e-7, divergence = 1.0e-8, pressure = 1.0e-3) : (;)
        sims = map((1, 1, 2)) do seed
            model = Model(spectral_grid; time_stepping = StochasticLeapfrog(spectral_grid; seed, kw...))
            simulation = seeded_initialize!(model)
            run!(simulation, steps = 10)
            simulation
        end
        @test same_state(sims[1], sims[2])                  # same seed
        @test !same_state(sims[1], sims[3])                 # other seed
        @test !any(isnan, sims[1].variables.prognostic.vorticity)

        sim_lf = seeded_initialize!(Model(spectral_grid; time_stepping = Leapfrog(spectral_grid)))
        run!(sim_lf, steps = 10)
        @test !same_state(sim_lf, sims[1])                  # the noise acts
    end

    # a few days of the PrimitiveWetModel with the default amplitudes
    spectral_grid = SpectralGrid(truncation = 31, nlayers = 8)
    model = PrimitiveWetModel(spectral_grid; time_stepping = StochasticLeapfrog(spectral_grid; seed = 123))
    simulation = initialize!(model)
    run!(simulation, period = Day(5))
    @test !simulation.model.feedback.nans_detected
end

@testset "StochasticLeapfrog: increment a/√Δt per variable, slot and layer" begin
    spectral_grid = SpectralGrid(truncation = 21, nlayers = 4)
    noise = AdditiveWienerNoise(
        spectral_grid; seed = 1, vorticity = 2, divergence = 3, temperature = 5,
        humidity = 7, pressure = 11, sea_surface_temperature = 13, tapering = σ -> σ,
    )
    L = StochasticLeapfrog(spectral_grid; noise)
    model = PrimitiveWetModel(spectral_grid; time_stepping = L)
    simulation = initialize!(model)
    SpeedyWeather.initialize!(simulation)       # scales vorticity, divergence by the radius
    vars = simulation.variables
    slots = SpeedyWeather.wiener_slots(model)
    nlayers = spectral_grid.nlayers
    @test slots.n == 4nlayers + 1
    @test size(vars.prognostic.wiener_pattern, 2) == slots.n
    @test noise.taper ≈ model.geometry.σ_levels_full

    for Δt in (Minute(20), Minute(40))
        set!(L, Δt)
        SpeedyWeather.wiener_draw!(vars, noise, model)  # a pattern to apply
        SpeedyWeather.reset_tendencies!(vars, L)
        SpeedyWeather.wiener_tendencies!(vars, noise, L, model)

        isdt = 1 / sqrt(L.Δt)
        R = vars.prognostic.scale[]
        spec = vars.prognostic.wiener_pattern
        grid = vars.grid.wiener_pattern
        tend(x) = get_step(x, 1).data      # the (only) tendency step, as plain array
        for k in 1:nlayers
            t = noise.taper[k]
            @test tend(vars.tendencies.vorticity)[:, k] == (noise.vorticity * isdt * R * t) .* spec.data[:, k]
            @test tend(vars.tendencies.divergence)[:, k] == (noise.divergence * isdt * R * t) .* spec.data[:, nlayers + k]
            @test tend(vars.tendencies.grid.temperature)[:, k] == (noise.temperature * isdt * t) .* grid.data[:, 2nlayers + k]
            @test tend(vars.tendencies.grid.humidity)[:, k] == (noise.humidity * isdt * t) .* grid.data[:, 3nlayers + k]
        end
        @test tend(vars.tendencies.grid.pressure) == (noise.pressure * isdt) .* grid.data[:, 4nlayers + 1]
        sst = tend(vars.tendencies.ocean.sea_surface_temperature)
        @test sst == (noise.sea_surface_temperature * isdt) .* vars.prognostic.wiener_pattern_grid.data
        @test any(!iszero, sst)
    end
end
