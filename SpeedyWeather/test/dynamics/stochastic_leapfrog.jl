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

const STOCHASTIC_TEST_MODELS = (BarotropicModel, ShallowWaterModel, PrimitiveWetModel)

@testset "StochasticLeapfrog is a Leapfrog with noise" begin
    spectral_grid = SpectralGrid(truncation = 21, nlayers = 4)
    L = StochasticLeapfrog(spectral_grid; seed = 1)
    @test L isa Leapfrog
    @test L isa StochasticLeapfrog
    @test L.noise isa AdditiveWienerNoise
    @test !(Leapfrog(spectral_grid) isa StochasticLeapfrog)
    @test Leapfrog(spectral_grid; noise = AdditiveWienerNoise(spectral_grid)) isa StochasticLeapfrog
    @test_throws ArgumentError StochasticLeapfrog(spectral_grid; noise = AdditiveWienerNoise(spectral_grid), seed = 1)
end

@testset "StochasticLeapfrog with zero amplitudes is Leapfrog" begin
    for Model in STOCHASTIC_TEST_MODELS
        nlayers = Model <: Union{BarotropicModel, ShallowWaterModel} ? 1 : 4
        spectral_grid = SpectralGrid(; truncation = 21, nlayers)
        zero_noise = AdditiveWienerNoise(
            spectral_grid; seed = 1, vorticity = 0, divergence = 0, temperature = 0,
            humidity = 0, sea_surface_temperature = 0,
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
        kw = Model <: BarotropicModel ? (; vorticity = 1.0e-7) :
            Model <: ShallowWaterModel ? (; vorticity = 1.0e-7, divergence = 1.0e-8, η = 1.0e-3) : (;)
        sims = map(((1, false), (1, false), (2, false), (1, true), (1, true))) do (seed, draw_on_cpu)
            model = Model(spectral_grid; time_stepping = StochasticLeapfrog(spectral_grid; seed, draw_on_cpu, kw...))
            simulation = seeded_initialize!(model)
            run!(simulation, steps = 10)
            simulation
        end
        @test same_state(sims[1], sims[2])                  # same seed
        @test !same_state(sims[1], sims[3])                 # other seed
        @test same_state(sims[4], sims[5])                  # same seed, drawn on the CPU
        @test same_state(sims[1], sims[4])                  # on CPU drawing on the architecture is drawing on the CPU
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

@testset "StochasticLeapfrog: noise variables from the tendencies of the model" begin
    for Model in STOCHASTIC_TEST_MODELS
        nlayers = Model <: Union{BarotropicModel, ShallowWaterModel} ? 1 : 4
        spectral_grid = SpectralGrid(; truncation = 21, nlayers)
        model = Model(spectral_grid; time_stepping = StochasticLeapfrog(spectral_grid; seed = 1))
        vars = SpeedyWeather.Variables(model)
        primitive = Model <: Union{PrimitiveDryModel, PrimitiveWetModel}

        # a slot per layer of every (atmospheric) prognostic variable
        n = sum(name -> size(get_step(vars.tendencies[name], 1).data, 2), SpeedyWeather.tendency_names(vars))
        @test size(vars.prognostic.wiener_pattern, 2) == n
        @test n == (Model <: BarotropicModel ? nlayers : Model <: ShallowWaterModel ? 2nlayers + 1 :
            Model <: PrimitiveDryModel ? 3nlayers + 1 : 4nlayers + 1)

        # all atmospheric variables are forced in spectral space; among ocean and land variables
        # only the sea surface temperature has a (default) amplitude
        @test !haskey(vars.grid, :wiener_pattern)
        @test haskey(vars.prognostic, :wiener_pattern_surface) == primitive
        primitive && @test size(vars.prognostic.wiener_pattern_surface, 2) == 1
    end

    # amplitudes for variables the model does not have are ignored, with a warning
    spectral_grid = SpectralGrid(truncation = 21, nlayers = 1)
    model = BarotropicModel(spectral_grid; time_stepping = StochasticLeapfrog(spectral_grid; seed = 1, η = 1))
    @test_logs (:warn, r"no tendency for \(:η,\)") initialize!(model)
end

@testset "StochasticLeapfrog: increment a/√Δt per variable, slot and layer" begin
    spectral_grid = SpectralGrid(truncation = 21, nlayers = 4)
    noise = AdditiveWienerNoise(
        spectral_grid; seed = 1, vorticity = 2, divergence = 3, temperature = 5, humidity = 7,
        pressure = 11, sea_surface_temperature = 13, soil_temperature = 17, tapering = σ -> σ,
    )
    L = StochasticLeapfrog(spectral_grid; noise)
    model = PrimitiveWetModel(spectral_grid; time_stepping = L)
    simulation = initialize!(model)
    SpeedyWeather.initialize!(simulation)       # scales vorticity, divergence by the radius
    vars = simulation.variables
    nlayers = spectral_grid.nlayers
    nsoil = size(get_step(vars.tendencies.land.soil_temperature, 1).data, 2)

    # variables with layers first (vorticity, divergence, temperature, humidity), then pressure;
    # surface: sea surface temperature, then soil temperature (the others have no amplitude)
    @test size(vars.prognostic.wiener_pattern, 2) == 4nlayers + 1
    @test size(vars.prognostic.wiener_pattern_surface, 2) == 1 + nsoil
    @test noise.taper ≈ model.geometry.σ_levels_full

    for Δt in (Minute(20), Minute(40))
        set!(L, Δt)
        SpeedyWeather.wiener_draw!(vars, noise)         # a pattern to apply
        SpeedyWeather.reset_tendencies!(vars, L)
        SpeedyWeather.wiener_tendencies!(vars, noise, L, model)

        isdt = 1 / sqrt(L.Δt)
        R = vars.prognostic.scale[]
        a = noise.amplitudes
        spec = vars.prognostic.wiener_pattern.data
        surface = vars.prognostic.wiener_pattern_surface.data
        tend(x) = get_step(x, 1).data      # the (only) tendency step, as plain array
        for k in 1:nlayers
            t = noise.taper[k]
            # spectral tendencies, radius-scaled as after the dynamical core (vorticity, divergence twice)
            @test tend(vars.tendencies.vorticity)[:, k] == (a.vorticity * isdt * R * R * t) .* spec[:, k]
            @test tend(vars.tendencies.divergence)[:, k] == (a.divergence * isdt * R * R * t) .* spec[:, nlayers + k]
            @test tend(vars.tendencies.temperature)[:, k] == (a.temperature * isdt * R * t) .* spec[:, 2nlayers + k]
            @test tend(vars.tendencies.humidity)[:, k] == (a.humidity * isdt * R * t) .* spec[:, 3nlayers + k]
        end
        @test tend(vars.tendencies.pressure) == (a.pressure * isdt * R) .* spec[:, 4nlayers + 1]
        @test all(iszero, tend(vars.tendencies.grid.temperature))   # nothing in grid space
        sst = tend(vars.tendencies.ocean.sea_surface_temperature)
        @test sst == (a.sea_surface_temperature * isdt) .* surface[:, 1]
        @test any(!iszero, sst)
        for k in 1:nsoil
            @test tend(vars.tendencies.land.soil_temperature)[:, k] == (a.soil_temperature * isdt) .* surface[:, 1 + k]
        end
        @test all(iszero, tend(vars.tendencies.ocean.sea_ice_concentration))
    end
end

@testset "StochasticLeapfrog: noise conserves the global means" begin
    # no noise on l = 0: with surface pressure noise the mass is conserved as without
    spectral_grid = SpectralGrid(truncation = 21, nlayers = 4)
    L = StochasticLeapfrog(spectral_grid; seed = 1, pressure = 1.0e-4)
    simulation = seeded_initialize!(PrimitiveWetModel(spectral_grid; time_stepping = L))
    mean_pressure() = Array(get_step(simulation.variables.prognostic.pressure, 1).data)[1]
    p0 = mean_pressure()
    run!(simulation, steps = 10)
    @test mean_pressure() == p0
    @test !simulation.model.feedback.nans_detected
end
