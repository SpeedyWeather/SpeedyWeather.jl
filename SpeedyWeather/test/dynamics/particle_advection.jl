@testset "Particle advection" begin
    @testset for Model in (
            BarotropicModel,
            ShallowWaterModel,
            PrimitiveDryModel,
            PrimitiveWetModel,
        )

        if Model <: PrimitiveEquation
            nlayers = 8
        else
            nlayers = 1
        end

        spectral_grid = SpectralGrid(truncation = 32, nlayers = nlayers)
        particle_advection = ParticleAdvection2D(spectral_grid, layer = 1, nparticles = 100)

        model = Model(spectral_grid; particle_advection)
        model.feedback.verbose = false

        tmp_tracker_path = mktempdir(pwd(), prefix = "tmp_tracker_")  # Cleaned up when the process exits
        add!(model.callbacks, ParticleTracker(spectral_grid, path = tmp_tracker_path, filename = "particles.nc"))

        simulation = initialize!(model)
        run!(simulation, period = Day(1))

        for particle in simulation.variables.prognostic.particles
            @test SpeedyWeather.ismod(particle)
            @test particle.σ == model.geometry.σ_levels_full[1]
        end
    end
end

@testset "find_vertical_bracket" begin
    σf = Float32[0.1, 0.3, 0.5, 0.7, 0.9]

    k_lo, k_hi, α = RingGrids.find_vertical_bracket(0.0f0, σf)
    @test k_lo == 1 && k_hi == 2 && α == 0.0f0       # below grid top → pin

    k_lo, k_hi, α = RingGrids.find_vertical_bracket(σf[1], σf)
    @test k_lo == 1 && k_hi == 2 && α == 0.0f0       # exactly at top layer

    k_lo, k_hi, α = RingGrids.find_vertical_bracket(0.4f0, σf)
    @test k_lo == 2 && k_hi == 3 && α ≈ 0.5f0        # midpoint between layers 2 and 3

    k_lo, k_hi, α = RingGrids.find_vertical_bracket(σf[end], σf)
    @test k_lo == 4 && k_hi == 5 && α == 1.0f0       # exactly at bottom layer

    k_lo, k_hi, α = RingGrids.find_vertical_bracket(1.0f0, σf)
    @test k_lo == 4 && k_hi == 5 && α == 1.0f0       # below grid bottom → pin

    @test all(0 ≤ α ≤ 1 for (_, _, α) in
        [RingGrids.find_vertical_bracket(σ, σf) for σ in range(0f0, 1f0, 50)])
end

@testset "ParticleAdvection3D" begin
    @testset for Model in (PrimitiveDryModel, PrimitiveWetModel)
        spectral_grid = SpectralGrid(truncation = 32, nlayers = 4)
        pa = ParticleAdvection3D(spectral_grid, nparticles = 20)
        @test pa.nparticles == 20
        @test pa.every_n_time_steps == 6

        model = Model(spectral_grid; particle_advection = pa)
        model.feedback.verbose = false
        simulation = initialize!(model)

        σ_initial = copy([p.σ for p in simulation.variables.prognostic.particles])

        run!(simulation, period = Day(1), output = false)

        for particle in simulation.variables.prognostic.particles
            @test SpeedyWeather.ismod(particle)
            @test 0 ≤ particle.σ ≤ 1
        end

        σ_final = [p.σ for p in simulation.variables.prognostic.particles]
        @test !all(σ_initial .≈ σ_final)   # vertical motion occurred

        # particles drift inside the column, not slammed onto σ=0/1. The checks above can't catch
        # that: `0 ≤ σ ≤ 1` holds by the clamp in mod(::Particle), and σ changing even more so.
        @test any(0 .< σ_final .< 1)
    end
end

@testset "ParticleAdvection3D vertical displacement" begin
    # Vertical drift must be gradual: σ̇ ~ O(1e-6) 1/s, so over a few advection steps particles
    # started mid-column stay mid-column. Deterministic (fixed σ), unlike the random init above.
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
    pa = ParticleAdvection3D(spectral_grid, nparticles = 10, every_n_time_steps = 2)
    model = PrimitiveWetModel(spectral_grid; particle_advection = pa)
    model.feedback.verbose = false
    simulation = initialize!(model)

    # place every particle mid-column, away from the boundaries where w vanishes by construction
    particles = simulation.variables.prognostic.particles
    particles .= [SpeedyWeather.set(p; σ = 0.5) for p in particles]
    SpeedyWeather.initialize!(simulation.variables, particles, model.particle_advection, model)

    run!(simulation, period = Hour(4), output = false)

    σ = [p.σ for p in particles]
    @test all(0 .< σ .< 1)             # no particle pinned to a boundary
    @test maximum(abs, σ .- 0.5) < 0.05 # drift is ~0.007 over 4h
end

# shortest signed longitude difference in (-180, 180], so displacements survive the 0/360 wrap
wrapped_dlon(a, b) = mod(a - b + 180, 360) - 180

@testset "Particle advection backwards" begin
    # `backwards` flips the sign of the time step, so identical particles in identical flows move
    # in opposite directions. Short window only, the symmetry decays as trajectories diverge.
    lons = collect(range(0, 350, 10))
    lats = collect(range(-60, 60, 10))
    k = 4

    @testset for Advection in (ParticleAdvection2D, ParticleAdvection3D)
        advect(spectral_grid, backwards) = Advection === ParticleAdvection2D ?
            ParticleAdvection2D(spectral_grid; nparticles = 10, layer = k, every_n_time_steps = 2, backwards) :
            ParticleAdvection3D(spectral_grid; nparticles = 10, every_n_time_steps = 2, backwards)

        simulations = map((false, true)) do backwards
            spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
            model = PrimitiveWetModel(spectral_grid; particle_advection = advect(spectral_grid, backwards))
            model.feedback.verbose = false
            simulation = initialize!(model)

            # identical deterministic seeding in both runs, on layer k's σ level
            P = eltype(simulation.variables.prognostic.particles)
            σk = model.geometry.σ_levels_full[k]
            particles = simulation.variables.prognostic.particles
            particles .= [P(true, lons[i], lats[i], σk) for i in eachindex(particles)]
            SpeedyWeather.initialize!(simulation.variables, particles, model.particle_advection, model)

            run!(simulation, period = Hour(2), output = false)
            simulation
        end

        forward, backward = simulations
        particles_forward = forward.variables.prognostic.particles
        particles_backward = backward.variables.prognostic.particles

        dlon_forward = [wrapped_dlon(particles_forward[i].lon, lons[i]) for i in eachindex(lons)]
        dlon_backward = [wrapped_dlon(particles_backward[i].lon, lons[i]) for i in eachindex(lons)]

        @test all(dlon_forward .* dlon_backward .< 0)                 # opposite directions
        @test maximum(abs, dlon_forward .+ dlon_backward) < 0.05 * maximum(abs, dlon_forward)
    end
end

@testset "ParticleAdvection2D/3D horizontal consistency" begin
    # A 3D particle sitting on layer k's σ level sees the same winds as a 2D particle advected on
    # layer k, so their horizontal tracks should agree while the 3D particle stays near that level.
    k = 4
    lons = collect(range(0, 350, 10))
    lats = collect(range(-60, 60, 10))

    spectral_grid_3d = SpectralGrid(truncation = 32, nlayers = 8)
    model_3d = PrimitiveWetModel(
        spectral_grid_3d;
        particle_advection = ParticleAdvection3D(spectral_grid_3d, nparticles = 10, every_n_time_steps = 2),
    )
    model_3d.feedback.verbose = false
    simulation_3d = initialize!(model_3d)
    σk = model_3d.geometry.σ_levels_full[k]

    spectral_grid_2d = SpectralGrid(truncation = 32, nlayers = 8)
    model_2d = PrimitiveWetModel(
        spectral_grid_2d;
        particle_advection = ParticleAdvection2D(spectral_grid_2d, nparticles = 10, layer = k, every_n_time_steps = 2),
    )
    model_2d.feedback.verbose = false
    simulation_2d = initialize!(model_2d)

    for simulation in (simulation_3d, simulation_2d)
        particles = simulation.variables.prognostic.particles
        P = eltype(particles)
        particles .= [P(true, lons[i], lats[i], σk) for i in eachindex(particles)]
    end
    SpeedyWeather.initialize!(simulation_3d.variables, simulation_3d.variables.prognostic.particles, model_3d.particle_advection, model_3d)
    SpeedyWeather.initialize!(simulation_2d.variables, simulation_2d.variables.prognostic.particles, model_2d.particle_advection, model_2d)

    run!(simulation_3d, period = Hour(6), output = false)
    run!(simulation_2d, period = Hour(6), output = false)

    particles_3d = simulation_3d.variables.prognostic.particles
    particles_2d = simulation_2d.variables.prognostic.particles
    dlon = maximum(abs(wrapped_dlon(particles_3d[i].lon, particles_2d[i].lon)) for i in eachindex(lons))
    dlat = maximum(abs(particles_3d[i].lat - particles_2d[i].lat) for i in eachindex(lons))

    # tracks agree to ~0.03˚ while the 3D particles stay near layer k; they diverged by ~4˚ when
    # the radius scaling of w flung them onto σ=0/1, into the top/bottom layer winds.
    @test dlon < 0.5
    @test dlat < 0.2
end

@testset "ParticleAdvection3D no boundary pile-up" begin
    # statistical counterpart to the above: particles seeded uniformly through the column stay
    # spread through it after a day, rather than collecting on σ=0/1.
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
    pa = ParticleAdvection3D(spectral_grid, nparticles = 50)
    model = PrimitiveWetModel(spectral_grid; particle_advection = pa)
    model.feedback.verbose = false
    simulation = initialize!(model)

    run!(simulation, period = Day(1), output = false)

    σ = [p.σ for p in simulation.variables.prognostic.particles]
    interior_fraction = count(s -> 0 < s < 1, σ) / length(σ)
    @test interior_fraction > 0.8   # 1.0 in practice; was 0.0, all 50 pinned, before the fix
end

