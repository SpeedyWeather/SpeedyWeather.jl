@testset "GPU 3D interpolation of random data" begin
    arch = SpeedyWeather.GPU()
    NF = Float32
    nlayers = 8
    nparticles = 100

    grid = FullGaussianGrid(24, arch)
    A = rand(NF, grid, nlayers)
    geometry = RingGrids.GridGeometry(A)
    locator = RingGrids.AnvilLocator(NF, nparticles, nlayers; architecture = arch)

    particles = on_architecture(arch, rand(Particle{NF}, nparticles))
    lons = on_architecture(arch, [p.lon for p in Array(particles)])
    lats = on_architecture(arch, [p.lat for p in Array(particles)])
    RingGrids.update_locator!(locator, geometry, lons, lats)

    σ_levels_full = on_architecture(arch, collect(range(NF(0.05), NF(0.95), length = nlayers)))

    Aout = on_architecture(arch, zeros(NF, nparticles))
    RingGrids.interpolate_3D!(Aout, A, locator, geometry, particles, σ_levels_full)

    @test all(isfinite, Array(Aout))
end

@testset "GPU 3D particle advection in isolation" begin
    arch = SpeedyWeather.GPU()
    spectral_grid = SpectralGrid(truncation = 31, nlayers = 4, architecture = arch)
    particle_advection = ParticleAdvection3D(spectral_grid, nparticles = 20)

    model = PrimitiveDryModel(spectral_grid; particle_advection)
    model.feedback.verbose = false
    simulation = initialize!(model)
    vars = simulation.variables

    # call particle_advection! directly to test the advection on its own, without the rest of
    # the time step; the full simulation is tested below
    for step in 1:10
        vars.prognostic.clock.time_step_counter = step
        SpeedyWeather.particle_advection!(vars, model)
    end

    σ_final = [p.σ for p in Array(vars.prognostic.particles)]
    @test all(0 .<= σ_final .<= 1)
end

# largest horizontal displacement of any particle between two snapshots
max_displacement(before, after) = maximum(max(abs(a.lon - b.lon), abs(a.lat - b.lat)) for (a, b) in zip(after, before))

# every_n_time_steps = 1 so that particles are advected within the two time steps
@testset "GPU 3D particle advection full simulation" begin
    arch = SpeedyWeather.GPU()
    spectral_grid = SpectralGrid(truncation = 31, nlayers = 4, architecture = arch)
    particle_advection = ParticleAdvection3D(spectral_grid, nparticles = 20, every_n_time_steps = 1)

    model = PrimitiveDryModel(spectral_grid; particle_advection)
    model.feedback.verbose = false
    simulation = initialize!(model)
    particles_initial = Array(simulation.variables.prognostic.particles)

    @test_nowarn run!(simulation, steps = 2, output = false)
    particles_final = Array(simulation.variables.prognostic.particles)
    @test max_displacement(particles_initial, particles_final) > 0
    @test all(0 <= p.σ <= 1 for p in particles_final)
end

@testset "GPU 2D particle advection full simulation" begin
    arch = SpeedyWeather.GPU()
    spectral_grid = SpectralGrid(truncation = 31, nlayers = 1, architecture = arch)
    particle_advection = ParticleAdvection2D(spectral_grid, nparticles = 20, every_n_time_steps = 1)

    model = BarotropicModel(spectral_grid; particle_advection)
    model.feedback.verbose = false
    simulation = initialize!(model)
    particles_initial = Array(simulation.variables.prognostic.particles)

    @test_nowarn run!(simulation, steps = 2, output = false)
    particles_final = Array(simulation.variables.prognostic.particles)
    @test max_displacement(particles_initial, particles_final) > 0
    @test all(p -> isfinite(p.lon) && isfinite(p.lat), particles_final)
end
