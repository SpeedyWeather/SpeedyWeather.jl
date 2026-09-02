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
    end
end

