@testset "Particle advection" begin
    for Model in (BarotropicModel,
                    ShallowWaterModel,
                    PrimitiveDryModel,
                    PrimitiveWetModel)

        if Model <: PrimitiveEquation
            nlev = 8
        else
            nlev = 1
        end

        spectral_grid = SpectralGrid(trunc=31, nlev=1, nparticles=100)
        particle_advection = ParticleAdvection2D(spectral_grid)

        model = Model(;spectral_grid,particle_advection)
        add!(model.callbacks, ParticleTracker(spectral_grid))

        simulation = initialize!(model)
        run!(simulation, period=Day(1))

        for particle in simulation.prognostic_variables.particles
            @test SpeedyWeather.ismod(particle)
        end
    end
end

