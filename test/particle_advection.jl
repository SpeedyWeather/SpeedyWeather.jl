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

        spectral_grid = SpectralGrid(trunc=31, nlev=1, n_particles=100)
        particle_advection = ParticleAdvection2D(spectral_grid)

        model = Model(;spectral_grid,particle_advection)
        simulation = initialize!(model)
        run!(simulation, period=Day(1))

        for particle in simulation.prognostic_variables.particles
            @test SpeedyWeather.ismod(particle)
        end
    end
end