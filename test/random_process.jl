using Statistics

@testset "Random process standard deviation" begin
    for trunc in (31, 42, 63, 85)
        for wavenumber in (12, 16, 20, 24)
            for σ in (1/16, 1/4, 1, 4, 16)

                # use equal-area HEALPixGrid, for unweighted mean and standard deviation later
                spectral_grid = SpectralGrid(; trunc, nlayers=1, Grid=HEALPixGrid)

                random_process = SpectralAR1Process(spectral_grid, wavenumber=wavenumber, standard_deviation=σ)
                model = BarotropicModel(spectral_grid; random_process)
                simulation = initialize!(model)

                # instead of running the model, just run the random process
                nsteps = 100
                for _ in 1:nsteps
                    SpeedyWeather.random_process!(simulation.prognostic_variables, random_process)
                end

                # transform to grid space (don't clamp)
                grid = simulation.diagnostic_variables.grid.random_pattern
                spec = simulation.prognostic_variables.random_pattern
                transform!(grid, spec, model.spectral_transform)

                @test mean(grid) ≈ 0 atol=5e-3
                @test std(grid) ≈ σ rtol=1e-1
            end
        end
    end
end

@testset "Random process seed" begin
    for seed in (123, 1234, 12345)

        spectral_grid = SpectralGrid(trunc=31, nlayers=1)
        random_process = SpectralAR1Process(spectral_grid; seed)

        model = BarotropicModel(spectral_grid; random_process)
        simulation = initialize!(model)

        # instead of running the model, just run the random process
        nsteps = 100
        for _ in 1:nsteps
            SpeedyWeather.random_process!(simulation.prognostic_variables, random_process)
        end

        spec1 = copy(simulation.prognostic_variables.random_pattern)

        # reinitialize the model with the same seed
        simulation = initialize!(model)
        nsteps = 100
        for _ in 1:nsteps
            SpeedyWeather.random_process!(simulation.prognostic_variables, random_process)
        end

        @test spec1 == simulation.prognostic_variables.random_pattern
    end
end