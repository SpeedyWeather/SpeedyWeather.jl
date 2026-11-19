using Statistics

@testset "Random process standard deviation" begin
    for truncation in (32, 43, 64, 86)
        for wavenumber in (12, 16, 20, 24)
            for σ in (1 / 16, 1 / 4, 1, 4, 16)

                # use equal-area HEALPixGrid, for unweighted mean and standard deviation later
                spectral_grid = SpectralGrid(; truncation, nlayers = 1, Grid = HEALPixGrid)

                seed = 123
                random_process = SpectralAR1Process(spectral_grid; seed, wavenumber, standard_deviation = σ)
                model = BarotropicModel(spectral_grid; random_process)
                simulation = initialize!(model)

                # instead of running the model, just run the random process
                nsteps = 100
                for _ in 1:nsteps
                    SpeedyWeather.update_random!(simulation.variables, random_process)
                end

                # transform to grid space (don't clamp)
                grid = simulation.variables.grid.random_pattern
                spec = simulation.variables.prognostic.random_pattern
                @test size(spec, 2) == size(grid, 2) == 1   # single pattern: trailing singleton dimension
                transform!(grid, spec, model.spectral_transform)

                @test mean(grid) ≈ 0 atol = 1.0e-2
                @test std(grid) ≈ σ rtol = 1.0e-1
            end
        end
    end
end

@testset "Random process seed" begin
    for seed in (0, 123, 1234, 12345)

        spectral_grid = SpectralGrid(truncation = 32, nlayers = 1)
        random_process = SpectralAR1Process(spectral_grid; seed)

        model = BarotropicModel(spectral_grid; random_process)
        simulation = initialize!(model)

        # instead of running the model, just run the random process
        nsteps = 100
        for _ in 1:nsteps
            SpeedyWeather.update_random!(simulation.variables, random_process)
        end

        spec1 = copy(simulation.variables.prognostic.random_pattern)

        # reinitialize the model with the same seed
        simulation = initialize!(model)
        nsteps = 100
        for _ in 1:nsteps
            SpeedyWeather.update_random!(simulation.variables, random_process)
        end

        if seed != 0    # non-zero seed is actual seed
            @test spec1 == simulation.variables.prognostic.random_pattern
        else            # but zero seed pulls random seed from global RNG
            @test spec1 != simulation.variables.prognostic.random_pattern
        end
    end
end

@testset "Random process independent patterns" begin
    # equal-area HEALPixGrid for unweighted standard deviation
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 1, Grid = HEALPixGrid)
    σ = 1 / 4
    random_process = SpectralAR1Process(spectral_grid; seed = 123, n_patterns = 3, standard_deviation = σ)
    model = BarotropicModel(spectral_grid; random_process)
    simulation = initialize!(model)
    spec = simulation.variables.prognostic.random_pattern
    grid = simulation.variables.grid.random_pattern
    @test size(spec, 2) == 3
    @test size(grid, 2) == 3

    for _ in 1:100
        SpeedyWeather.update_random!(simulation.variables, random_process)
    end
    transform!(grid, spec, model.spectral_transform)

    # every pattern is an AR1 process with the same standard deviation, independent of the others
    c(a, b) = real(sum(a .* conj.(b))) / sqrt(sum(abs2, a) * sum(abs2, b))
    for k in 1:3
        @test std(grid.data[:, k]) ≈ σ rtol = 1.0e-1
    end
    @test abs(c(spec.data[:, 1], spec.data[:, 2])) < 0.2
    @test abs(c(spec.data[:, 2], spec.data[:, 3])) < 0.2
end

@testset "AR1 process: spectral structure" begin
    # Gaussian grid: spectral -> grid -> spectral is exact to rounding for band-limited fields
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 1)
    for wavenumber in (0, 12)
        random_process = SpectralAR1Process(spectral_grid; seed = 123, wavenumber, clamp = (-Inf, Inf))
        model = BarotropicModel(spectral_grid; random_process)
        simulation = initialize!(model)

        # no noise on the global mean (l = 0) and on the last degree l = lmax, flat for wavenumber = 0
        F = random_process.noise_factors
        @test F[1] == 0
        @test F[end] == 0
        @test all(F[2:(end - 1)] .> 0)
        wavenumber == 0 && @test all(F[2:(end - 1)] .== F[2])

        for _ in 1:10
            SpeedyWeather.update_random!(simulation.variables, random_process)
        end

        # m = 0 coefficients are real
        spec = simulation.variables.prognostic.random_pattern
        m0 = spec.spectrum.m_indices .== 1
        @test all(imag.(spec.data[m0, 1]) .== 0)

        # band-limited: a grid round trip returns the pattern
        grid = simulation.variables.grid.random_pattern
        transform!(grid, spec, model.spectral_transform)
        spec2 = zero(spec)
        transform!(spec2, grid, model.spectral_transform)
        @test maximum(abs, spec2.data .- spec.data) / maximum(abs, spec.data) < 1.0e-5
    end
end

@testset "White noise is an AR1 process without memory" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 1)
    white = SpectralWhiteNoiseProcess(spectral_grid; seed = 123, wavenumber = 12)
    @test white isa SpectralAR1Process
    @test white.time_scale == Second(0)
    @test white.wavenumber == 12                    # keyword arguments override the defaults
    @test white.clamp == (-Inf, Inf)
    @test white.random_numbers === nothing          # drawn directly into the pattern

    # a time scale much shorter than the time step: a = exp(-Δt/τ) underflows to 0,
    # the AR1 process with its buffer draws the same patterns
    ar1 = SpectralAR1Process(spectral_grid; seed = 123, wavenumber = 12, time_scale = Second(1))
    @test ar1.random_numbers !== nothing
    patterns = map((white, ar1)) do random_process
        simulation = initialize!(BarotropicModel(spectral_grid; random_process))
        SpeedyWeather.update_random!(simulation.variables, random_process)
        SpeedyWeather.update_random!(simulation.variables, random_process)
        copy(simulation.variables.prognostic.random_pattern)
    end
    @test ar1.autoregressive_factor[] == 0
    @test ar1.noise_factors == white.noise_factors
    @test patterns[1] == patterns[2]

    # negative time scales are rejected
    random_process = SpectralAR1Process(spectral_grid; time_scale = Second(-1))
    @test_throws ArgumentError initialize!(BarotropicModel(spectral_grid; random_process))
end
