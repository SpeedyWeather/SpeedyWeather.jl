using Statistics, Random

@testset "Spectral white noise process: standard deviation, structure, round trip" begin
    for truncation in (32, 64)
        for wavenumber in (0, 12)
            for σ in (1 / 4, 1)
                # equal-area HEALPixGrid for unweighted mean and standard deviation
                spectral_grid = SpectralGrid(; truncation, nlayers = 1, Grid = HEALPixGrid)
                random_process = SpectralWhiteNoiseProcess(spectral_grid; seed = 123, wavenumber, standard_deviation = σ)
                model = BarotropicModel(spectral_grid; random_process)
                simulation = initialize!(model)

                # no noise on the global mean (l = 0) and on the last degree l = lmax
                @test random_process.noise_factors[1] == 0
                @test random_process.noise_factors[end] == 0
                @test all(random_process.noise_factors[2:(end - 1)] .> 0)

                spec = simulation.variables.prognostic.random_pattern
                grid = simulation.variables.grid.random_pattern
                stds = zeros(20)
                for i in eachindex(stds)
                    SpeedyWeather.update_random!(simulation.variables, random_process)
                    transform!(grid, spec, model.spectral_transform)
                    stds[i] = std(grid)
                    @test mean(grid) ≈ 0 atol = 1.0e-2
                end
                @test mean(stds) ≈ σ rtol = 5.0e-2

                # m = 0 coefficients are real
                m0 = spec.spectrum.m_indices .== 1
                @test all(imag.(spec.data[m0]) .== 0)
            end
        end
    end
end

@testset "Spectral white noise process: band-limited, grid round trip" begin
    # Gaussian grid: spectral -> grid -> spectral is exact to rounding for band-limited fields
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 1)
    random_process = SpectralWhiteNoiseProcess(spectral_grid; seed = 123)
    model = BarotropicModel(spectral_grid; random_process)
    simulation = initialize!(model)
    SpeedyWeather.update_random!(simulation.variables, random_process)
    spec = simulation.variables.prognostic.random_pattern
    grid = simulation.variables.grid.random_pattern
    transform!(grid, spec, model.spectral_transform)        # no clamp by default anyway
    spec2 = zero(spec)
    transform!(spec2, grid, model.spectral_transform)
    @test maximum(abs, spec2.data .- spec.data) / maximum(abs, spec.data) < 1.0e-5
end

@testset "Spectral white noise process: memoryless, seed, independent patterns" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 1)
    for seed in (0, 123)
        random_process = SpectralWhiteNoiseProcess(spectral_grid; seed, n_patterns = 3)
        model = BarotropicModel(spectral_grid; random_process)
        simulation = initialize!(model)
        spec = simulation.variables.prognostic.random_pattern
        @test size(spec, 2) == 3

        SpeedyWeather.update_random!(simulation.variables, random_process)
        spec1 = copy(spec)
        SpeedyWeather.update_random!(simulation.variables, random_process)

        # memoryless: consecutive draws and the patterns among themselves are uncorrelated
        c(a, b) = real(sum(a .* conj.(b))) / sqrt(sum(abs2, a) * sum(abs2, b))
        @test abs(c(spec1.data[:, 1], spec.data[:, 1])) < 0.2
        @test abs(c(spec.data[:, 1], spec.data[:, 2])) < 0.2
        @test abs(c(spec.data[:, 2], spec.data[:, 3])) < 0.2

        # reinitialize: same seed, same draws; seed = 0 takes a random seed
        simulation = initialize!(model)
        SpeedyWeather.update_random!(simulation.variables, random_process)
        if seed != 0
            @test spec1 == simulation.variables.prognostic.random_pattern
        else
            @test spec1 != simulation.variables.prognostic.random_pattern
        end
    end
end

@testset "Grid white noise process" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 1)
    σ = 0.5
    random_process = GridWhiteNoiseProcess(spectral_grid; seed = 123, standard_deviation = σ, n_patterns = 2)
    model = BarotropicModel(spectral_grid; random_process)
    simulation = initialize!(model)
    SpeedyWeather.update_random!(simulation.variables, random_process)
    pattern = copy(simulation.variables.grid.random_pattern)
    @test size(pattern, 2) == 2
    @test std(pattern.data) ≈ σ rtol = 5.0e-2
    @test mean(pattern.data) ≈ 0 atol = 2.0e-2

    # the next draw is new, a reinitialized process repeats the first
    SpeedyWeather.update_random!(simulation.variables, random_process)
    @test pattern != simulation.variables.grid.random_pattern
    simulation = initialize!(model)
    SpeedyWeather.update_random!(simulation.variables, random_process)
    @test pattern == simulation.variables.grid.random_pattern
end

@testset "Random number generator per architecture" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 1)
    @test SpeedyWeather.architecture_rng(SpeedyWeather.CPU(), 1) isa Random.Xoshiro
    @test SpeedyWeather.architecture_rng(SpeedyWeather.CPU(), 1; draw_on_cpu = true) isa Random.Xoshiro
    for Process in (SpectralAR1Process, SpectralWhiteNoiseProcess, GridWhiteNoiseProcess)
        @test Process(spectral_grid; seed = 1).random_number_generator isa Random.Xoshiro
        @test Process(spectral_grid; seed = 1, draw_on_cpu = true).random_number_generator isa Random.Xoshiro
    end
end
