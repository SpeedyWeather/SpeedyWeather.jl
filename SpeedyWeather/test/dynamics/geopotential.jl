@testset "Geopotential reasonable" begin
    for NF in (Float32, Float64)
        nlayers = 8
        spectral_grid = SpectralGrid(; NF, nlayers, Grid = FullGaussianGrid)
        model = PrimitiveWetModel(spectral_grid)
        simulation = initialize!(model)
        progn = simulation.prognostic_variables
        diagn = simulation.diagnostic_variables
        temp = get_step(progn.temp, 1)
        humid = get_step(progn.humid, 1)

        # give every layer some constant temperature by setting the l=m=0 mode (index 1) for all k
        temp0 = 280      # in Kelvin
        temp[1, :] .= temp0 * model.spectral_transform.norm_sphere
        humid .= 0

        lf = 1      # leapfrog time step
        SpeedyWeather.transform!(diagn, progn, lf, model)
        SpeedyWeather.linear_virtual_temperature!(diagn, progn, lf, model)
        SpeedyWeather.geopotential!(diagn, model.geopotential, model.orography)
        geopot_grid = transform(diagn.dynamics.geopotential, model.spectral_transform)

        # approximate heights [m] for this setup
        heights = [27000, 18000, 13000, 9000, 6000, 3700, 1800, 700]

        height_over_ocean = geopot_grid[2304, :] / model.planet.gravity       # somwhere in the tropics
        for k in eachmatrix(temp)
            @test heights[k] ≈ height_over_ocean[k] rtol = 0.5                # very large error allowed
        end
    end
end

@testset "Add geopotential and kinetic energy, compute -∇²B term, no errors" begin
    for NF in (Float32, Float64)
        spectral_grid = SpectralGrid(; NF, nlayers = 8, Grid = FullGaussianGrid)
        model = PrimitiveWetModel(spectral_grid)
        simulation = initialize!(model)
        progn = simulation.prognostic_variables
        diagn = simulation.diagnostic_variables
        temp = get_step(progn.temp, 1)
        humid = get_step(progn.humid, 1)

        # give every layer some constant temperature by setting the l=m=0 mode (index 1) for all k
        temp0 = 280      # in Kelvin
        temp[1, :] .= temp0 * model.spectral_transform.norm_sphere
        humid .= 0

        lf = 1      # leapfrog time step
        SpeedyWeather.transform!(diagn, progn, lf, model)
        SpeedyWeather.linear_virtual_temperature!(diagn, progn, lf, model)
        SpeedyWeather.geopotential!(diagn, model.geopotential, model.orography)
        SpeedyWeather.bernoulli_potential!(diagn, model.spectral_transform)
    end
end

@testset "Virtual temperature calculation" begin
    for NF in (Float32, Float64)
        spectral_grid = SpectralGrid(; NF)
        atmosphere = EarthAtmosphere(spectral_grid)

        for T in NF[250.0, 280.0, 310.0]

            # zero humidity
            Tv = SpeedyWeather.virtual_temperature(T, 0, atmosphere)
            @test Tv == T
            @test typeof(Tv) == typeof(T) == NF

            Tv = SpeedyWeather.linear_virtual_temperature(T, 0, T, atmosphere)
            @test Tv == T
            @test typeof(Tv) == typeof(T) == NF

            for q in NF[0.01, 0.02, 0.03, 0.04]
                Tv = SpeedyWeather.virtual_temperature(T, q, atmosphere)
                @test 1.03 * T > Tv > T       # within a few percent
                @test typeof(Tv) == typeof(T) == NF

                Tv = SpeedyWeather.linear_virtual_temperature(T, q, T, atmosphere)
                @test 1.03 * T > Tv > T
                @test typeof(Tv) == typeof(T) == NF
            end
        end

        atmosphere = EarthDryAtmosphere(spectral_grid)
        for T in NF[250.0, 280.0, 310.0]

            # zero humidity
            Tv = SpeedyWeather.virtual_temperature(T, 0, atmosphere)
            @test Tv == T
            @test typeof(Tv) == typeof(T) == NF

            Tv = SpeedyWeather.linear_virtual_temperature(T, 0, T, atmosphere)
            @test Tv == T
            @test typeof(Tv) == typeof(T) == NF

            for q in NF[0.01, 0.02, 0.03, 0.04]
                Tv = SpeedyWeather.virtual_temperature(T, q, atmosphere)
                @test Tv == T
                @test typeof(Tv) == typeof(T) == NF

                Tv = SpeedyWeather.linear_virtual_temperature(T, q, T, atmosphere)
                @test Tv == T
                @test typeof(Tv) == typeof(T) == NF
            end
        end
    end
end
