@testset "Virtual temperature calculation" begin
    for NF in (Float32, Float64)
        spectral_grid = SpectralGrid(; NF)
        atmosphere = EarthAtmosphere(spectral_grid)

        for T in NF[250.0, 280.0, 310.0]

            # zero humidity: virtual temperature equals absolute temperature
            Tv = SpeedyWeather.virtual_temperature(T, 0, atmosphere)
            @test Tv == T
            @test typeof(Tv) == NF

            Tv = SpeedyWeather.linear_virtual_temperature(T, 0, T, atmosphere)
            @test Tv == T
            @test typeof(Tv) == NF

            for q in NF[0.01, 0.02, 0.03, 0.04]
                # moist air is lighter, so virtual temperature > absolute temperature
                Tv = SpeedyWeather.virtual_temperature(T, q, atmosphere)
                @test 1.03 * T > Tv > T
                @test typeof(Tv) == NF

                Tv = SpeedyWeather.linear_virtual_temperature(T, q, T, atmosphere)
                @test 1.03 * T > Tv > T
                @test typeof(Tv) == NF
            end
        end

        atmosphere = EarthDryAtmosphere(spectral_grid)
        for T in NF[250.0, 280.0, 310.0]
            # dry atmosphere: virtual temperature always equals absolute temperature
            Tv = SpeedyWeather.virtual_temperature(T, 0, atmosphere)
            @test Tv == T

            for q in NF[0.01, 0.02, 0.03, 0.04]
                Tv = SpeedyWeather.virtual_temperature(T, q, atmosphere)
                @test Tv == T
                Tv = SpeedyWeather.linear_virtual_temperature(T, q, T, atmosphere)
                @test Tv == T
            end
        end
    end
end

@testset "Spectral geopotential: no errors" begin
    for NF in (Float32, Float64)
        spectral_grid = SpectralGrid(; NF, nlayers = 5)
        model = PrimitiveWetModel(spectral_grid)
        simulation = initialize!(model)
        vars = simulation.variables

        temp = get_step(vars.prognostic.temperature, 1)
        humid = get_step(vars.prognostic.humidity, 1)
        temp[1, :] .= 280 * model.spectral_transform.norm_sphere
        humid .= 0

        SpeedyWeather.transform!(vars, model)
        SpeedyWeather.linear_virtual_temperature!(vars, model)
        SpeedyWeather.geopotential!(vars, model.geopotential, model.orography)
        SpeedyWeather.bernoulli_potential!(vars, model.spectral_transform, model.time_stepping)

        Φ = transform(vars.dynamics.spectral_geopotential, model.spectral_transform)
        # geopotential decreases from top to surface (layer 1 is highest, last is lowest)
        for k in 1:spectral_grid.nlayers-1
            @test all(Φ[:, k] .> Φ[:, k+1])
        end
    end
end

@testset "Grid-space geopotential: analytical check" begin
    # With constant temperature T₀, zero orography, and uniform surface pressure pₛ
    # the hydrostatic equation gives Φ(k) = R * T₀ * log(pₛ / p_k) at every grid point.
    # Non-equispaced half levels: finer spacing near the surface and top.
    σ_half = [0.0, 0.05, 0.15, 0.35, 0.65, 0.85, 0.95, 1.0]
    nlayers = length(σ_half) - 1

    for NF in (Float32, Float64)
        spectral_grid = SpectralGrid(; NF, nlayers)

        coords = (
            SigmaCoordinates(spectral_grid, σ_half),
            SigmaPressureCoordinates(spectral_grid; σ_half, transition = _ -> 1.0),
        )

        for coord in coords
            model = PrimitiveDryModel(spectral_grid;
                orography = NoOrography(spectral_grid),
                geometry  = Geometry(spectral_grid; vertical_coordinates = coord),
            )
            simulation = initialize!(model)
            vars = simulation.variables

            T₀ = NF(280)
            p_surf = NF(1e5)
            vars.grid.temperature .= T₀
            vars.parameterizations.surface_pressure .= p_surf

            SpeedyWeather.geopotential!(vars, model)

            R = model.atmosphere.R_dry
            Φ = Array(vars.dynamics.geopotential)

            for k in 1:nlayers
                pₖ = SpeedyWeather.pressure(k, p_surf, coord)
                expected = R * T₀ * log(p_surf / pₖ)
                @test all(Φ[:, k] .≈ expected)
            end

            # geopotential is strictly decreasing from top (k=1) to near-surface (k=nlayers)
            @test all(Φ[:, 1] .> Φ[:, end])
            for k in 1:(nlayers - 1)
                @test all(Φ[:, k] .> Φ[:, k + 1])
            end
        end
    end
end

@testset "Grid-space geopotential: orography raises geopotential at all levels" begin
    spectral_grid = SpectralGrid(nlayers = 4)
    model_flat = PrimitiveDryModel(spectral_grid; orography = NoOrography(spectral_grid))
    model_oro  = PrimitiveDryModel(spectral_grid; orography = NoOrography(spectral_grid))

    sim_flat = initialize!(model_flat)
    sim_oro  = initialize!(model_oro)

    T₀ = 280f0
    p_surf = 1f5
    for (vars, model) in ((sim_flat.variables, model_flat), (sim_oro.variables, model_oro))
        vars.grid.temperature .= T₀
        vars.parameterizations.surface_pressure .= p_surf
    end

    # add a uniform orography of 1 km to the second model
    z_oro = 1000f0
    sim_oro.variables.dynamics.geopotential  # touch to ensure allocated
    model_oro.orography.orography .= z_oro

    SpeedyWeather.geopotential!(sim_flat.variables, model_flat)
    SpeedyWeather.geopotential!(sim_oro.variables, model_oro)

    Φ_flat = Array(sim_flat.variables.dynamics.geopotential)
    Φ_oro  = Array(sim_oro.variables.dynamics.geopotential)
    g = model_flat.planet.gravity

    # uniform orography shifts geopotential at every level and grid point by g * z_oro
    @test Φ_oro ≈ Φ_flat .+ g * z_oro rtol = 1e-5
end

@testset "Grid-space geopotential: warmer atmosphere gives higher geopotential" begin
    # In sigma coordinates with constant T, Φ(k) = R*T*log(1/σ_k), which is
    # independent of surface pressure but directly proportional to temperature.
    spectral_grid = SpectralGrid(nlayers = 4)
    model = PrimitiveDryModel(spectral_grid; orography = NoOrography(spectral_grid))
    simulation = initialize!(model)
    vars = simulation.variables
    vars.parameterizations.surface_pressure .= 1f5

    vars.grid.temperature .= 250f0
    SpeedyWeather.geopotential!(vars, model)
    Φ_cold = copy(Array(vars.dynamics.geopotential))

    vars.grid.temperature .= 310f0
    SpeedyWeather.geopotential!(vars, model)
    Φ_warm = Array(vars.dynamics.geopotential)

    # warmer column expands → higher geopotential at every level
    @test all(Φ_warm .> Φ_cold)

    # proportionality: same log factor, so ratio equals temperature ratio
    @test Φ_warm ≈ Φ_cold .* (310f0 / 250f0) rtol = 1e-5
end

@testset "Grid-space geopotential: SigmaPressureCoordinates with A=0 matches SigmaCoordinates" begin
    # SigmaPressureCoordinates with transition = 1 everywhere (A=0, B=σ) is
    # mathematically identical to SigmaCoordinates; the two kernels must agree.
    spectral_grid = SpectralGrid(nlayers = 12)
    σ = SigmaCoordinates(spectral_grid)
    S = SigmaPressureCoordinates(spectral_grid; transition = _ -> 1.0)

    model_σ = PrimitiveDryModel(spectral_grid;
        geometry  = Geometry(spectral_grid; vertical_coordinates = σ),
    )
    model_S = PrimitiveDryModel(spectral_grid;
        geometry  = Geometry(spectral_grid; vertical_coordinates = S),
    )

    sim_σ = initialize!(model_σ)
    sim_S = initialize!(model_S)

    initialize!(sim_σ)
    initialize!(sim_S)

    # above executes this:
    # SpeedyWeather.geopotential!(sim_σ.variables, model_σ)
    # SpeedyWeather.geopotential!(sim_S.variables, model_S)

    Φ_σ = sim_σ.variables.dynamics.geopotential
    Φ_S = sim_S.variables.dynamics.geopotential

    @test Φ_σ ≈ Φ_S rtol = 1e-4
end
