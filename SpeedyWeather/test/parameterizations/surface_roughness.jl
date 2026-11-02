@testset "ConstantSurfaceRoughness" begin
    @testset for Model in (PrimitiveWetModel, PrimitiveDryModel)
        spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)

        # use non-default roughness lengths to make sure they are actually used
        z₀M_land = 0.3
        z₀H_land = 0.2
        z₀Q_land = 0.1
        z₀M_ocean = 2.0e-4
        z₀H_ocean = 3.0e-4
        z₀Q_ocean = 4.0e-4

        surface_roughness = ConstantSurfaceRoughness(
            spectral_grid;
            momentum_roughness_land = z₀M_land,
            heat_roughness_land = z₀H_land,
            humidity_roughness_land = z₀Q_land,
            momentum_roughness_ocean = z₀M_ocean,
            heat_roughness_ocean = z₀H_ocean,
            humidity_roughness_ocean = z₀Q_ocean,
        )
        boundary_layer = BoundaryLayer(spectral_grid; surface_roughness)

        # AQUA PLANET: land_fraction == 0 everywhere -> pure ocean roughness
        model = Model(spectral_grid; boundary_layer, land_sea_mask = AquaPlanetMask(spectral_grid))
        initialize!(model.land_sea_mask, model)
        initialize!(boundary_layer, model)
        vars = Variables(model)
        SpeedyWeather.parameterization_tendencies!(vars, model)

        # unpack as ==(x) requires type equality too
        sr = boundary_layer.surface_roughness
        z₀M_ocean = sr.momentum_roughness_ocean
        z₀H_ocean = sr.heat_roughness_ocean
        z₀Q_ocean = sr.humidity_roughness_ocean
        z₀M_land = sr.momentum_roughness_land
        z₀H_land = sr.heat_roughness_land
        z₀Q_land = sr.humidity_roughness_land

        @test all(==(z₀M_ocean), vars.parameterizations.momentum_roughness)
        @test all(==(z₀M_ocean), vars.parameterizations.ocean.momentum_roughness)
        @test all(==(0), vars.parameterizations.land.momentum_roughness)

        @test all(==(z₀H_ocean), vars.parameterizations.heat_roughness)
        @test all(==(z₀H_ocean), vars.parameterizations.ocean.heat_roughness)
        @test all(==(0), vars.parameterizations.land.heat_roughness)

        @test all(==(z₀Q_ocean), vars.parameterizations.humidity_roughness)
        @test all(==(z₀Q_ocean), vars.parameterizations.ocean.humidity_roughness)
        @test all(==(0), vars.parameterizations.land.humidity_roughness)

        # ROCKY PLANET: land_fraction == 1 everywhere -> pure land roughness
        model = Model(spectral_grid; boundary_layer, land_sea_mask = RockyPlanetMask(spectral_grid))
        initialize!(model.land_sea_mask, model)
        vars = Variables(model)
        SpeedyWeather.parameterization_tendencies!(vars, model)

        @test all(==(z₀M_land), vars.parameterizations.momentum_roughness)
        @test all(==(z₀M_land), vars.parameterizations.land.momentum_roughness)
        @test all(==(0), vars.parameterizations.ocean.momentum_roughness)

        @test all(==(z₀H_land), vars.parameterizations.heat_roughness)
        @test all(==(z₀H_land), vars.parameterizations.land.heat_roughness)
        @test all(==(0), vars.parameterizations.ocean.heat_roughness)

        @test all(==(z₀Q_land), vars.parameterizations.humidity_roughness)
        @test all(==(z₀Q_land), vars.parameterizations.land.humidity_roughness)
        @test all(==(0), vars.parameterizations.ocean.humidity_roughness)

        # EARTH: fractional land-sea mask -> area-weighted blend of land and ocean
        model = Model(spectral_grid; boundary_layer)     # default EarthLandSeaMask
        initialize!(model.land_sea_mask, model)
        vars = Variables(model)
        SpeedyWeather.parameterization_tendencies!(vars, model)

        srM = vars.parameterizations.momentum_roughness
        lsrM = vars.parameterizations.land.momentum_roughness
        osrM = vars.parameterizations.ocean.momentum_roughness

        srH = vars.parameterizations.heat_roughness
        lsrH = vars.parameterizations.land.heat_roughness
        osrH = vars.parameterizations.ocean.heat_roughness

        srQ = vars.parameterizations.humidity_roughness
        lsrQ = vars.parameterizations.land.humidity_roughness
        osrQ = vars.parameterizations.ocean.humidity_roughness

        land_fraction = model.land_sea_mask.land_fraction

        # the combined roughness is the area-weighted blend of land and ocean roughness
        @test srM ≈ @. land_fraction * lsrM + (1 - land_fraction) * osrM
        @test srH ≈ @. land_fraction * lsrH + (1 - land_fraction) * osrH
        @test srQ ≈ @. land_fraction * lsrQ + (1 - land_fraction) * osrQ

        # all roughness lengths are finite and the blend stays strictly positive
        for (sr, lsr, osr) in (
                (srM, lsrM, osrM),
                (srH, lsrH, osrH),
                (srQ, lsrQ, osrQ),
            )
            @test all(isfinite, sr)
            @test all(>(0), sr)
            @test all(>=(0), lsr)
            @test all(>=(0), osr)
        end
    end
end

@testset "Learned ocean momentum roughness" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
    scheme = LearnedOceanRoughness(spectral_grid)

    for uₙ in (1.0, 5.0, 10.0, 20.0, 40.0, 60.0)
        α = log(scheme.charnock(uₙ))

        A = uₙ^scheme.m_c1
        B = (scheme.m_c3 * uₙ + α + exp(scheme.m_c2)) / (exp(scheme.m_c4 * uₙ - scheme.m_c5) + scheme.m_c6)
        C = scheme.m_c7 * α - scheme.m_c8

        z₀_base = exp(A - B + C)
        z₀_max = scheme.max_ocean_momentum_roughness(uₙ)
        z₀_expected = min(z₀_base, z₀_max)

        z₀_actual = SpeedyWeather.ocean_momentum_roughness(scheme, uₙ, α, 0.0)

        @test z₀_actual ≈ z₀_expected
        @test isfinite(z₀_actual)
        @test z₀_actual > 0
        @test z₀_actual <= z₀_max
    end
end

@testset "No neutral wind speed" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)

    # without neutral wind speed the learned roughness uses the surface wind speed instead
    @testset for surface_roughness in (
            ConstantSurfaceRoughness(spectral_grid),
            LearnedSurfaceRoughness(spectral_grid),
        )
        boundary_layer = BoundaryLayer(spectral_grid; surface_roughness, neutral_wind_speed = nothing)
        model = PrimitiveWetModel(spectral_grid; boundary_layer)
        simulation = initialize!(model)
        run!(simulation, steps = 2)

        (; parameterizations) = simulation.variables
        @test !haskey(parameterizations, :neutral_wind_speed)
        @test all(isfinite, parameterizations.boundary_layer_drag_momentum)
        @test all(>(0), parameterizations.boundary_layer_drag_momentum)
    end
end
