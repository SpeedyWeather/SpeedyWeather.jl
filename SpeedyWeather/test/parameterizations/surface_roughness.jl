@testset "ConstantSurfaceRoughness" begin
    @testset for Model in (PrimitiveWetModel, PrimitiveDryModel)
        spectral_grid = SpectralGrid(trunc = 31, nlayers = 8)
        NF = spectral_grid.NF

        # use non-default roughness lengths to make sure they are actually used
        z₀_land = NF(0.3)
        z₀_ocean = NF(2.0e-4)
        surface_roughness = ConstantSurfaceRoughness{NF}(
            roughness_length_land = z₀_land,
            roughness_length_ocean = z₀_ocean,
        )

        # AQUA PLANET: land_fraction == 0 everywhere -> pure ocean roughness
        model = Model(spectral_grid; surface_roughness, land_sea_mask = AquaPlanetMask(spectral_grid))
        initialize!(model.land_sea_mask, model)
        initialize!(surface_roughness, model)
        vars = Variables(model)
        SpeedyWeather.parameterization_tendencies!(vars, model)

        @test all(==(z₀_ocean), vars.parameterizations.surface_roughness)
        @test all(==(z₀_ocean), vars.parameterizations.ocean.surface_roughness)
        @test all(==(zero(NF)), vars.parameterizations.land.surface_roughness)

        # ROCKY PLANET: land_fraction == 1 everywhere -> pure land roughness
        model = Model(spectral_grid; surface_roughness, land_sea_mask = RockyPlanetMask(spectral_grid))
        initialize!(model.land_sea_mask, model)
        vars = Variables(model)
        SpeedyWeather.parameterization_tendencies!(vars, model)

        @test all(==(z₀_land), vars.parameterizations.surface_roughness)
        @test all(==(z₀_land), vars.parameterizations.land.surface_roughness)
        @test all(==(zero(NF)), vars.parameterizations.ocean.surface_roughness)

        # EARTH: fractional land-sea mask -> area-weighted blend of land and ocean
        model = Model(spectral_grid; surface_roughness)     # default EarthLandSeaMask
        initialize!(model.land_sea_mask, model)
        vars = Variables(model)
        SpeedyWeather.parameterization_tendencies!(vars, model)

        sr = vars.parameterizations.surface_roughness
        lsr = vars.parameterizations.land.surface_roughness
        osr = vars.parameterizations.ocean.surface_roughness
        mask = model.land_sea_mask.mask

        # the combined roughness is the area-weighted blend of land and ocean roughness
        @test sr ≈ @. mask * lsr + (1 - mask) * osr

        # all roughness lengths are finite and the blend stays strictly positive
        @test all(isfinite, sr)
        @test all(>(0), sr)
        @test all(>=(0), lsr)
        @test all(>=(0), osr)
    end
end