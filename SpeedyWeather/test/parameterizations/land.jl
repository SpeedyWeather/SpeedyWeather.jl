@testset "Dry land models" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)

    for Temperature in (SeasonalLandTemperature, ConstantLandTemperature, LandBucketTemperature)
        for Model in (PrimitiveDryModel, PrimitiveWetModel)

            geometry = LandGeometry(spectral_grid)
            temperature = Temperature(spectral_grid, geometry)
            land = DryLandModel(spectral_grid; temperature)

            model = Model(spectral_grid; land)
            initialize!(model.land, model)

            vars = Variables(model)
            SpeedyWeather.land_timestep!(vars, model)
        end
    end
end

@testset "Wet land models" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)

    for Temperature in (SeasonalLandTemperature, ConstantLandTemperature, LandBucketTemperature)
        for SoilMoisture in (Nothing, SeasonalSoilMoisture, LandBucketMoisture)
            for Vegetation in (NoVegetation, VegetationClimatology)
                for Model in (PrimitiveDryModel, PrimitiveWetModel)

                    # constructors construct without geometry
                    temperature = Temperature(spectral_grid)
                    soil_moisture = SoilMoisture(spectral_grid)
                    vegetation = Vegetation(spectral_grid)

                    # ... and with geometry
                    geometry = LandGeometry(spectral_grid)
                    temperature = Temperature(spectral_grid, geometry)
                    soil_moisture = SoilMoisture(spectral_grid, geometry)
                    vegetation = Vegetation(spectral_grid, geometry)
                    land = LandModel(spectral_grid; geometry, temperature, soil_moisture, vegetation)
                    model = Model(spectral_grid; land)

                    # just test that no errors are thrown
                    initialize!(land, model)
                    vars = Variables(model)
                    SpeedyWeather.land_timestep!(vars, model)
                end
            end
        end
    end
end

@testset "With or without snow" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)

    for Snow in (Nothing, SnowModel)
        for Model in (PrimitiveDryModel, PrimitiveWetModel)

            geometry = LandGeometry(spectral_grid)
            snow = Snow(spectral_grid) # redundant but tests constructor
            snow = Snow(spectral_grid, geometry)
            land = LandModel(spectral_grid; geometry, snow)
            model = Model(spectral_grid; land)

            # just test that no errors are thrown
            initialize!(land, model)
            vars = Variables(model)
            SpeedyWeather.land_timestep!(vars, model)

            if snow isa SnowModel
                @test all(isfinite.(vars.prognostic.land.snow_depth))
            else
                @test !haskey(vars.prognostic.land, :snow_depth)
            end
        end
    end
end

@testset "Initial snow depth from soil temperature" begin
    spectral_grid = SpectralGrid(truncation = 31, nlayers = 8)
    model = PrimitiveWetModel(spectral_grid)
    simulation = initialize!(model)
    snow_depth = simulation.variables.prognostic.land.snow_depth[:, 1].data          # first step
    soil_temperature = simulation.variables.prognostic.land.soil_temperature[:, 1, 1].data
    (; melting_threshold, initial_snow_depth_max) = model.land.snow
    land = model.land_sea_mask.land_fraction.data .> 0

    @test all(0 .<= snow_depth .<= initial_snow_depth_max)
    @test all(snow_depth[.!land] .== 0)                                     # no snow over ocean
    @test all(snow_depth[land .& (soil_temperature .>= melting_threshold)] .== 0)
    @test any(snow_depth[land] .> 0)                                        # some snow over cold land
end

@testset "LandGeometry default constructor" begin
    SG = SpectralGrid(truncation = 22, nlayers = 2)
    geom = LandGeometry(SG)
    @test geom.layer_thickness isa Vector{<:AbstractFloat}

    SG = SpectralGrid(truncation = 22, nlayers = 5)
    geom = LandGeometry(SG)
    @test geom.layer_thickness isa Vector{<:AbstractFloat}
end
