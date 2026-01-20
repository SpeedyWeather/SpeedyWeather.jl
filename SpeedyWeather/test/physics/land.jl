@testset "Dry land models" begin
    spectral_grid = SpectralGrid(trunc = 31, nlayers = 8)

    for Temperature in (SeasonalLandTemperature, ConstantLandTemperature, LandBucketTemperature)
        for Model in (PrimitiveDryModel, PrimitiveWetModel)

            geometry = LandGeometry(spectral_grid)
            temperature = Temperature(spectral_grid, geometry)
            land = DryLandModel(spectral_grid; temperature)

            model = Model(spectral_grid; land)
            initialize!(model.land, model)

            progn = PrognosticVariables(model)
            diagn = DiagnosticVariables(model)
            SpeedyWeather.land_timestep!(progn, diagn, model)
        end
    end
end

@testset "Wet land models" begin
    spectral_grid = SpectralGrid(trunc = 31, nlayers = 8)

    for Temperature in (SeasonalLandTemperature, ConstantLandTemperature, LandBucketTemperature)
        for SoilMoisture in (Nothing, SeasonalSoilMoisture, LandBucketMoisture)
            for Vegetation in (NoVegetation, VegetationClimatology)
                for Model in (PrimitiveDryModel, PrimitiveWetModel)
                    
                    geometry = LandGeometry(spectral_grid)
                    temperature = Temperature(spectral_grid, geometry)
                    soil_moisture = SoilMoisture(spectral_grid, geometry)
                    vegetation = Vegetation(spectral_grid, geometry)
                    land = LandModel(spectral_grid; geometry, temperature, soil_moisture, vegetation)
                    model = Model(spectral_grid; land)

                    # just test that no errors are thrown
                    initialize!(land, model)
                    progn = PrognosticVariables(model)
                    diagn = DiagnosticVariables(model)
                    SpeedyWeather.land_timestep!(progn, diagn, model)
                end
            end
        end
    end
end

@testset "With or without snow" begin
    spectral_grid = SpectralGrid(trunc = 31, nlayers = 8)

    for Snow in (Nothing, SnowModel)
        for Model in (PrimitiveDryModel, PrimitiveWetModel)

            geometry = LandGeometry(spectral_grid)
            snow = Snow(spectral_grid, geometry)
            land = LandModel(spectral_grid; geometry, snow)
            model = Model(spectral_grid; land)

            # just test that no errors are thrown
            initialize!(land, model)
            progn = PrognosticVariables(model)
            diagn = DiagnosticVariables(model)
            SpeedyWeather.land_timestep!(progn, diagn, model)

            if snow isa SnowModel
                @test all(isfinite.(progn.land.snow_depth))
            else
                @test !haskey(progn.land, :snow_depth)
            end
        end
    end
end

@testset "LandGeometry default constructor" begin
    SG = SpectralGrid(trunc = 21, nlayers = 2)
    geom = LandGeometry(SG)
    @test geom.layer_thickness isa Vector{<:AbstractFloat}

    SG = SpectralGrid(trunc = 21, nlayers = 5)
    geom = LandGeometry(SG)
    @test geom.layer_thickness isa Vector{<:AbstractFloat}
end
