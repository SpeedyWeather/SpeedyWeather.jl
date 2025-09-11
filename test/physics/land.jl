@testset "Dry land models" begin
    spectral_grid = SpectralGrid(trunc=31, nlayers=8)

    for Temperature in (SeasonalLandTemperature, ConstantLandTemperature, LandBucketTemperature)
        for Model in (PrimitiveDryModel, PrimitiveWetModel)

            temperature = Temperature(spectral_grid)
            land = DryLandModel(spectral_grid; temperature)

            model = Model(spectral_grid; land)
            simulation = initialize!(model)
            
            progn = simulation.prognostic_variables
            diagn = simulation.diagnostic_variables
            SpeedyWeather.land_timestep!(progn, diagn, model)
        end
    end
end

@testset "Wet land models" begin
    spectral_grid = SpectralGrid(trunc=31, nlayers=8)

    for Temperature in (SeasonalLandTemperature, ConstantLandTemperature, LandBucketTemperature)
        for SoilMoisture in (Nothing, SeasonalSoilMoisture, LandBucketMoisture)
            for Vegetation in (NoVegetation, VegetationClimatology)
                for Model in (PrimitiveDryModel, PrimitiveWetModel)

                    temperature = Temperature(spectral_grid)
                    soil_moisture = SoilMoisture(spectral_grid)
                    vegetation = Vegetation(spectral_grid)
                    land = LandModel(spectral_grid; temperature, soil_moisture, vegetation)
                    model = Model(spectral_grid; land)

                    # just test that no errors are thrown
                    initialize!(land, model)
                    progn = PrognosticVariables(spectral_grid)
                    diagn = DiagnosticVariables(spectral_grid, model)
                    SpeedyWeather.land_timestep!(progn, diagn, model)
                end
            end
        end
    end
end

@testset "With or without snow" begin
    spectral_grid = SpectralGrid(trunc=31, nlayers=8)

    for Snow in (Nothing, SnowModel)
        for Model in (PrimitiveDryModel, PrimitiveWetModel)
    
            snow = Snow(spectral_grid)
            land = LandModel(spectral_grid; snow)
            model = Model(spectral_grid; land)

            # just test that no errors are thrown
            initialize!(land, model)
            progn = PrognosticVariables(spectral_grid)
            diagn = DiagnosticVariables(spectral_grid, model)
            SpeedyWeather.land_timestep!(progn, diagn, model)

            @test all(isfinite.(progn.land.snow_depth))
        end
    end
end