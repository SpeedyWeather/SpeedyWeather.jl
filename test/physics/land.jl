@testset "Dry land models" begin
    spectral_grid = SpectralGrid(trunc=31, nlayers=8)

    for Temperature in (SeasonalLandTemperature, ConstantLandTemperature, LandBucketTemperature)
        for Model in (PrimitiveDryModel, PrimitiveWetModel)

            temperature = Temperature(spectral_grid)
            land = DryLandModel(spectral_grid; temperature)

            model = Model(spectral_grid; land)
            simulation = initialize!(model)
            run!(simulation, period=Hour(2))
        end
    end
end

@testset "Wet land models" begin
    spectral_grid = SpectralGrid(trunc=31, nlayers=8)

    for Temperature in (SeasonalLandTemperature, ConstantLandTemperature, LandBucketTemperature)
        for SoilMoisture in (NoSoilMoisture, SeasonalSoilMoisture, LandBucketMoisture)
            for Vegetation in (NoVegetation, VegetationClimatology)
                for Model in (PrimitiveDryModel, PrimitiveWetModel)

                    temperature = Temperature(spectral_grid)
                    soil_moisture = SoilMoisture(spectral_grid)
                    vegetation = Vegetation(spectral_grid)
                    land = LandModel(spectral_grid; temperature, soil_moisture, vegetation)
    
                    model = Model(spectral_grid; land)
                    simulation = initialize!(model)
                    run!(simulation, period=Hour(2))
                end
            end
        end
    end
end