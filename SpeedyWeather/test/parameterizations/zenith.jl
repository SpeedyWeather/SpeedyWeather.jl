@testset "SolarZenith" begin
    @testset for Z in (SolarZenith, SolarZenithSeason)
        spectral_grid = SpectralGrid()

        # create directly
        zenith = Z(spectral_grid)
        model = PrimitiveDryModel(spectral_grid; solar_zenith = zenith)
        simulation = initialize!(model, time = DateTime(2000, 6, 21))

        @test zenith.initial_time[] == DateTime(2000, 6, 21)

        # create via planet
        planet = Earth(spectral_grid; daily_cycle = true, seasonal_cycle = true)
        model = PrimitiveDryModel(spectral_grid; planet)
        @test model.solar_zenith isa SolarZenith

        planet = Earth(spectral_grid; daily_cycle = false, seasonal_cycle = true)
        model = PrimitiveDryModel(spectral_grid; planet)
        @test model.solar_zenith isa SolarZenithSeason
    end
end

@testset "cos_zenith!" begin
    @testset for Z in (SolarZenith, SolarZenithSeason)
        spectral_grid = SpectralGrid()
        zenith = Z(spectral_grid)
        model = PrimitiveDryModel(spectral_grid; solar_zenith = zenith)
        initialize!(model, time = DateTime(2000, 6, 21))

        cos_zenith = zeros(Float32, spectral_grid.grid)

        # June solstice: exercises polar day/night paths and acos clamp
        SpeedyWeather.cos_zenith!(cos_zenith, model.solar_zenith, DateTime(2000, 6, 21), model.geometry)
        @test all(cos_zenith .>= 0)
        @test all(cos_zenith .<= 1)
        june = copy(cos_zenith)

        # December solstice: sun in southern hemisphere, spatial pattern should differ
        SpeedyWeather.cos_zenith!(cos_zenith, model.solar_zenith, DateTime(2000, 12, 21), model.geometry)
        @test all(cos_zenith .>= 0)
        @test all(cos_zenith .<= 1)
        @test cos_zenith != june
    end
end
