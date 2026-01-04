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
