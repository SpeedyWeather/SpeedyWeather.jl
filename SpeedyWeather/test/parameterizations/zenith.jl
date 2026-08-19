using Dates

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

@testset "year_angle" begin
    @testset for NF in (Float32, Float64)
        # 2000 is a leap year (366 days), 2001 a common year (365 days)
        @testset for year in (2000, 2001)
            year_angle(time) = SpeedyWeather.year_angle(NF, time)

            # sawtooth from 0 at Jan-01T00:00 to 2π at the end of Dec-31
            @test year_angle(DateTime(year, 1, 1)) ≈ 0 atol = 1.0e-5
            @test year_angle(DateTime(year, 12, 31, 23, 59, 59)) ≈ 2π rtol = 1.0e-5

            # no jump across the New Year, only one second of a year (regardless of leap year)
            g_before = year_angle(DateTime(year, 12, 31, 23, 59, 59))
            g_after = year_angle(DateTime(year + 1, 1, 1))
            @test mod(g_after - g_before, 2π) < 1.0e-4

            # monotonically increasing throughout the year
            g = [year_angle(DateTime(year, month, 1)) for month in 1:12]
            @test issorted(g)
            @test all(0 .<= g .< 2π)

            # half way through the year is half way to 2π (Hour to stay exact for 365 days)
            mid_year = DateTime(year, 1, 1) + Hour(24 * Dates.daysinyear(year)) ÷ 2
            @test year_angle(mid_year) ≈ π rtol = 1.0e-5
        end
    end
end

@testset "solar_hour_angle" begin
    @testset for NF in (Float32, Float64)
        solar_hour_angle(time, λ = 0) = SpeedyWeather.solar_hour_angle(NF, time, λ)

        # sawtooth from -π at midnight through 0 at noon to π at the end of the day
        @test solar_hour_angle(DateTime(2000, 1, 1)) ≈ -π rtol = 1.0e-5
        @test solar_hour_angle(DateTime(2000, 1, 1, 12)) ≈ 0 atol = 1.0e-5
        # 23:59:59 is one second short of π, hence the looser tolerance
        @test solar_hour_angle(DateTime(2000, 1, 1, 23, 59, 59)) ≈ π rtol = 1.0e-4

        # no jump across midnight, only one second of a day
        h_before = solar_hour_angle(DateTime(2000, 1, 1, 23, 59, 59))
        h_after = solar_hour_angle(DateTime(2000, 1, 2))
        @test mod(h_after - h_before, 2π) < 1.0e-4

        # monotonically increasing throughout the day
        @test issorted([solar_hour_angle(DateTime(2000, 1, 1, hour)) for hour in 0:23])

        # longitude offsets the hour angle directly
        @test solar_hour_angle(DateTime(2000, 1, 1, 12), NF(1)) ≈ 1 rtol = 1.0e-5
    end
end

@testset "cos_zenith!" begin
    @testset for Z in (SolarZenith, SolarZenithSeason)
        spectral_grid = SpectralGrid()
        model = PrimitiveDryModel(spectral_grid; solar_zenith = Z(spectral_grid))

        cos_zenith = zeros(spectral_grid.NF, spectral_grid.grid)

        # June solstice: exercises polar day/night paths and acos clamp
        june, december = DateTime(2000, 6, 21), DateTime(2000, 12, 21)
        SpeedyWeather.cos_zenith!(cos_zenith, model.solar_zenith, june, june, model.geometry)
        @test all(0 .<= cos_zenith .<= 1)
        @test any(cos_zenith .> 0)     # sun shines somewhere
        cos_zenith_june = copy(cos_zenith)

        # December solstice: sun in southern hemisphere
        SpeedyWeather.cos_zenith!(cos_zenith, model.solar_zenith, december, december, model.geometry)
        @test all(0 .<= cos_zenith .<= 1)
        @test any(cos_zenith .> 0)     # sun shines somewhere
        @test cos_zenith != cos_zenith_june     # seasonal cycle changes the insolation

        # rotation and orbit time are used independently: same orbit time (season) but a
        # different rotation time (time of day) only changes the daily cycle. SolarZenithSeason
        # is a daily average and ignores rotation time, so it's back to the June insolation
        SpeedyWeather.cos_zenith!(cos_zenith, model.solar_zenith, june + Hour(6), june, model.geometry)
        @test all(0 .<= cos_zenith .<= 1)
        @test (cos_zenith == cos_zenith_june) == (Z == SolarZenithSeason)
    end
end

@testset "Clock rotation and orbit time" begin
    # set! synchronizes rotation and orbit time with the model time
    clock = SpeedyWeather.Clock()
    clock.rotation_time = DateTime(2010, 1, 1)      # deliberately out of sync
    SpeedyWeather.set!(clock, time = DateTime(2000, 6, 1), start = DateTime(2000, 6, 1))
    @test clock.time == clock.rotation_time == clock.orbit_time == DateTime(2000, 6, 1)

    spectral_grid = SpectralGrid(truncation = 31, nlayers = 4)

    # length of day and year, and the resulting rotation and orbit dilation
    cases = (
        (Hour(24), Day(365) + Hour(6), 1, 1),        # Earth
        (Hour(12), Day(100), 2, 3.6525),             # faster day and year
        (Hour(-24), Day(365) + Hour(6), -1, 1),      # backwards rotation
    )

    @testset for (length_of_day, length_of_year, rotation_dilation, orbit_dilation) in cases
        planet = Earth(spectral_grid; length_of_day, length_of_year)
        model = PrimitiveDryModel(spectral_grid; planet)
        model.feedback.verbose = false
        simulation = initialize!(model, time = DateTime(2000, 1, 1))
        clock = simulation.variables.prognostic.clock

        # all times are synchronized initially
        @test clock.time == clock.rotation_time == clock.orbit_time == DateTime(2000, 1, 1)

        run!(simulation, steps = 5)

        @test clock.rotation_dilation ≈ rotation_dilation rtol = 1.0e-5
        @test clock.orbit_dilation ≈ orbit_dilation rtol = 1.0e-5

        # rotation and orbit time advance at their dilated rate relative to the model time
        Δtime = Millisecond(clock.time - clock.start).value
        @test Millisecond(clock.rotation_time - clock.start).value ≈ rotation_dilation * Δtime rtol = 1.0e-5
        @test Millisecond(clock.orbit_time - clock.start).value ≈ orbit_dilation * Δtime rtol = 1.0e-5
    end
end
