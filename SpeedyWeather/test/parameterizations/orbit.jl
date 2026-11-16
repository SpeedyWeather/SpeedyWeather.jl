@testset "Kepler equation" begin
    # residual of Kepler's equation E - e sin(E) = M after the fixed Newton iterations
    for e in 0:0.05:0.8, M in range(-π, 3π, length = 41)
        E = SpeedyWeather.eccentric_anomaly(M, e)
        @test abs(E - e * sin(E) - M) < 1.0e-10
    end
end

@testset "Orbit position" begin
    spectral_grid = SpectralGrid()
    zenith = SolarZenith(spectral_grid)
    T = Float64

    # one year of uniformly spaced times (2001 is not a leap year)
    times = DateTime(2001, 1, 1) .+ Minute.(round.(Int, range(0, 365 * 24 * 60, length = 3651)[1:(end - 1)]))

    @testset "CircularOrbit" begin
        planet = Earth(spectral_grid)
        @test planet.orbit isa CircularOrbit
        for time in times[1:100:end]
            declination, distance_factor = SpeedyWeather.orbit_position(T, planet.orbit, zenith, planet, time)
            @test declination == SpeedyWeather.solar_declination(T, zenith, time)
            @test distance_factor == 1
        end
    end

    @testset "KeplerOrbit with zero eccentricity" begin
        planet = Earth(spectral_grid; orbit = KeplerOrbit(eccentricity = 0.0, longitude_of_perihelion = 100.0))
        tilt = deg2rad(T(planet.axial_tilt))
        g_equinox = SpeedyWeather.equinox_angle(T, planet.equinox)
        for time in times[1:100:end]
            declination, distance_factor = SpeedyWeather.orbit_position(T, planet.orbit, zenith, planet, time)
            g = SpeedyWeather.year_angle(T, zenith, time)
            # circular orbit: the true solar longitude advances uniformly from the vernal equinox
            @test declination ≈ asin(sin(tilt) * sin(g - g_equinox)) atol = 1.0e-10
            @test distance_factor ≈ 1
        end
    end

    @testset for e in (0.0167, 0.1, 0.3), ϖ in (0.0, 102.9, 282.9)
        planet = Earth(spectral_grid; orbit = KeplerOrbit(eccentricity = e, longitude_of_perihelion = ϖ))
        tilt = deg2rad(T(planet.axial_tilt))

        # vernal equinox: true solar longitude 0, so no declination. Up to the leap-year offset
        # (≲ half a day) as the equinox phase uses the mean year like SinSolarDeclination does,
        # an offset that grows near perihelion where the solar longitude advances fastest
        declination, _ = SpeedyWeather.orbit_position(T, planet.orbit, zenith, planet, planet.equinox)
        @test declination ≈ 0 atol = 5.0e-3

        positions = [SpeedyWeather.orbit_position(T, planet.orbit, zenith, planet, time) for time in times]
        declinations = first.(positions)
        distance_factors = last.(positions)

        # time-mean of (a/r)² over an orbit is 1/√(1 - e²)
        @test sum(distance_factors) / length(distance_factors) ≈ 1 / sqrt(1 - e^2) rtol = 1.0e-4

        # distance factor between aphelion and perihelion
        @test maximum(distance_factors) ≈ 1 / (1 - e)^2 rtol = 1.0e-4
        @test minimum(distance_factors) ≈ 1 / (1 + e)^2 rtol = 1.0e-4

        # at perihelion the true solar longitude is the longitude of perihelion
        perihelion = argmax(distance_factors)
        @test declinations[perihelion] ≈ asin(sin(tilt) * sind(ϖ)) atol = 2.0e-3

        # declination stays within ± tilt
        @test all(abs.(declinations) .<= tilt + 1.0e-10)
    end

    @testset "Seasonal ozone follows the orbit" begin
        planet = Earth(spectral_grid; orbit = KeplerOrbit(eccentricity = 0.3, longitude_of_perihelion = 90.0))
        time = DateTime(2001, 5, 1)
        @test SpeedyWeather.solar_declination(T, zenith, time, planet) ==
            SpeedyWeather.orbit_position(T, planet.orbit, zenith, planet, time)[1]
    end
end

@testset "Simulation with an eccentric orbit" begin
    spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)

    # default: circular orbit, distance factor 1
    model = PrimitiveWetModel(spectral_grid)
    simulation = initialize!(model)
    run!(simulation, period = Day(1))
    @test simulation.variables.parameterizations.solar_distance_factor[] == 1

    for e in (0.1, 0.3)
        planet = Earth(spectral_grid; orbit = KeplerOrbit(eccentricity = e))
        model = PrimitiveWetModel(spectral_grid; planet)
        simulation = initialize!(model)
        run!(simulation, period = Day(5))
        distance_factor = simulation.variables.parameterizations.solar_distance_factor[]
        @test 1 / (1 + e)^2 - 1.0e-4 <= distance_factor <= 1 / (1 - e)^2 + 1.0e-4
        @test distance_factor != 1
        @test !simulation.model.feedback.nans_detected
    end
end
