import Parameters: @unpack
using Distributions

@testset "humidity.jl" begin
    @testset "get_saturation_vapour_pressure!" begin
        _, diag, model = SpeedyWeather.initialize_speedy()
        @unpack nlon, nlat, nlev = model.geospectral.geometry
        @unpack sat_vap_pressure = diag.parametrization_variables

        temp_grid = rand(Uniform(200, 350), nlon, nlat, nlev)  # Typical values in K

        SpeedyWeather.get_saturation_vapour_pressure!(sat_vap_pressure, temp_grid, model)

        @test all(sat_vap_pressure .> 0.0) && all(sat_vap_pressure .< 500.0)
    end

    @testset "get_saturation_specific_humidity!" begin
        _, diag, model = SpeedyWeather.initialize_speedy()
        @unpack nlon, nlat, nlev = model.geospectral.geometry
        @unpack sat_vap_pressure, sat_spec_humidity = diag.parametrization_variables

        temp_grid = rand(Uniform(200, 350), nlon, nlat, nlev)       # Typical values in K
        pres_grid = rand(Uniform(log(300), log(2000)), nlon, nlat)  # Typical values in log(hPa)

        SpeedyWeather.get_saturation_specific_humidity!(sat_spec_humidity, sat_vap_pressure, temp_grid, pres_grid, model)

        @test all(isfinite.(sat_spec_humidity))
        @test !any(iszero.(sat_spec_humidity))
    end
end
