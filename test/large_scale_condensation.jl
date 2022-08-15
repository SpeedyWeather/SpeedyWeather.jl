@testset "humidity.jl" begin
    @testset "get_saturation_vapour_pressure!" begin
        NF = Float32
        _, diag, model = SpeedyWeather.initialize_speedy(NF)
        (;nlon, nlat, nlev) = model.geospectral.geometry
        (;sat_vap_pressure) = diag.parametrization_variables

        temp_grid = 200 .+ 150 * rand(NF, nlon, nlat, nlev)  # Typical values between 200-350K

        SpeedyWeather.get_saturation_vapour_pressure!(sat_vap_pressure, temp_grid, model)

        @test all(sat_vap_pressure .> 0.0) && all(sat_vap_pressure .< 500.0)
    end

    @testset "get_saturation_specific_humidity!" begin
        NF = Float32
        _, diag, model = SpeedyWeather.initialize_speedy(NF)
        (;nlon, nlat, nlev) = model.geospectral.geometry
        (;sat_vap_pressure, sat_spec_humidity) = diag.parametrization_variables

        temp_grid = 200 .+ 150 * rand(NF, nlon, nlat, nlev)  # Typical values between 200-350 K
        pres_grid = 300 .* 1700 * rand(NF, nlon, nlat)       # Typical values between 300-2000 hPa

        SpeedyWeather.get_saturation_specific_humidity!(sat_spec_humidity, sat_vap_pressure, temp_grid, pres_grid, model)

        @test all(isfinite.(sat_spec_humidity))
        @test !any(iszero.(sat_spec_humidity))
    end
end


@testset "large_scale_condensation.jl" begin
    @testset "get_large_scale_condensation_tendencies!" begin
        _, diag, model = SpeedyWeather.initialize_speedy()

        # For now, just check that it runs without errors
        SpeedyWeather.get_large_scale_condensation_tendencies!(diag, model)
    end
end
