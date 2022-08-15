@testset "Parametrization: Saturation vapour pressure" begin
    @testset for NF in (Float32,Float64)
        _, diagn, model = SpeedyWeather.initialize_speedy(NF,model=:primitive)
        (;nlev) = model.geometry

        column = ColumnVariables{NF}(;nlev)
        column.temp .= 200 .+ 150 * rand(NF, nlev)  # Typical values between 200-350K

        SpeedyWeather.get_saturation_vapour_pressure!(column, model)
        (;sat_vap_pres) = column

        @test all(sat_vap_pres .> 0.0)
        @test all(sat_vap_pres .< 500.0)
    end
end

@testset "Parametrization: Saturation specific humidity" begin
    @testset for NF in (Float32,Float64)
        _, diag, model = SpeedyWeather.initialize_speedy(NF,model=:primitive)
        (;nlev) = model.geometry

        column = ColumnVariables{NF}(;nlev)
        (;sat_vap_pres, sat_spec_humid) = column

        column.temp = 200 .+ 150*rand(NF, nlev)     # Typical values between 200-350 K
        column.pres = 300 + 1700*rand(NF)           # Typical values between 300-2000 hPa

        SpeedyWeather.get_saturation_vapour_pressure!(column, model)
        SpeedyWeather.get_saturation_specific_humidity!(column, model)

        @test all(isfinite.(sat_spec_humid))
        @test !any(iszero.(sat_spec_humid))
    end
end

@testset "Parametrization: large scale condensation" begin
    @testset for NF in (Float32,Float64)
        _, diagn, model = SpeedyWeather.initialize_speedy(NF,model=:primitive)

        column = ColumnVariables{NF}(nlev=diagn.nlev)

        for ij in SpeedyWeather.eachgridpoint(diagn)
            SpeedyWeather.reset_column!(column)
            SpeedyWeather.get_column!(column,diagn,ij,model.geometry)
            SpeedyWeather.large_scale_condensation!(column, model)
        end
    end
end
