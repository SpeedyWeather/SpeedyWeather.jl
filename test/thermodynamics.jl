@testset "Thermodynamics" begin
    @testset "get_thermodynamics!" begin
        @testset for NF in (Float32, Float64)
            _, _, model = SpeedyWeather.initialize_speedy(NF, model = :primitive)
            (; nlev) = model.geometry

            column = ColumnVariables{NF}(; nlev)
            column.temp .= 200 .+ 150 * rand(NF, nlev)  # Typical values between 200-350K
            column.humid .= rand(NF, nlev)              # What are typical values for humidity?
            column.pres = 300 + 1700 * rand(NF)         # Typical values between 300-2000 hPa
            column.geopot .= rand(NF, nlev)             # What are typical values for geopotential?

            # For now, test that it runs with no errors
            SpeedyWeather.get_thermodynamics!(column, model)
        end
    end

    @testset "interpolate!" begin
        @testset for NF in (Float32, Float64)
            _, _, model = SpeedyWeather.initialize_speedy(NF, model = :primitive)
            (; nlev) = model.geometry

            column = ColumnVariables{NF}(; nlev)

            A_full_level = rand(NF, nlev)
            A_half_level = zeros(NF, nlev)

            SpeedyWeather.interpolate!(A_full_level, A_half_level, column, model)

            # Test that the half-level values lie between the enclosing full level values.
            @test all(
                (A_full_level[2:nlev] .> A_half_level[1:nlev-1] .> A_full_level[1:nlev-1])
                .||
                (A_full_level[2:nlev] .< A_half_level[1:nlev-1] .< A_full_level[1:nlev-1])
                )
        end
    end

    @testset "Saturation vapour pressure" begin
        @testset for NF in (Float32, Float64)
            _, diagn, model = SpeedyWeather.initialize_speedy(NF, model = :primitive)
            (; nlev) = model.geometry

            column = ColumnVariables{NF}(; nlev)
            column.temp .= 200 .+ 150 * rand(NF, nlev)  # Typical values between 200-350K

            SpeedyWeather.saturation_vapour_pressure!(column, model)
            (; sat_vap_pres) = column

            @test all(sat_vap_pres .> 0.0)
            @test all(sat_vap_pres .< 500.0)
        end
    end

    @testset "Saturation specific humidity" begin
        @testset for NF in (Float32, Float64)
            _, diag, model = SpeedyWeather.initialize_speedy(NF, model = :primitive)
            (; nlev) = model.geometry

            column = ColumnVariables{NF}(; nlev)
            (; sat_vap_pres, sat_humid) = column

            column.temp = 200 .+ 150 * rand(NF, nlev)     # Typical values between 200-350 K
            column.pres = 300 + 1700 * rand(NF)           # Typical values between 300-2000 hPa

            SpeedyWeather.saturation_vapour_pressure!(column, model)
            SpeedyWeather.saturation_specific_humidity!(column, model)

            @test all(isfinite.(sat_humid))
            @test !any(iszero.(sat_humid))
        end
    end

    @testset "Dry static energy" begin
        @testset for NF in (Float32, Float64)
            _, _, model = SpeedyWeather.initialize_speedy(NF, model = :primitive)
            (; nlev) = model.geometry

            column = ColumnVariables{NF}(; nlev)
            column.temp .= 200 .+ 150 * rand(NF, nlev)  # Typical values between 200-350K
            column.geopot .= rand(NF, nlev)

            SpeedyWeather.dry_static_energy!(column, model)
            (; dry_static_energy) = column

            @test all(isfinite.(dry_static_energy))
            @test !any(iszero.(dry_static_energy))
        end
    end

    @testset "Moist static energy" begin
        @testset for NF in (Float32, Float64)
            _, _, model = SpeedyWeather.initialize_speedy(NF, model = :primitive)
            (; nlev) = model.geometry

            column = ColumnVariables{NF}(; nlev)
            column.dry_static_energy .= rand(NF, nlev)
            column.humid .= rand(NF, nlev)

            SpeedyWeather.moist_static_energy!(column, model)
            (; moist_static_energy) = column

            @test all(isfinite.(moist_static_energy))
            @test !any(iszero.(moist_static_energy))
        end
    end

    @testset "Saturation moist static energy" begin
        @testset for NF in (Float32, Float64)
            _, _, model = SpeedyWeather.initialize_speedy(NF, model = :primitive)
            (; nlev) = model.geometry

            column = ColumnVariables{NF}(; nlev)
            column.dry_static_energy .= rand(NF, nlev)
            column.sat_humid .= rand(NF, nlev)

            SpeedyWeather.saturation_moist_static_energy!(column, model)
            (; sat_moist_static_energy) = column

            @test all(isfinite.(sat_moist_static_energy))
            @test !any(iszero.(sat_moist_static_energy))
        end
    end
end
