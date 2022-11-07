@testset "Parametrization: convection" begin
    @testset "diagnose_convection!" begin
        @testset for NF in (Float32, Float64)
            _, diagn, model = SpeedyWeather.initialize_speedy(NF, model=PrimitiveEquation)
            nlev = diagn.nlev
            pres_thresh_cnv = model.constants.pres_thresh_cnv

            column = ColumnVariables{NF}(nlev = nlev)
            column.humid .= 0 .+ 50 * rand(NF, nlev)         # Typical values between 0-50 g/kg
            column.pres = pres_thresh_cnv + 0.1              # Normalised pressure above threshold for convection
            column.sat_humid .= 0 .+ 50 * rand(NF, nlev)     # Typical values between 0-50 g/kg
            column.dry_static_energy .= rand(NF, nlev)
            column.moist_static_energy .= rand(NF, nlev)
            column.sat_moist_static_energy .= rand(NF, nlev)
            column.sat_moist_static_energy_half .= rand(NF, nlev)

            SpeedyWeather.diagnose_convection!(column, model)

            @test isfinite(column.excess_humidity)
            @test column.excess_humidity >= 0
        end
    end

    @testset "convection!" begin
        @testset for NF in (Float32, Float64)
            _, diagn, model = SpeedyWeather.initialize_speedy(NF, model=PrimitiveEquation)
            nlev = diagn.nlev
            pres_thresh_cnv = model.constants.pres_thresh_cnv

            column = ColumnVariables{NF}(nlev = nlev)
            column.humid .= 0 .+ 50 * rand(NF, nlev)         # Typical values between 0-50 g/kg
            column.pres = pres_thresh_cnv + 0.1              # Normalised pressure above threshold for convection
            column.sat_humid .= 0 .+ 50 * rand(NF, nlev)     # Typical values between 0-50 g/kg
            column.dry_static_energy .= rand(NF, nlev)
            column.moist_static_energy .= rand(NF, nlev)
            column.sat_moist_static_energy .= rand(NF, nlev)
            column.sat_moist_static_energy_half .= rand(NF, nlev)

            SpeedyWeather.convection!(column, model)

            @test isfinite(column.cloud_base_mass_flux)
            @test isfinite(column.precip_convection)
            @test all(isfinite.(column.net_flux_humid))
            @test all(isfinite.(column.net_flux_dry_static_energy))
        end
    end
end
