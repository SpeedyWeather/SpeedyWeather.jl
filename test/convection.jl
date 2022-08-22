@testset "Parametrization: convection" begin
    @testset "conditional_instability!" begin
        @testset for NF in (Float32, Float64)
            _, diagn, model = SpeedyWeather.initialize_speedy(NF, model = :primitive)
            nlev = diagn.nlev

            column = ColumnVariables{NF}(nlev = nlev)
            column.humid .= rand(NF, nlev)
            column.pres = rand(NF)
            column.sat_humid .= rand(NF, nlev)
            column.dry_static_energy .= rand(NF, nlev)
            column.moist_static_energy .= rand(NF, nlev)
            column.sat_moist_static_energy .= rand(NF, nlev)
            column.sat_moist_static_energy_half .= rand(NF, nlev)

            SpeedyWeather.conditional_instability!(column, model)
        end
    end

    @testset "convection!" begin
        @testset for NF in (Float32, Float64)
            _, diagn, model = SpeedyWeather.initialize_speedy(NF, model = :primitive)
            nlev = diagn.nlev

            column = ColumnVariables{NF}(nlev = nlev)
            column.humid .= rand(NF, nlev)
            column.pres = rand(NF)
            column.sat_humid .= rand(NF, nlev)
            column.dry_static_energy .= rand(NF, nlev)
            column.moist_static_energy .= rand(NF, nlev)
            column.sat_moist_static_energy .= rand(NF, nlev)
            column.sat_moist_static_energy_half .= rand(NF, nlev)

            SpeedyWeather.convection!(column, model)
        end
    end
end
