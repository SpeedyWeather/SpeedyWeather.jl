@testset "ColumnVariables initialisation" begin

    # no number format provided
    # column = ColumnVariables(nlev=8)
    # @test eltype(column.temp) == Float64

    @testset for NF in (Float16,Float32,Float64)
        column = ColumnVariables{NF}(nlev=8)

        @test eltype(column.temp) == NF
        SpeedyWeather.reset_column!(column)

        @test all(column.temp_tend .=== zero(NF))
        @test all(column.humid_tend .=== zero(NF))
        @test all(column.u_tend .=== zero(NF))
        @test all(column.v_tend .=== zero(NF))

        # Convection
        @test column.cloud_top === column.nlev+1
        @test column.conditional_instability === false
        @test column.activate_convection === false
        @test column.excess_humidity === zero(NF)
        @test column.precip_convection === zero(NF)
        @test column.cloud_base_mass_flux === zero(NF)
        @test all(column.net_flux_humid .=== zero(NF))
        @test all(column.net_flux_dry_static_energy .=== zero(NF))

        # Large-scale condensation
        @test column.precip_large_scale === zero(NF)
    end
end

@testset "ColumnVariables initialisation" begin
    @testset for NF in (Float32,Float64)

        nlev = 8
        _,diagn,model = initialize_speedy(NF,nlev=nlev,model=PrimitiveEquationModel)
        column = ColumnVariables{NF}(;nlev)

        SpeedyWeather.reset_column!(column)
        SpeedyWeather.get_column!(column,diagn,1,model.geometry)

        # set a tendency to something
        humid_tend = rand(NF,nlev)
        column.humid_tend .= humid_tend

        # copy into diagn
        SpeedyWeather.write_column_tendencies!(diagn,column,1)

        # and check that that worked
        for (k,layer) in enumerate(diagn.layers)
            @test humid_tend[k] == layer.tendencies.humid_grid_tend[1]
        end
    end
end
