@testset "ColumnVariables initialisation" begin
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

        # Large-scale condensation
        @test column.precip_large_scale === zero(NF)
        @test column.precip_convection === zero(NF)
    end
end

@testset "ColumnVariables initialisation" begin
    @testset for NF in (Float32,Float64)

        nlev = 8
        spectral_grid = SpectralGrid(NF;nlev)
        model = PrimitiveDryModel(;spectral_grid)
        simulation = initialize!(model)
        diagn = simulation.diagnostic_variables
        progn = simulation.prognostic_variables

        column = ColumnVariables{NF}(;nlev)

        SpeedyWeather.reset_column!(column)
        SpeedyWeather.get_column!(column,diagn,progn,1,model)

        # set a tendency to something
        humid_tend = rand(NF,nlev)
        column.humid_tend .= humid_tend

        # copy into diagn
        SpeedyWeather.write_column_tendencies!(diagn,column,model.planet,1)

        # and check that that worked
        for (k,layer) in enumerate(diagn.layers)
            @test humid_tend[k] == layer.tendencies.humid_tend_grid[1]
        end
    end
end
