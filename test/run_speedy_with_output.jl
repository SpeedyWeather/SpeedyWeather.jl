@testset "Output on various grids" begin
    p = run_speedy(Float64,output=true)
    @test all(isfinite.(p.layers[1].leapfrog[1].vor))

    p = run_speedy(Float32,output=true)
    @test all(isfinite.(p.layers[1].leapfrog[1].vor))

    p = run_speedy(Float64,Grid=FullClenshawGrid,output=true)
    @test all(isfinite.(p.layers[1].leapfrog[1].vor))

    p = run_speedy(Float64,Grid=OctahedralGaussianGrid,output=true)
    @test all(isfinite.(p.layers[1].leapfrog[1].vor))

    p = run_speedy(Float64,Grid=OctahedralClenshawGrid,output=true)
    @test all(isfinite.(p.layers[1].leapfrog[1].vor))

    p = run_speedy(Float64,Grid=OctahedralClenshawGrid,output_matrix=true,output=true)
    @test all(isfinite.(p.layers[1].leapfrog[1].vor))

    p = run_speedy(Float64,Grid=OctahedralClenshawGrid,output_matrix=true,output_NF=Float32,output=true)
    @test all(isfinite.(p.layers[1].leapfrog[1].vor))
end

@testset "Restart from output file" begin
    
    p1, d1, m1 = initialize_speedy(Float32, ShallowWater, output=true, run_id="restart-test")
    run_speedy!(p1, d1, m1)
 
    p2, d2, m2 = initialize_speedy(Float32, ShallowWaterModel, initial_conditions=StartFromFile, restart_id="restart-test")

    for varname in propertynames(p1.layers[1].leapfrog[1])
        if SpeedyWeather.has(p1, varname)
            for (var_new, var_old) in zip(SpeedyWeather.get_var(p1, varname), SpeedyWeather.get_var(p2, varname))
                @test all(var_new .== var_old)
            end
        end
    end
    @test all(SpeedyWeather.get_pressure(p1) .== SpeedyWeather.get_pressure(p2))
end 