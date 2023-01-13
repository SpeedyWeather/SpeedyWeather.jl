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
    
    p = run_speedy(Float32, initial_conditions=StartFromFile, restart_id=1)

    p, d, m = initialize_speedy(Float32, ShallowWater, output=true, run_id="restart-test")
    SpeedyWeather.time_stepping!(p, d, m)
 
    progn, diagn, model = initialize_speedy(Float32, ShallowWaterModel, initial_conditions=StartFromFile, restart_id="restart-test")

    for varname in propertynames(progn.layers[1].leapfrog[1])
        if SpeedyWeather.has(progn, varname)
            for (var_new, var_old) in zip(SpeedyWeather.get_var(p, varname), SpeedyWeather.get_var(progn, varname))
                @test all(var_new .== var_old)
            end
        end
    end
    @test all(SpeedyWeather.get_pressure(p) .== SpeedyWeather.get_pressure(progn))
end 