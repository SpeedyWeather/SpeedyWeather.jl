@testset "run_speedy no errors, no blowup" begin
    p = run_speedy(Float32)
    @test all(isfinite.(p.layers[1].timesteps[1].vor))
    
    p = run_speedy(Float64)
    @test all(isfinite.(p.layers[1].timesteps[1].vor))
    
    p = run_speedy(Float32,Grid=OctahedralClenshawGrid)
    @test all(isfinite.(p.layers[1].timesteps[1].vor))
    
    p = run_speedy(Float64,Grid=HEALPixGrid)
    @test all(isfinite.(p.layers[1].timesteps[1].vor))
    
    p = run_speedy(Float64,Grid=OctahedralGaussianGrid)
    @test all(isfinite.(p.layers[1].timesteps[1].vor))

    p,d,m = initialize_speedy(Float32)
    run_speedy!(p,d,m)
    @test all(isfinite.(p.layers[1].timesteps[1].vor)) 
end