@testset "run_speedy no errors, no blowup" begin
    p = run_speedy(Barotropic)
    @test all(isfinite.(p.layers[1].timesteps[1].vor))
    
    p = run_speedy(ShallowWater)
    @test all(isfinite.(p.layers[1].timesteps[1].vor))
    
    p = run_speedy(PrimitiveDryCore)
    @test all(isfinite.(p.layers[1].timesteps[1].vor))
    
    p = run_speedy(PrimitiveWetCore)
    @test all(isfinite.(p.layers[1].timesteps[1].vor))
    
    p,d,m = initialize_speedy()
    run_speedy!(p,d,m)
    @test all(isfinite.(p.layers[1].timesteps[1].vor)) 
end