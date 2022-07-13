@testset "run_speedy no errors, no blowup" begin
    p = run_speedy(Float32)
    @test all(isfinite.(p.vor))
    
    p = run_speedy(Float32,output=true)
    @test all(isfinite.(p.vor))
    
    p = run_speedy(Float64)
    @test all(isfinite.(p.vor))

    p = run_speedy(Float64,output=true)
    @test all(isfinite.(p.vor))
end