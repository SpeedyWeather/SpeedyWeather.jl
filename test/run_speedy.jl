@testset "run_speedy no errors" begin
    run_speedy(Float32)
    run_speedy(Float32,output=true)
    
    run_speedy(Float64)
    run_speedy(Float64,output=true)
end