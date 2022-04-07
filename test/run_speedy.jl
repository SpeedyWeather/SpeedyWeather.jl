@testset "Call run_speedy" begin
    run_speedy(Float64)     # just check that no error is triggered
    run_speedy(Float32)
end