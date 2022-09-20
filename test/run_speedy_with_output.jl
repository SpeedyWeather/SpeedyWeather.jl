@testset "Output on various grids" begin
    p = run_speedy(Float64,output=true)
    @test all(isfinite.(p.layers[1].leapfrog[1].vor))

    p = run_speedy(Float32,output=true)
    @test all(isfinite.(p.layers[1].leapfrog[1].vor))

    p = run_speedy(Float64,output=true, output_dense=true)
    @test all(isfinite.(p.layers[1].leapfrog[1].vor))

    p = run_speedy(Float64,Grid=FullClenshawGrid,output=true)
    @test all(isfinite.(p.layers[1].leapfrog[1].vor))

    p = run_speedy(Float64,Grid=OctahedralGaussianGrid,output=true)
    @test all(isfinite.(p.layers[1].leapfrog[1].vor))

    p = run_speedy(Float64,Grid=OctahedralClenshawGrid,output=true)
    @test all(isfinite.(p.layers[1].leapfrog[1].vor))

    p = run_speedy(Float64,Grid=OctahedralClenshawGrid,output_grid=:matrix,output=true)
    @test all(isfinite.(p.layers[1].leapfrog[1].vor))

    p = run_speedy(Float64,Grid=OctahedralClenshawGrid,output_grid=:matrix,output_NF=Float32,output=true)
    @test all(isfinite.(p.layers[1].leapfrog[1].vor))
end