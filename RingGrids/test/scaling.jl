@testset "Scale, unscale coslat" begin
    @testset for NF in (Float32, Float64)
        for Grid in (
                FullGaussianGrid,
                FullClenshawGrid,
                OctahedralGaussianGrid,
                OctahedralClenshawGrid,
                HEALPixGrid,
            )

            grid = Grid(24)
            A = randn(NF, grid)
            B = copy(A)

            RingGrids.scale_coslat⁻¹!(A)
            RingGrids.scale_coslat!(A)

            @test all(isapprox.(A, B, rtol = 10 * eps(NF)))

            RingGrids.scale_coslat²!(A)
            RingGrids.scale_coslat⁻²!(A)

            @test all(isapprox.(A, B, rtol = 10 * eps(NF)))
        end
    end
end
