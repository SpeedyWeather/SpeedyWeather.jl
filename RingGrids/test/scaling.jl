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

@testset "Scale coslat n-dimensional fields" begin
    @testset for Grid in (FullGaussianGrid, OctahedralGaussianGrid, HEALPixGrid)
        grid = Grid(4)
        @testset for trailing_dims in ((2,), (2, 3), (2, 3, 2))
            field = randn(Float64, grid, trailing_dims...)
            scaled = RingGrids.scale_coslat(field)

            # scaling all layers at once == scaling every 2D layer independently
            for c in CartesianIndices(trailing_dims)
                layer = Field(field.data[:, c], grid)
                @test RingGrids.scale_coslat(layer) == Field(scaled.data[:, c], grid)
            end
        end
    end
end
