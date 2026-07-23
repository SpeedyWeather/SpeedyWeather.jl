@testset "Rotate grids" begin
    @testset for Grid in (
            FullGaussianGrid,
            FullClenshawGrid,
            FullHEALPixGrid,
            FullOctaHEALPixGrid,
            OctahedralGaussianGrid,
            OctahedralClenshawGrid,
            HEALPixGrid,
            OctaHEALPixGrid,
            OctaminimalGaussianGrid,
        )

        @testset for nlat_half in (4, 6, 8)
            @testset for k in 0:3
                if k == 0
                    grid = rand(Grid, nlat_half)
                else
                    grid = rand(Grid, nlat_half, k)
                end

                grid2 = deepcopy(grid)
                @test rotate!(grid2, 0) == grid
                @test rotate!(rotate!(rotate!(rotate!(grid2, 90), 90), 90), 90) == grid
                @test rotate!(rotate!(grid2, 180), 180) == grid
                @test rotate!(rotate!(grid2, 270), 90) == grid

                @test rotate!(rotate!(grid2, 90), -90) == grid
                @test rotate!(rotate!(grid2, -90), 90) == grid
            end
        end
    end
end

@testset "Rotate n-dimensional fields" begin
    @testset for Grid in (FullGaussianGrid, OctahedralGaussianGrid, HEALPixGrid)
        grid = Grid(4)
        @testset for trailing_dims in ((2,), (2, 3), (2, 3, 2))
            field = rand(Float32, grid, trailing_dims...)
            @testset for degree in (90, 180, 270)
                rotated = rotate(field, degree)

                # rotating all layers at once == rotating every 2D layer independently
                for c in CartesianIndices(trailing_dims)
                    layer = Field(field.data[:, c], grid)
                    @test rotate(layer, degree) == Field(rotated.data[:, c], grid)
                end
            end
        end
    end
end
