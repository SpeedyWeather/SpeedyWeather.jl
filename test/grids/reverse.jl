@testset "Reverse grids" begin
    @testset for Grid in (
        FullGaussianGrid,
        FullClenshawGrid,
        FullHEALPixGrid,
        FullOctaHEALPixGrid,
        OctahedralGaussianGrid,
        OctahedralClenshawGrid,
        HEALPixGrid,
        OctaHEALPixGrid,
        OctaminimalGaussianGrid
        )

        @testset for nlat_half in (4, 6, 8)
            @testset for k in 0:3
                if k == 0
                    grid = rand(Grid, nlat_half)
                else
                    grid = rand(Grid, nlat_half, k)
                end

                grid2 = deepcopy(grid)
                reverse!(grid2)
                @test reverse!(grid2) == grid
                @test reverse(reverse(grid)) == grid

                reverse!(grid2, dims=:lat)
                @test reverse!(grid2, dims=:latitude) == grid
                @test reverse(reverse(grid, dims=:lat), dims=:lat) == grid

                reverse!(grid2, dims=:lon)
                @test reverse!(grid2, dims=:longitude) == grid
                @test reverse(reverse(grid, dims=:lon), dims=:lon) == grid

                if k < 2
                    @test reverse(reverse(grid, dims=:lat), dims=:lon) == reverse(grid)
                    @test reverse(reverse(grid, dims=:lon), dims=:lat) == reverse(grid)
                end
            end
        end
    end
end