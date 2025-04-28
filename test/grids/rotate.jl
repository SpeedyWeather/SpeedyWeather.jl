import SpeedyWeather.RingGrids: rotate!, rotate

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