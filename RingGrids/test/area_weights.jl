@testset "Area weights" begin
    for Grid in (
            FullGaussianGrid, FullClenshawGrid, OctahedralGaussianGrid, OctahedralClenshawGrid,
            HEALPixGrid, OctaHEALPixGrid,
        )
        nlat_half = 8
        grid = Grid(nlat_half)
        weights = get_area_weights(grid)
        @test weights == get_area_weights(Grid, nlat_half)
        @test length(weights) == get_npoints(grid)
        @test sum(weights) ≈ 1
        @test all(>(0), weights)

        # global mean of a constant is the constant, of an antisymmetric field zero
        field = ones(grid)
        @test sum(weights .* field.data) ≈ 1
        lat = get_latd(grid)
        field.data .= [lat[j] for j in RingGrids.whichring(grid)]
        @test abs(sum(weights .* field.data)) < 1.0e-10

        # equal area grids have equal weights
        if Grid <: Union{HEALPixGrid, OctaHEALPixGrid}
            @test all(weights .≈ 1 / get_npoints(grid))
        end
    end
end
