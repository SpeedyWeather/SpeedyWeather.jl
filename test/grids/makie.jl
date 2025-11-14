using SpeedyWeather.RingGrids
using CairoMakie, GeoMakie
using Test

RingGrids.globe(FullGaussianGrid, 12, interactive = true)
plt = RingGrids.globe(randn(FullGaussianGrid(12)))

const grid_types = (
    FullClenshawGrid,
    FullGaussianGrid,
    OctahedralGaussianGrid,
    OctahedralClenshawGrid,
    OctaminimalGaussianGrid,
    HEALPixGrid,
    OctaHEALPixGrid,
    FullHEALPixGrid,
    FullOctaHEALPixGrid,
)


# Baisc tests for grid/field plotting functions to check for errors
@testset "globe (Grid)" begin
    for G in grid_types
        grid = G(12)
        @test isa(RingGrids.globe(grid), Figure)
    end
end

@testset "globe (Field)" begin
    for G in grid_types
        grid = G(12)
        @test isa(RingGrids.globe(randn(grid)), Figure)
    end
end

@testset "heatmap (Field)" begin
    for G in grid_types
        grid = G(12)
        @test isa(heatmap(randn(grid)), Figure)
    end
end
