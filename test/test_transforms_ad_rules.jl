using SpeedyWeather
using EnzymeTestUtils, Enzyme
using FiniteDifferences
import FiniteDifferences: j′vp, central_fdm

grid_types = [FullGaussianGrid, OctahedralClenshawGrid] # one full and one reduced grid 

@testset "SpeedyTransforms: AD Rules" begin
    @testset "_fourier! Enzyme rules" begin
        
        @testset "reverse rule" begin
            for grid_type in grid_types,
                Tf_n in (Const,),
                Tf_s in (Const,),
                Tf_grids in (Duplicated,),
                Tf_S in (Const,),
                fun in (SpeedyWeather.SpeedyTransforms._fourier!, )

                spectral_grid = SpectralGrid(Grid=grid_type)
                S = SpectralTransform(spectral_grid)
                grid = rand(spectral_grid.Grid{spectral_grid.NF}, spectral_grid.nlat_half, spectral_grid.nlayers)
                f_north = S.scratch_memory_north
                f_south = S.scratch_memory_south

                test_reverse(fun, Const, (f_north, Tf_n), (f_south, Tf_s), (grid, Tf_grids), (S, Tf_S))
            end
        end
    end 

    @testset "Complete Transform Enzyme" begin
        # make a high level finite difference test of the whole transform
        # can't use Enzyme or ChainRule Test tools for tests for that
        for grid_type in grids 

            spectral_grid = SpectralGrid(Grid=grid_type)
            S = SpectralTransform(spectral_grid)
            grid = rand(spectral_grid.Grid{spectral_grid.NF}, spectral_grid.nlat_half, spectral_grid.nlayers)
            dgrid = zero(grid)
            specs = rand(LowerTriangularArray{Complex{spectral_grid.NF}}, spectral_grid.trunc+2, spectral_grid.trunc+1, spectral_grid.nlayers)

            autodiff(Reverse, transform!, Const, Const(specs), Duplicated(grid, dgrid), Const(S))

            # finite difference comparision

            fd_jvpi = j′vp(central_fdm(5,1), x -> transform(x, S), specs, grid)
 
        end 
    end 

    @testset "Complete Transform ChainRules" begin 
        # WIP
    end 

end 