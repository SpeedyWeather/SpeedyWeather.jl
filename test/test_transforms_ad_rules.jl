using EnzymeTestUtils, Enzyme, FiniteDifferences

grids = [FullGaussianGrid, OctahedralClenshawGrid] # one full and one reduced grid 

@testset "SpeedyTransforms: AD Rules" begin


    @testset "_fourier! Enzyme rules" begin
        
        @testset "reverse rule" begin
            for grid in grids,
                Tf_n in (Duplicated,),
                Tf_s in (Duplicated,),
                Tf_grids in (Duplicated,),
                Tf_S in (Const,),
                fun in (SpeedyWeather.SpeedyTransforms._fourier!, )

                spectral_grid = SpectralGrid(Grid=grid)
                S = SpectralTransform(spectral_grid)
                grid = rand(spectral_grid.Grid{spectral_grid.NF}, spectral_grid.nlat_half, spectral_grid.nlayers)
                f_north = S.scratch_memory_north
                f_south = S.scratch_memory_south

                test_reverse(fun, Const, (f_north, Tf_n), (f_south, Tf_s), (grid, Tf_grids), (S, Const))
            end
        end
    end 

    @testset "_fourier! ChainRules" begin 
        # WIP
    end 


    @testset "Complete Transform" begin
        # make a high level finite difference test of the whole transform
        # can't use Enzyme or ChainRule Test tools for tests for that
        for grid in grids 

            spectral_grid = SpectralGrid(Grid=grid)
            S = SpectralTransform(spectral_grid)
            grid = rand(spectral_grid.Grid{spectral_grid.NF}, spectral_grid.nlat_half, spectral_grid.nlayers)
            dgrid = zero(grid)
            specs = rand(LowerTriangularArray{Complex{spectral_grid.NF}}, spectral_grid.trunc+2, spectral_grid.trunc+1, spectral_grid.nlayers)
            dspecs = zero(specs)

            autodiff(Reverse, transform!, Const, Duplicated(specs, dspecs), Duplicated(grid, dgrid), Const(S))

            # finite difference comparision 

    end 
end 