using SpeedyWeather
using EnzymeTestUtils, Enzyme
using FiniteDifferences
import FiniteDifferences: j′vp, grad, central_fdm

grid_types = [FullGaussianGrid, OctahedralClenshawGrid] # one full and one reduced grid 
grid_type = grid_types[2]
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

                function speedy_fourier!(f, grid, S)
                    nothing
                end 

                #dgrid = zero(grid)

                #autodiff(Reverse, _fourier!, Const, Const(f_north), Const(f_south), Duplicated(grid, dgrid) )
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
            specs = zeros(LowerTriangularArray{Complex{spectral_grid.NF}}, spectral_grid.trunc+2, spectral_grid.trunc+1, spectral_grid.nlayers)
            dspecs = zero(specs)
            fill!(dspecs, 1+1im)
            dS = deepcopy(S)

            autodiff(Reverse, transform!, Const, Duplicated(specs, dspecs), Duplicated(grid, dgrid), Duplicated(S, dS))

            dspecs2 = zero(specs)
            fill!(dspecs2, 1+1im)

            # finite difference comparision, seeded with a one adjoint to get the direct gradient
            fd_grad = FiniteDifferences.j′vp(central_fdm(5,1), x -> transform(x, S), dspecs2, grid)
            @test isapprox(dgrid, fd_jvp[1])
        end 
    end 

    @testset "Complete Transform ChainRules" begin 
        # WIP
    end 

end 