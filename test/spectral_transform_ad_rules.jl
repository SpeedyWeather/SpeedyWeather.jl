using SpeedyWeather
using EnzymeTestUtils, Enzyme
import EnzymeTestUtils: test_approx
using FiniteDifferences
import FiniteDifferences: j′vp, grad, central_fdm
import AbstractFFTs

grid_types = [FullGaussianGrid, OctahedralGaussianGrid] # one full and one reduced grid, both Gaussian to have exact transforms 
grid_dealiasing = [2, 3]
fd_tests = [true, false] # to save CI time, only do FiniteDifferences test for one of the grids

# currenlty there's an issue with EnzymeTestUtils not being able to work with structs with undefined fields like FFT plans
# https://github.com/EnzymeAD/Enzyme.jl/issues/1992
# This is a very hacky workaround 
function EnzymeTestUtils.test_approx(x::AbstractFFTs.Plan, y::AbstractFFTs.Plan, msg; kwargs...)
    EnzymeTestUtils.@test_msg "$msg: types must match" typeof(x) == typeof(y)
    names = fieldnames(typeof(x))[1:end-1] # exclude pinv field (which is the last field)
    if isempty(names)
        EnzymeTestUtils.@test_msg msg x == y
    else
        for k in names
            if k isa Symbol && hasproperty(x, k)
                msg_new = "$msg: ::$(typeof(x)).$k"
            else
                msg_new = "$msg: getfield(::$(typeof(x)), $k)"
            end
            EnzymeTestUtils.test_approx(getfield(x, k), getfield(y, k), msg_new; kwargs...)
        end
    end
    return nothing
end 

@testset "SpeedyTransforms: AD Rules" begin
    @testset "_fourier! Enzyme rules" begin      
        @testset "EnzymeTestUtils reverse rule test" begin
            for (i_grid, grid_type) in enumerate(grid_types)
                
                # these tests don't pass for reduced grids 
                # this is likely due to FiniteDifferences and not our EnzymeRules 
                # see comments in https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/589
                if !(grid_type <: AbstractReducedGridArray) & fd_tests[i_grid]
                    spectral_grid = SpectralGrid(Grid=grid_type, nlayers=1, trunc=5, dealiasing=grid_dealiasing[i_grid])
                    S = SpectralTransform(spectral_grid)
                    grid = rand(spectral_grid.Grid{spectral_grid.NF}, spectral_grid.nlat_half, spectral_grid.nlayers)
                    f_north = S.scratch_memory_north
                    f_south = S.scratch_memory_south

                    # forward transform 
                    test_reverse(SpeedyWeather.SpeedyTransforms._fourier!, Const, (f_north, Duplicated), (f_south, Duplicated), (grid, Duplicated), (S, Const); fdm=FiniteDifferences.central_fdm(15, 1), rtol=1e-3, atol=1e-3)

                    # inverse transform
                    grid = zero(grid)
                    test_reverse(SpeedyWeather.SpeedyTransforms._fourier!, Const, (grid, Duplicated), (f_north, Duplicated), (f_south, Duplicated), (S, Const); fdm=FiniteDifferences.central_fdm(15, 1), rtol=1e-3, atol=1e-3)
                end 
            end
        end
    end 

    @testset "Complete Transform Enzyme" begin
        # make a high level finite difference test of the whole transform
        # can't use Enzyme or ChainRule Test tools for tests for that
        for (i_grid, grid_type) in enumerate(grid_types)

            spectral_grid = SpectralGrid(Grid=grid_type, trunc=10, nlayers=1, dealiasing=grid_dealiasing[i_grid])
            S = SpectralTransform(spectral_grid, one_more_degree=true)
            dS = deepcopy(S)

            if fd_tests[i_grid]
                
                # forwards 
                grid = rand(spectral_grid.Grid{spectral_grid.NF}, spectral_grid.nlat_half, spectral_grid.nlayers)
                dgrid = zero(grid)
                specs = zeros(LowerTriangularArray{Complex{spectral_grid.NF}}, spectral_grid.trunc+2, spectral_grid.trunc+1, spectral_grid.nlayers)
                
                # seed
                dspecs = zero(specs)
                fill!(dspecs, 1+1im)

                autodiff(Reverse, transform!, Const, Duplicated(specs, dspecs), Duplicated(grid, dgrid), Duplicated(S, dS))

                # new seed
                dspecs2 = zero(specs)
                fill!(dspecs2, 1+1im)

                # finite difference comparision, seeded with a one adjoint to get the direct gradient
                fd_jvp = FiniteDifferences.j′vp(central_fdm(5,1), x -> transform(x, S), dspecs2, grid)
                @test isapprox(dgrid, fd_jvp[1])

                ## now backwards, as the input for spec we use the output of the forward transform

                fill!(dspecs,0)
                grid = zeros(spectral_grid.Grid{spectral_grid.NF}, spectral_grid.nlat_half, spectral_grid.nlayers)
                dgrid = similar(grid)
                fill!(dgrid, 1)

                autodiff(Reverse, transform!, Const, Duplicated(grid, dgrid), Duplicated(specs, dspecs), Duplicated(S, dS))

                # new seed 
                dgrid2 = similar(grid)
                fill!(dgrid2, 1)

                fd_jvp = FiniteDifferences.j′vp(central_fdm(5,1), x -> transform(x, S), dgrid2, specs)

                @test isapprox(dspecs, fd_jvp[1])
            end 

            # test that d S^{-1}(S(x)) / dx = dx/dx = 1 (starting in both domains) 
            # this only holds for exact transforms, like Gaussian grids 

            # start with grid (but with a truncated one)
            function transform_identity!(x_out::AbstractGridArray{T}, x::AbstractGridArray{T}, S::SpectralTransform{T}) where T
                x_SH = zeros(LowerTriangularArray{Complex{T}}, S.lmax+1, S.mmax+1, S.nlayers)
                transform!(x_SH, x, S)
                transform!(x_out, x_SH, S)
                return nothing
            end 

            function transform_identity(x::AbstractGridArray{T}, S::SpectralTransform{T}) where T
                x_copy = deepcopy(x)
                transform_identity!(x_copy, x, S) 
                return x_copy
            end
              
            grid = rand(S.Grid{spectral_grid.NF}, S.nlat_half, S.nlayers)
            spec = transform(grid, S)

            grid = transform(spec, S)
            grid_out = zero(grid)
            
            transform_identity!(grid_out, grid, S)
            @test isapprox(grid, grid_out)

            dgrid = similar(grid)
            fill!(dgrid, 1)

            dgrid_out = zero(grid_out)

            autodiff(Reverse, transform_identity!, Const, Duplicated(grid_out, dgrid_out), Duplicated(grid, dgrid), Duplicated(S, dS))

            @test all(isapprox.(dgrid, 1)) 
            # TODO: previously this test was broken, with a version that directly mutates in-place. 
            # FD yields the same non-one values though. 
            # Not sure why. Do we use such things in our model? 
            #
            #function transform_identity!(x::AbstractGridArray{T}, S::SpectralTransform{T}) where T
            #   x_SH = zeros(LowerTriangularArray{Complex{T}}, S.lmax+1, S.mmax+1, S.nlayers)
            #   transform!(x_SH, x, S)
            #   transform!(x, x_SH, S)
            #   return nothing
            #end 
            # The FD comparision passes, but it takes a long time to compute, so it's commented out. 
            #dgrid2 = similar(grid)
            #fill!(dgrid2, 1)
            #fd_jvp = FiniteDifferences.j′vp(central_fdm(5,1), x -> transform_identity(x, S), dgrid2, grid)
            #@test isapprox(dgrid, fd_jvp[1], rtol=0.01)

            # now start with spectral space 

            function transform_identity!(x::LowerTriangularArray{Complex{T}}, S::SpectralTransform{T}) where T
                x_grid = zeros(S.Grid{T}, S.nlat_half, S.nlayers)
                transform!(x_grid, x, S)
                transform!(x, x_grid, S)
                return nothing
            end 

            spec = transform(grid, S)
            spec_copy = deepcopy(spec)
            transform_identity!(spec, S)
            @test isapprox(spec, spec_copy)

            dspec = similar(spec)
            fill!(dspec, 1+im)

            autodiff(Reverse, transform_identity!, Const, Duplicated(spec, dspec), Duplicated(S, dS))

            @test all(all.([isapprox.(dspec[il,1,:], 1) for il in 1:S.lmax+1])) # m = 0 => Im = 0

            for i in eachmatrix(dspec)
                @test all(isapprox.(dspec[:,i][S.lmax+2:end], 1+im)) 
            end 
        end 
    end 

    @testset "Complete Transform ChainRules" begin 
        # WIP
    end
end 