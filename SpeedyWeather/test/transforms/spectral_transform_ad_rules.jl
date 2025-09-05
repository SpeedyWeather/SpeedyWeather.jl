using EnzymeTestUtils, Enzyme
import EnzymeTestUtils: test_approx 
import AbstractFFTs
using FiniteDifferences

grid_types = [FullGaussianGrid, OctahedralGaussianGrid] # one full and one reduced grid, both Gaussian to have exact transforms 
grid_dealiasing = [2.0, 3.0]
fd_tests = [true, true] 

# currently there's an issue with EnzymeTestUtils not being able to work with structs with undefined fields like FFT plans
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
                if !(grid_type <: AbstractReducedGrid) & fd_tests[i_grid]
                    spectral_grid = SpectralGrid(Grid=grid_type, nlayers=1, trunc=5, dealiasing=grid_dealiasing[i_grid])
                    S = SpectralTransform(spectral_grid)
                    field = rand(spectral_grid.NF, spectral_grid.grid, spectral_grid.nlayers)
                    f_north = S.scratch_memory.north
                    f_south = S.scratch_memory.south

                    # forward transform 
                    test_reverse(SpeedyWeather.SpeedyTransforms._fourier!, Const, (f_north, Duplicated), (f_south, Duplicated), (field, Duplicated), (S, Const); fdm=FiniteDifferences.central_fdm(5, 1), rtol=1e-2, atol=1e-2)

                    # inverse transform
                    field = zero(field)
                    test_reverse(SpeedyWeather.SpeedyTransforms._fourier!, Const, (field, Duplicated), (f_north, Duplicated), (f_south, Duplicated), (S, Const); fdm=FiniteDifferences.central_fdm(5, 1), rtol=1e-2, atol=1e-2)
                end 
            end
        end
    end 
    @testset "Complete Transform ChainRules" begin 
        # WIP
    end
end 
