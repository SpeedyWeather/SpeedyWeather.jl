using EnzymeTestUtils, Enzyme
import EnzymeTestUtils: test_approx 
import AbstractFFTs
using FiniteDifferences

grid_types = [FullGaussianGrid, OctahedralGaussianGrid] # one full and one reduced grid, both Gaussian to have exact transforms 
grid_dealiasing = [2, 3]
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
                if !(grid_type <: AbstractReducedGridArray) & fd_tests[i_grid]
                    spectral_grid = SpectralGrid(Grid=grid_type, nlayers=1, trunc=5, dealiasing=grid_dealiasing[i_grid])
                    S = SpectralTransform(spectral_grid)
                    grid = rand(spectral_grid.Grid{spectral_grid.NF}, spectral_grid.nlat_half, spectral_grid.nlayers)
                    f_north = S.scratch_memory.north
                    f_south = S.scratch_memory.south
                    # forward transform 
                    test_reverse(SpeedyWeather.SpeedyTransforms._fourier!, Const, (f_north, Duplicated), (f_south, Duplicated), (grid, Duplicated), (S, Const); fdm=FiniteDifferences.central_fdm(15, 1), rtol=1e-3, atol=1e-3)

                    # inverse transform
                    grid = zero(grid)
                    test_reverse(SpeedyWeather.SpeedyTransforms._fourier!, Const, (grid, Duplicated), (f_north, Duplicated), (f_south, Duplicated), (S, Const); fdm=FiniteDifferences.central_fdm(15, 1), rtol=1e-3, atol=1e-3)
                end 
            end
        end
    end 
    @testset "Complete Transform ChainRules" begin 
        # WIP
    end
end 

# Enzyme and Julia 1.11 still has some problems, and the test below is broken
# in Julia 1.11
if VERSION <= v"1.11.0"
    @testset "Complete Differentiability" begin 
        # We do extensive correctness checks and tests on the differentiability 
        # in a seperate test set. But we do want to ensure in the regular CI that 
        # we don't commit some kind of problem for the Enzyme differentiability
        # so, we test here if we get a non-zero gradient from the timestepping.  
        spectral_grid = SpectralGrid(trunc=8, nlayers=1)          # define resolution
        model = PrimitiveWetModel(; spectral_grid)   # construct model
        simulation = initialize!(model)  
        initialize!(simulation)
        run!(simulation, period=Day(1))
        
        (; prognostic_variables, diagnostic_variables, model) = simulation
        (; Δt, Δt_millisec) = model.time_stepping
        dt = 2Δt

        progn = prognostic_variables
        diagn = diagnostic_variables

        diagn_copy = deepcopy(diagn)
        progn_copy = deepcopy(progn)

        d_progn = zero(progn)
        d_diag = make_zero(diagn)
        d_model = make_zero(model)

        progn_new = zero(progn)
        dprogn_new = one(progn) # seed 

        function timestep_oop!(progn_new::PrognosticVariables, progn_old::PrognosticVariables, diagn, dt, model, lf1=2, lf2=2)
            copy!(progn_new, progn_old)
            SpeedyWeather.timestep!(progn_new, diagn, dt, model, lf1, lf2)
            return nothing
        end 

        autodiff(Reverse, timestep_oop!, Const, Duplicated(progn_new, dprogn_new), Duplicated(progn, d_progn), Duplicated(diagn, d_diag), Const(dt), Duplicated(model, d_model))
        @test sum(to_vec(d_progn)[1]) != 0

        # with Const(model)
        autodiff(set_runtime_activity(Reverse), timestep_oop!, Const, Duplicated(progn_new, dprogn_new), Duplicated(progn, d_progn), Duplicated(diagn, d_diag), Const(dt), Const(model))
        @test sum(to_vec(d_progn)[1]) != 0
        @test progn != d_progn
    end 
else 
    @testset "Complete Differentiability" begin
        @test_broken false # we report a broken test here on v1.11, just to indicate that this (properly) doens't work yet
    end 
end 