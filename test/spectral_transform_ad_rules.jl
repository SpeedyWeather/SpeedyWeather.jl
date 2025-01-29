using EnzymeTestUtils, Enzyme, FiniteDifferences
import EnzymeTestUtils: test_approx 
import FiniteDifferences: j′vp, grad, central_fdm
import AbstractFFTs

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
    @testset "Complete Transform ChainRules" begin 
        # WIP
    end
end 

@testset "Complete Differentiability" begin 
    # We do extensive correctness checks and tests on the differentiability 
    # in a seperate test set. But we do want to ensure in the regular CI that 
    # we don't commit some kind of problem for the Enzyme differentiability
    # so, we test here if we get a non-zero gradient from the timestepping.  
    spectral_grid = SpectralGrid(trunc=15, nlayers=3)          # define resolution
    model = PrimitiveWetModel(; spectral_grid)   # construct model
    simulation = initialize!(model)  
    initialize!(simulation)

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

    autodiff(Reverse, timestep!, Const, Duplicated(progn, d_progn), Duplicated(diagn, d_diag), Const(dt), Duplicated(model, d_model))
    @test sum(to_vec(d_progn)[1]) != 0
end 