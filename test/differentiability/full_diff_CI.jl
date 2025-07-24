
# Enzyme and Julia 1.11 still has some problems, and the test below is broken
# in Julia 1.11
using Enzyme, FiniteDifferences

if VERSION <= v"1.11.0"
    @testset "Complete Differentiability" begin 
        # We do extensive correctness checks and tests on the differentiability 
        # in a seperate test set. But we do want to ensure in the regular CI that 
        # we don't commit some kind of problem for the Enzyme differentiability
        # so, we test here if we get a non-zero gradient from the timestepping.  
        spectral_grid = SpectralGrid(trunc=5, nlayers=1)          # define resolution
        model = PrimitiveWetModel(; spectral_grid)   # construct model
        simulation = initialize!(model)  
        initialize!(simulation)
        run!(simulation, period=Hour(6))
        
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
        # currently not activated to keep the CI fast 
        #autodiff(set_runtime_activity(Reverse), timestep_oop!, Const, Duplicated(progn_new, dprogn_new), Duplicated(progn, d_progn), Duplicated(diagn, d_diag), Const(dt), Const(model))
        #@test sum(to_vec(d_progn)[1]) != 0
        #@test progn != d_progn
    end 
else 
    @testset "Complete Differentiability" begin
        @test_broken false # we report a broken test here on v1.11, just to indicate that this (properly) doesn't work yet
    end 
end 