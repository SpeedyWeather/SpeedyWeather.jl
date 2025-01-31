# High-level tests whether time stepping in all models work

@testset "Differentiability: Timestepping" begin 
    # T15 still yields somewhat sensible dynamics, that's why it's chosen here
    model_types = [ShallowWaterModel, PrimitiveDryModel, PrimitiveWetModel]

    for model_type in model_types 
        
        nlayer = model_type == ShallowWaterModel ? 1 : 1

        # FiniteDifferences struggles with the NaN when we have a land-sea mask, so we have to test on AquaPlanets 
        spectral_grid = SpectralGrid(trunc=8, nlayers=nlayer)          # define resolution
        model = model_type(; spectral_grid)   # construct model
        simulation = initialize!(model)  
        initialize!(simulation)
        run!(simulation, period=Day(3))

        (; prognostic_variables, diagnostic_variables, model) = simulation
        (; Δt, Δt_millisec) = model.time_stepping
        dt = 2Δt

        progn = prognostic_variables
        diagn = diagnostic_variables

        # TO-DO: The first time we execute this, the gradient is different. Why?
        timestep_oop!(make_zero(progn), progn, diagn, dt, model)

        diagn_copy = deepcopy(diagn)
        progn_copy = deepcopy(progn)

        d_progn = zero(progn)
        d_diag = make_zero(diagn)

        progn_new = zero(progn)
        dprogn_new = one(progn) # seed 

        # test that we can differentiate wrt an IC 
        autodiff(Reverse, timestep_oop!, Const, Duplicated(progn_new, dprogn_new), Duplicated(progn, d_progn), Duplicated(diagn, d_diag), Const(dt), Duplicated(model, make_zero(model)))

        # nonzero gradient
        @test sum(to_vec(d_progn)[1]) != 0

        # FD comparison 
        dprogn_2 = one(progn) # seed 

        # this takes a long time 
        # with the FD comparision we have to go to quite low tolerences for the full time step 
        fd_vjp = FiniteDifferences.j′vp(central_fdm(5,1), x -> timestep_oop(x, diagn_copy, dt, model), dprogn_2, progn_copy)

        @test isapprox(to_vec(fd_vjp[1])[1], to_vec(d_progn)[1])
    end 
end 