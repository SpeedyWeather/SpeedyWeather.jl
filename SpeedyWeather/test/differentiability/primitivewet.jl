### Experiments going a bit deeper into the timestepping of the primtive wet model
# this script / these tests were mainly written for debugging, we might exclude it in future 
# tests because it is quite maintanance heavy code 
@testset "Differentiability: Primitive Wet Model Components" begin 

    spectral_grid = SpectralGrid(trunc=8, nlayers=1)          # define resolution

    model = PrimitiveWetModel(; spectral_grid)  # construct model
    simulation = initialize!(model)  
    initialize!(simulation)
    run!(simulation, period=Day(5)) # spin-up to get nonzero values for all fields
    initialize!(simulation; period=Day(1))
    
    (; prognostic_variables, diagnostic_variables, model) = simulation
    (; Δt, Δt_millisec) = model.time_stepping
    dt = 2Δt

    progn = prognostic_variables
    diagn = diagnostic_variables

    # TO-DO: The first time we execute this, the gradient is different. Why?
    timestep_oop!(make_zero(progn), progn, diagn, dt, model)

    #
    # We go individually through all components of the time stepping and check 
    # correctness
    #

    fill!(diagn.tendencies, 0, PrimitiveWetModel)
    (; time) = progn.clock 

    #
    # model physics 
    # 

    progn_copy = deepcopy(progn)
    dprogn = one(progn)
    ddiagn = one(diagn)
    ddiagn_copy = deepcopy(ddiagn)
    diagn_copy = deepcopy(diagn)

    autodiff(Reverse, SpeedyWeather.parameterization_tendencies!, Const, Duplicated(diagn, ddiagn), Duplicated(progn, dprogn), Const(progn.clock.time), Const(model))

    function parameterization_tendencies(diagn, progn, time, model)
        diagn_new = deepcopy(diagn)
        SpeedyWeather.parameterization_tendencies!(diagn_new, deepcopy(progn), deepcopy(time), deepcopy(model))
        return diagn_new
    end 

    fd_vjp = FiniteDifferences.j′vp(central_fdm(11,1), x -> parameterization_tendencies(diagn_copy, x, progn.clock.time, model), ddiagn_copy, progn_copy)
    
    # TO-DO this test is broken, they gradients don't line up
    #@test all(isapprox.(to_vec(fd_vjp[1])[1], to_vec(dprogn)[1],rtol=1e-4,atol=1e-1))
    
    #
    # ocean 
    # 

    progn_copy = deepcopy(progn)
    dprogn = one(progn)
    ddiagn = make_zero(diagn)
    dprogn_copy = deepcopy(dprogn)
    diagn_copy = deepcopy(diagn)

    autodiff(Reverse, SpeedyWeather.ocean_timestep!, Const, Duplicated(progn, dprogn), Duplicated(diagn, ddiagn), Const(model))

    function ocean_timestep(progn, diagn, model)
        progn_new = deepcopy(progn)
        SpeedyWeather.ocean_timestep!(progn_new, deepcopy(diagn), deepcopy(model))
        return progn_new
    end 

    fd_vjp = FiniteDifferences.j′vp(central_fdm(5,1), x -> ocean_timestep(progn_copy, x, model), dprogn_copy, diagn_copy)

    # pass
    @test all(isapprox.(to_vec(fd_vjp[1])[1], to_vec(ddiagn)[1],rtol=1e-4,atol=1e-2))

    #
    # land 
    # 

    progn_copy = deepcopy(progn)
    dprogn = one(progn)
    ddiagn = make_zero(diagn)
    dprogn_copy = deepcopy(dprogn)
    diagn_copy = deepcopy(diagn)

    autodiff(Reverse, SpeedyWeather.land_timestep!, Const, Duplicated(progn, dprogn), Duplicated(diagn, ddiagn), Const(model))

    function land_timestep(progn, diagn, model)
        progn_new = deepcopy(progn)
        SpeedyWeather.ocean_timestep!(progn_new, deepcopy(diagn), deepcopy(model))
        return progn_new
    end 

    fd_vjp = FiniteDifferences.j′vp(central_fdm(5,1), x -> land_timestep(progn_copy, x, model), dprogn_copy, diagn_copy)

    # pass
    @test all(isapprox.(to_vec(fd_vjp[1])[1], to_vec(ddiagn)[1],rtol=1e-4,atol=1e-2))

    #####
    # DYNAMICS 
    lf2 = 2 
   
    # 
    # dynamics_tendencies!
    #

    diagn_copy = deepcopy(diagn)
    ddiag = one(diagn_copy)
    ddiag_copy = deepcopy(ddiag)
    progn_copy = deepcopy(progn)
    dprogn = make_zero(progn)
   
    autodiff(Reverse, SpeedyWeather.dynamics_tendencies!, Const, Duplicated(diagn, ddiag), Duplicated(progn, dprogn), Const(lf2), Const(model))

    function dynamics_tendencies(diagn, progn, lf, model)
        diagn_new = deepcopy(diagn)
        SpeedyWeather.dynamics_tendencies!(diagn_new, deepcopy(progn), lf, deepcopy(model))
        return diagn_new
    end 

    fd_vjp = FiniteDifferences.j′vp(central_fdm(5,1), x -> dynamics_tendencies(diagn_copy, x, lf2, model), ddiag_copy, progn_copy)

    # there are some NaNs in the FD, that's why this test is currently broken
    @test all(isapprox.(to_vec(fd_vjp[1])[1], to_vec(dprogn)[1],rtol=1e-4,atol=1e-1))
    
    #
    # Implicit correction 
    #
    # continue here
    diagn_copy = deepcopy(diagn)
    ddiag = make_zero(diagn_copy)
    ddiag_copy = deepcopy(ddiag)
    progn_copy = deepcopy(progn)
    dprogn = one(progn)

    autodiff(Reverse, SpeedyWeather.implicit_correction!, Const, Duplicated(diagn, ddiag), Const(model.implicit), Duplicated(progn, dprogn))

    function implicit_correction(diagn, implicit, progn)
        diagn_new = deepcopy(diagn)
        SpeedyWeather.implicit_correction!(diagn_new, deepcopy(implicit), deepcopy(progn))
        return diagn_new
    end 

    fd_vjp = FiniteDifferences.j′vp(central_fdm(9,1), x -> implicit_correction(diagn_copy, model.implicit, x), ddiag_copy, progn_copy)

    @test all(isapprox.(to_vec(fd_vjp[1])[1], to_vec(dprogn)[1],rtol=1e-4,atol=1e-1))

    #
    # transform!(diagn, progn, lf2, model)
    #

    diag_copy = deepcopy(diagn)
    
    ddiag = one(diagn)
    ddiag_copy = deepcopy(ddiag)

    progn_copy = deepcopy(progn)
    dprogn = make_zero(progn)

    autodiff(Reverse, SpeedyWeather.transform!, Const, Duplicated(diagn, ddiag), Duplicated(progn, dprogn), Const(lf2), Duplicated(model, make_zero(model)))
    #autodiff(Reverse, SpeedyWeather.transform!, Const, Duplicated(diagn, ddiag), Duplicated(progn, dprogn), Const(lf2), Const(model))

    function transform_diagn(diag, progn, lf2, model)
        diag_copy = deepcopy(diag)
        transform!(diag_copy, deepcopy(progn), lf2, deepcopy(model))
        return diag_copy
    end 

    fd_vjp = FiniteDifferences.j′vp(central_fdm(5,1), x -> transform_diagn(diag_copy, x, lf2, model), ddiag_copy, progn_copy)
    
    @test all(isapprox.(to_vec(fd_vjp[1])[1], to_vec(dprogn)[1],rtol=1e-3,atol=1e-3))
end


