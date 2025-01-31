### Experiments going a bit deeper into the timestepping of the primtive wet model
# this script / tests was mainly written for debugging, we might exclude it in future 
# tests because it is quite maintanance heavy code 
@testset "Differentiability: Primitive Wet Model Components" begin 

    spectral_grid = SpectralGrid(trunc=8, nlayers=1)          # define resolution

    # FiniteDifferences struggles with the NaN when we have a land-sea mask, so we have to test on AquaPlanets 
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

    autodiff(Reverse, SpeedyWeather.parameterization_tendencies!, Const, Duplicated(diagn, ddiagn), Duplicated(progn, dprogn), Const(progn.clock.time), Duplicated(model, make_zero(model)))

    function parameterization_tendencies(diagn, progn, time, model)
        diagn_new = deepcopy(diagn)
        SpeedyWeather.parameterization_tendencies!(diagn_new, progn, time, model)
        return diagn_new
    end 

    fd_vjp = FiniteDifferences.j′vp(central_fdm(5,1), x -> parameterization_tendencies(diagn_copy, x, progn.clock.time, model), ddiagn_copy, progn_copy)
    
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

    autodiff(Reverse, SpeedyWeather.ocean_timestep!, Const, Duplicated(progn, dprogn), Duplicated(diagn, ddiagn), Duplicated(model, make_zero(model)))

    function ocean_timestep(progn, diagn, model)
        progn_new = deepcopy(progn)
        SpeedyWeather.ocean_timestep!(progn_new, diagn, model)
        return progn_new
    end 

    fd_vjp = FiniteDifferences.j′vp(central_fdm(5,1), x -> ocean_timestep(progn_copy, x, model), dprogn_copy, diagn_copy)

    # also broken? 
    @test all(isapprox.(to_vec(fd_vjp[1])[1], to_vec(ddiagn)[1],rtol=1e-4,atol=1e-1))

    #
    # land 
    # 

    progn_copy = deepcopy(progn)
    dprogn = one(progn)
    ddiagn = make_zero(diagn)
    dprogn_copy = deepcopy(dprogn)
    diagn_copy = deepcopy(diagn)

    autodiff(Reverse, SpeedyWeather.land_timestep!, Const, Duplicated(progn, dprogn), Duplicated(diagn, ddiagn), Duplicated(model, make_zero(model)))

    function land_timestep(progn, diagn, model)
        progn_new = deepcopy(progn)
        SpeedyWeather.ocean_timestep!(progn_new, diagn, model)
        return progn_new
    end 

    fd_vjp = FiniteDifferences.j′vp(central_fdm(5,1), x -> land_timestep(progn_copy, x, model), dprogn_copy, diagn_copy)

    # also broken currently
    @test all(isapprox.(to_vec(fd_vjp[1])[1], to_vec(ddiagn)[1],rtol=1e-4,atol=1e-1))

    #####
    # DYNAMICS 
    lf2 = 2 
    #
    # drag! 
    #

    diagn_copy = deepcopy(diagn)
    ddiag = one(diagn_copy)
    ddiag_copy = deepcopy(ddiag)
    progn_copy = deepcopy(progn)
    dprogn = make_zero(progn)

    autodiff(Reverse, SpeedyWeather.drag!, Const, Duplicated(diagn, ddiag), Duplicated(progn, dprogn), Const(lf2), Duplicated(model, make_zero(model)))

    function drag(diagn, progn, lf, model)
        diagn_new = deepcopy(diagn)
        SpeedyWeather.drag!(diagn_new, progn, lf, model)
        return diagn_new
    end 

    fd_vjp = FiniteDifferences.j′vp(central_fdm(9,1), x -> drag(diagn_copy, x, lf2, model), ddiag_copy, progn_copy)

    @test all(isapprox.(to_vec(fd_vjp[1])[1], to_vec(dprogn)[1],rtol=1e-4,atol=1e-1))

    # in the default configuration without forcing or drag, the barotropic model's don't dependent on the previous prognostic state 
    @test sum(to_vec(dprogn)[1]) ≈ 0 

    # 
    # dynamics_tendencies!
    #

    diagn_copy = deepcopy(diagn)
    ddiag = one(diagn_copy)
    ddiag_copy = deepcopy(ddiag)
    progn_copy = deepcopy(progn)
    dprogn = make_zero(progn)

    autodiff(Reverse, SpeedyWeather.dynamics_tendencies!, Const, Duplicated(diagn, ddiag), Duplicated(progn, dprogn), Const(lf2), Duplicated(model, make_zero(model)))

    function dynamics_tendencies(diagn, progn, lf, model)
        diagn_new = deepcopy(diagn)
        SpeedyWeather.dynamics_tendencies!(diagn_new, progn, lf, model)
        return diagn_new
    end 

    fd_vjp = FiniteDifferences.j′vp(central_fdm(9,1), x -> dynamics_tendencies(diagn_copy, x, lf2, model), ddiag_copy, progn_copy)

    @test all(isapprox.(to_vec(fd_vjp)[1], to_vec(dprogn)[1],rtol=1e-4,atol=1e-1))

    # in the default configuration without forcing or drag, the barotropic model's don't dependent on the previous prognostic state 
    @test sum(to_vec(dprogn)[1]) ≈ 0 
    
    #
    # Implicit correction 
    #

    diagn_copy = deepcopy(diagn)
    ddiag = make_zero(diagn_copy)
    ddiag_copy = deepcopy(ddiag)
    progn_copy = deepcopy(progn)
    dprogn = one(progn)

    autodiff(Reverse, SpeedyWeather.implicit_correction!, Const, Duplicated(diagn, ddiag), Duplicated(model.implicit, make_zero(model.implicit)), Duplicated(progn, dprogn))

    function implicit_correction(diagn, implicit, progn)
        diagn_new = deepcopy(diagn)
        SpeedyWeather.implicit_correction(diagn_new, implicit, progn)
        return diagn_new
    end 

    fd_vjp = FiniteDifferences.j′vp(central_fdm(9,1), x -> implicit_correction(diagn_copy, model.implicit, x), ddiag_copy, progn_copy)

    @test all(isapprox.(to_vec(fd_vjp[1])[1], to_vec(dprogn)[1],rtol=1e-4,atol=1e-1))

    #
    # horizontal_diffusion!
    #

    lf1 = 1
    diagn_copy = deepcopy(diagn)
    ddiag = one(diagn_copy)
    ddiag_copy = deepcopy(ddiag)

    progn_copy = deepcopy(progn)
    dprogn = make_zero(progn)

    autodiff(Reverse, SpeedyWeather.horizontal_diffusion!, Const, Duplicated(diagn, ddiag), Duplicated(progn, dprogn), Const(model.horizontal_diffusion), Duplicated(model, make_zero(model)), Const(lf1))

    # FD comparision not necessary, we have the exact values 
    #function horizontal_diffusion(diagn, progn, diffusion, model, lf)
    #    diagn_new = deepcopy(diagn)
    #    SpeedyWeather.horizontal_diffusion!(diagn_new, progn, diffusion, model, lf)
    #    return diagn_new
    #end 

    #fd_vjp = FiniteDifferences.j′vp(central_fdm(5,1), x -> horizontal_diffusion(diagn_copy, x, model.horizontal_diffusion, model, lf1), ddiag_copy, progn_copy)
    #@test all(isapprox.(to_vec(fd_vjp[1])[1], to_vec(dprogn)[1],rtol=1e-4,atol=1e-2))

    # ∂(progn)
    # should be row-wise `model.horizontal_diffusion.impl .* model.horizontal_diffusion.expl`
    # for all variables that are diffused 
    diff_coefficient = model.horizontal_diffusion.impl .* model.horizontal_diffusion.expl
    l_indices = [(1:l) for l=1:progn.vor[1].n]
    for (i,il) in enumerate(l_indices)
        @test all(real.(Matrix(dprogn.vor[lf1][:,1])[i, il]) .≈ diff_coefficient[i])
    end 

    # ∂(tend_old)
    # should be row-wise `model.horizontal_diffusion.impl` 
    for (i,il) in enumerate(l_indices)
        @test all(real.(Matrix(ddiag.tendencies.vor_tend[:,1])[i, il]) .≈ model.horizontal_diffusion.impl[i])
    end 

    #
    # Test the leapfrog 
    # 

    lf1 = 2 
    lf2 = 2 

    progn_copy = deepcopy(progn)
    dprogn = one(progn_copy)
    dprogn_copy = one(progn_copy)

    tend = diagn.tendencies
    tend_copy = deepcopy(tend)
    dtend = make_zero(tend)
    dmodel = make_zero(model)

    autodiff(Reverse, SpeedyWeather.leapfrog!, Const, Duplicated(progn, dprogn), Duplicated(tend, dtend), Const(dt), Const(lf1), Const(model))

    function leapfrog_step(progn_new::PrognosticVariables, progn::PrognosticVariables, tend, dt, lf, model)
        copy!(progn_new, progn)
        SpeedyWeather.leapfrog!(progn_new, tend, dt, lf, model)
        return progn_new
    end 

    prog_new = zero(progn_copy)

    fd_vjp = FiniteDifferences.j′vp(central_fdm(5,1), x -> leapfrog_step(prog_new, progn_copy, x, dt, lf1, model), dprogn_copy, tend_copy)

    @test all(isapprox.(to_vec(fd_vjp[1])[1], to_vec(dtend)[1],rtol=1e-5,atol=1e-5))

    #
    # transform!(diagn, progn, lf2, model)
    #

    fill!(diagn.tendencies, 0, PrimitiveWetModel)
    diag_copy = deepcopy(diagn)
    
    ddiag = one(diagn)
    ddiag_copy = deepcopy(ddiag)

    progn_copy = deepcopy(progn)
    dprogn = make_zero(progn)

    autodiff(Reverse, SpeedyWeather.transform!, Const, Duplicated(diagn, ddiag), Duplicated(progn, dprogn), Const(lf2), Duplicated(model, make_zero(model)))

    function transform_diagn(diag, progn, lf2, model)
        diag_copy = deepcopy(diag)
        transform!(diag_copy, progn, lf2, model)
        return diag_copy
    end 

    fd_vjp = FiniteDifferences.j′vp(central_fdm(5,1), x -> transform_diagn(diag_copy, x, lf2, model), ddiag_copy, progn_copy)
    
    @test all(isapprox.(to_vec(fd_vjp[1])[1], to_vec(dprogn)[1],rtol=1e-3,atol=1e-3))
end


