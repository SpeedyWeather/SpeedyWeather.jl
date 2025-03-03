### Experiments going a bit deeper into the timestepping of the barotropic model
@testset "Differentiability: Barotropic Model Components" begin 
    # T15 still yields somewhat sensible dynamics, that's why it's chosen here
    spectral_grid = SpectralGrid(trunc=15, nlayers=1)          # define resolution
    model = BarotropicModel(; spectral_grid)   # construct model
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

    diagn_copy = deepcopy(diagn) # first copy for FD comparison 
    progn_copy = deepcopy(progn)

    diagn_copy_2 = deepcopy(diagn) # 2nd copy for Const activity test 
    progn_copy_2 = deepcopy(progn)

    d_progn = zero(progn)
    d_diag = make_zero(diagn)

    progn_new = zero(progn)
    dprogn_new = one(progn) # seed 

    # test if differentiation works wrt copy! (there were some problems with it before)
    autodiff(Reverse, copy!, Const, Duplicated(progn_new, dprogn_new), Duplicated(progn, d_progn))

    progn_new = zero(progn)
    dprogn_new = one(progn) # seed 

    # test that we can differentiate wrt to everything  
    autodiff(Reverse, timestep_oop!, Const, Duplicated(progn_new, dprogn_new), Duplicated(progn, d_progn), Duplicated(diagn, d_diag), Const(dt), Duplicated(model, make_zero(model)))

    # nonzero gradient
    @test sum(to_vec(d_progn)[1]) != 0

    # FD comparison 
    dprogn_2 = one(progn) # seed 

    # for the full timestep, we need a bit higher precision 
    fd_vjp = FiniteDifferences.j′vp(central_fdm(15,1), x -> timestep_oop(x, diagn_copy, dt, model), dprogn_2, progn_copy)

    @test isapprox(to_vec(fd_vjp[1])[1], to_vec(d_progn)[1], rtol=0.05) # we have to go really quite low with the tolerances here
    @test mean(abs.(to_vec(fd_vjp[1])[1] - to_vec(d_progn)[1])) < 0.002 # so we check a few extra statistics
    @test maximum(to_vec(fd_vjp[1].vor)[1] - to_vec(d_progn.vor)[1]) < 0.05
    
    # test that we can differentiante with Const(Model) only wrt to the state
    d_progn = zero(progn)
    d_diag = make_zero(diagn)

    progn_new = zero(progn)
    dprogn_new = one(progn) # seed 

    autodiff(set_runtime_activity(Reverse), timestep_oop!, Const, Duplicated(progn_new, dprogn_new), Duplicated(progn_copy_2, d_progn), Duplicated(diagn_copy_2, d_diag), Const(dt), Const(model))

    # use the same FD comparision 

    @test isapprox(to_vec(fd_vjp[1])[1], to_vec(d_progn)[1], rtol=0.05) # we have to go really quite low with the tolerances here
    @test mean(abs.(to_vec(fd_vjp[1])[1] - to_vec(d_progn)[1])) < 0.002 # so we check a few extra statistics
    @test maximum(to_vec(fd_vjp[1].vor)[1] - to_vec(d_progn.vor)[1]) < 0.05
    
    #
    # We go individually through all components of the time stepping and check 
    # correctness
    #

    fill!(diagn.tendencies, 0, Barotropic)

    #
    # dynamics_tendencies!
    #

    lf2 = 2 

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

    fd_vjp = FiniteDifferences.j′vp(central_fdm(9,1), x -> dynamics_tendencies(diagn_copy, x, lf2, model), ddiag_copy, progn_copy)

    @test all(isapprox.(to_vec(fd_vjp[1])[1], to_vec(dprogn)[1],rtol=1e-4,atol=1e-1))

    # in the default configuration without forcing or drag, the barotropic model's don't dependent on the previous prognostic state 
    @test sum(to_vec(dprogn)[1]) ≈ 0 

    #
    # horizontal_diffusion!
    #

    lf1 = 1
    diagn_copy = deepcopy(diagn)
    ddiag = one(diagn_copy)
    ddiag_copy = deepcopy(ddiag)

    progn_copy = deepcopy(progn)
    dprogn = make_zero(progn)

    autodiff(Reverse, SpeedyWeather.horizontal_diffusion!, Const, Duplicated(diagn, ddiag), Duplicated(progn, dprogn), Const(model.horizontal_diffusion), Const(model), Const(lf1))

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

    autodiff(Reverse, SpeedyWeather.leapfrog!, Const, Duplicated(progn, dprogn), Duplicated(tend, dtend), Const(dt), Const(lf1), Const(model))

    function leapfrog_step(progn_new::PrognosticVariables, progn::PrognosticVariables, tend, dt, lf, model)
        copy!(progn_new, progn)
        SpeedyWeather.leapfrog!(progn_new, tend, dt, lf, model)
        return progn_new
    end 

    prog_new = zero(progn_copy)

    fd_vjp = FiniteDifferences.j′vp(central_fdm(5,1), x -> leapfrog_step(prog_new, progn_copy, x, dt, lf1, model), dprogn_copy, tend_copy)

    @test all(isapprox.(to_vec(fd_vjp[1])[1], to_vec(dtend)[1],rtol=1e-4,atol=1e-4))

    # 
    # single variable leapfrog step 
    # 

    A_old = progn.vor[1]
    A_old_copy = copy(A_old)
    dA_old = one(A_old)

    A_new = progn.vor[2]
    A_new_copy = copy(A_new)
    dA_new = one(A_new)

    tendency = diagn.tendencies.vor_tend
    tendency_copy = copy(tendency)
    dtendency = zero(tendency)

    L = model.time_stepping

    autodiff(Reverse, SpeedyWeather.leapfrog!, Const, Duplicated(A_old, dA_old), Duplicated(A_new, dA_new), Duplicated(tendency, dtendency), Const(dt), Const(lf1), Const(L))

    w1 = L.robert_filter*L.williams_filter/2   
    w2 = L.robert_filter*(1-L.williams_filter)/2   
    @test all(dtendency .≈ dt*(1+w1-w2))
    # ∂(tend) needs to be: dt* ( 1 + w1 - w2) (for every coefficient)

    #
    # transform!(diagn, progn, lf2, model)
    #

    fill!(diagn.tendencies, 0, PrimitiveWetModel)
    diag_copy = deepcopy(diagn)
    
    ddiag = one(diagn)
    ddiag_copy = deepcopy(ddiag)

    progn_copy = deepcopy(progn)
    dprogn = make_zero(progn)

    autodiff(Reverse, SpeedyWeather.transform!, Const, Duplicated(diagn, ddiag), Duplicated(progn, dprogn), Const(lf2), Const(model))
    autodiff(Reverse, SpeedyWeather.transform!, Const, Duplicated(diagn, ddiag), Duplicated(progn, dprogn), Const(lf2), Duplicated(model, make_zero(model)))

    function transform_diagn(diag, progn, lf2, model)
        diag_copy = deepcopy(diag)
        transform!(diag_copy, progn, lf2, model)
        return diag_copy
    end 

    fd_vjp = FiniteDifferences.j′vp(central_fdm(5,1), x -> transform_diagn(diag_copy, x, lf2, model), ddiag_copy, progn_copy)
    
    @test all(isapprox.(to_vec(fd_vjp[1])[1], to_vec(dprogn)[1],rtol=1e-3,atol=1e-3))
end




(; vor_grid, u_grid, v_grid ) = diagn.grid
(; scratch_memory) = diagn.dynamics  
U = diagn.dynamics.a            # reuse work arrays for velocities in spectral
V = diagn.dynamics.b            # reuse work arrays for velocities in spectral
   
                            # U = u*coslat, V=v*coslat
vor = progn.vor[2]             # relative vorticity at leapfrog step lf
S = model.spectral_transform

dvor_grid = similar(vor_grid)
dvor_grid .= 1 



autodiff(Reverse, transform!, Const, Duplicated(vor_grid, dvor_grid), Duplicated(vor, make_zero(vor)), Duplicated(scratch_memory, make_zero(scratch_memory)), Const(S))






function transform_test(diagn, progn, lf, scratch_memory, model)

    (; vor_grid, u_grid, v_grid ) = diagn.grid
    #(; scratch_memory) = diagn.dynamics  
    U = diagn.dynamics.a            # reuse work arrays for velocities in spectral
    V = diagn.dynamics.b            # reuse work arrays for velocities in spectral
    
                                # U = u*coslat, V=v*coslat
    vor = progn.vor[2]             # relative vorticity at leapfrog step lf
    S = model.spectral_transform

    transform!(vor_grid, vor, scratch_memory, S)    # get vorticity on grid from spectral vor

    SpeedyWeather.UV_from_vor!(U, V, vor, S)

    return nothing 
end 
    
ddiag = one(diagn)
dprogn = make_zero(progn)

autodiff(Reverse, transform_test, Const, Duplicated(diagn, ddiag), Duplicated(progn, dprogn), Const(lf2), Const(model))
autodiff(Reverse, transform_test, Const, Duplicated(diagn, ddiag), Duplicated(progn, dprogn), Const(lf2), Duplicated(diagn.dynamics.scratch_memory, make_zero(diagn.dynamics.scratch_memory)),Const(model))

U = diagn.dynamics.a 
V = diagn.dynamics.b            # reuse work arrays for velocities in spectral

dU = one(U)
dV = one(V)

vor = progn.vor[2]             # relative vorticity at leapfrog step lf
dvor = similar(vor)
dvor .= 1 
S = model.spectral_transform


autodiff(Reverse, SpeedyWeather.UV_from_vor!, Const, Duplicated(U, dU), Duplicated(V, dV), Duplicated(vor, dvor), Const(S))


Base.@propagate_inbounds call_transform(diagn, progn, lf, model) = transform_test(diagn, progn, lf, diagn.dynamics.scratch_memory, model)

autodiff(Reverse, call_transform, Const, Duplicated(diagn, ddiag), Duplicated(progn, dprogn), Const(lf2), Const(model))
