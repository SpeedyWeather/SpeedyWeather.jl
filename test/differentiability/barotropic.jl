### Experiments going a bit deeper into the timestepping of the barotropic model
@testset "Differentiability: Barotropic Model Timestepping" begin 
    # T15 still yields somewhat sensible dynamics, that's why it's chosen here
    spectral_grid = SpectralGrid(trunc=15, nlayers=1)          # define resolution
    model = BarotropicModel(; spectral_grid)   # construct model
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

    # test if differentiation works wrt copy! (there were some problems with it before)
    autodiff(Reverse, copy!, Const, Duplicated(progn_new, dprogn_new), Duplicated(progn, d_progn))

    # it's identity: 
    @test all(flatten(d_progn)[1] .== 1)
    @test all(flatten(d_progn)[2] .== 1)

    progn_new = zero(progn)
    dprogn_new = one(progn) # seed 

    # test that we can differentiate wrt an IC 
    autodiff(Reverse, timestep_oop!, Const, Duplicated(progn_new, dprogn_new), Duplicated(progn, d_progn), Duplicated(diagn, d_diag), Const(dt), Duplicated(model, d_model))

    # nonzero gradient
    @test sum(to_vec(d_progn)[1]) != 0

    # FD comparison 
    dprogn_2 = one(progn) # seed 

    fd_jvp = FiniteDifferences.j′vp(central_fdm(5,1), x -> timestep_oop(x, diagn_copy, dt, model), dprogn_2, progn_copy)

    @test isapprox(to_vec(fd_jvp[1])[1], to_vec(d_progn)[1])

    ### Test the leapfrog 

    lf1 = 2 
    lf2 = 2 

    # set the tendencies back to zero for accumulation
    fill!(diagn.tendencies, 0, Barotropic)

    # TENDENCIES, DIFFUSION, LEAPFROGGING AND TRANSFORM SPECTRAL STATE TO GRID
    SpeedyWeather.dynamics_tendencies!(diagn, progn, lf2, model)
    SpeedyWeather.horizontal_diffusion!(diagn, progn, model.horizontal_diffusion, model)


    progn_copy = deepcopy(progn)
    dprogn = one(progn_copy)
    dprogn_copy = one(progn_copy)

    tend = diagn.tendencies
    tend_copy = deepcopy(tend)
    dtend = make_zero(tend)
    dmodel = make_zero(model)

    #leapfrog!(progn, diagn.tendencies, dt, lf1, model)
    autodiff(Reverse, SpeedyWeather.leapfrog!, Const, Duplicated(progn, dprogn), Duplicated(tend, dtend), Const(dt), Const(lf1), Const(model))

    function leapfrog_step(progn::PrognosticVariables, tend, dt, lf, model)
        SpeedyWeather.leapfrog!(progn, tend, dt, lf, model)
        return progn
    end 

    tend_one = deepcopy(tend_copy)
    fill!(tend_one, 1, Barotropic)

    fd_jvp = FiniteDifferences.j'vp(central_fdm(5,1), x -> leapfrog_step(progn, x, dt, lf1, model),(tend_copy, tend_one))

    # single variable leapfrog step 

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
    @test all(dtendency .== dt*(1+w1-w2))
    # d tend needs to be: dt* ( 1 + w1 - w2) (for every coefficient)
    # enzyme shows it is 
    
    function leapfrog(A_old, A_new, tendency, dt, lf, L)
        A_old_new = copy(A_old)
        A_new_new = copy(A_new)
        SpeedyWeather.leapfrog!(A_old_new, A_new_new, tendency, dt, lf, L)
        return A_new_new
    end 

    fd_jvp = FiniteDifferences.jvp(central_fdm(5,1), (x) ->leapfrog(A_old_copy, A_old_copy, x, dt, lf1, L), (tendency_copy, one(A_new_copy)))
    # in this case it need to be dt*(1 - w2) (no contributation from A_old in FD)

    
    # differnetiate wrt parameter 
    # write this as function (model, progn, diagn, 2\Delta t) -> progn_new
end


