
@testset "Differentiability: Barotropic Model Timestepping" begin 
    spectral_grid = SpectralGrid(trunc=31, nlayers=1)          # define resolution
    model = BarotropicModel(; spectral_grid)   # construct model
    simulation = initialize!(model)  
    initialize!(simulation)

    (; prognostic_variables, diagnostic_variables, model) = simulation
    (; Δt, Δt_millisec) = model.time_stepping

    progn = prognostic_variables
    diagn = diagnostic_variables

    diagn_copy = deepcopy(diagn)
    progn_copy = deepcopy(progn)

    d_progn = zero(progn)
    d_diag = DiagnosticVariables(spectral_grid, model)
    d_model = deepcopy(model)

    
    #SpeedyWeather.timestep!(progn, diagn, 2Δt, model) # calculate tendencies and leapfrog forward


    progn_new = zero(progn)
    dprogn_new = one(progn) # seed 

    # this is a test of the full timestep of the model. We do this just to see that it causes no errors
    # and the gradient is non-zero. we don't have something to compare it against
    autodiff(Reverse, timestep_oop!, Const, Duplicated(progn_new, dprogn_new), Duplicated(progn, d_progn), Duplicated(diagn, d_diag), Const(2Δt), Duplicated(model, d_model))

    @test sum(d_progn.vor[1]) != 0
    #dprogn_new.vor[1]

    # differnetiate vorticity wrt initial conditions / previous state

    # new seeed 
    vor = progn.vor[1] 
    vor_copy = deepcopy(vor)
    dvor = zero(vor)

    vor_new = zero(vor)
    dvor_new = zero(vor)
    fill!(dvor_new, 1+1im)

    #somethings not right

    diagn_copy = deepcopy(diagn)
    progn_copy = deepcopy(progn)
    d_progn = zero(progn)
    d_diagn = DiagnosticVariables(spectral_grid, model)

    # step_vorticity!(vor_new, vor, progn, diagn, 2Δt, model)
    autodiff(Reverse, step_vorticity!, Const, Duplicated(vor_new, dvor_new), Duplicated(vor, dvor), Duplicated(progn, d_progn), Duplicated(diagn, d_diag), Const(2Δt), Duplicated(model, d_model))

    # To-DO: there's a problem with copy!(LTA), do something about it 

    # write this as functions (progn_old, diagn_old, 2\Delta t, model -> progn, diagn)
    dvor_2 = one(progn.vor[1])
    fd_jvp = FiniteDifferences.j′vp(central_fdm(5,1), x -> step_vorticity(x, progn_copy, diagn_copy, 2Δt, model), dvor_2, vor_copy)



    # differnetiate wrt parameter 
    # write this as function (model, progn, diagn, 2\Delta t) -> progn_new

end