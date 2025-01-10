
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

    autodiff(Reverse, timestep_oop!, Const, Duplicated(progn_new, dprogn_new), Duplicated(progn, d_progn), Duplicated(diagn, d_diag), Const(2Δt), Duplicated(model, d_model))

    # differnetiate wrt initial conditions / previous state

    # new seeed 
    dprogn_new_2 = one(progn)

    fd_jvp = FiniteDifferences.j′vp(central_fdm(5,1), x -> timestep_oop(x, diagn_copy, 2Δt, model), dprogn_new_2, progn_copy )
    

    # write this as functions (progn_old, diagn_old, 2\Delta t, model -> progn, diagn)
    vor = progn.vor[1]

    fd_jvp = FiniteDifferences.j′vp(central_fdm(5,1), x -> step_vorticity(x, progn_copy, diagn_copy, 2Δt, model), dprogn_new_2, progn_copy )


    # differnetiate wrt parameter 
    # write this as function (model, progn, diagn, 2\Delta t) -> progn_new

end