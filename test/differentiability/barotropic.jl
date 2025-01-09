
@testset "Differentiability: Barotropic Model Timestepping" begin 
    spectral_grid = SpectralGrid(trunc=31, nlayers=1)          # define resolution
    model = BarotropicModel(; spectral_grid)   # construct model
    simulation = initialize!(model)  
    initialize!(simulation)

    (; prognostic_variables, diagnostic_variables, model) = simulation
    (; Δt, Δt_millisec) = model.time_stepping

    progn = prognostic_variables
    diagn = diagnostic_variables

    d_progn = PrognosticVariables(spectral_grid)
    d_diag = DiagnosticVariables(spectral_grid)
    d_model = deepcopy(model)

    #SpeedyWeather.timestep!(progn, diagn, 2Δt, model) # calculate tendencies and leapfrog forward

    autodiff(Reverse, SpeedyWeather.timestep!, Const, Duplicated(progn, d_progn), Duplicated(diagn, d_diag), Const(2Δt), Duplicated(model, d_model))

    # differnetiate wrt initial conditions / previous state
    # write this as functions (progn_old, diagn_old, 2\Delta t, model -> progn, diagn)

    # differnetiate wrt parameter 
    # write this as function (model, progn, diagn, 2\Delta t) -> progn_new

end