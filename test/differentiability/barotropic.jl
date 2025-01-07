
@testset "Differentiability: Barotropic Model Timestepping" begin 
    spectral_grid = SpectralGrid(trunc=31, nlayers=1)          # define resolution
    model = BarotropicModel(; spectral_grid)   # construct model
    simulation = initialize!(model)  

    (; prognostic_variables, diagnostic_variables, model) = simulation
    (; clock) = prognostic_variables

    period = Day(10)
    progn = prognostic_variables
    diagn = diagnostic_variables

    # CLOCK
    SpeedyWeather.set_period!(clock, period)              # set how long to integrate for
    SpeedyWeather.initialize!(clock, model.time_stepping) # store the start date, reset counter

    (; Δt, Δt_millisec) = model.time_stepping

    # SCALING: we use vorticity*radius, divergence*radius in the dynamical core
    SpeedyWeather.scale!(progn, diagn, model.spectral_grid.radius)

    # OUTPUT INITIALISATION AND STORING INITIAL CONDITIONS + FEEDBACK
    # propagate spectral state to grid variables for initial condition output
    (; output, feedback) = model
    lf = 1                                  # use first leapfrog index
    SpeedyWeather.transform!(diagn, progn, lf, model, initialize=true)
    SpeedyWeather.initialize!(progn.particles, progn, diagn, model.particle_advection)
    SpeedyWeather.initialize!(output, feedback, progn, diagn, model)
    SpeedyWeather.initialize!(model.callbacks, progn, diagn, model)

    # FIRST TIMESTEPS: EULER FORWARD THEN 1x LEAPFROG
    # considered part of the model initialisation
    SpeedyWeather.first_timesteps!(progn, diagn, model)

    # only now initialise feedback for benchmark accuracy
    SpeedyWeather.initialize!(feedback, clock, model)

    d_progn = PrognosticVariables(spectral_grid)
    d_diag = DiagnosticVariables(spectral_grid)

    d_model = deepcopy(model)


    SpeedyWeather.timestep!(progn, diagn, 2Δt, model) # calculate tendencies and leapfrog forward

    autodiff(Reverse, SpeedyWeather.timestep!, Const, Duplicated(progn, d_progn), Duplicated(diagn, d_diag), Const(2Δt), Duplicated(model, d_model))

    # differnetiate wrt initial conditions / previous state
    # write this as functions (progn_old, diagn_old, 2\Delta t, model -> progn, diagn)

    # differnetiate wrt parameter 
    # write this as function (model, progn, diagn, 2\Delta t) -> progn_new

end