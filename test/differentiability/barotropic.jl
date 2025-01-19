
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
    d_diag = DiagnosticVariables(spectral_grid, model)
    d_model = deepcopy(model)

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

    # differnetiate wrt parameter 
    # write this as function (model, progn, diagn, 2\Delta t) -> progn_new

end