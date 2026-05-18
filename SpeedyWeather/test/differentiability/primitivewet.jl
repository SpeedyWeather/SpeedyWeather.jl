### Experiments going a bit deeper into the timestepping of the primtive wet model
# this script / these tests were mainly written for debugging, we might exclude it in future
# tests because it is quite maintanance heavy code
@testset "Differentiability: Primitive Wet Model Components" begin

    spectral_grid = SpectralGrid(trunc = 8, nlayers = 1)          # define resolution

    model = PrimitiveWetModel(; spectral_grid)  # construct model
    simulation = initialize!(model)
    initialize!(simulation)
    run!(simulation, period = Day(5)) # spin-up to get nonzero values for all fields
    initialize!(simulation; period = Day(1))


    adsim = ADSimulation(simulation)
    (; vars, model) = simulation
    (; Δt, Δt_millisec) = model.time_stepping
    dt = 2Δt

    # TO-DO: The first time we execute this, the gradient is different. Why?
    timestep_oop!(make_zero(vars), vars, dt, model)

    #
    # We go individually through all components of the time stepping and check
    # correctness
    #

    fill!(vars.tendencies, 0, PrimitiveWetModel)
    (; time) = progn.clock

    #
    # model physics
    #
    vars, dvars = ADseed(adsim, :tendencies)
    vars_copy = deepcopy(vars)
    dvars_copy = deepcopy(dvars)

    autodiff(Reverse, SpeedyWeather.parameterization_tendencies!, Const, Duplicated(vars, dvars), Const(model))

    function parameterization_tendencies(vars, model)
        vars_new = deepcopy(vars)
        SpeedyWeather.parameterization_tendencies!(vars_new, deepcopy(model))
        return vars_new
    end

    fd_vjp = FiniteDifferences.j′vp(central_fdm(11, 1), x -> parameterization_tendencies(x, model), dvars_copy, vars_copy)

    # TO-DO this test is broken, they gradients don't line up
    # old test checked dprogn
    @test all(isapprox.(to_vec(fd_vjp[1])[1], to_vec(dvars)[1], rtol = 1.0e-4, atol = 1.0e-1))

    #
    # ocean
    #
    vars, dvars = ADseed(adsim, :prognostic)
    vars_copy = deepcopy(vars)
    dvars_copy = deepcopy(dvars)

    autodiff(Reverse, SpeedyWeather.ocean_timestep!, Const, Duplicated(vars, dvars), Const(model))

    function ocean_timestep(vars, model)
        vars_new = deepcopy(vars)
        SpeedyWeather.ocean_timestep!(vars_new, deepcopy(model))
        return progn_new
    end

    fd_vjp = FiniteDifferences.j′vp(central_fdm(5, 1), x -> ocean_timestep(x, model), dvars_copy, vars_copy)

    # pass
    # old test checked ddiags
    @test all(isapprox.(to_vec(fd_vjp[1])[1], to_vec(dvars)[1], rtol = 1.0e-4, atol = 1.0e-2))

    #
    # land
    #

    vars, dvars = ADSeed(adsim, :prognostic)
    vars_copy = deepcopy(vars)
    dvars_copy = deepcopy(dvars)

    autodiff(Reverse, SpeedyWeather.land_timestep!, Const, Duplicated(vars, dvars), Const(model))

    function land_timestep(vars, model)
        vars_new = deepcopy(vars)
        SpeedyWeather.ocean_timestep!(vars_new, deepcopy(model))
        return vars_new
    end

    fd_vjp = FiniteDifferences.j′vp(central_fdm(5, 1), x -> land_timestep(x, model), dvars_copy, vars_copy)

    # pass
    # old test checked ddiagn
    @test all(isapprox.(to_vec(fd_vjp[1])[1], to_vec(dvars)[1], rtol = 1.0e-4, atol = 1.0e-2))

    #####
    # DYNAMICS
    lf2 = 2

    #
    # dynamics_tendencies!
    #

    vars, dvars = ADSeed(adsim, :tendencies)
    vars_copy = deepcopy(vars)
    dvars_copy = deepcopy(dvars)

    autodiff(Reverse, SpeedyWeather.dynamics_tendencies!, Const, Duplicated(vars, dvars), Const(lf2), Const(model))

    function dynamics_tendencies(vars, lf, model)
        vars_new = deepcopy(vars)
        SpeedyWeather.dynamics_tendencies!(vars_new, lf, deepcopy(model))
        return vars_new
    end

    fd_vjp = FiniteDifferences.j′vp(central_fdm(5, 1), x -> dynamics_tendencies(x, lf2, model), dvars_copy, vars_copy)

    # there are some NaNs in the FD, that's why this test is currently broken
    # old test checked dprogn
    @test all(isapprox.(to_vec(fd_vjp[1])[1], to_vec(dvars)[1], rtol = 1.0e-4, atol = 1.0e-1))

    #
    # Implicit correction
    #
    # continue here
    vars, dvars = ADSeed(adsim, :tendencies)
    vars_copy = deepcopy(vars)
    dvars_copy = deepcopy(dvars)

    autodiff(Reverse, SpeedyWeather.implicit_correction!, Const, Duplicated(vars, dvars), Const(model.implicit), Const(model))

    function implicit_correction(vars, implicit, model)
        vars_new = deepcopy(vars)
        SpeedyWeather.implicit_correction!(vars_new, deepcopy(implicit), deepcopy(model))
        return vars_new
    end

    fd_vjp = FiniteDifferences.j′vp(central_fdm(9, 1), x -> implicit_correction(x, model.implicit, model), dvars_copy, vars_copy)

    # old test checked dprogn
    @test all(isapprox.(to_vec(fd_vjp[1])[1], to_vec(dprogn)[1], rtol = 1.0e-4, atol = 1.0e-1))

    #
    # transform!(diagn, progn, lf2, model)
    #
    vars, dvars = ADSeed(adsim, :tendencies)
    vars_copy = deepcopy(vars)
    dvars_copy = deepcopy(dvars)

    autodiff(Reverse, SpeedyWeather.transform!, Const, Duplicated(vars, dvars), Const(lf2), Duplicated(model, make_zero(model)))

    function transform_diagn(vars, lf2, model)
        vars_copy = deepcopy(vars)
        transform!(vars_copy, lf2, deepcopy(model))
        return vars_copy
    end

    fd_vjp = FiniteDifferences.j′vp(central_fdm(5, 1), x -> transform_diagn(x, lf2, model), dvars_copy, vars_copy)

    # old test checked dprogn
    @test all(isapprox.(to_vec(fd_vjp[1])[1], to_vec(dvars)[1], rtol = 1.0e-3, atol = 1.0e-3))
end
