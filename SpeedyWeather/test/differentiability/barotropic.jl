### Experiments going a bit deeper into the timestepping of the barotropic model

#
# We go individually through all components of the time stepping and check
# correctness
#

#
# dynamics_tendencies!
#
@testset "Differentiability: dynamics_tendencies! on Barotropic model" begin
    # T15 still yields somewhat sensible dynamics, that's why it's chosen here
    spectral_grid = SpectralGrid(trunc = 9, nlayers = 1) # define resolution
    model = BarotropicModel(; spectral_grid) # construct model
    simulation = initialize_with_spinup!(model)
    lf2 = 2

    adsim = ADSimulation(simulation)
    vars, dvars = diagnosticseed(adsim)

    @info "Running reverse-mode AD"
    @time autodiff(Reverse, SpeedyWeather.dynamics_tendencies!, Const, Duplicated(vars, dvars), Const(lf2), Const(model))

    # basic sanity checks for VJP
    dvec, _ = to_vec(dvars)
    @test all(isfinite.(dvec))
    @test any(abs.(dvec) .> 0)

    adsim2 = ADSimulation(simulation)
    vars2, dvars2 = prognosticseed(adsim2)

    function dynamics_tendencies(vars, lf, model)
        vars_new = deepcopy(vars)
        SpeedyWeather.dynamics_tendencies!(vars_new, lf, deepcopy(model))
        return vars_new
    end

    fdsim = ADSimulation(simulation)
    @info "Running finite differences"
    fd_vjp = @time FiniteDifferences.j′vp(central_fdm(15, 1), x -> dynamics_tendencies(x, lf2, model), make_zero(adsim.dvars), fdsim.vars)

    # this is currently failing, possibly due to problems with finite diff?
    @test_broken all(isapprox.(to_vec(fd_vjp[1])[1], to_vec(dvars)[1], rtol = 1.0e-1, atol = 1.0e-1))
end

#
# horizontal_diffusion!
#
@testset "Differentiability: horizontal_diffusion! on Barotropic model" begin
    spectral_grid = SpectralGrid(trunc = 9, nlayers = 1)
    model = BarotropicModel(; spectral_grid)
    simulation = initialize_with_spinup!(model)

    lf1 = 1
    adsim = ADSimulation(simulation)
    vars, dvars = diagnosticseed(adsim)

    @time autodiff(Reverse, SpeedyWeather.horizontal_diffusion!, Const, Duplicated(vars, dvars), Const(model.horizontal_diffusion), Const(model), Const(lf1))

    # ∂(progn)
    # should be row-wise `model.horizontal_diffusion.impl .* model.horizontal_diffusion.expl`
    # for all variables that are diffused
    diff_coefficient = model.horizontal_diffusion.impl .* model.horizontal_diffusion.expl
    l_indices = [(1:l) for l in 1:spectral_grid.spectrum.mmax]
    for (i, il) in enumerate(l_indices)
        @test all(real.(Matrix(dvars.prognostic.vor[:, 1, lf1])[i, il]) .≈ diff_coefficient[i])
    end

    # ∂(tend_old)
    # should be row-wise `model.horizontal_diffusion.impl`
    for (i, il) in enumerate(l_indices)
        @test all(real.(Matrix(dvars.tendencies.vor[:, 1])[i, il]) .≈ model.horizontal_diffusion.impl[i])
    end
end

#
# Test the leapfrog
#
@testset "Differentiability: leapfrog! on Barotropic model" begin
    spectral_grid = SpectralGrid(trunc = 9, nlayers = 1)
    model = BarotropicModel(; spectral_grid)
    simulation = initialize_with_spinup!(model)

    (; Δt, Δt_millisec) = simulation.model.time_stepping
    dt = 2Δt
    lf1 = 2
    lf2 = 2

    adsim = ADSimulation(simulation)
    vars, dvars = prognosticseed(adsim)

    tend = adsim.vars.tendencies
    tend_copy = deepcopy(tend)
    dtend = make_zero(tend)

    @info "Running reverse-mode AD"
    @time autodiff(Reverse, SpeedyWeather.leapfrog!, Const, Duplicated(vars, dvars), Const(dt), Const(lf1), Const(model))

    function leapfrog_step(vars_new::Variables, vars_in::Variables, dt, lf, model)
        copy!(vars_new, vars_in)
        SpeedyWeather.leapfrog!(vars_new, dt, lf, model)
        return vars_new
    end

    fdsim = ADSimulation(simulation)
    vars_new = deepcopy(fdsim.vars)

    @info "Running finite differences"
    fd_vjp = @time FiniteDifferences.j′vp(central_fdm(5, 1), x -> leapfrog_step(vars_new, x, dt, lf1, model), make_zero(dvars), tend_copy)

    @test all(isapprox.(to_vec(fd_vjp[1])[1], to_vec(dvars)[1], rtol = 1.0e-3, atol = 1.0e-3))

    #
    # single variable leapfrog step
    #

    A_old = vars.prognostic.vor[:, :, 1]
    A_old_copy = deepcopy(A_old)
    dA_old = one(A_old)

    A_new = vars.prognostic.vor[:, :, 2]
    A_new_copy = deepcopy(A_new)
    dA_new = one(A_new)

    tendency = adsim.vars.tendencies.vor
    tendency_copy = deepcopy(tendency)
    dtendency = make_zero(tendency)

    L = model.time_stepping

    @info "Running reverse-mode AD"
    @time autodiff(Reverse, SpeedyWeather.leapfrog!, Const, Duplicated(A_old, dA_old), Duplicated(A_new, dA_new), Duplicated(tendency, dtendency), Const(dt), Const(lf1), Const(L))

    w1 = L.robert_filter * L.williams_filter / 2
    w2 = L.robert_filter * (1 - L.williams_filter) / 2
    @test all(dtendency .≈ dt * (1 + w1 - w2))
    # ∂(tend) needs to be: dt* ( 1 + w1 - w2) (for every coefficient)
end

@testset "Differentiability: transform!(::Diagnostic, ::Prognostic)" begin
    spectral_grid = SpectralGrid(trunc = 9, nlayers = 1)
    model = BarotropicModel(; spectral_grid)
    simulation = initialize_with_spinup!(model)

    (; Δt, Δt_millisec) = simulation.model.time_stepping
    dt = 2Δt
    lf1 = 2
    lf2 = 2

    adsim = ADSimulation(simulation)
    vars, dvars = prognosticseed(adsim)

    @info "Running reverse-mode AD"
    @time autodiff(Reverse, SpeedyWeather.transform!, Const, Duplicated(vars, dvars), Const(lf2), Const(model))

    function transform_step(vars_in::Variables, lf, model)
        vars_new = deepcopy(vars_in)
        SpeedyWeather.transform!(vars_new, lf, model)
        return vars_new
    end

    fdsim = ADSimulation(simulation)

    @info "Running finite differences"
    fd_vjp = @time FiniteDifferences.j′vp(central_fdm(5, 1), x -> transform_step(x, lf2, model), make_zero(dvars), fdsim.vars)

    @test all(isapprox.(to_vec(fd_vjp[1])[1], to_vec(dvars)[1], rtol = 1.0e-3, atol = 1.0e-3))

end

@testset "Differentiability: timestep! on Barotropic model" begin
    # T15 still yields somewhat sensible dynamics, that's why it's chosen here
    spectral_grid = SpectralGrid(trunc = 9, nlayers = 1)
    model = BarotropicModel(; spectral_grid)
    simulation = initialize_with_spinup!(model)
    (; Δt, Δt_millisec) = simulation.model.time_stepping
    dt = 2Δt

    fdsim = ADSimulation(simulation)
    @info "Running finite differences"
    # for the full timestep, we need a bit higher precision
    fd_vjp = @time FiniteDifferences.j′vp(
        central_fdm(15, 1),
        x -> timestep_oop(x, dt, deepcopy(fdsim.model)),
        make_zero(fdsim.vars),
        deepcopy(fdsim.vars)
    )

    adsim = ADSimulation(simulation)
    vars_new, dvars_new = prognosticseed(adsim)

    @info "Running reverse-mode AD"
    # test that we can differentiate wrt to everything
    @time autodiff(
        Reverse,
        timestep_oop!,
        Const,
        Duplicated(vars_new, dvars_new),
        Duplicated(adsim.vars, adsim.dvars),
        Const(dt),
        Duplicated(model, make_zero(model))
    )

    # nonzero gradient
    dvars = adsim.dvars
    @test sum(to_vec(dvars)[1]) != 0

    @test_broken isapprox(to_vec(fd_vjp[1])[1], to_vec(dvars)[1], rtol = 0.05) # we have to go really quite high with the tolerances here
    @test mean(abs.(to_vec(fd_vjp[1])[1] - to_vec(dvars)[1])) < 0.002 # so we check a few extra statistics
    @test maximum(to_vec(fd_vjp[1].prognostic.vor)[1] - to_vec(dvars.prognostic.vor)[1]) < 0.05

    # test that we can differentiate with Const(Model) only wrt to the state
    adsim2 = ADSimulation(simulation)
    vars_new, dvars_new = prognosticseed(adsim2)

    @time autodiff(
        set_runtime_activity(Reverse),
        timestep_oop!,
        Const,
        Duplicated(vars_new, dvars_new),
        Duplicated(adsim2.vars, adsim2.dvars),
        Const(dt),
        Const(model)
    )

    d_vars = adsim2.dvars
    @test_broken isapprox(to_vec(fd_vjp[1])[1], to_vec(d_vars)[1], rtol = 0.05) # we have to go really quite high with the tolerances here
    @test mean(abs.(to_vec(fd_vjp[1])[1] - to_vec(d_vars)[1])) < 0.002 # so we check a few extra statistics
    @test maximum(to_vec(fd_vjp[1].prognostic.vor)[1] - to_vec(d_vars.prognostic.vor)[1]) < 0.05
end

@testset "Differentiability: Barotropic model parameters" begin
    # T15 still yields somewhat sensible dynamics, that's why it's chosen here
    spectral_grid = SpectralGrid(trunc = 9, nlayers = 1)          # define resolution
    model = BarotropicModel(; spectral_grid)   # construct model
    simulation = initialize_with_spinup!(model)
    (; Δt, Δt_sec) = simulation.model.time_stepping
    dt = Δt
    ps = parameters(model)
    pvec = vec(ps)
    adsim = ADSimulation(simulation)
    dp = zero(pvec)
    vars_new, dvars_new = deepcopy(adsim.vars), make_zero(adsim.dvars)
    @time autodiff(
        Reverse,
        timestep_oop!,
        Const,
        Duplicated(vars_new, dvars_new),
        Duplicated(adsim.vars, adsim.dvars),
        Const(dt),
        Duplicated(adsim.model, make_zero(adsim.model)),
        Duplicated(pvec, dp),
    )

    fdsim = ADSimulation(simulation)
    fd_vjp = @time FiniteDifferences.j′vp(
        central_fdm(10, 1),
        x -> timestep_oop(deepcopy(fdsim.vars), dt, deepcopy(fdsim.model), x),
        make_zero(fdsim.vars),
        copy(pvec),
    )
    @test all(isapprox.(dp, fd_vjp[1], atol = 1.0e-5, rtol = 1.0e-3))
end
