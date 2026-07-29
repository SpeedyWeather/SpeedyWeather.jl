### Experiments going a bit deeper into the timestepping of the barotropic model
#
# We go individually through all components of the time stepping and check
# correctness.
#
# Seed the global RNG: the BarotropicModel's RandomVelocity process uses `seed = 0`,
# i.e. it reseeds itself from Julia's global RNG on every `initialize!`. Without a
# fixed global seed the dynamics — and hence the finite-difference comparisons below —
# are nondeterministic and can spuriously pass/fail between runs.
import Random
Random.seed!(123)

#
# dynamics_tendencies!
#
@testset "Differentiability: dynamics_tendencies! on Barotropic model" begin
    # T9 still yields somewhat sensible dynamics, that's why it's chosen here
    spectral_grid = SpectralGrid(trunc = 9, nlayers = 1) # define resolution
    model = BarotropicModel(; spectral_grid, time_stepping = Leapfrog(spectral_grid)) # construct model
    simulation = initialize_with_spinup!(model)

    adsim = ADSimulation(simulation)

    vars, dvars = ADseed(adsim, :tendencies)

    @info "Running reverse-mode AD"
    @time autodiff(set_runtime_activity(Reverse), SpeedyWeather.dynamics_tendencies!, Const, Duplicated(vars, dvars), Const(model))

    # basic sanity checks for VJP: gradient should be finite and nonzero somewhere
    dvec, _ = to_vec(dvars)
    @test all(isfinite.(dvec))
    @test any(abs.(dvec) .> 0)

    vars2, dvars2 = ADseed(adsim, :tendencies)

    function dynamics_tendencies(vars, model)
        vars_new = deepcopy(vars)
        SpeedyWeather.dynamics_tendencies!(vars_new, deepcopy(model))
        return vars_new
    end

    @info "Running finite differences"
    fd_vjp = @time FiniteDifferences.j′vp(central_fdm(15, 1), x -> dynamics_tendencies(x, model), dvars2, vars2)

    # this is currently failing, possibly due to problems with finite diff?
    @test_broken all(isapprox.(to_vec(fd_vjp[1].grid)[1], to_vec(dvars.grid)[1], rtol = 1.0e-1, atol = 1.0e-1))
end

#
# horizontal_diffusion!
#
@testset "Differentiability: horizontal_diffusion! on Barotropic model" begin
    spectral_grid = SpectralGrid(trunc = 9, nlayers = 1)
    model = BarotropicModel(; spectral_grid, time_stepping = Leapfrog(spectral_grid))
    simulation = initialize_with_spinup!(model)

    adsim = ADSimulation(simulation)
    vars, dvars = ADseed(adsim, :tendencies)

    @time autodiff(set_runtime_activity(Reverse), SpeedyWeather.horizontal_diffusion!, Const, Duplicated(vars, dvars), Const(model.horizontal_diffusion), Const(model))

    # the barotropic diffusion reads prognostic vorticity step 1 (which_prognostic_step
    # for AbstractHorizontalDiffusion == 1), tendency step 1, and applies (per degree l)
    #   tendency[lm] = (tendency[lm] + expl[l]*vor[lm]) * impl[l]
    # so   ∂tendency_out/∂vor          = expl[l]*impl[l]
    # and  ∂tendency_out/∂tendency_in  = impl[l]
    # with the seed (1+im) on the (mutated) output tendency the real part recovers these.
    expl = model.horizontal_diffusion.expl
    impl = model.horizontal_diffusion.impl

    dvor = dvars.prognostic.vorticity[:, 1, 1]      # LowerTriangularMatrix, step 1 (the step diffusion reads)
    dtend = dvars.tendencies.vorticity[:, 1, 1]     # LowerTriangularMatrix, tendency step 1
    l_idx = dvor.spectrum.l_indices                 # maps packed index lm -> degree l

    # ∂(progn): row-wise `expl .* impl`
    @test all(real(dvor[lm]) ≈ expl[l_idx[lm], 1] * impl[l_idx[lm], 1] for lm in eachindex(dvor))
    # ∂(tend_old): row-wise `impl`
    @test all(real(dtend[lm]) ≈ impl[l_idx[lm], 1] for lm in eachindex(dtend))
end

#
# Test the leapfrog time step (now `update_prognostic!`)
#
@testset "Differentiability: update_prognostic! (leapfrog) on Barotropic model" begin
    spectral_grid = SpectralGrid(trunc = 9, nlayers = 1)
    model = BarotropicModel(; spectral_grid, time_stepping = Leapfrog(spectral_grid))
    simulation = initialize_with_spinup!(model)

    adsim = ADSimulation(simulation)

    vars, dvars = ADseed(adsim, :prognostic)

    @info "Running reverse-mode AD"
    @time autodiff(set_runtime_activity(Reverse), SpeedyWeather.update_prognostic!, Const, Duplicated(vars, dvars), Const(model))

    function update_prognostic_step(vars_in::Variables, model)
        vars_new = deepcopy(vars_in)
        SpeedyWeather.update_prognostic!(vars_new, deepcopy(model))
        return vars_new
    end

    vars_new, dvars_new = ADseed(adsim, :prognostic)

    @info "Running finite differences"
    fd_vjp = @time FiniteDifferences.j′vp(central_fdm(5, 1), x -> update_prognostic_step(x, model), dvars_new, vars_new)

    @test all(isapprox.(to_vec(fd_vjp[1])[1], to_vec(dvars)[1], rtol = 1.0e-3, atol = 1.0e-3))

    #
    # single variable leapfrog step via the per-variable update_prognostic!
    # (the direct analogue of the old single-variable `leapfrog!`)
    #
    vars2, _ = ADseed(adsim, :prognostic)
    L = model.time_stepping
    clock = vars2.prognostic.clock
    clock.step_counter = 2  # force the filtered leapfrog step (lf == 2) to exercise Robert+Williams filters

    vor = vars2.prognostic.vorticity            # full LowerTriangularArray with the 2 leapfrog steps
    dvor = one(vor)                              # seed both steps with 1

    tendency = deepcopy(vars2.tendencies.vorticity)
    dtendency = make_zero(tendency)

    implicit = model.implicit                   # `nothing` for the BarotropicModel

    @info "Running reverse-mode AD"
    @time autodiff(
        set_runtime_activity(Reverse), SpeedyWeather.update_prognostic!, Const,
        Duplicated(vor, dvor), Duplicated(tendency, dtendency),
        Const(clock), Const(L), Const(implicit), Const(model),
    )

    lf = SpeedyWeather.prognostic_step(L, clock)
    Δt_eff = SpeedyWeather.time_step(L, clock)
    w1 = (lf - 1) * L.robert_filter * L.williams_filter / 2
    w2 = (lf - 1) * L.robert_filter * (1 - L.williams_filter) / 2
    # ∂(tend) needs to be: Δt * (1 + w1 - w2) (for every coefficient)
    @test all(dtendency .≈ Δt_eff * (1 + w1 - w2))
end

@testset "Differentiability: transform!(::Variables)" begin
    spectral_grid = SpectralGrid(trunc = 9, nlayers = 1)
    model = BarotropicModel(; spectral_grid, time_stepping = Leapfrog(spectral_grid))
    simulation = initialize_with_spinup!(model)

    adsim = ADSimulation(simulation)
    vars, dvars = ADseed(adsim, :grid)

    @info "Running reverse-mode AD"
    @time autodiff(set_runtime_activity(Reverse), SpeedyWeather.transform!, Const, Duplicated(vars, dvars), Const(model))

    function transform_step(vars_in::Variables, model)
        vars_new = deepcopy(vars_in)
        SpeedyWeather.transform!(vars_new, deepcopy(model))
        return vars_new
    end

    vars_new, dvars_new = ADseed(adsim, :grid)

    @info "Running finite differences"
    fd_vjp = @time FiniteDifferences.j′vp(central_fdm(11, 1), x -> transform_step(x, model), dvars_new, vars_new)

    @test all(isapprox.(to_vec(fd_vjp[1].prognostic)[1], to_vec(dvars.prognostic)[1], rtol = 1.0e-3, atol = 1.0e-3))
end

@testset "Differentiability: time_step! on Barotropic model" begin
    # T9 still yields somewhat sensible dynamics, that's why it's chosen here
    spectral_grid = SpectralGrid(trunc = 9, nlayers = 1)
    model = BarotropicModel(; spectral_grid, time_stepping = Leapfrog(spectral_grid))
    simulation = initialize_with_spinup!(model)

    adsim = ADSimulation(simulation)
    vars_fd, dvars_fd = ADseed(adsim, :prognostic)

    @info "Running finite differences"
    # for the full timestep, we need a bit higher precision
    fd_vjp = @time FiniteDifferences.j′vp(
        central_fdm(21, 1),
        x -> timestep_oop(x, deepcopy(adsim.model)),
        dvars_fd,
        vars_fd
    )

    vars_new, dvars_new = ADseed(adsim, :prognostic)

    @info "Running reverse-mode AD"
    # test that we can differentiate wrt to everything
    @time autodiff(
        set_runtime_activity(Reverse),
        timestep_oop!,
        Const,
        Duplicated(vars_new, dvars_new),
        Duplicated(adsim.vars, adsim.dvars),
        Duplicated(model, make_zero(model))
    )

    # nonzero gradient of the prognostic variables
    dvars = adsim.dvars
    @test sum(to_vec(dvars.prognostic)[1]) != 0


    @test isapprox(to_vec(fd_vjp[1])[1], to_vec(dvars)[1], rtol = 0.05) # we have to go really quite high with the tolerances here
    @test mean(abs.(to_vec(fd_vjp[1])[1] - to_vec(dvars)[1])) < 0.002 # so we check a few extra statistics
    @test maximum(to_vec(fd_vjp[1].prognostic.vorticity)[1] - to_vec(dvars.prognostic.vorticity)[1]) < 0.05

    # test that we can differentiate with Const(Model) only wrt to the state
    vars_new, dvars_new = ADseed(adsim, :prognostic)

    @time autodiff(
        set_runtime_activity(Reverse),
        timestep_oop!,
        Const,
        Duplicated(vars_new, dvars_new),
        Duplicated(adsim.vars, adsim.dvars),
        Const(model)
    )

    d_vars = adsim.dvars
    @test isapprox(to_vec(fd_vjp[1])[1], to_vec(d_vars)[1], rtol = 0.05) # we have to go really quite high with the tolerances here
    @test mean(abs.(to_vec(fd_vjp[1])[1] - to_vec(d_vars)[1])) < 0.002 # so we check a few extra statistics
    @test maximum(to_vec(fd_vjp[1].prognostic.vorticity)[1] - to_vec(d_vars.prognostic.vorticity)[1]) < 0.05
end

@testset "Differentiability: Barotropic model parameters" begin
    # T9 still yields somewhat sensible dynamics, that's why it's chosen here
    spectral_grid = SpectralGrid(trunc = 9, nlayers = 1)          # define resolution
    model = BarotropicModel(; spectral_grid, time_stepping = Leapfrog(spectral_grid))   # construct model
    simulation = initialize_with_spinup!(model)
    ps = parameters(model)
    pvec = vec(ps)
    adsim = ADSimulation(simulation)
    dp = zero(pvec)
    vars_new, dvars_new = ADseed(adsim, :prognostic)
    @time autodiff(
        set_runtime_activity(Reverse),
        timestep_oop!,
        Const,
        Duplicated(vars_new, dvars_new),
        Duplicated(adsim.vars, adsim.dvars),
        Duplicated(adsim.model, make_zero(adsim.model)),
        Duplicated(pvec, dp),
    )

    fdsim = ADSimulation(simulation)
    fd_vjp = @time FiniteDifferences.j′vp(
        central_fdm(10, 1),
        x -> timestep_oop(deepcopy(fdsim.vars), deepcopy(fdsim.model), x),
        make_zero(fdsim.vars),
        copy(pvec),
    )
    @test all(isapprox.(dp, fd_vjp[1], atol = 1.0e-5, rtol = 1.0e-3))
end

#
# Same full time step, but with the (default) NCycleLorenz time stepper instead of Leapfrog.
# Here we only check that Enzyme differentiates without error and yields a finite, nonzero
# gradient (plus a loose finite-difference sanity statistic).
#
@testset "Differentiability: time_step! on Barotropic model with NCycleLorenz" begin
    spectral_grid = SpectralGrid(trunc = 9, nlayers = 1)
    model = BarotropicModel(; spectral_grid, time_stepping = SpeedyWeather.NCycleLorenz(spectral_grid))
    simulation = initialize_with_spinup!(model)

    @test model.time_stepping isa SpeedyWeather.NCycleLorenz

    adsim = ADSimulation(simulation)
    vars_new, dvars_new = ADseed(adsim, :prognostic)

    @info "Running reverse-mode AD (NCycleLorenz)"
    @time autodiff(
        set_runtime_activity(Reverse),
        timestep_oop!,
        Const,
        Duplicated(vars_new, dvars_new),
        Duplicated(adsim.vars, adsim.dvars),
        Duplicated(model, make_zero(model)),
    )

    dvars = adsim.dvars
    @test all(isfinite.(to_vec(dvars.prognostic)[1]))
    @test sum(to_vec(dvars.prognostic)[1]) != 0

    vars_fd, dvars_fd = ADseed(adsim, :prognostic)

    @info "Running finite differences (NCycleLorenz)"
    fd_vjp = @time FiniteDifferences.j′vp(
        central_fdm(11, 1),
        x -> timestep_oop(x, deepcopy(adsim.model)),
        dvars_fd,
        vars_fd,
    )

    @test isapprox(to_vec(fd_vjp[1])[1], to_vec(dvars)[1], rtol = 0.05)
    @test mean(abs.(to_vec(fd_vjp[1])[1] - to_vec(dvars)[1])) < 0.01
end
