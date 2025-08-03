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
    spectral_grid = SpectralGrid(trunc=9, nlayers=1) # define resolution
    model = BarotropicModel(; spectral_grid) # construct model
    simulation = initialize_with_spinup!(model)
    lf2 = 2 

    adsim = ADSimulation(simulation)
    diagn, ddiagn = diagnosticseed(adsim)

    @info "Running reverse-mode AD"
    @time autodiff(Reverse, SpeedyWeather.dynamics_tendencies!, Const, Duplicated(diagn, ddiagn), Duplicated(adsim.progvars, adsim.dprogvars), Const(lf2), Const(model))

    # basic sanity checks for VJP
    dprogvec, _ = to_vec(adsim.dprogvars)
    @test all(isfinite.(dprogvec))
    @test any(abs.(dprogvec) .> 0)

    adsim2 = ADSimulation(simulation)
    progn, dprogn = prognosticseed(adsim2)
    # doesn't work currently, missing rules
    # @time autodiff(Forward, SpeedyWeather.dynamics_tendencies!, Const, Duplicated(adsim2.diagvars, adsim2.ddiagvars), Duplicated(progn, dprogn), Const(lf2), Duplicated(model, make_zero(model)))

    function dynamics_tendencies(diagn, progn, lf, model)
        diagn_new = deepcopy(diagn)
        SpeedyWeather.dynamics_tendencies!(diagn_new, deepcopy(progn), lf, deepcopy(model))
        return diagn_new
    end
    
    fdsim = ADSimulation(simulation)
    @info "Running finite differences"
    fd_vjp = @time FiniteDifferences.j′vp(central_fdm(15,1), x -> dynamics_tendencies(fdsim.diagvars, x, lf2, model), one(ddiagn), fdsim.progvars)
    
    # this is currently failing, possibly due to problems with finite diff?
    @test_broken all(isapprox.(to_vec(fd_vjp[1])[1], to_vec(adsim.dprogvars)[1], rtol=1e-1, atol=1e-1))
end

#
# horizontal_diffusion!
#
@testset "Differentiability: horizontal_diffusion! on Barotropic model" begin
    spectral_grid = SpectralGrid(trunc=9, nlayers=1)
    model = BarotropicModel(; spectral_grid)
    simulation = initialize_with_spinup!(model)

    lf1 = 1
    adsim = ADSimulation(simulation)
    dprogn = adsim.dprogvars
    diagn, ddiagn = diagnosticseed(adsim)

    @time autodiff(Reverse, SpeedyWeather.horizontal_diffusion!, Const, Duplicated(diagn, ddiagn), Duplicated(adsim.progvars, adsim.dprogvars), Const(model.horizontal_diffusion), Const(model), Const(lf1))

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
    l_indices = [(1:l) for l=1:spectral_grid.spectrum.mmax]
    for (i,il) in enumerate(l_indices)
        @test all(real.(Matrix(dprogn.vor[:,1,lf1])[i, il]) .≈ diff_coefficient[i])
    end 

    # ∂(tend_old)
    # should be row-wise `model.horizontal_diffusion.impl` 
    for (i,il) in enumerate(l_indices)
        @test all(real.(Matrix(ddiagn.tendencies.vor_tend[:,1])[i, il]) .≈ model.horizontal_diffusion.impl[i])
    end
end

#
# Test the leapfrog 
# 
@testset "Differentiability: leapfrog! on Barotropic model" begin
    spectral_grid = SpectralGrid(trunc=9, nlayers=1)
    model = BarotropicModel(; spectral_grid)
    simulation = initialize_with_spinup!(model)

    lf1 = 2
    lf2 = 2 

    adsim = ADSimulation(simulation)
    progn, dprogn = prognosticseed(adsim)

    tend = adsim.diagvars.tendencies
    tend_copy = deepcopy(tend)
    dtend = make_zero(tend)

    @info "Running reverse-mode AD"
    @time autodiff(Reverse, SpeedyWeather.leapfrog!, Const, Duplicated(progn, dprogn), Duplicated(tend, dtend), Const(dt), Const(lf1), Const(model))

    function leapfrog_step(progn_new::PrognosticVariables, progn::PrognosticVariables, tend, dt, lf, model)
        copy!(progn_new, progn)
        SpeedyWeather.leapfrog!(progn_new, tend, dt, lf, model)
        return progn_new
    end 

    fdsim = ADSimulation(simulation)
    progn_new = deepcopy(fdsim.progvars)

    @info "Running finite differences"
    fd_vjp = @time FiniteDifferences.j′vp(central_fdm(5,1), x -> leapfrog_step(progn_new, fdsim.progvars, x, dt, lf1, model), one(dprogn), tend_copy)

    @test all(isapprox.(to_vec(fd_vjp[1])[1], to_vec(dtend)[1], rtol=1e-3, atol=1e-3))

    # 
    # single variable leapfrog step 
    # 

    A_old = progn.vor[:,:,1]
    A_old_copy = deepcopy(A_old)
    dA_old = one(A_old)

    A_new = progn.vor[:,:,2]
    A_new_copy = deepcopy(A_new)
    dA_new = one(A_new)

    tendency = diagn.tendencies.vor_tend
    tendency_copy = deepcopy(tendency)
    dtendency = make_zero(tendency)

    L = model.time_stepping

    @info "Running reverse-mode AD"
    @time autodiff(Reverse, SpeedyWeather.leapfrog!, Const, Duplicated(A_old, dA_old), Duplicated(A_new, dA_new), Duplicated(tendency, dtendency), Const(dt), Const(lf1), Const(L))

    w1 = L.robert_filter*L.williams_filter/2   
    w2 = L.robert_filter*(1-L.williams_filter)/2   
    @test all(dtendency .≈ dt*(1+w1-w2))
    # ∂(tend) needs to be: dt* ( 1 + w1 - w2) (for every coefficient)
end

@testset "Differentiability: timestep! on Barotropic model" begin
    # T15 still yields somewhat sensible dynamics, that's why it's chosen here
    spectral_grid = SpectralGrid(trunc=9, nlayers=1)
    model = BarotropicModel(; spectral_grid)
    simulation = initialize_with_spinup!(model)
    (; Δt, Δt_millisec) = simulation.model.time_stepping
    dt = 2Δt

    fdsim = ADSimulation(simulation)
    @info "Running finite differences"
    # for the full timestep, we need a bit higher precision 
    fd_vjp = @time FiniteDifferences.j′vp(
        central_fdm(15,1),
        x -> timestep_oop(x, deepcopy(fdsim.diagvars), dt, deepcopy(fdsim.model)),
        one(fdsim.progvars),
        deepcopy(fdsim.progvars)
    )

    adsim = ADSimulation(simulation)
    progn_new, dprogn_new = prognosticseed(adsim)
    
    @info "Running reverse-mode AD"
    # test that we can differentiate wrt to everything
    @time autodiff(
        Reverse,
        timestep_oop!,
        Const,
        Duplicated(progn_new, dprogn_new),
        Duplicated(adsim.progvars, adsim.dprogvars),
        Duplicated(adsim.diagvars, adsim.ddiagvars),
        Const(dt),
        Duplicated(model, make_zero(model))
    )

    # nonzero gradient
    dprogn = adsim.dprogvars
    @test sum(to_vec(dprogn)[1]) != 0

    @test_broken isapprox(to_vec(fd_vjp[1])[1], to_vec(dprogn)[1], rtol=0.05) # we have to go really quite high with the tolerances here
    @test mean(abs.(to_vec(fd_vjp[1])[1] - to_vec(dprogn)[1])) < 0.002 # so we check a few extra statistics
    @test maximum(to_vec(fd_vjp[1].vor)[1] - to_vec(dprogn.vor)[1]) < 0.05
    
    # test that we can differentiante with Const(Model) only wrt to the state
    adsim2 = ADSimulation(simulation)
    progn_new, dprogn_new = prognosticseed(adsim2)

    @time autodiff(
        set_runtime_activity(Reverse),
        timestep_oop!,
        Const,
        Duplicated(progn_new, dprogn_new),
        Duplicated(adsim2.progvars, adsim2.dprogvars),
        Duplicated(adsim2.diagvars, adsim2.ddiagvars),
        Const(dt),
        Const(model)
    )

    # use the same FD comparision 

    @test_broken isapprox(to_vec(fd_vjp[1])[1], to_vec(d_progn)[1], rtol=0.05) # we have to go really quite high with the tolerances here
    @test mean(abs.(to_vec(fd_vjp[1])[1] - to_vec(d_progn)[1])) < 0.002 # so we check a few extra statistics
    @test maximum(to_vec(fd_vjp[1].vor)[1] - to_vec(d_progn.vor)[1]) < 0.05
end

# finite differences on transform! has issues with numerical accuracy
# @testset "VJP: transform!" begin
#     # T15 still yields somewhat sensible dynamics, that's why it's chosen here
#     spectral_grid = SpectralGrid(trunc=9, nlayers=1)
#     model = BarotropicModel(; spectral_grid)
#     simulation = initialize_with_spinup!(model)

#     diag_copy = deepcopy(diagn)

#     ddiag = one(diagn)
#     ddiag_copy = deepcopy(ddiag)

#     progn_copy = deepcopy(progn)
#     dprogn = make_zero(progn)

#     autodiff(Reverse, SpeedyWeather.transform!, Const, Duplicated(diagn, ddiag), Duplicated(progn, dprogn), Const(lf2), Const(model))
#     autodiff(Reverse, SpeedyWeather.transform!, Const, Duplicated(diagn, ddiag), Duplicated(progn, dprogn), Const(lf2), Duplicated(model, make_zero(model)))

#     function transform_diagn(diag, progn, lf2, model)
#         diag_copy = deepcopy(diag)
#         transform!(diag_copy, progn, lf2, model)
#         return diag_copy
#     end 

#     fd_vjp = FiniteDifferences.j′vp(central_fdm(12,1), x -> transform_diagn(diag_copy, x, lf2, model), ddiag_copy, progn_copy)

#     @test all(isapprox.(to_vec(fd_vjp[1])[1], to_vec(dprogn)[1],rtol=1e-2,atol=1e-2))
# end

@testset "Differentiability: Barotropic model parameters" begin
    # T15 still yields somewhat sensible dynamics, that's why it's chosen here
    spectral_grid = SpectralGrid(trunc=9, nlayers=1)          # define resolution
    model = BarotropicModel(; spectral_grid)   # construct model
    simulation = initialize_with_spinup!(model)
    (; Δt, Δt_sec) = simulation.model.time_stepping
    dt = Δt
    ps = parameters(model)
    pvec = vec(ps)
    adsim = ADSimulation(simulation)
    # timestep_oop!(deepcopy(adsim.progvars), adsim.progvars, adsim.diagvars, dt, adsim.model, pvec)
    dp = zero(pvec)
    progn_new, dprogn_new = deepcopy(adsim.progvars), one(adsim.dprogvars)
    @time autodiff(
        Reverse,
        timestep_oop!,
        Const,
        Duplicated(progn_new, dprogn_new),
        Duplicated(adsim.progvars, adsim.dprogvars),
        Duplicated(adsim.diagvars, adsim.ddiagvars),
        Const(dt),
        Duplicated(adsim.model, make_zero(adsim.model)),
        Duplicated(pvec, dp),
    )

    fdsim = ADSimulation(simulation)
    fd_vjp = @time FiniteDifferences.j′vp(
        central_fdm(10,1),
        x -> timestep_oop(deepcopy(fdsim.progvars), deepcopy(fdsim.diagvars), dt, deepcopy(fdsim.model), x),
        one(fdsim.progvars),
        copy(pvec),
    )
    @test all(isapprox.(dp, fd_vjp[1], atol=1e-5, rtol=1e-3))
end
