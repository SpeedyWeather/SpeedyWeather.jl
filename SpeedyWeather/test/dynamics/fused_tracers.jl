using Random

# two simulations of `Model` with `ntracers` tracers from the same initial conditions, one processing the tracers
# together (fused), the other one by one
function fused_and_per_tracer_simulations(Model; nlayers, ntracers = 3, truncation = 16)
    spectral_grid = SpectralGrid(; truncation, nlayers)
    names = ntuple(i -> Symbol(:tracer, i), ntracers)
    simulations = map((true, false)) do fuse_tracers
        tracers = NamedTuple{names}(ntuple(_ -> Tracer(), ntracers))
        model = Model(spectral_grid; tracers, fuse_tracers)
        model.feedback.verbose = false
        return initialize!(model)
    end
    fused, per_tracer = simulations

    # a different initial condition per tracer, varying in all directions
    set!(
        fused; namespace = :tracers,
        (name => ((λ, φ, σ) -> exp(-(λ - 36i)^2 / 30^2) * (1 + i * σ) + cosd(φ)) for (i, name) in enumerate(names))...
    )
    copy!(per_tracer.variables, fused.variables)
    return fused, per_tracer, names
end

@testset "Fused tracers: same as one by one" begin
    for (Model, nlayers) in ((BarotropicModel, 1), (ShallowWaterModel, 1), (PrimitiveDryModel, 4), (PrimitiveWetModel, 4))
        fused, per_tracer, names = fused_and_per_tracer_simulations(Model; nlayers)
        @test SpeedyWeather.fused_tracers(fused.variables, fused.model)
        @test !SpeedyWeather.fused_tracers(per_tracer.variables, per_tracer.model)

        # the tracers are stacked in the fused parents, the variables of each tracer are views of them
        fused_parents = fused.variables.fused
        @test keys(fused_parents.tracers_spectral.slot_map) == names
        @test keys(fused_parents.tracer_tendencies_grid.slot_map) ==
            (names..., map(name -> Symbol(name, :_u_flux), names)..., map(name -> Symbol(name, :_v_flux), names)...)

        tracer0 = copy(fused.variables.prognostic.tracers.tracer1.data)
        run!(fused, steps = 10)
        run!(per_tracer, steps = 10)

        for name in names
            for group in (:prognostic, :grid)
                tracer_fused = getfield(getfield(fused.variables, group).tracers, name)
                tracer_per_tracer = getfield(getfield(per_tracer.variables, group).tracers, name)
                @test tracer_fused == tracer_per_tracer     # same chunks of the transforms on CPU, so bitwise
            end
        end
        @test fused.variables.prognostic.tracers.tracer1.data != tracer0     # advected
        @test fused.variables.prognostic.vorticity == per_tracer.variables.prognostic.vorticity
    end
end

@testset "Fused tracers: more tracers than fit into one transform" begin
    # with nlayers = 1 the largest transform batch is 10 layers, i.e. the tendencies and fluxes of 4 tracers (12) and
    # with nlayers = 2 the tracers (22) and their tendencies and fluxes (66 layers) of 11 tracers are split into chunks
    for (Model, nlayers, ntracers) in ((ShallowWaterModel, 1, 4), (PrimitiveDryModel, 2, 11))
        fused, per_tracer, names = fused_and_per_tracer_simulations(Model; nlayers, ntracers)
        S = fused.model.spectral_transform
        @test 3 * ntracers * nlayers > S.nlayers

        run!(fused, steps = 4)
        run!(per_tracer, steps = 4)
        for name in names
            @test getfield(fused.variables.prognostic.tracers, name) ==
                getfield(per_tracer.variables.prognostic.tracers, name)
            @test getfield(fused.variables.grid.tracers, name) == getfield(per_tracer.variables.grid.tracers, name)
        end
    end
end

@testset "Fused tracers: a deactivated tracer switches to one by one" begin
    for (Model, nlayers) in ((ShallowWaterModel, 1), (PrimitiveWetModel, 4))
        fused, per_tracer, names = fused_and_per_tracer_simulations(Model; nlayers)
        run!(fused, steps = 2)

        deactivate!(fused, :tracer2)
        @test !SpeedyWeather.fused_tracers(fused.variables, fused.model)

        # the current step on the grid (the previous one of all tracers is moved back in primitive equation models)
        (; time_stepping, spectral_transform) = fused.model
        current_grid(name) = SpeedyWeather.get_prognostic_step(
            getfield(fused.variables.grid.tracers, name), time_stepping, spectral_transform, fused.model
        )
        (; prognostic) = fused.variables
        prognostic0 = Dict(name => copy(getfield(prognostic.tracers, name).data) for name in names)
        grid0 = Dict(name => copy(current_grid(name).data) for name in names)
        run!(fused, steps = 4)

        @test prognostic.tracers.tracer2.data == prognostic0[:tracer2]   # frozen
        @test current_grid(:tracer2).data == grid0[:tracer2]
        for name in (:tracer1, :tracer3)
            @test getfield(prognostic.tracers, name).data != prognostic0[name]
            @test current_grid(name).data != grid0[name]
        end

        activate!(fused, :tracer2)
        @test SpeedyWeather.fused_tracers(fused.variables, fused.model)
    end
end

@testset "Fused tracers: stacked kernels" begin
    spectral_grid = SpectralGrid(truncation = 16, nlayers = 4)
    model = PrimitiveDryModel(spectral_grid; tracers = (; a = Tracer(), b = Tracer()))
    simulation = initialize!(model)
    run!(simulation, steps = 2)      # some vertical velocity
    (; variables) = simulation
    L = spectral_grid.nlayers

    # vertical advection: the stacked kernel as the kernel for each tracer, the stencil doesn't cross tracers
    for advection_scheme in (
            CenteredVerticalAdvection(spectral_grid), UpwindVerticalAdvection(spectral_grid),
            WENOVerticalAdvection(spectral_grid),
        )
        tracers = parent(variables.fused.tracers_grid)
        rand!(tracers.data)
        tendencies = parent(variables.fused.tracer_tendencies_grid)
        fill!(tendencies, 0)
        (; w) = variables.dynamics
        Δσ = model.geometry.σ_levels_thick
        SpeedyWeather.vertical_advection_tracers!(variables, model, w, Δσ, advection_scheme)
        stacked = copy(tendencies.data[:, 1:(2L), :])

        fill!(tendencies, 0)
        SpeedyWeather.foreach_active_tracer(SpeedyWeather.vertical_advection_tracer!, variables, model, w, Δσ, advection_scheme)
        @test stacked == tendencies.data[:, 1:(2L), :]
        @test any(!iszero, stacked)
    end

    # horizontal diffusion, explicit and implicit
    diffusion = model.horizontal_diffusion
    for impl in (nothing, diffusion.impl)
        tracers = get_step(parent(variables.fused.tracers_spectral), 1)
        tendencies = get_step(parent(variables.fused.tracer_tendencies_spectral), 1)
        tendencies_tracers = SpeedyWeather.wrapped_view(tendencies, :, 1:(2L))
        rand!(tracers.data)
        rand!(tendencies.data)
        expected = copy(tendencies_tracers)
        for n in 1:2
            layers = ((n - 1) * L + 1):(n * L)
            SpeedyWeather.horizontal_diffusion!(
                SpeedyWeather.wrapped_view(expected, :, layers), SpeedyWeather.wrapped_view(tracers, :, layers),
                diffusion.expl, impl
            )
        end
        SpeedyWeather.horizontal_diffusion_stacked!(tendencies_tracers, tracers, diffusion.expl, impl, L)
        @test tendencies_tracers == expected
    end
end
