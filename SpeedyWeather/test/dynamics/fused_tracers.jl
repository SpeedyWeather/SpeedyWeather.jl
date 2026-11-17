using Random

# two simulations of `Model` with `ntracers` tracers from the same initial conditions, one processing the tracers
# together (fused), the other one by one. Same types for both (`fuse_tracers` is a runtime switch), compiled once
function fused_and_per_tracer_simulations(Model; nlayers, ntracers, truncation = 16)
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

# (model, nlayers, number of tracers). The largest transform batch is 9nlayers + 1 layers: with nlayers = 1 (10) the
# tendencies and fluxes of 4 tracers (12 layers) and with nlayers = 2 (19) all transforms of 11 tracers (22 and 66
# layers) are split into chunks of whole tracers
const FUSED_TRACER_CASES = (
    (BarotropicModel, 1, 3), (ShallowWaterModel, 1, 4), (PrimitiveDryModel, 2, 11), (PrimitiveWetModel, 4, 3),
)

@testset "Fused tracers: same as one by one" begin
    for (Model, nlayers, ntracers) in FUSED_TRACER_CASES
        fused, per_tracer, names = fused_and_per_tracer_simulations(Model; nlayers, ntracers)
        @test SpeedyWeather.fused_tracers(fused.variables, fused.model)
        @test !SpeedyWeather.fused_tracers(per_tracer.variables, per_tracer.model)

        # the tracers are stacked in the fused parents, the variables of each tracer are views of them
        fused_parents = fused.variables.fused
        @test keys(fused_parents.tracers_spectral.slot_map) == names
        @test keys(fused_parents.tracer_tendencies_grid.slot_map) ==
            (names..., map(name -> Symbol(name, :_u_flux), names)..., map(name -> Symbol(name, :_v_flux), names)...)
        Model in (ShallowWaterModel, PrimitiveDryModel) && @test 3ntracers * nlayers > fused.model.spectral_transform.nlayers

        tracer0 = copy(fused.variables.prognostic.tracers.tracer1.data)
        run!(fused, steps = 6)
        run!(per_tracer, steps = 6)

        for name in names, group in (:prognostic, :grid)
            tracer_fused = getfield(getfield(fused.variables, group).tracers, name)
            tracer_per_tracer = getfield(getfield(per_tracer.variables, group).tracers, name)
            @test tracer_fused == tracer_per_tracer     # same chunks of the transforms on CPU, so bitwise
        end
        @test fused.variables.prognostic.tracers.tracer1.data != tracer0     # advected
        @test fused.variables.prognostic.vorticity == per_tracer.variables.prognostic.vorticity
    end
end

@testset "Fused tracers: a deactivated tracer switches to one by one" begin
    for (Model, nlayers, ntracers) in FUSED_TRACER_CASES[2:2:4]     # ShallowWater, PrimitiveWet
        fused, _, names = fused_and_per_tracer_simulations(Model; nlayers, ntracers)

        run!(fused, steps = 2)              # all tracers together
        deactivate!(fused, :tracer2)
        @test !SpeedyWeather.fused_tracers(fused.variables, fused.model)

        # the 1st run! after deactivating a tracer still copies its 1st leapfrog step to the 2nd (with or
        # without fused tracers), it is frozen from then on
        run!(fused, steps = 1)

        (; prognostic, grid) = fused.variables
        prognostic0 = Dict(name => copy(getfield(prognostic.tracers, name).data) for name in names)
        grid0 = Dict(name => copy(getfield(grid.tracers, name).data) for name in names)
        run!(fused, steps = 4)

        @test prognostic.tracers.tracer2.data == prognostic0[:tracer2]   # frozen
        @test grid.tracers.tracer2.data == grid0[:tracer2]
        for name in (:tracer1, :tracer3)
            @test getfield(prognostic.tracers, name).data != prognostic0[name]
            @test getfield(grid.tracers, name).data != grid0[name]
        end

        # all active again: processed together again
        activate!(fused, :tracer2)
        @test SpeedyWeather.fused_tracers(fused.variables, fused.model)
        run!(fused, steps = 2)
        @test prognostic.tracers.tracer2.data != prognostic0[:tracer2]
    end
end

@testset "Fused tracers: stacked kernels" begin
    Model, nlayers, ntracers = FUSED_TRACER_CASES[3]     # PrimitiveDry
    simulation, _, _ = fused_and_per_tracer_simulations(Model; nlayers, ntracers)
    run!(simulation, steps = 2)      # some vertical velocity
    (; variables, model) = simulation
    (; spectral_grid) = model
    L = nlayers
    NL = ntracers * L

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
        stacked = copy(tendencies.data[:, 1:NL, :])

        fill!(tendencies, 0)
        SpeedyWeather.foreach_active_tracer(SpeedyWeather.vertical_advection_tracer!, variables, model, w, Δσ, advection_scheme)
        @test stacked == tendencies.data[:, 1:NL, :]
        @test any(!iszero, stacked)
    end

    # horizontal diffusion, explicit and implicit
    diffusion = model.horizontal_diffusion
    for impl in (nothing, diffusion.impl)
        tracers = get_step(parent(variables.fused.tracers_spectral), 1)
        tendencies = get_step(parent(variables.fused.tracer_tendencies_spectral), 1)
        tendencies_tracers = SpeedyWeather.wrapped_view(tendencies, :, 1:NL)
        rand!(tracers.data)
        rand!(tendencies.data)
        expected = copy(tendencies_tracers)
        for n in 1:ntracers
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
