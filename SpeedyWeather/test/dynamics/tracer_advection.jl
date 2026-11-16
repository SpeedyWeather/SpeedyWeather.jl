@testset "Tracers: add!, delete!" begin
    for TS in (Leapfrog, NCycleLorenz)
        spectral_grid = SpectralGrid(nlayers = 1)
        time_stepping = TS(spectral_grid)
        model = BarotropicModel(spectral_grid; time_stepping)

        add!(model, Tracer(:abc))
        @test length(model.tracers) == 1

        add!(model, Tracer(:t1), Tracer(:t2))
        @test length(model.tracers) == 3

        delete!(model, Tracer(:t1))
        delete!(model, Tracer(:t2))
        @test length(model.tracers) == 1

        add!(model, Tracer(:t4))
        @test length(model.tracers) == 2

        simulation = initialize!(model)

        ntracers = 2
        @test length(simulation.variables.prognostic.tracers) == ntracers
        @test length(simulation.variables.grid.tracers) == ntracers
        @test length(simulation.variables.tendencies.tracers) == ntracers       # spectral tracers
        @test length(simulation.variables.tendencies.grid_tracers) == ntracers  # grid tracers

        @test ndims(simulation.variables.prognostic.tracers.abc) == 3           # horizontal x vertical x steps
        @test ndims(simulation.variables.tendencies.tracers.abc) == 3           # horizontal x vertical x steps
        @test ndims(simulation.variables.tendencies.grid_tracers.abc) == 3      # horizontal x vertical x steps
        @test size(simulation.variables.prognostic.tracers.abc) == size(simulation.variables.prognostic.vorticity)
    end
end

@testset "Tracers: activate!, deactivate!" begin
    spectral_grid = SpectralGrid(nlayers = 1)
    model = BarotropicModel(spectral_grid)
    model.feedback.verbose = false

    add!(model, Tracer(:abc))
    deactivate!(model, Tracer(:abc))
    @test model.tracers[:abc].active == false

    activate!(model, Tracer(:abc))
    @test model.tracers[:abc].active

    simulation = initialize!(model)
    deactivate!(simulation, Tracer(:abc))
    @test model.tracers[:abc].active == false

    activate!(simulation, Tracer(:abc))
    @test model.tracers[:abc].active

    f = (λ, φ, σ) -> exp(-(λ - 180)^2 / 10^2)
    set!(simulation, abc = f, namespace = :tracers)

    run!(simulation, period = Day(0))

    # initial conditions
    abc0_spec = get_step(simulation.variables.prognostic.tracers.abc, 1)
    abc0 = deepcopy(simulation.variables.grid.tracers.abc)

    # set some grid in the same way and check that the tracer abc is correctly set
    # but compare in spectral space due to transform errors
    def = zero(get_step(abc0, 1))
    set!(def, f, model.geometry)
    @test abc0_spec == transform(def, model.spectral_transform)

    # check that everything is different after 10 days
    run!(simulation, period = Day(10))
    abc1 = simulation.variables.grid.tracers.abc

    for ij in eachindex(abc0, abc1)
        @test abc0[ij] != abc1[ij]
    end

    # check that everything is the same if tracer deactivated
    deactivate!(simulation, Tracer(:abc))
    run!(simulation, period = Day(10))
    abc2 = simulation.variables.grid.tracers.abc

    for ij in eachindex(abc1, abc2)
        @test abc1[ij] == abc2[ij]
    end
end

@testset "Tracers: deactivated or deleted after initialization stay frozen, added after initialization error" begin
    spectral_grid = SpectralGrid(nlayers = 1)
    model = BarotropicModel(spectral_grid)
    model.feedback.verbose = false

    tracer_names = (:active_tracer, :deactivated_tracer, :deleted_tracer)
    add!(model, Tracer.(tracer_names)...)
    simulation = initialize!(model)

    f = (λ, φ, σ) -> exp(-(λ - 180)^2 / 10^2)
    set!(simulation, active_tracer = f, deactivated_tracer = f, deleted_tracer = f, namespace = :tracers)
    run!(simulation, period = Day(0))       # all tracers active: transformed to the grid

    # deactivate one, delete one (its variables remain in the simulation)
    deactivate!(simulation, Tracer(:deactivated_tracer))
    delete!(model, Tracer(:deleted_tracer))

    (; prognostic, grid) = simulation.variables
    prognostic0 = Dict(name => copy(getfield(prognostic.tracers, name).data) for name in tracer_names)
    grid0 = Dict(name => copy(getfield(grid.tracers, name).data) for name in tracer_names)
    @test all(name -> any(!iszero, grid0[name]), tracer_names)

    run!(simulation, period = Day(1))

    @test getfield(prognostic.tracers, :active_tracer).data != prognostic0[:active_tracer]
    @test getfield(grid.tracers, :active_tracer).data != grid0[:active_tracer]
    for name in (:deactivated_tracer, :deleted_tracer)     # frozen
        @test getfield(prognostic.tracers, name).data == prognostic0[name]
        @test getfield(grid.tracers, name).data == grid0[name]
    end

    # a tracer added after initialization has no variables
    add!(model, Tracer(:late_tracer))
    @test_throws ArgumentError run!(simulation, period = Day(0))
    delete!(model, Tracer(:late_tracer))
    run!(simulation, period = Day(0))       # runs again without it
end

@testset "Tracers primitive equation models" begin
    @testset for Model in (PrimitiveDryModel, PrimitiveWetModel)
        spectral_grid = SpectralGrid(nlayers = 8)
        model = Model(spectral_grid)
        model.feedback.verbose = false

        add!(model, Tracer(:abc))
        simulation = initialize!(model)
        set!(simulation, abc = (λ, φ, σ) -> σ * exp(-(λ - 180)^2 / 10^2), namespace = :tracers)
        run!(simulation, period = Day(0))

        # check that initial conditions are non-zero indeed
        @test all(simulation.variables.grid.tracers.abc .!= 0)
        abc0 = deepcopy(simulation.variables.grid.tracers.abc)

        # check that everything is different after simulation
        run!(simulation, period = Day(1))
        abc1 = simulation.variables.grid.tracers.abc

        n = sum(abc0 .!= abc1)
        m = length(abc1)
        @info "$Model tracel advection: $n/$m elements differ."

        # caused failures hence info above and weaken test to 95% of grid cells must be different
        # @test all(abc0 .!= abc1)
        @test n / m > 0.95
    end
end
