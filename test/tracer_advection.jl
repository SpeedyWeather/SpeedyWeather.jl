@testset "Tracers: add!, delete!" begin
    spectral_grid = SpectralGrid(nlayers=1)
    model = BarotropicModel(spectral_grid)

    add!(model, Tracer(:abc))
    @test length(model.tracers) == 1

    add!(model, Tracer(:t1), Tracer(:t2))
    @test length(model.tracers) == 3

    delete!(model, Tracer(:t1))
    delete!(model, Tracer(:t2))
    @test length(model.tracers) == 1

    simulation = initialize!(model)

    add!(simulation, Tracer(:co2))
    @test length(model.tracers) == 2
    @test length(simulation.prognostic_variables.tracers) == 2

    delete!(simulation, Tracer(:abc))
    @test length(model.tracers) == 1
    @test length(simulation.prognostic_variables.tracers) == 1
    @test length(simulation.diagnostic_variables.grid.tracers_grid) == 1
    @test length(simulation.diagnostic_variables.grid.tracers_grid_prev) == 1
    @test length(simulation.diagnostic_variables.tendencies.tracers_tend) == 1
    @test length(simulation.diagnostic_variables.tendencies.tracers_tend_grid) == 1
end

@testset "Tracers: activate!, deactivate!" begin
    spectral_grid = SpectralGrid(nlayers=1)
    model = BarotropicModel(spectral_grid)

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

    set!(simulation, abc = (λ, φ, σ) -> exp(-(λ-180)^2/10^2))
    run!(simulation, period=Day(0))

    # initial conditions
    abc0 = simulation.diagnostic_variables.grid.tracers_grid[:abc][:, 1]
    
    # check that everything is different after 10 days
    run!(simulation, period=Day(10))
    abc1 = simulation.diagnostic_variables.grid.tracers_grid[:abc][:, 1]

    for ij in eachindex(abc0, abc1)
        @test abc0[ij] != abc1[ij]
    end

    # check that everything is the same if tracer deactivated
    deactivate!(simulation, Tracer(:abc))
    run!(simulation, period=Day(10))
    abc2 = simulation.diagnostic_variables.grid.tracers_grid[:abc][:, 1]

    for ij in eachindex(abc1, abc2)
        @test abc1[ij] == abc2[ij]
    end
end

@testset "Tracers primitive equation models" begin
    @testset for Model in (PrimitiveDryModel, PrimitiveWetModel)
        spectral_grid = SpectralGrid(nlayers=8)
        model = Model(spectral_grid)

        add!(model, Tracer(:abc))
        simulation = initialize!(model)
        set!(simulation, abc = (λ, φ, σ) -> σ*exp(-(λ-180)^2/10^2))
        run!(simulation, period=Day(0))
        abc0 = simulation.diagnostic_variables.grid.tracers_grid[:abc][:, :]

        # check that everything is different after simulation
        run!(simulation, period=Day(1))
        abc1 = simulation.diagnostic_variables.grid.tracers_grid[:abc][:, :]

        for ij in eachindex(abc0, abc1)
            @test abc0[ij] != abc1[ij]
        end
    end
end