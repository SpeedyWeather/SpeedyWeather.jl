@testset "SPPT" begin
    spectral_grid = SpectralGrid(trunc=31, nlayers=8)

    @testset for S in (Nothing, StochasticallyPerturbedParameterizationTendencies)
        s = S(spectral_grid)
        random_process = SpectralAR1Process(spectral_grid, seed=1)

        # construct model
        model = PrimitiveWetModel(spectral_grid; random_process, stochastic_physics=s)

        initialize!(model.stochastic_physics, model)

        progn = PrognosticVariables(model)
        diagn = DiagnosticVariables(model)

        SpeedyWeather.random_process!(progn, model.random_process)
        r_grid = diagn.grid.random_pattern
        r_spec = progn.random_pattern
        transform!(r_grid, r_spec, diagn.dynamics.scratch_memory, model.spectral_transform)

        # execute only one parameterization
        ij = rand(1:model.spectral_grid.npoints)
        SpeedyWeather.parameterization!(ij, diagn, progn, model.longwave_radiation, model)
        
        dTdt = diagn.tendencies.temp_tend_grid[ij, end]
        @test dTdt != 0
        
        # now stochastic perturbation and check it's not the same
        SpeedyWeather.parameterization!(ij, diagn, progn, model.stochastic_physics, model)

        dTdt2 = diagn.tendencies.temp_tend_grid[ij, end]
        @test dTdt2 != dTdt != 0
    end
end