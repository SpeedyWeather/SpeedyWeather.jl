@testset "SPPT" begin
    spectral_grid = SpectralGrid(trunc = 31, nlayers = 8)

    @testset for S in (Nothing, StochasticallyPerturbedParameterizationTendencies)
        s = S(spectral_grid)
        random_process = SpectralAR1Process(spectral_grid, seed = 1)

        # construct model
        model = PrimitiveWetModel(spectral_grid; random_process, stochastic_physics = s)

        initialize!(model.stochastic_physics, model)

        vars = Variables(model)

        SpeedyWeather.random_process!(vars, model.random_process)

        # execute only one parameterization
        ij = rand(1:model.geometry.npoints)
        SpeedyWeather.parameterization!(ij, vars, model.longwave_radiation, model)

        dTdt = vars.tendencies.grid.temperature[ij, end]
        @test dTdt != 0

        # now stochastic perturbation and check it's not the same
        SpeedyWeather.parameterization!(ij, vars, model.stochastic_physics, model)

        dTdt2 = vars.tendencies.grid.temperature[ij, end]
        @test dTdt2 != dTdt != 0
    end
end
