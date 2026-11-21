using Terrarium: Terrarium

@testset "Differentiability: Terrarium land coupling" begin
    # Regression test for the SpeedyWeather-Terrarium coupling (SpeedyWeatherTerrariumExt): one
    # reverse-mode time step of a PrimitiveWetModel with a Terrarium land model must compile and give
    # a finite, nonzero gradient that reaches the Terrarium state. Kept as small as possible (lowest
    # resolution, one atmospheric layer, two soil layers, all land, no vegetation), since the cost is
    # dominated by compilation rather than by the size of the model.
    spectral_grid = SpectralGrid(truncation = 6, nlayers = 1)
    column_grid = Terrarium.ColumnRingGrid(
        Terrarium.CPU(), Float32, Terrarium.UniformSpacing(Δz = 0.5, N = 2), spectral_grid.grid,
    )
    terrarium_model = Terrarium.LandModel(column_grid; vegetation = nothing)
    land = SpeedyWeather.LandModel(spectral_grid, terrarium_model)
    model = PrimitiveWetModel(
        spectral_grid;
        land,
        land_sea_mask = RockyPlanetMask(spectral_grid),
        surface_heat_flux = SurfaceHeatFlux(spectral_grid, land = PrescribedLandHeatFlux()),
        surface_humidity_flux = SurfaceHumidityFlux(spectral_grid, land = PrescribedLandHumidityFlux()),
    )
    simulation = initialize!(model)
    run!(simulation, period = Hour(1))

    (; variables, model) = simulation
    vars = variables
    dvars = make_zero(vars)
    dmodel = make_zero(model)

    # seed the atmospheric temperature and humidity after the step
    dvars.prognostic.temperature .= 1 + im
    dvars.prognostic.humidity .= 1 + im

    autodiff(
        set_runtime_activity(Reverse), SpeedyWeather.time_step!, Const,
        Duplicated(vars, dvars),
        Duplicated(model.time_stepping, dmodel.time_stepping),
        Duplicated(model, dmodel),
    )

    @test all(isfinite, to_vec(dvars)[1])
    @test sum(to_vec(dvars)[1]) != 0
    # the gradient must reach the land: the soil temperature of the Terrarium state
    dstate = dvars.prognostic.land.terrarium.state
    @test all(isfinite, Terrarium.interior(dstate.temperature))
    @test any(!iszero, Terrarium.interior(dstate.temperature))
end
