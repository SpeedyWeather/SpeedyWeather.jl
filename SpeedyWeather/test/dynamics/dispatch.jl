using JET

# Regression guard: the barotropic `time_step!` must contain NO runtime dynamic dispatch.

@testset "barotropic time_step! free of runtime dispatch (JET)" begin
    spectral_grid = SpectralGrid(truncation = 6, nlayers = 1)
    model = BarotropicModel(; spectral_grid)
    model.feedback.verbose = false
    simulation = initialize!(model)
    (; variables, model) = simulation

    initialize!(simulation)
    # `ignored_modules = (Base,)` filters dispatches inside Base itself (e.g. generic show/printing
    # machinery) that this package cannot fix; everything reachable in SpeedyWeather & friends counts.
    @test_opt ignored_modules = (Base,) SpeedyWeather.time_step!(variables, model.time_stepping, model)
end

@testset "PrimitiveWet time_step! free of runtime dispatch (JET)" begin
    spectral_grid = SpectralGrid(truncation = 6, nlayers = 8)
    model = PrimitiveWetModel(; spectral_grid)
    model.feedback.verbose = false
    simulation = initialize!(model)
    (; variables, model) = simulation
    initialize!(simulation)

    @test_opt ignored_modules = (Base,) SpeedyWeather.time_step!(variables, model.time_stepping, model)
end

@testset "StochasticLeapfrog time_step! free of runtime dispatch (JET)" begin
    # Barotropic: no ocean and land pattern; PrimitiveWet: noise on variables with and without
    # layers, ocean (sea surface temperature, default amplitude) and land (soil temperature)
    for (Model, nlayers, amplitudes) in (
            (BarotropicModel, 1, (; vorticity = 1.0e-7)),
            (PrimitiveWetModel, 8, (; pressure = 1.0e-5, soil_temperature = 1.0e-4)),
        )
        spectral_grid = SpectralGrid(; truncation = 6, nlayers)
        time_stepping = StochasticLeapfrog(spectral_grid; seed = 1, amplitudes...)
        model = Model(spectral_grid; time_stepping)
        model.feedback.verbose = false
        simulation = initialize!(model)
        (; variables, model) = simulation
        initialize!(simulation)

        @test_opt ignored_modules = (Base,) SpeedyWeather.time_step!(variables, model.time_stepping, model)
    end
end
