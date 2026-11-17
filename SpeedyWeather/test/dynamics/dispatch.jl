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

@testset "barotropic time_step! with a tracer free of runtime dispatch (JET)" begin
    spectral_grid = SpectralGrid(truncation = 6, nlayers = 1)
    model = BarotropicModel(; spectral_grid, tracers = (; abc = Tracer()))
    model.feedback.verbose = false
    simulation = initialize!(model)
    (; variables, model) = simulation
    initialize!(simulation)

    # the tracer loops are unrolled over the tracer names in the type of `variables`
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
