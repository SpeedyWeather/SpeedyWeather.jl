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

@testset "time_step! with fused tracers free of runtime dispatch (JET)" begin
    for (Model, nlayers) in ((ShallowWaterModel, 1), (PrimitiveWetModel, 8))
        spectral_grid = SpectralGrid(truncation = 6; nlayers)
        tracers = (; a = Tracer(), b = Tracer(), c = Tracer())
        model = Model(; spectral_grid, tracers, fuse_tracers = true)
        model.feedback.verbose = false
        simulation = initialize!(model)
        (; variables, model) = simulation
        initialize!(simulation)
        @test SpeedyWeather.fused_tracers(variables, model)

        # all tracers processed together, through their fused parents
        @test_opt ignored_modules = (Base,) SpeedyWeather.time_step!(variables, model.time_stepping, model)
    end
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
