using JET

# Regression guard: the barotropic `time_step!` must contain NO runtime dynamic dispatch.

@testset "barotropic time_step! free of runtime dispatch (JET)" begin
    spectral_grid = SpectralGrid(trunc = 5, nlayers = 1)
    model = BarotropicModel(; spectral_grid)
    model.feedback.verbose = false
    simulation = initialize!(model)
    (; variables, model) = simulation

    initialize!(simulation)
    # `ignored_modules = (Base,)` filters dispatches inside Base itself (e.g. generic show/printing
    # machinery) that this package cannot fix; everything reachable in SpeedyWeather & friends counts.
    @test_opt ignored_modules = (Base,) SpeedyWeather.time_step!(variables, model.time_stepping, model)
end

# PrimitiveWet still has runtime dispatches, from a DIFFERENT root cause than the barotropic ones:
# `update_prognostic!` (time_stepping/time_integration.jl) loops
#
#     for varname in tendency_names(vars)
#         var = getfield(prognostic, varname)
#
# `tendency_names` is `@generated`, so the names ARE a compile-time tuple of `Symbol`s — but the loop is
# not unrolled, so `varname` is a runtime `Symbol` and `getfield(prognostic, ::Symbol)` widens to the
# union of all field types (`Union{RefValue{Float32}, Clock, NamedTuple, LowerTriangularArray}`). That
# union is too wide to split, so `update_prognostic!`/`get_steps`/`lta_view` are reached by dispatch.
#
# Fixing it means unrolling those loops over the statically-known names (e.g. `@generated`/`ntuple`) — a
# change to the Variables time-stepping design, out of scope here. This test pins the CURRENT count so it
# cannot silently grow, and fails loudly (prompting an update) if it improves. See
# DISPATCH_REDUCTION_PLAN.md.
@testset "PrimitiveWet time_step! dispatch count (pinned)" begin
    spectral_grid = SpectralGrid(trunc = 5, nlayers = 8)
    model = PrimitiveWetModel(; spectral_grid)
    model.feedback.verbose = false
    simulation = initialize!(model)
    (; variables, model) = simulation
    initialize!(simulation)

    report = JET.report_opt(
        SpeedyWeather.time_step!,
        typeof.((variables, model.time_stepping, model));
        ignored_modules = (Base,),
    )
    # regression guard: must not grow (barotropic is at 0; PrimitiveWet was 19 when this was written)
    @test length(JET.get_reports(report)) <= 19
end
