# Enzyme and Julia 1.11 still has some problems, and the test below is broken
# in Julia 1.11
import Pkg
Pkg.activate(@__DIR__)
using SpeedyWeather, Enzyme, FiniteDifferences, Test

@testset "Complete Differentiability" begin
    # We do extensive correctness checks and tests on the differentiability
    # in a seperate test set. But we do want to ensure in the regular CI that
    # we don't commit some kind of problem for the Enzyme differentiability
    # so, we test here if we get a non-zero gradient from the timestepping.
    spectral_grid = SpectralGrid(trunc = 5, nlayers = 1)          # define resolution
    model = PrimitiveWetModel(; spectral_grid)   # construct model
    simulation = initialize!(model)
    initialize!(simulation)
    run!(simulation, period = Hour(6))

    (; variables, model) = simulation

    vars = variables
    dvars = make_zero(vars)
    dmodel = make_zero(model)

    # set a seeed for the prognostic variables
    dvars.prognostic.vorticity .= 1 + im
    dvars.prognostic.divergence .= 1 + im
    dvars.prognostic.humidity .= 1 + im
    dvars.prognostic.temperature .= 1 + im
    dvars.prognostic.pressure .= 1 + im

    # differentiate time_step!(vars, time_stepping, model), the inner time step
    # without clock/output/feedback; pass the time_stepping of model (and its shadow)
    # explicitly to keep the aliasing with the model (and its shadow) consistent
    autodiff(
        set_runtime_activity(Reverse), SpeedyWeather.time_step!, Const,
        Duplicated(vars, dvars),
        Duplicated(model.time_stepping, dmodel.time_stepping),
        Duplicated(model, dmodel),
    )

    @test sum(to_vec(dvars)[1]) != 0
end
