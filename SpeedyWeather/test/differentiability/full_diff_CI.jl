# Enzyme and Julia 1.11 still has some problems, and the test below is broken
# in Julia 1.11
import Pkg
Pkg.activate(@__DIR__)
using SpeedyWeather, Enzyme, FiniteDifferences, Test

if VERSION <= v"1.11.0"
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
        (; Δt, Δt_millisec) = model.time_stepping
        dt = 2Δt
        lf1 = 1
        lf2 = 2

        vars = variables
        dvars = make_zero(vars)

        # set a seeed for the prognostic variables
        dvars.prognostic.vor .= 1 + im
        dvars.prognostic.div .= 1 + im
        dvars.prognostic.humid .= 1 + im
        dvars.prognostic.temp .= 1 + im
        dvars.prognostic.pres .= 1 + im

        dmodel = make_zero(model)

        autodiff(Reverse, SpeedyWeather.timestep!, Const, Duplicated(vars, dvars), Const(dt), Duplicated(model, dmodel), Const(lf1), Const(lf2))
        @test sum(to_vec(dvars)[1]) != 0

    end
else
    @testset "Complete Differentiability" begin
        @test_broken false
    end
end
