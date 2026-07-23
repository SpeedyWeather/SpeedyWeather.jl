# full_diff_CI content but with maxtypeoffset!(4096) set explicitly, to test whether the offset
# fixes the 1.12 parameter-AD (Duplicated(model)) failure/churn of the full time_step!.
# Run: julia --project=SpeedyWeather/test/differentiability --check-bounds=yes <file>
using SpeedyWeather, Enzyme, FiniteDifferences, Test

println("Julia ", VERSION, ", Enzyme ", pkgversion(Enzyme))
offset = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 4096
Enzyme.API.maxtypeoffset!(offset)
println("maxtypeoffset set to ", offset)

@testset "Complete Differentiability (offset $offset)" begin
    spectral_grid = SpectralGrid(trunc = 5, nlayers = 1)
    model = PrimitiveWetModel(; spectral_grid)
    simulation = initialize!(model)
    initialize!(simulation)
    run!(simulation, period = Hour(6))

    (; variables, model) = simulation
    vars = variables
    dvars = make_zero(vars)
    dmodel = make_zero(model)

    dvars.prognostic.vorticity .= 1 + im
    dvars.prognostic.divergence .= 1 + im
    dvars.prognostic.humidity .= 1 + im
    dvars.prognostic.temperature .= 1 + im
    dvars.prognostic.pressure .= 1 + im

    println(">>> entering autodiff time_step! (Duplicated model)")
    flush(stdout)
    autodiff(
        set_runtime_activity(Reverse), SpeedyWeather.time_step!, Const,
        Duplicated(vars, dvars),
        Duplicated(model.time_stepping, dmodel.time_stepping),
        Duplicated(model, dmodel),
    )
    println(">>> autodiff returned")
    @test sum(to_vec(dvars)[1]) != 0
end
