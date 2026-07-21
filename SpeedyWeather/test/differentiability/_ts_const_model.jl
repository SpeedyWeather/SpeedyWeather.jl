using SpeedyWeather, Enzyme, Test
println("Julia ", VERSION)
spectral_grid = SpectralGrid(trunc = 5, nlayers = 1)
model = PrimitiveWetModel(; spectral_grid)
simulation = initialize!(model); initialize!(simulation)
run!(simulation, period = Hour(6))
(; variables, model) = simulation
vars = variables; dvars = make_zero(vars)
dvars.prognostic.vorticity .= 1 + im
dvars.prognostic.divergence .= 1 + im
dvars.prognostic.humidity .= 1 + im
dvars.prognostic.temperature .= 1 + im
dvars.prognostic.pressure .= 1 + im
println(">>> entering autodiff time_step! with Const(model) [STATE AD]"); flush(stdout)
autodiff(set_runtime_activity(Reverse), SpeedyWeather.time_step!, Const,
    Duplicated(vars, dvars), Const(model.time_stepping), Const(model))
println(">>> autodiff returned; sum(dvars)!=0 : ", sum(to_vec(dvars)[1]) != 0)
