# full time_step! but dynamics_only=true (NO physics parameterizations) — isolate whether the
# 1.12 Enzyme compile hang is in the physics or the dyn-core time_step! chain.
using SpeedyWeather, Enzyme, Test
println("Julia ", VERSION)
spectral_grid = SpectralGrid(trunc = 5, nlayers = 1)
model = PrimitiveWetModel(; spectral_grid, dynamics_only = true)
simulation = initialize!(model); initialize!(simulation)
run!(simulation, period = Hour(6))
(; variables, model) = simulation
vars = variables; dvars = make_zero(vars)
dvars.prognostic.vorticity .= 1 + im
dvars.prognostic.divergence .= 1 + im
dvars.prognostic.humidity .= 1 + im
dvars.prognostic.temperature .= 1 + im
dvars.prognostic.pressure .= 1 + im
dmodel = make_zero(model)
println(">>> entering autodiff time_step! dynamics_only, Duplicated(model)"); flush(stdout)
autodiff(set_runtime_activity(Reverse), SpeedyWeather.time_step!, Const,
    Duplicated(vars, dvars), Duplicated(model.time_stepping, dmodel.time_stepping), Duplicated(model, dmodel))
println(">>> autodiff returned; nonzero: ", sum(to_vec(dvars)[1]) != 0)
