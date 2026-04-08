using SpeedyWeather
using Dates

spectral_grid = SpectralGrid(trunc=63, nlayers=1)

time_stepping = LorenzNCycle(spectral_grid; N=4, version=2)
implicit      = ImplicitShallowWater(spectral_grid; α=0.5)

model = ShallowWaterModel(spectral_grid;
    time_stepping,
    implicit,
    initial_conditions = ZonalJet(spectral_grid),
)

simulation = initialize!(model)
run!(simulation, period=Day(6))

progn = simulation.prognostic_variables
vor_ok  = all(isfinite, SpeedyWeather.get_step(progn.vor,  1).data)
div_ok  = all(isfinite, SpeedyWeather.get_step(progn.div,  1).data)
pres_ok = all(isfinite, SpeedyWeather.get_step(progn.pres, 1).data)

if vor_ok && div_ok && pres_ok
    println("6-day ZonalJet integration: OK (all fields finite)")
else
    println("6-day ZonalJet integration: FAILED — NaNs in: ",
        join([v for (v, ok) in [("vor", vor_ok), ("div", div_ok), ("pres", pres_ok)] if !ok], ", "))
end
