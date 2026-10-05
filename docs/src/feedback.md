# [Feedback](@id feedback)

`Feedback` is the model component that prints a progress line of a running simulation to the
REPL and checks the prognostic variables for NaNs. Its options are described in its docstring
(`?Feedback`), this page explains how to customize the progress line.

## Progress line

While running, `Feedback` prints a progress line like
```
 23% ━━━━━━     ETA: 0:00:07 (2000-01-03, 153.56 years/day,  44 m/s, [ -84,   35] ˚C)
```
with the simulation date, speed, maximum zonal wind and temperature range. The line is built from
ProgressMeter.jl elements, `default_elements(feedback)` returns the default ones (controlled by the
options `showspeed`, `show_time`, `show_umax`, `show_temperature_range`). Pass your own tuple with
`elements` to change it, e.g. to also show the maximum vertical Courant number
```julia
feedback = Feedback(elements = (default_elements(Feedback())..., VerticalCourantNumber()))
model = PrimitiveWetModel(spectral_grid; feedback)
```
SpeedyWeather's elements (`SimulationTime`, `SimulationSpeed`, `MaximumWindSpeed`, `TemperatureRange`,
`VerticalCourantNumber`) are created without arguments and bound to the simulation when it starts.
Diagnostics are only computed when the line is redrawn (every `feedback_dt` seconds).
