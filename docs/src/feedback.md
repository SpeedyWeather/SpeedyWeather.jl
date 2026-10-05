# [Feedback](@id feedback)

`Feedback` is the model component that prints a progress line of a running simulation to the
REPL and checks the prognostic variables for NaNs. Its options are described in its docstring
(`?Feedback`), this page explains how to customize the progress line.

## Progress line

While running, `Feedback` prints a progress line like
```
 23% ━━━━━━     ETA: 0:00:07, 2000-01-03, 153.56 years/day,  44 m/s, [ -84,   35] ˚C
```
with the simulation date, speed, maximum zonal wind and temperature range. The line is built from
*elements*, which are collected in the module `ProgressElements`: the generic ones from ProgressMeter.jl
(`Description`, `Percentage`, `Bar`, `ETA`, `Speed`, `ElapsedTime`, `Counter`, `Colored`) and the ones
specific to SpeedyWeather (`SimulationTime`, `SimulationSpeed`, `MaximumWindSpeed`, `TemperatureRange`,
`VerticalCourantNumber`). `ProgressElements.default_elements()` returns the default ones. Pass your own
tuple with `elements` to change the line, e.g. to also show the maximum vertical Courant number
```julia
feedback = Feedback(elements = (ProgressElements.default_elements()..., ProgressElements.VerticalCourantNumber()))
model = PrimitiveWetModel(spectral_grid; feedback)
```
which appends `, Cᵥ = 0.04` to the line. SpeedyWeather's elements only print their value, `Feedback`
puts its `separator` (default `", "`) in front of each of them unless nothing but the description
comes before. ProgressMeter.jl's elements and strings (e.g. `" "`) are printed as they are, so
```julia
feedback = Feedback(
    elements = (ProgressElements.Percentage(), ProgressElements.SimulationTime(), ProgressElements.MaximumWindSpeed()),
    separator = " | ",
)
```
prints ` 23% | 2000-01-03 |  44 m/s`. Elements that do not apply to a model are left out, e.g.
`TemperatureRange` and `VerticalCourantNumber` in a `BarotropicModel`. SpeedyWeather's elements are
created without arguments and bound to the simulation when it starts. Diagnostics are only computed
when the line is redrawn (every `feedback_dt` seconds).

A custom element is a subtype of `ProgressElements.AbstractProgressElement` that extends
`ProgressElements.print_element(element, progress)` to return the text it shows, see the
documentation of ProgressMeter.jl. Subtype `ProgressElements.AbstractSimulationElement` instead to
have it separated like SpeedyWeather's elements.
