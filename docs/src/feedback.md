# [Feedback](@id feedback)

`Feedback` is the model component that prints a progress line of a running simulation to the
REPL and checks the prognostic variables for NaNs. Its options are described in its docstring
(`?Feedback`), this page explains how to customize the progress line.

## Progress line

While running, `Feedback` prints a progress line like
```
 23%, ETA: 0:00:07, 2000-01-03, 153.56 years/day,  44 m/s, [ -84,   35] ˚C
```
with the percentage of time steps done, the estimated remaining time, the simulation date, speed,
maximum zonal wind and temperature range. The line is built from *elements*, which are collected in
the module `ProgressElements`: the generic ones (`Percentage`, `Bar`, `ETA`, `ElapsedTime`, `Counter`)
and the ones specific to SpeedyWeather (`SimulationTime`, `SimulationSpeed`, `MaximumWindSpeed`,
`TemperatureRange`, `VerticalCourantNumber`). `ProgressElements.default_elements()` returns the
default ones, which are used if `elements` is not set (with Reactant the maximum wind speed and the
temperature range are left out by default as computing them is slow there, list them explicitly
to show them). Pass your own tuple with `elements` to change the line, e.g. to also show the maximum
vertical Courant number
```julia
using SpeedyWeather
spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
feedback = Feedback(elements = (ProgressElements.default_elements()..., ProgressElements.VerticalCourantNumber()))
model = PrimitiveWetModel(spectral_grid; feedback)
```
which appends `, Cᵥ = 0.04` to the line.
Elements print only their value, `Feedback` puts its `separator` (default `", "`) between them, and a
single space next to a bar. `Bar(length = nothing)` fills the remaining width of the terminal,
`Bar(length = 20, fill = '#', empty = '.')` changes its length and characters. So
```julia
(; Percentage, Bar, SimulationTime, MaximumWindSpeed) = ProgressElements
feedback = Feedback(elements = (Percentage(), Bar(length = 10), SimulationTime(), MaximumWindSpeed()), separator = " | ")
```
prints ` 23% ━━         2000-01-03 |  44 m/s`. Elements that do not apply to a model are left out, e.g.
`TemperatureRange` and `VerticalCourantNumber` in a `BarotropicModel`. SpeedyWeather's elements are
created without arguments and bound to the simulation when it starts. Diagnostics are only computed
when the line is redrawn (every `feedback_dt` seconds).

## Custom elements

A custom element is a subtype of `ProgressElements.AbstractProgressElement` that extends
`ProgressElements.print_element(element, progress)` to return the text it shows. It is called on every
redraw with the progress meter `progress` (a `SpeedyWeatherInternals.ProgressMeters.Progress`), whose
fields `counter`, `n`, `tinit` (`time()` when the simulation started), `tcurrent` (`time()` of this redraw)
and `finished` (`true` for the last line) can be used, e.g.
```julia
struct StepsLeft <: ProgressElements.AbstractProgressElement end
ProgressElements.print_element(::StepsLeft, p) = "$(p.n - p.counter) steps left"

feedback = Feedback(elements = (ProgressElements.default_elements()..., StepsLeft()))
```
An element returning an empty string is left out of the line. To show something of the simulation,
hold a reference to it in the element (e.g. a `Ref` that a callback updates).

## NaN detection

With `Feedback(debug = true)` (default) the simulation is checked for NaNs (or Infs), which stops it.
NaNs propagate through the spectral transform to all coefficients, so whether the 1-0 mode of the
surface vorticity is finite is recorded on every time step. As copying it to the CPU would synchronize
a GPU, the record is only checked every `interval` time steps (default 50) and at the end, so the
simulation may run for up to `interval` time steps after the blow-up, but the warning reports the
time step and simulation time at which it happened, e.g.
```
┌ Warning: NaN or Inf detected at time step 1234 (2000-02-01T12:00:00)
```
which are also stored in `simulation.model.feedback.nans_detected_step` and `nans_detected_time`.

## Acknowledgement

The progress meter of SpeedyWeather (`SpeedyWeatherInternals.ProgressMeters`) is based on
[ProgressMeter.jl](https://github.com/timholy/ProgressMeter.jl) by Tim Holy and contributors (MIT
License), which SpeedyWeather used until v0.24. As SpeedyWeather needs a progress line built from its
own elements, which ProgressMeter.jl only allowed through type piracy, and the corresponding proposal
for ProgressMeter.jl ([#369](https://github.com/timholy/ProgressMeter.jl/pull/369)) has not been
merged, the small part of ProgressMeter.jl that SpeedyWeather uses was reimplemented with a simpler
element interface. Thanks to the ProgressMeter.jl authors.

