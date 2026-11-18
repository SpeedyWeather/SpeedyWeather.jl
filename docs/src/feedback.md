# [Feedback](@id feedback)

`Feedback` is the model component that prints a progress line of a running simulation to the
REPL and checks the prognostic variables for NaNs. Its options are described in its docstring
(`?Feedback`), this page explains how to customize the progress line.

## Progress line

While running, `Feedback` prints two progress lines like
```
 23% ━━━━━━━━━━━━━━━━━━━━━━━                                   ETA: 0:00:07, 2000-01-03, 153.56 years/day
     44 m/s, [-84, 35] ˚C, Cᵥ = 0.05, precip 2.91 mm/day, TOA +1.2 W/m²
```
with the time and performance information on the first line: percentage of time steps done, a bar, the
estimated remaining time, the simulation date and speed. The second line shows diagnostics of the
simulation: the maximum zonal wind, the temperature range, the maximum vertical Courant number
(`|σ̇|Δt/Δσ`, a stability indicator), the global mean precipitation and the global mean energy
imbalance at the top of the atmosphere (incoming shortwave `S₀/4` minus outgoing shortwave minus
outgoing longwave radiation). Elements that do not apply to a model are left out, e.g. the temperature
range, Courant number, precipitation and energy imbalance in a `BarotropicModel`.

The lines are built from *elements*, which are collected in the module `ProgressElements`: the generic
ones (`Percentage`, `Bar`, `ETA`, `ElapsedTime`, `Counter`, `NewLine`) and the ones specific to SpeedyWeather
(`SimulationTime`, `SimulationSpeed`, `MaximumWindSpeed`, `TemperatureRange`, `VerticalCourantNumber`,
`Precipitation`, `EnergyImbalance`). `ProgressElements.default_elements()` returns the default ones,
```julia
(Percentage(), Bar(), ETA(), SimulationTime(), SimulationSpeed(), NewLine(indent = 5),
    MaximumWindSpeed(), TemperatureRange(), VerticalCourantNumber(), Precipitation(), EnergyImbalance())
```
which are used if `elements` is not set (with Reactant only the first line, as computing the
diagnostics is slow there, list them explicitly to show them). Pass your own tuple with `elements` to
change the lines, e.g. for a single line without the bar
```julia
using SpeedyWeather
spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
(; Percentage, ETA, SimulationTime, MaximumWindSpeed, EnergyImbalance) = ProgressElements
feedback = Feedback(elements = (Percentage(), ETA(), SimulationTime(), MaximumWindSpeed(), EnergyImbalance()))
model = PrimitiveWetModel(spectral_grid; feedback)
```
Elements print only their value, `Feedback` puts its `separator` (default `", "`) between them, and a
single space next to a bar. `NewLine(indent = 5)` starts a new line, indented by 5 spaces (plus the
width of the description). `Bar(length = nothing)` fills the remaining width of the terminal,
`Bar(length = 20, fill = '#', empty = '.')` changes its length and characters. So
```julia
(; Percentage, Bar, SimulationTime, MaximumWindSpeed) = ProgressElements
feedback = Feedback(elements = (Percentage(), Bar(length = 10), SimulationTime(), MaximumWindSpeed()), separator = " | ")
```
prints ` 23% ━━         2000-01-03 | 44 m/s`. Lines are cut to the width of the terminal so that they
don't wrap. SpeedyWeather's elements are created without arguments and bound to the simulation when it
starts. Diagnostics are only computed when the lines are redrawn (every `feedback_dt` seconds).

## Custom elements

A custom element is a subtype of `ProgressElements.AbstractProgressElement` that extends
`ProgressElements.print_element(element, progress)` to return the text it shows. It is called on every
redraw with the progress meter `progress` (a `SpeedyWeatherInternals.ProgressMeters.Progress`), whose
fields `counter`, `n`, `tinit` (`time()` when the simulation started), `tcurrent` (`time()` of this redraw)
and `finished` (`true` for the last line) can be used, e.g.
```julia
struct StepsLeft <: ProgressElements.AbstractProgressElement end
ProgressElements.print_element(::StepsLeft, p) = "$(p.n - p.counter) steps left"

feedback = Feedback(elements = (ProgressElements.default_elements()..., StepsLeft()))   # appended to the 2nd line
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

