# Diagnostics

Diagnostics compute variables that are only used for output, e.g. accumulated precipitation
or the daily maximum temperature. They are calculated on the fly every time step and stored
in `simulation.variables.output`. Nothing in the dynamics or parameterizations reads from
`variables.output`, so diagnostics cannot change the simulation, they only add work to the
time step. The primitive equation models hold them in `model.diagnostics`, a `NamedTuple`
of `AbstractDiagnostic`s.

## Default diagnostics

The `PrimitiveWetModel` accumulates large-scale and convective rain and snow as well as the
river runoff, all in meters. The parameterizations only compute rates (in m/s, in
`variables.parameterizations`) and the [`Accumulate`](@ref) diagnostic sums them up over time

```@example diagnostics
using SpeedyWeather
spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
model = PrimitiveWetModel(spectral_grid)
model.diagnostics
```

The diagnostics define their own variables in `variables.output`, so after running a simulation

```@example diagnostics
simulation = initialize!(model)
run!(simulation, period = Day(2))
simulation.variables.output
```

Accumulated variables are also written to restart files and restored with
[`StartFromFile`](@ref) so that accumulations continue across restarts.

## Adding diagnostics

To add a diagnostic, merge it into the default diagnostics, e.g. the daily maximum
surface air temperature with [`PeriodicMaximum`](@ref)

```@example diagnostics
diagnostics = (; PrecipitationDiagnostics()..., temperature_max = PeriodicMaximum())
model = PrimitiveWetModel(spectral_grid; diagnostics)
simulation = initialize!(model)
run!(simulation, period = Day(2))
simulation.variables.output.surface_air_temperature_maximum
```

`PeriodicMaximum` keeps a running maximum, which it copies into `surface_air_temperature_maximum`
at the end of every period (a day by default), and resets. Before the first period has passed the
maximum is `NaN`. Pass `schedule` for another period, e.g. `PeriodicSchedule(every = Hour(6))`,
and `variable`, `group` and `units` to take the maximum of another 2D grid variable.
Without `PrecipitationDiagnostics()` in `diagnostics` the precipitation is not accumulated, and
precipitation output variables like `ConvectiveRain` are not available.

## Custom diagnostics

A diagnostic is a subtype of `AbstractDiagnostic` which

1. defines the variables it needs via `variables(diagnostic)` returning `OutputVariable`s
   (and possibly the variables it reads from) so that they are allocated in `variables.output`,
2. optionally initializes these with `initialize!(variables, diagnostic, model)`,
3. and calculates them in `diagnostic!(variables, diagnostic, model)`, called once per time
   step after the time step, before the output is written.

E.g. a diagnostic for the time-mean surface air temperature, summing it up and counting the
time steps

```@example diagnostics
struct MeanSurfaceAirTemperature <: SpeedyWeather.AbstractDiagnostic end

SpeedyWeather.variables(::MeanSurfaceAirTemperature) = (
    OutputVariable(:surface_air_temperature_sum, SpeedyWeather.Grid2D(), units = "K"),
    OutputVariable(:surface_air_temperature_count, SpeedyWeather.ScalarDim()),
)

function SpeedyWeather.diagnostic!(vars, ::MeanSurfaceAirTemperature, model)
    vars.output.surface_air_temperature_sum .+= vars.parameterizations.surface_air_temperature
    vars.output.surface_air_temperature_count[] += 1
    return nothing
end

diagnostics = (; PrecipitationDiagnostics()..., temperature_mean = MeanSurfaceAirTemperature())
model = PrimitiveWetModel(spectral_grid; diagnostics)
simulation = initialize!(model)
run!(simulation, period = Day(1))

(; surface_air_temperature_sum, surface_air_temperature_count) = simulation.variables.output
surface_air_temperature_mean = surface_air_temperature_sum / surface_air_temperature_count[]
```

To write diagnostics into the netCDF output, define an output variable whose
`path` points into `simulation.variables.output`, see [Customizing netCDF output](@ref).

```@docs; canonical=false
Accumulate
PeriodicMaximum
```
