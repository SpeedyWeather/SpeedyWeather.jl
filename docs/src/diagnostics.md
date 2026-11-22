# Diagnostics

Diagnostics compute variables that are only used for output, e.g. accumulated precipitation
or the daily maximum temperature. They are calculated on the fly every time step and stored
in `simulation.variables.output`. Nothing in the dynamics or parameterizations reads from
`variables.output`, so diagnostics cannot change the simulation, they only add work to the
time step. The primitive equation models hold them in `model.diagnostics`, a `NamedTuple`
of `AbstractDiagnostic`s.

The diagnostics are responsible for the accumulation (sum, maximum, mean, ...) every time step,
the output writer resets them after every output. That way the output interval determines
the period over which a diagnostic is taken: precipitation is accumulated since the last
output, with daily output [`Maximum`](@ref) is the daily maximum and with weekly output
[`Mean`](@ref) the weekly mean. Without output, diagnostics are never reset and taken over
the whole simulation.

## Default diagnostics

The `PrimitiveWetModel` accumulates large-scale and convective rain and snow as well as the
river runoff, all in meters. The parameterizations only compute rates (in m/s, in
`variables.parameterizations`) and the [`Accumulate`](@ref) diagnostic integrates them
over time

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

The output variables are also written to restart files and restored with
[`StartFromFile`](@ref) so that, e.g., accumulations continue across restarts.

## Adding diagnostics

To add a diagnostic, merge it into the default diagnostics, e.g. the maximum and mean
surface air temperature with [`Maximum`](@ref) and [`Mean`](@ref). With daily output
they are the daily maximum and daily mean

```@example diagnostics
diagnostics = (;
    PrecipitationDiagnostics()...,
    temperature_max = Maximum(),
    temperature_mean = Mean(),
)
output = ArrayOutput(interval = Day(1), groups = (:output,))
model = PrimitiveWetModel(spectral_grid; diagnostics, output)
simulation = initialize!(model)
run!(simulation, period = Day(2), output = true)
model.output.output[end].output.surface_air_temperature_maximum
```

Pass `variable`, `group`, `namespace` and `units` to take the maximum or mean of another
2D grid variable. Without `PrecipitationDiagnostics()` in `diagnostics` the precipitation is
not accumulated, and precipitation output variables like `ConvectiveRain` are not available.

## Custom diagnostics

A diagnostic is a subtype of `AbstractDiagnostic` which

1. defines the variables it needs via `variables(diagnostic)` returning `OutputVariable`s
   so that they are allocated in `variables.output`,
2. optionally initializes these with `initialize!(variables, diagnostic, model)`,
3. calculates them in `diagnostic!(variables, diagnostic, model)`, called once per time
   step after the time step, before the output is written,
4. and resets them in `reset_diagnostic!(variables, diagnostic, model)`, called by the
   output writer after every output.

E.g. a diagnostic for the minimum surface air temperature

```@example diagnostics
struct MinimumSurfaceAirTemperature <: SpeedyWeather.AbstractDiagnostic end

SpeedyWeather.variables(::MinimumSurfaceAirTemperature) =
    (OutputVariable(:surface_air_temperature_minimum, SpeedyWeather.Grid2D(), units = "K"),)

SpeedyWeather.initialize!(vars, diagnostic::MinimumSurfaceAirTemperature, model) =
    SpeedyWeather.reset_diagnostic!(vars, diagnostic, model)

SpeedyWeather.reset_diagnostic!(vars, ::MinimumSurfaceAirTemperature, model) =
    fill!(vars.output.surface_air_temperature_minimum, Inf)

function SpeedyWeather.diagnostic!(vars, ::MinimumSurfaceAirTemperature, model)
    (; surface_air_temperature_minimum) = vars.output
    surface_air_temperature_minimum .= min.(surface_air_temperature_minimum, vars.parameterizations.surface_air_temperature)
    return nothing
end

diagnostics = (; PrecipitationDiagnostics()..., temperature_min = MinimumSurfaceAirTemperature())
model = PrimitiveWetModel(spectral_grid; diagnostics)
simulation = initialize!(model)
run!(simulation, period = Day(1))
simulation.variables.output.surface_air_temperature_minimum
```

To write diagnostics into the netCDF output, define an output variable whose
`path` points into `simulation.variables.output`, see [Customizing netCDF output](@ref).

```@docs; canonical=false
Accumulate
Maximum
Mean
```
