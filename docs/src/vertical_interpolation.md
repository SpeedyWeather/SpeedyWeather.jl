# [Vertical interpolation onto pressure levels](@id vertical_interpolation)

SpeedyWeather.jl integrates on terrain-following vertical coordinates, either
[Sigma coordinates](@ref sigma_coordinates_usage) or
[Hybrid sigma-pressure coordinates](@ref hybrid_sigma_pressure_usage), see
[Vertical coordinates](@ref vertical_coordinates_page). Model levels therefore sit at a
different pressure in every column, bending around the orography. For analysis and for
comparison with observations or reanalyses one usually wants the data on fixed pressure
levels instead, e.g. temperature at 850 hPa everywhere. This page describes how to
interpolate from the model's vertical levels onto pressure levels.

## Interpolating a field onto pressure levels

`SpeedyWeather.interpolate_pressure_levels!` interpolates a field on model levels onto pressure
levels, writing into an output field that you allocate yourself

```julia
SpeedyWeather.interpolate_pressure_levels!(
    out_field,          # OUTPUT: (horizontal, npressure)
    in_field,           # INPUT: (horizontal, nlayers) on model levels
    surface_pressure,   # INPUT: (horizontal,) surface pressure [Pa]
    p,                  # pressure levels [Pa]
    coordinates,        # model.geometry.vertical_coordinates
    interpolation,      # optional, LinearInLogPressure() by default
    extrapolation,      # optional, ConstantExtrapolation() by default
)
```

The pressure of every model level in every column follows from the surface pressure and
the vertical `coordinates`, so the same function covers sigma and hybrid sigma-pressure
coordinates. The interpolation is column-local and therefore launched as a kernel over
horizontal grid points and pressure levels, running on CPU and GPU. Consequently the
pressure levels `p` have to be on the same architecture (see [GPU and Architectures](@ref)) as the
fields, and ideally of the same number format. They do not have to be sorted.

The number of pressure levels is independent of the number of model layers, it is the
second dimension of `out_field` that decides. Note that `in_field` must not have a time
step dimension, so for prognostic variables select the step first, see
[Step dimension](@ref).

## Interpolation methods

Interpolating a variable ``\xi`` between two model levels at pressures ``p_1 \leq p \leq p_2``
uses the weight ``w`` of the lower level

```math
\xi(p) = (1 - w) \xi_1 + w \xi_2
```

with the two available choices for ``w`` being linear in pressure

```math
w = \frac{p - p_1}{p_2 - p_1} \qquad \text{(\texttt{LinearInPressure})}
```

or linear in the logarithm of pressure

```math
w = \frac{\log(p/p_1)}{\log(p_2/p_1)} \qquad \text{(\texttt{LinearInLogPressure})}
```

`SpeedyWeather.LinearInLogPressure` is the default as most variables vary more linearly with
``\log p`` than with ``p``.

## Extrapolation beyond the model levels

A pressure level can also lie outside the range spanned by the full model levels: above
the top-most level, or below the bottom-most level. The latter splits into two cases,
because the lowest full model level is not the surface: a level below the lowest full
level can still be above ground (`p < pₛ`), or genuinely below ground (`p > pₛ`).
The following are available

- `ConstantExtrapolation` (default) holds the outer-most model level constant.
- `DryAdiabaticExtrapolation` descends dry-adiabatically below the lowest model
  level, ``T(p) = T_\text{bottom} (p/p_\text{bottom})^\kappa``, which is what you want for
  a temperature at, say, 1000 hPa below a lowest model level that sits above it. It is the
  same adiabat that the mean sea-level pressure output uses.
- `SubsurfaceMask` masks everything below the surface with a `missing_value`
  (`NaN` by default) and uses another extrapolation between the lowest model level and the
  surface.

Above the model top every method currently holds the top-most level constant.

## Example

```@example vertical_interpolation
using SpeedyWeather
spectral_grid = SpectralGrid(truncation = 31, nlayers = 8)
model = PrimitiveDryModel(spectral_grid)
simulation = initialize!(model)
run!(simulation, period = Hour(6))
nothing # hide
```

Take the temperature on model levels (selecting the time step to interpolate) and the
surface pressure

```@example vertical_interpolation
(; variables) = simulation
temp = SpeedyWeather.get_prognostic_step(variables.grid.temperature, model.time_stepping, model.output)
pₛ = variables.parameterizations.surface_pressure
summary(temp)
```

allocate an output field for the pressure levels of interest, and interpolate

```@example vertical_interpolation
p = spectral_grid.NF[850e2, 500e2, 200e2]       # pressure levels in Pa
temp_p = zeros(spectral_grid.NF, spectral_grid.grid, length(p))

SpeedyWeather.interpolate_pressure_levels!(temp_p, temp, pₛ, p, model.geometry.vertical_coordinates)
[sum(temp_p[:, k]) / length(pₛ) for k in eachindex(p)]   # mean temperature [K] per level
```

A 1000 hPa level is below ground over much of the globe. Extrapolating
dry-adiabatically down to it while masking the points that are genuinely below the
surface

```@example vertical_interpolation
p = spectral_grid.NF[1000e2]
temp_1000 = zeros(spectral_grid.NF, spectral_grid.grid, length(p))

SpeedyWeather.interpolate_pressure_levels!(
    temp_1000, temp, pₛ, p, model.geometry.vertical_coordinates,
    SpeedyWeather.LinearInLogPressure(),
    SpeedyWeather.SubsurfaceMask(
        above_surface = SpeedyWeather.DryAdiabaticExtrapolation(model.atmosphere.κ),
    ),
)
count(isnan, temp_1000), length(pₛ)     # masked points of total
```

## Output on pressure levels

You normally don't need to call any of this yourself: every output writer takes a `levels`
keyword argument to write its 3D variables on pressure levels instead of model levels

```julia
output = NetCDFOutput(spectral_grid, PrimitiveWet, levels = PressureLevels([850, 500, 200] .* 100))
```

See [Output levels](@ref output_levels) for the details, including how a variable chooses
its extrapolation below the lowest model level.
