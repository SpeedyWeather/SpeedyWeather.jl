# [Vertical interpolation onto pressure layers](@id vertical_interpolation)

SpeedyWeather.jl integrates on terrain-following vertical coordinates, either
[Sigma coordinates](@ref sigma_coordinates_usage) or
[Hybrid sigma-pressure coordinates](@ref hybrid_sigma_pressure_usage), see
[Vertical coordinates](@ref vertical_coordinates_page). Model layers therefore sit at a
different pressure in every column, bending around the orography. For analysis and for
comparison with observations or reanalyses one usually wants the data on fixed pressure
layers instead, e.g. temperature at 850 hPa everywhere. This page describes how to
interpolate from the model's vertical layers onto pressure layers.

## Interpolating a field onto pressure layers

`SpeedyWeather.interpolate_pressure_layers!` interpolates a field on model layers onto pressure
layers, writing into an output field that you allocate yourself

```julia
SpeedyWeather.interpolate_pressure_layers!(
    out_field,          # OUTPUT: (horizontal, npressure)
    in_field,           # INPUT: (horizontal, nlayers) on model layers
    surface_pressure,   # INPUT: (horizontal,) surface pressure [Pa]
    p,                  # pressure layers [Pa]
    coordinates,        # model.geometry.vertical_coordinates
    interpolation,      # optional, LinearInLogPressure() by default
    extrapolation,      # optional, ConstantExtrapolation() by default
)
```

The pressure of every model layer in every column follows from the surface pressure and
the vertical `coordinates`, so the same function covers sigma and hybrid sigma-pressure
coordinates. The interpolation is column-local and therefore launched as a kernel over
horizontal grid points and pressure layers, running on CPU and GPU. Consequently the
pressure layers `p` have to be on the same architecture (see [GPU and Architectures](@ref)) as the
fields, and ideally of the same number format. They do not have to be sorted.

The number of pressure layers is independent of the number of model layers, it is the
second dimension of `out_field` that decides. Note that `in_field` must not have a time
step dimension, so for prognostic variables select the step first, see
[Step dimension](@ref).

## Interpolation methods

Interpolating a variable ``\xi`` between two model layers at pressures ``p_1 \leq p \leq p_2``
uses the weight ``w`` of the lower layer

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

## Extrapolation beyond the model layers

A pressure layer can also lie outside the range spanned by the full model layers: above
the top-most layer, or below the bottom-most layer. The latter splits into two cases,
because the lowest full model layer is not the surface: a layer below the lowest full
layer can still be above ground (`p < pₛ`), or genuinely below ground (`p > pₛ`).
The following are available

- `ConstantExtrapolation` (default) holds the outer-most model layer constant.
- `DryAdiabaticExtrapolation` descends dry-adiabatically below the lowest model
  layer, ``T(p) = T_\text{bottom} (p/p_\text{bottom})^\kappa``, which is what you want for
  a temperature at, say, 1000 hPa below a lowest model layer that sits above it. It is the
  same adiabat that the mean sea-level pressure output uses.
- `SubsurfaceMask` masks everything below the surface with a `missing_value`
  (`NaN` by default) and uses another extrapolation between the lowest model layer and the
  surface.

Above the model top every method currently holds the top-most layer constant.

## Example

```@example vertical_interpolation
using SpeedyWeather
spectral_grid = SpectralGrid(truncation = 31, nlayers = 8)
model = PrimitiveDryModel(spectral_grid)
simulation = initialize!(model)
run!(simulation, period = Hour(6))
nothing # hide
```

Take the temperature on model layers (selecting the time step to interpolate) and the
surface pressure

```@example vertical_interpolation
(; variables) = simulation
temp = SpeedyWeather.get_prognostic_step(variables.grid.temperature, model.time_stepping, model.output)
pₛ = variables.parameterizations.surface_pressure
summary(temp)
```

allocate an output field for the pressure layers of interest, and interpolate

```@example vertical_interpolation
p = spectral_grid.NF[850e2, 500e2, 200e2]       # pressure layers in Pa
temp_p = zeros(spectral_grid.NF, spectral_grid.grid, length(p))

SpeedyWeather.interpolate_pressure_layers!(temp_p, temp, pₛ, p, model.geometry.vertical_coordinates)
[sum(temp_p[:, k]) / length(pₛ) for k in eachindex(p)]   # mean temperature [K] per layer
```

A 1000 hPa layer is below ground over much of the globe. Extrapolating
dry-adiabatically down to it while masking the points that are genuinely below the
surface

```@example vertical_interpolation
p = spectral_grid.NF[1000e2]
temp_1000 = zeros(spectral_grid.NF, spectral_grid.grid, length(p))

SpeedyWeather.interpolate_pressure_layers!(
    temp_1000, temp, pₛ, p, model.geometry.vertical_coordinates,
    SpeedyWeather.LinearInLogPressure(),
    SpeedyWeather.SubsurfaceMask(
        above_surface = SpeedyWeather.DryAdiabaticExtrapolation(model.atmosphere.κ),
    ),
)
count(isnan, temp_1000), length(pₛ)     # masked points of total
```

## Output on pressure layers

You normally don't need to call any of this yourself: every output writer takes a `layers`
keyword argument to write its 3D variables on pressure layers instead of model layers

```julia
output = NetCDFOutput(spectral_grid, PrimitiveWet, layers = PressureLayers(spectral_grid, [850, 500, 200] .* 100))
```

See [Output layers](@ref output_layers) for the details, including how a variable chooses
its extrapolation below the lowest model layer.

## [Vertical interpolation onto isentropic layers](@id isentropic_interpolation)

Isentropic layers are surfaces of constant dry potential temperature
``\theta = T (p_0/p)^\kappa``, with ``p_0`` the atmosphere's `reference_pressure`
(1000 hPa) and ``\kappa = R_\text{dry}/c_p``. `SpeedyWeather.interpolate_isentropic_layers!`
interpolates a field on model layers onto them

```julia
SpeedyWeather.interpolate_isentropic_layers!(
    out_field,          # OUTPUT: (horizontal, nθ)
    in_field,           # INPUT: (horizontal, nlayers) on model layers
    temperature,        # INPUT: (horizontal, nlayers) absolute temperature [K]
    surface_pressure,   # INPUT: (horizontal,) surface pressure [Pa]
    θ,                  # isentropic layers [K]
    coordinates,        # model.geometry.vertical_coordinates
    κ, p₀,              # model.atmosphere.κ, model.atmosphere.reference_pressure
    extrapolation,      # optional, ConstantExtrapolation() by default
)
```

The potential temperature of every model layer follows from the temperature and the
pressure of that layer, which, as for pressure layers, follows from the surface pressure
and the vertical coordinates. Like `interpolate_pressure_layers!` this runs as a kernel
over horizontal grid points and isentropic layers on CPU and GPU, and needs at least two
model layers.

### Interpolation

Between two model layers at pressures ``p_1 < p_2`` with potential temperatures
``\theta_1, \theta_2``, ``\theta`` is taken to be linear in ``\ln p``. The pressure ``p^*``
of the isentrope ``\theta`` then follows in closed form

```math
\ln p^* = \ln p_1 + w (\ln p_2 - \ln p_1), \qquad w = \frac{\theta - \theta_1}{\theta_2 - \theta_1}
```

and every variable is interpolated linearly in ``\ln p`` at ``p^*``, as
`LinearInLogPressure` does for pressure layers. That is the same weight ``w``,

```math
\xi(\theta) = (1 - w) \xi_1 + w \xi_2.
```

In a statically unstable column ``\theta`` is not monotonic and an isentrope can cross
it more than once. Every column is searched from the top, so the first crossing from the
top is used.

### Extrapolation

Isentropes often lie below the lowest model layer, and below the ground: e.g. 290 K is
underground in the tropics. As for pressure layers, the `extrapolation` methods above
decide, at the isentrope's pressure ``p^*``. Outside the column, ``p^*`` is extrapolated with
``\theta`` linear in ``\ln p`` from the two outer-most model layers. So
`DryAdiabaticExtrapolation` descends from the lowest model layer to ``p^*``.
An isentrope below the ground, ``p^* > p_s``, follows the Lorenz condition: it is put at
the surface, ``p^* = p_s``, as for a weakly stratified lowest layer pair the extrapolated
``p^*`` can be arbitrarily large. `SubsurfaceMask` masks these isentropes instead. If the
two outer-most layers are not stably stratified, no ``p^*`` can be extrapolated and the
Lorenz condition is used as well, ``p^* = p_s`` below and ``p^*`` = the top-most layer's
pressure above.

### Example

Continuing the example above, temperature and zonal wind on the 300 K and 330 K isentropes

```@example vertical_interpolation
u = SpeedyWeather.get_prognostic_step(variables.grid.u, model.time_stepping, model.output)
θ = spectral_grid.NF[300, 330]                  # isentropic layers in K
temp_θ = zeros(spectral_grid.NF, spectral_grid.grid, length(θ))
u_θ = zeros(spectral_grid.NF, spectral_grid.grid, length(θ))

(; κ, reference_pressure) = model.atmosphere
coordinates = model.geometry.vertical_coordinates
SpeedyWeather.interpolate_isentropic_layers!(
    temp_θ, temp, temp, pₛ, θ, coordinates, κ, reference_pressure,
    SpeedyWeather.DryAdiabaticExtrapolation(κ),
)
SpeedyWeather.interpolate_isentropic_layers!(u_θ, u, temp, pₛ, θ, coordinates, κ, reference_pressure)
[sum(temp_θ[:, k]) / length(pₛ) for k in eachindex(θ)]   # mean temperature [K] per isentrope
```

## Output on isentropic layers

Output writers take `IsentropicLayers` for their `layers` keyword argument, for models
with temperature (`PrimitiveDry`, `PrimitiveWet`)

```julia
output = NetCDFOutput(spectral_grid, PrimitiveDry, layers = IsentropicLayers(spectral_grid, [300, 315, 330, 350]))
```

See [Output layers](@ref output_layers).
