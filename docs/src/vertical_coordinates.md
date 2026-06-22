# [Vertical coordinates](@id vertical_coordinates_page)

SpeedyWeather.jl supports two families of terrain-following vertical coordinates.
Both are implemented as subtypes of `AbstractVerticalCoordinate` and are interchangeable
in the model constructor. For the mathematical background on coordinate transformations
see [Vertical coordinates](@ref) in the Primitive equation model documentation.

## [Sigma coordinates](@id sigma_coordinates_usage)

Sigma coordinates use the fraction of surface pressure as the vertical coordinate

```math
\sigma = \frac{p}{p_s}, \qquad \sigma \in [0, 1]
```

with ``\sigma = 0`` at the model top and ``\sigma = 1`` at the surface.
Sigma levels are terrain-following because ``\sigma = 1`` is always at surface pressure,
bending around every mountain. One specifies the half levels ``\sigma_{k+\tfrac{1}{2}}`` and
the full levels are obtained as midpoints

```math
\sigma_k = \frac{\sigma_{k+\tfrac{1}{2}} + \sigma_{k-\tfrac{1}{2}}}{2}, \qquad
p_k = \sigma_k p_s, \qquad
\Delta p_k = \Delta\sigma_k \, p_s
```

### Creating sigma coordinates

```@example vertical_coordinates
using SpeedyWeather
spectral_grid = SpectralGrid(nlayers = 8)
SigmaCoordinates(spectral_grid)
```

From a custom vector of half levels (`nlayers` is inferred from the length):

```@example vertical_coordinates
SigmaCoordinates([0.0, 0.1, 0.25, 0.45, 0.62, 0.78, 0.9, 0.97, 1.0])
```

From a range:

```@example vertical_coordinates
SigmaCoordinates(0:0.2:1)
```

Frierson (2006) spacing with finer resolution near the surface and stratosphere:

```@example vertical_coordinates
FriersonSigmaCoordinates(spectral_grid)
```

### Using sigma coordinates in a model

Pass the coordinate to `Geometry` and then to the model constructor:

```@example vertical_coordinates
σ = FriersonSigmaCoordinates(spectral_grid)
geometry = Geometry(spectral_grid; vertical_coordinates = σ)
model = PrimitiveDryModel(spectral_grid; geometry)
simulation = initialize!(model)
nothing # hide
```

## [Hybrid sigma-pressure coordinates](@id hybrid_sigma_pressure_usage)

!!! warning "Work in progress"
    Hybrid sigma-pressure coordinates are implemented for the coordinate geometry,
    the `pressure`, `pressure_thickness`, and `sigma` functions, and all
    parameterizations, but **the dynamical core still uses the pure sigma formulation**
    throughout (vertical advection, surface pressure tendency, geopotential, and adiabatic
    conversion). See the TODO comments in those source files for details. Using
    `SigmaPressureCoordinates` with a transition close to pure sigma (e.g. the default
    ``f(\sigma) = \sigma``) therefore remains consistent with the dynamics.

Pure sigma coordinates tilt sharply with orography even at high altitudes where
terrain-following levels are not needed. Hybrid sigma-pressure coordinates solve this
by blending constant-pressure surfaces near the model top with terrain-following sigma
surfaces near the surface.

The pressure at layer ``k`` is

```math
p_k = A_k \, p_{\mathrm{ref}} + B_k \, p_s
```

where ``p_{\mathrm{ref}}`` is a constant reference pressure and ``A_k``, ``B_k`` are
layer coefficients. Given a nominal sigma level ``\sigma_k`` and a transition function
``f(\sigma) \in [0, 1]`` the coefficients are

```math
A_k = \sigma_k \bigl(1 - f(\sigma_k)\bigr), \qquad B_k = \sigma_k f(\sigma_k)
```

so that ``A_k + B_k = \sigma_k`` always holds. When ``f = 0`` the level is a pure
constant-pressure surface; when ``f = 1`` it is a pure sigma surface. The layer
thickness in pressure is

```math
\Delta p_k = \Delta A_k \, p_{\mathrm{ref}} + \Delta B_k \, p_s
```

with ``\Delta A_k = A_{k+\tfrac{1}{2}} - A_{k-\tfrac{1}{2}}`` and equivalently for ``B``.

### Creating hybrid sigma-pressure coordinates

The default transition is ``f(\sigma) = \sigma``, giving a linear blend:

```@example vertical_coordinates
SigmaPressureCoordinates(spectral_grid)
```

Pure pressure levels everywhere (``f = 0``):

```@example vertical_coordinates
SigmaPressureCoordinates(spectral_grid; transition = _ -> 0.0)
```

Pure sigma levels everywhere (``f = 1``, equivalent to `SigmaCoordinates`):

```@example vertical_coordinates
SigmaPressureCoordinates(spectral_grid; transition = _ -> 1.0)
```

A ready-made `CubicSigmaPressureCoordinates` uses a cubic smoothstep transition:
pure pressure for ``\sigma \le \sigma_{\mathrm{low}}``, pure sigma for
``\sigma \ge \sigma_{\mathrm{high}}``, and a smooth cubic interpolation in between.
The bar on each full level shows the A/B split (█ = sigma fraction, ░ = pressure fraction):

```@example vertical_coordinates
CubicSigmaPressureCoordinates(spectral_grid)
```

Custom thresholds can be set:

```@example vertical_coordinates
CubicSigmaPressureCoordinates(spectral_grid; pressure_only_above = 0.1, σ_only_below = 0.9)
```

### Using hybrid coordinates in a model

```@example vertical_coordinates
S = CubicSigmaPressureCoordinates(spectral_grid)
geometry = Geometry(spectral_grid; vertical_coordinates = S)
model = PrimitiveDryModel(spectral_grid; geometry)
simulation = initialize!(model)
nothing # hide
```

!!! info "Reference pressure and atmosphere"
    If the `reference_pressure` of the `SigmaPressureCoordinates` differs from the
    `reference_pressure` of the atmosphere (default 1e5 Pa), a warning is issued during
    `initialize!`. Make sure they match for physical consistency.

## Coordinate functions

Three functions provide a coordinate-agnostic interface for evaluating the vertical
coordinate at a given full level ``k``, working for both `SigmaCoordinates` and
`SigmaPressureCoordinates`:

| Function | Returns |
|---|---|
| `SpeedyWeather.pressure(k, pₛ, coord)` | Pressure [Pa] at full level ``k`` |
| `SpeedyWeather.pressure_thickness(k, pₛ, coord)` | Pressure thickness [Pa] of layer ``k`` |
| `SpeedyWeather.sigma(k, coord)` | Nominal sigma level ``\sigma_k`` (surface-pressure independent) |

For `SigmaCoordinates` these reduce to ``\sigma_k p_s``, ``\Delta\sigma_k p_s``, and
``\sigma_k`` respectively. For `SigmaPressureCoordinates` they use the full
``A_k p_{\mathrm{ref}} + B_k p_s`` formula. Note that `sigma` always equals
``A_k + B_k = \sigma_k`` regardless of the transition.
