# [Vertical coordinates](@id vertical_coordinates_page)

SpeedyWeather.jl supports two families of terrain-following vertical coordinates.
Both are implemented as subtypes of `AbstractVerticalCoordinates` and are interchangeable
in the model constructor. For more mathematical background on coordinate transformations
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
pure pressure above (in the vertical sense) a given sigma value called `pressure_only_above`,
typically in the upper atmosphere, pure sigma below `sigma_only_below` typically between the 
boundary layer and the mid-troposphere, and a smooth cubic interpolation in between.
The bar on each full level shows the A/B split (█ = sigma fraction, ░ = pressure fraction):

```@example vertical_coordinates
CubicSigmaPressureCoordinates(spectral_grid)
```

Custom thresholds can be set:

```@example vertical_coordinates
CubicSigmaPressureCoordinates(spectral_grid; pressure_only_above = 0.1, σ_only_below = 0.9)
```

Or custom sigma spacing too, to have, regardless sigma or sigma-pressure coordinates
more layers near the surface and the top of the atmosphere:

```@example vertical_coordinates
(; σ_half) = FriersonSigmaCoordinates(spectral_grid)
CubicSigmaPressureCoordinates(spectral_grid;  σ_half)
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

Functions provide a coordinate-agnostic interface for evaluating the vertical
coordinate at a given full level ``k``, working for both `SigmaCoordinates` and
`SigmaPressureCoordinates`:

| Function | Returns |
|---|---|
| `SpeedyWeather.pressure(k, pₛ, coord)` | Pressure [Pa] at full level ``k`` |
| `SpeedyWeather.pressure_thickness(k, pₛ, coord)` | Pressure thickness [Pa] of layer ``k`` |
| `SpeedyWeather.sigma(k, coord)` | Nominal sigma level ``\sigma_k`` (surface-pressure independent) |
| `sigma_half(k, coord)` | Nominal sigma level ``\sigma_{k-\tfrac{1}{2}}`` at half level ``k`` (surface-pressure independent) |

For `SigmaCoordinates` these reduce to ``\sigma_k p_s``, ``\Delta\sigma_k p_s``, and
``\sigma_k`` respectively. For `SigmaPressureCoordinates` they use the full
``A_k p_{\mathrm{ref}} + B_k p_s`` formula. Note that `sigma` and `sigma_half` always equal
``A_k + B_k = \sigma_k`` (respectively at half levels) regardless of the transition —
this is what `Geometry` stores as `σ_levels_full`/`σ_levels_half`/`σ_levels_thick`
(nominal sigma, ``p_s``-independent, **not** the ``B`` coefficients alone).

A few further internal (documented but not exported) accessors return the
``p_s``-dependent/-sensitivity quantities the dynamical core needs; see the next section.

## What the dynamical core does with hybrid coordinates

The guiding design principle is that `SigmaCoordinates` code paths stay byte-identical to
before hybrid coordinates were introduced (same arrays, same kernels), while
`SigmaPressureCoordinates` is reached through separate dispatch. This section explains what
changes for the hybrid path, why, and at what cost. See Simmons and Burridge (1981) for the
underlying continuous equations, generalised here to hybrid coordinates; the pure sigma case
is `SigmaCoordinates`' existing implementation.

### Three layer weights that collapse to ``\Delta\sigma_k``

The dynamical core carries ``\ln p_s`` as the prognostic variable and therefore works with
several per-layer quantities normalised by ``p_s``. For `SigmaCoordinates` all three below
are the same constant ``\Delta \sigma_k``, which is why the (unmodified) sigma code path can
use a single precomputed array throughout. For `SigmaPressureCoordinates` they are three
genuinely different quantities:

| weight | meaning | sigma limit | ``p_s``-dependent? |
|---|---|---|---|
| ``\delta_k(p_s) = \Delta p_k / p_s = \Delta A_k (p_{\mathrm{ref}}/p_s) + \Delta B_k`` | mass of layer ``k`` per unit ``p_s`` | ``\Delta\sigma_k`` | **yes** |
| ``\Delta B_k = \partial \Delta p_k / \partial p_s`` | sensitivity of layer mass to ``p_s`` | ``\Delta\sigma_k`` | no |
| ``B_{k+\tfrac{1}{2}} = \partial p_{k+\tfrac{1}{2}} / \partial p_s`` | sensitivity of interface pressure to ``p_s`` | ``\sigma_{k+\tfrac{1}{2}}`` | no |

`SpeedyWeather.pressure_thickness_ratio(k, pₛ, coord)` returns ``\delta_k(p_s)``,
`SpeedyWeather.pressure_thickness_sensitivity(k, coord)` returns ``\Delta B_k`` and
`SpeedyWeather.pressure_sensitivity_half(k, coord)` returns ``B_{k-\tfrac{1}{2}}``
(and `SpeedyWeather.pressure_sensitivity(k, coord)` returns ``B_k``). All four are `@inline`,
allocation-free and reduce exactly to ``\Delta\sigma_k``/``\sigma_k`` for `SigmaCoordinates`.

### Continuity, vertical velocity and adiabatic conversion in hybrid form

Using ``\mathcal{D}_k`` for the layer divergence, the exact mass-flux divergence normalised
by ``p_s`` is ``c_k = \delta_k(p_s)\mathcal{D}_k + \Delta B_k (\mathbf{u}_k \cdot \nabla \ln
p_s)``, so that ``\partial_t \ln p_s = -\sum_k c_k``. The (spectral, ``\Delta\sigma_k``-weighted)
vertical integral used for the semi-implicit terms differs from the exact, spatially-varying
``\delta_k``-weighted one by ``\hat{C} = \sum_k (\delta_k(p_s) - \Delta\sigma_k)\mathcal{D}_k``,
which SpeedyWeather.jl computes explicitly on the grid (`vars.dynamics.div_mean_correction`,
identically zero for `SigmaCoordinates`) and adds to the surface pressure tendency. The
vertical velocity `w` (what the sigma path calls ``\dot\sigma``) and the adiabatic conversion
coefficients generalise analogously — see the corresponding notes in
[Surface pressure tendency](@ref), [Vertical velocity](@ref) and [Temperature equation](@ref)
of the primitive equation model documentation for the exact formulas. All of these reduce
term by term to the existing sigma expressions when `transition = _ -> 1`.

### What stays linearised, and why that is fine

Two parts of the dynamical core are *not* generalised to the exact hybrid pressures, by
deliberate scope decision:

* **Geopotential integration constants.** The precomputed spectral constants
  ``\Delta\Phi_k = R \ln(\sigma_{k+\tfrac{1}{2}}/\sigma_k)`` are the sigma form of
  ``R\ln(p_{k+\tfrac{1}{2}}/p_k)``. In the pure-pressure limit (``B=0``) the ``p_{\mathrm{ref}}``
  factors cancel in the ratio and the two expressions are *identical*; in the pure-sigma limit
  (``A=0``) the ``p_s`` factors cancel and they are again *identical*. Only in the blend zone
  is there a discrepancy, of order ``\Delta A \Delta B / \sigma^2 \cdot (p_s/p_{\mathrm{ref}} -
  1)``. The exact, ``p_s``-dependent grid-space kernel is used for the grid geopotential (read
  by the parameterizations), so this only affects the spectral geopotential used inside the
  dynamical core itself.
* **Semi-implicit operators** (``\mathbf{R}, \mathbf{L}, \mathbf{U}, \mathbf{W}``, see
  [Semi-implicit time stepping](@ref implicit_primitive)) are the linearisation of the system
  about a resting reference state. Linearising the hybrid weights about ``p_s = p_{\mathrm{ref}}``
  gives ``\delta_k(p_{\mathrm{ref}}) = \Delta A_k + \Delta B_k = \Delta\sigma_k`` and
  ``p_k(p_{\mathrm{ref}}) = \sigma_k p_{\mathrm{ref}}``. So with nominal ``\sigma`` in
  `σ_levels_*` the existing sigma-coordinate operators *are* the correct linearisation of the
  hybrid system — no reformulation is needed. As with every other nonlinearity in the
  semi-implicit scheme, the residual between the linearised and the exact hybrid operator is
  handled explicitly, though the correction is less well conditioned the further ``p_s``
  departs from ``p_{\mathrm{ref}}``.

### Performance

The hybrid path costs a handful of extra `log`/`exp` evaluations per grid point and layer
(vertical velocity, adiabatic conversion, the grid geopotential) plus one extra grid-space
correction term for the surface pressure tendency; `SigmaCoordinates` performance is
unchanged by construction, since its code paths are untouched (same arrays, same kernels, no
extra transcendental function calls).
