# Surface fluxes

The surfaces fluxes in SpeedyWeather represent the exchange of momentum, heat,
and humidity/moisture between ocean and land as surface into or out of
the lowermost atmospheric layer (see [Vertical coordinates](@ref vertical_coordinates_physics) for how vertical
layers are indexed). Surface fluxes of momentum represent
a drag that the boundary layer wind experiences due to friction over more
or less rough ground on land or over sea. Surface fluxes of heat represent
a sensible heat flux from  a warmer or colder ocean or land into or out of
the surface layer of the atmosphere. Surface fluxes of humidity represent
evaporation of sea water or condensation of humidity onto the ocean surface.
Similarly, humidity fluxes over land with a given soil moisture
and vegetation's evapotranspiration.

## Surface flux implementations

Currently implemented surface fluxes of momentum are

```@example surface_fluxes
using InteractiveUtils # hide
using SpeedyWeather
subtypes(SpeedyWeather.AbstractSurfaceMomentumFlux)
```

!!! note "Interdependence of surface flux computations"
    `model.boundary_layer.surface_condition` extrapolates atmospheric variables (wind, temperature, humidity)
    to the surface. `model.boundary_layer` computes centrally a drag coefficient for surface
    momentum, heat and humidity fluxes. Setting `model.boundary_layer.surface_condition = nothing` or
    `model.boundary_layer = nothing` will therefore disable all other surface fluxes unless
    they have `use_boundary_layer_drag = false` which lets them use independently their
    own drag coefficient.

with more explanation below. The surface heat fluxes currently implemented are

```@example surface_fluxes
subtypes(SpeedyWeather.AbstractSurfaceHeatFlux)
```

and the surface humidity fluxes (this does not include [Convection](@ref) or
[Large-scale condensation](@ref) those are precipitation) implemented are

```@example surface_fluxes
subtypes(SpeedyWeather.AbstractSurfaceHumidityFlux)
```

The calculation of thermodynamic quantities at the surface (air density, temperature, humidity) 
and the computation of drag coefficients (which is used by default for the surface fluxes above)
are handled through the `model.boundary_layer` where currently implemented are

```@example surface_fluxes
subtypes(SpeedyWeather.AbstractBoundaryLayer)
```

See [Parameterizations](@ref) for the general parameterization interface these implement.

## Fluxes to tendencies

In SpeedyWeather.jl, parameterizations can be defined either in terms of tendencies for a given
layer or as fluxes between two layers including the surface flux and a top-of-the-atmosphere flux.
The upward flux ``F^\uparrow`` out of layer ``k+1`` into layer ``k`` (vertical indexing ``k`` increases
downwards) is ``F^\uparrow_{k+h}`` (``h = \frac{1}{2}`` for half, as the flux sits on the cell face,
the half-layer in between ``k`` and ``k+1``) and similarly ``F^\downarrow_{k+h}`` is the downward flux.
For clarity we may define fluxes as either upward or downward depending on the process although an
upward flux can always be regarded as a negative downward flux.
The absorbed flux in layer ``k`` is

```math
\Delta F_k = (F^\uparrow_{k+h} - F^\uparrow_{k-h}) + (F^\downarrow_{k-h} - F^\downarrow_{k+h})
```

A quick overview of the units

| Quantity | Variable | Unit | Flux unit |
| -------- | -------- | ---- | --------- |
| Momentum | Velocity ``u, v`` | ``m/s``  | ``Pa = kg/m/s^2`` |
| Heat | Temperature ``T`` | ``m/s`` | ``W/m^2 = kg/s^3`` |
| Moisture | Specific humidity ``q`` | ``kg/kg`` | ``kg/m^2/s`` |

The time-stepping in SpeedyWeather.jl (see [Time integration](@ref leapfrog)) eventually requires
tendencies which are calculated from the absorbed fluxes of momentum ``u`` or ``v`` and moisture ``q`` as

```math
\frac{\partial u_k}{\partial t} = \frac{g \Delta F_k}{\Delta p_k}
```
with gravity ``g`` and layer-thickness ``\Delta p_k`` (see [Sigma coordinates](@ref sigma_coordinates_physics)) so that
the right-hand side divides the absorbed flux by the mass of layer ``k`` (per unit area).
Tendencies for ``v, q`` equivalently with their respective absorbed fluxes.

The temperature tendency is calculated from the absorbed heat flux as

```math
\frac{\partial T_k}{\partial t} = \frac{g \Delta F_k}{c_p \Delta p_k}
```

with heat capacity of air at constant pressure ``c_p`` with a typical value of ``1004~J/kg/K``.
Because we define the heat flux as having units of ``W/m^2`` the conversion includes
the division by the heat capacity to convert to a rate of temperature change.

## [Time integration of surface fluxes](@id surface_flux_time_integration)

### Surface fluxes as a relaxation

The bulk surface fluxes of momentum, heat and humidity below all relax a variable ``x``
of the lowermost layer ``N`` (``u, v, T`` or ``q``) towards a surface value ``x_s``
(``0`` for momentum, the skin or sea surface temperature for heat, the saturation
humidity at the surface for humidity). Each of them can be written as

```math
F^\uparrow = -a (x - x_s), \qquad a = \rho_s C V_0
```

with the conductance ``a = -\partial F^\uparrow / \partial x`` in ``kg/m^2/s``.
``a`` is a mass exchange rate: the mass of air per square metre and second that effectively
exchanges its momentum, heat or moisture with the surface. ``a`` contains the surface air
density ``\rho_s``, the (dimensionless) drag coefficient ``C`` and the surface wind speed ``V_0``,
it is additionally scaled by ``f_w`` for momentum, and divided by the sea-ice or snow insulation
factors for heat and humidity. All fluxes are positive upward.

| Flux | ``x`` | ``x_s`` | ``F^\uparrow`` units | ``a`` |
| ---- | ----- | ------- | -------------------- | ----- |
| Momentum | ``u, v`` [``m/s``] | ``0`` | ``kg/m/s^2 = Pa`` | ``\rho_s C V_0 f_w`` |
| Heat | ``T`` [``K``] | ``T_{skin}`` or SST | ``K\,kg/m^2/s``, see below | ``\rho_s C V_0`` / insulation |
| Humidity | ``q`` [``kg/kg``] | ``\alpha q^*(T_{skin})`` or ``q^*(SST)`` | ``kg/m^2/s`` | ``\rho_s C V_0`` / insulation |

with ``q^*`` the saturation humidity and ``\alpha`` the soil moisture availability over land.
Note that ``\alpha`` only scales the surface value ``x_s``, not ``x``, so it does not enter ``a``.

A flux ``F^\uparrow`` through the bottom of layer ``N`` with pressure thickness ``\Delta p_N``
changes the content of that layer. The layer has a mass of ``\Delta p_N / g`` per square metre,
so its momentum or moisture content is ``(\Delta p_N / g) x`` per square metre and therefore

```math
\frac{\partial x}{\partial t} = \frac{g F^\uparrow}{\Delta p_N} = -\frac{g a}{\Delta p_N}(x - x_s) = -\frac{x - x_s}{\tau},
\qquad \tau = \frac{\Delta p_N}{g a}
```

The relaxation time scale ``\tau`` is the time it takes to exchange the mass of the lowermost
layer ``\Delta p_N/g`` with the surface at a rate ``a``, which is the same for momentum, heat
and humidity if they use the same drag coefficient. ``\tau`` is short for a thin lowermost layer
(small ``\Delta p_N``), strong winds or high drag. With `FriersonSigmaCoordinates` for example,
the lowermost layer is about 3 times thinner than with equi-distant ``\sigma`` levels,
which can reduce ``\tau`` below the time step.

### Heat fluxes and the heat capacity

The sensible heat flux is defined in ``W/m^2`` as ``F_T^\uparrow = -c_p a (T - T_s)``,
so that ``-\partial F_T^\uparrow / \partial T = c_p a`` in ``W/m^2/K``, not ``a``. However,
the heat content of the lowermost layer is ``c_p (\Delta p_N / g) T`` per square metre,
so the temperature tendency is

```math
\frac{\partial T}{\partial t} = \frac{g F_T^\uparrow}{c_p \Delta p_N} = -\frac{g a}{\Delta p_N}(T - T_s)
```

and ``c_p`` cancels: it appears both in the flux and in the heat content of the layer.
The relaxation time scale ``\tau = \Delta p_N/(g a)`` is therefore the same as above and
does not depend on ``c_p``. In the code the heat flux is accordingly calculated without ``c_p``,
i.e. as ``F^\uparrow = -a (T - T_s)`` in ``K\,kg/m^2/s``, for the time integration and
only multiplied by ``c_p`` when stored in ``W/m^2``.

### Exact solution and flux damping

Holding ``a`` and ``x_s`` fixed over a time step ``\Delta t`` (``2\Delta t`` for leapfrog,
as the tendency is applied from ``t - \Delta t`` to ``t + \Delta t``), the relaxation equation
has the exact solution of an exponential decay towards the surface value

```math
x(t + \Delta t) - x_s = e^{-z} (x(t) - x_s), \qquad z = \frac{\Delta t}{\tau}
```

We can write this exact change as an explicit step with a modified flux. With
``F^\uparrow`` the flux evaluated with ``x(t)`` at the beginning of the time step, the
explicit step changes ``x`` by ``-z (x - x_s)``, while the exact solution changes it by
``-(1 - e^{-z}) (x - x_s)``. Hence the exact solution is obtained by multiplying the
explicit flux with a damping factor ``\varphi(z)``

```math
F^\uparrow \to \varphi(z) F^\uparrow, \qquad x^{new} - x_s = (1 - z \varphi(z))(x - x_s), \qquad
\varphi(z) = \frac{1 - e^{-z}}{z}
```

This is exponential time differencing. ``\varphi`` only depends on ``z``, so the damped flux is
calculated directly from the explicit one without iterating or solving a system.
``\varphi = 1 - z/2 + z^2/6 - \dots`` for small ``z``, and ``\varphi \approx 1/z`` for large ``z``,
meaning ``x`` is relaxed fully to ``x_s`` within one time step.

However, evaluating ``\varphi`` requires an exponential for every grid point, every surface flux
and every time step. This is relatively expensive: at T32L8 using `ExponentialFluxes` for all
surface fluxes increases the time per time step by about 5%, compared to about 1% for
`ImplicitFluxes`. We therefore look for cheaper approximations of ``\varphi`` that keep its
essential properties: ``\varphi \approx 1`` for small ``z`` (the flux is unchanged if
the time step resolves the relaxation) and ``\varphi \approx 1/z`` for large ``z``
(``x`` relaxes to ``x_s`` but does not overshoot it). The other time integration schemes
correspond to such approximations of the exact ``\varphi``:

- **Explicit**, ``\varphi = 1``, only approximates the exact solution to first order in ``z``.
  ``1 - z`` replaces ``e^{-z}``, which overshoots the surface value for ``z > 1``
  and is unstable (``|1 - z| > 1``) for ``z > 2``.
- **Implicit** (backward Euler) evaluates the flux with ``x`` at the end of the time step,
  ``x^{new} - x_s = (x - x_s) - z (x^{new} - x_s)``, so that ``\varphi = 1/(1 + z)``.
  This is unconditionally stable and monotone, and has the correct ``1/z`` behaviour for large ``z``,
  but is only first-order accurate, ``\varphi \approx 1 - z``, so it overdamps for small ``z``.
- **Rational**, ``\varphi = (6 + z)/(6 + 4z + z^2)``, matches the expansion of the exact
  ``\varphi`` to second order, ``1 - z/2 + z^2/6``, and has the same ``1/z`` behaviour for large ``z``,
  at the cost of the implicit scheme (one division) instead of an `exp`. ``1 - z\varphi(z)``
  becomes negative for ``z > 3`` with a minimum of about ``-0.1`` near ``z \approx 8``,
  meaning that ``x`` overshoots the surface value by at most 10% of the initial deviation,
  then tends to zero for large ``z``.

| Scheme | ``\varphi(z)`` | ``1 - z\varphi(z)`` | Single step ``x \to x_s`` | Cost |
| ------ | -------------- | ------------------- | ------------------------- | ---- |
| `ExplicitFluxes()` | ``1`` | ``1 - z`` | overshoots for ``z > 1``, unstable for ``z > 2`` | none |
| `ImplicitFluxes()` | ``\frac{1}{1 + z}`` | ``\frac{1}{1 + z}`` | monotone, stable, first order | 1 division |
| `ExponentialFluxes()` | ``\frac{1 - e^{-z}}{z}`` | ``e^{-z}`` | monotone, stable, exact for fixed ``a, x_s`` | `exp`, ~5% of a time step at T32L8 |
| `RationalFluxes()` (default) | ``\frac{6 + z}{6 + 4z + z^2}`` | ``\frac{6 - 2z}{6 + 4z + z^2}`` | stable, ≤10% overshoot for ``z > 3`` | 1 division |

Because ``\varphi`` is applied to the flux (not only to the tendency) before it is stored, the fluxes
passed on to the ocean and land models and to the output are consistent with those the atmosphere
receives, conserving momentum, energy and water between the components. For the land and ocean
fluxes, ``a`` is weighted by the land or ocean fraction respectively, as only that fraction of the flux
enters the lowermost layer of the atmosphere. The drag coefficient ``C``, the wind speed ``V_0`` and
the surface value ``x_s`` are held fixed over the time step (linearisation), which makes ``a`` and
``\tau`` constant within a time step. Prescribed fluxes (e.g. `PrescribedOceanHeatFlux`) are applied as given.

### Choosing the time integration

The scheme is chosen per flux, for example

```julia
spectral_grid = SpectralGrid(truncation = 31, nlayers = 8)
surface_momentum_flux = SurfaceMomentumFlux(spectral_grid, time_integration = ExponentialFluxes())
surface_heat_flux = SurfaceHeatFlux(
    spectral_grid,
    ocean = SurfaceOceanHeatFlux(spectral_grid, time_integration = ImplicitFluxes()),
    land = SurfaceLandHeatFlux(spectral_grid, time_integration = ExplicitFluxes()),
)
model = PrimitiveWetModel(spectral_grid; surface_momentum_flux, surface_heat_flux)
```

and similarly for `SurfaceOceanHumidityFlux` and `SurfaceLandHumidityFlux`.

## Bulk Richardson-based drag coefficient

All surface fluxes depend on a dimensionless drag coefficient ``C`` which
we calculate as a function of the bulk Richardson number ``Ri``
(see also [Vertical diffusion](@ref) for a similar Richardson-number-based diffusion coefficient)
following [Frierson2006](@citet) with some simplification as outlined below.
We use the same drag coefficient for momentum, heat and moisture fluxes.
The bulk Richardson number at the lowermost model layer ``k = N`` of height ``z_N`` is

```math
Ri_N = \frac{gz_N \left[ \Theta_v(z_N) - \Theta_v(0) \right]}{|v(z_N)|^2 \Theta_v(0)}
```

with ``gz_N = \Phi_N`` the [Geopotential](@ref) at ``z = z_N``, ``\Theta =  c_pT_v + gz``
the virtual dry static energy and ``T_v`` the [Virtual temperature](@ref).
Then the drag coefficient ``C`` follows as

```math
C = \begin{cases}
        \kappa^2 \left[ \log(\frac{z_N}{z_0}) \right]^{-2} \quad &\text{for} \quad Ri_N  \leq 0\\
        \kappa^2 \left[ \log(\frac{z_N}{z_0}) \right]^{-2} \left(1 - \frac{Ri_N}{Ri_c}\right)^2 \quad &\text{for} \quad 0 < Ri_N < Ri_c \\
        0 \quad &\text{for} \quad Ri_N \geq Ri_c, \\
    \end{cases}
```

with ``\kappa = 0.4`` the von Kármán constant, ``z_0`` the
roughness length (we currently use different values over ocean and land).
There is a maximum drag ``C`` for negative bulk Richardson numbers ``Ri_N``
but the drag becomes 0 for bulk Richardson numbers being larger than a critical
``Ri_c = 10`` with a smooth transition in between. The height of the ``N``-th
model layer is ``z_N = \tfrac{\Phi_N - \Phi_0}{g}`` with the geopotential
```math
\Phi_N = \Phi_{0} + T_N R_d ( \log p_{N+h} - \log p_N)
```
which depends on the temperature ``T_N`` of that layer. This model layer
height above the surface ``z_N`` is recalculated on every time step for
every grid cell as it depends on temperature. The maximum drag coefficient
then follows as

```math
C_{max} = \left[\frac{\kappa}{\log(\frac{z}{z_0})} \right]^2
```

to simplify the drag coefficient calculation to

```math
C = C_{max} \left[1 - \frac{Ri_N^*}{Ri_c}\right]^2
```
with ``Ri_N^* = \max(0, \min(Ri_N, Ri_c))`` the clamped ``Ri_N`` which
is at least ``0`` and at most ``Ri_c``.

## Surface momentum fluxes

The surface momentum flux is calculated from the surface wind velocities
```math
u_s = f_w u_N, \quad v_s = f_w v_N
```

meaning it is scaled down by ``f_w = 0.95`` (Fortran SPEEDY default;
[KucharskiMolteniBracco2006](@citep))
from the lowermost layer wind velocities ``u_N, v_N``. A wind speed scale
accounting for gustiness with ``V_{gust} = 5~m/s`` (Fortran SPEEDY default;
[KucharskiMolteniBracco2006](@citep))
is then defined as

```math
V_0 = \sqrt{u_s^2 + v_s^2 + V_{gust}^2}
```

such that for low wind speeds the fluxes are somewhat higher because of
unresolved winds on smaller time and length scales. The momentum flux is then

```math
\begin{aligned}
F_u^\uparrow &= - \rho_s C V_0 u_s \\
F_v^\uparrow &= - \rho_s C V_0 v_s
\end{aligned}
```

with ``\rho_s = \frac{p_s}{R_d T_N}`` the surface air density calculated from
surface pressure ``p_s`` and lowermost layer temperature ``T_N``.
Better would be to extrapolate ``T_N`` to ``T_s`` a surface air temperature
assuming adiabatic descent but that is currently not implemented. How this flux is integrated in time is described in [Time integration of surface fluxes](@ref surface_flux_time_integration).

## Surface heat fluxes

The surface heat flux is proportional to the difference of the surface air temperature ``T_s``
and the land or ocean skin temperature ``T_{skin}``.
We currently just approximate as ``T_N`` the lowermost layer temperature although
an adiabatic descent from pressure ``p_N`` to surface pressure ``p_s`` would be more accurate.
We also use land and sea surface temperature to approximate ``T_{skin}`` although
future improvements should account for faster radiative effects on ``T_{skin}`` compared
to sea and land surface temperatures determined by a higher heat capacity of the relevant
land surface layer or the mixed layer in the ocean. We then compute

```math
F_T^\uparrow = \rho_s C V_0 c_p (T_{skin} - T_s)
```

in ``W/m^2``, positive up. How this flux is integrated in time is described in [Time integration of surface fluxes](@ref surface_flux_time_integration).

## Surface humidity fluxes

The surface humidity flux, i.e. evaporation/condensation of soil moisture over land and
evaporation of sea water over the ocean or condensation of humidity on the ocean surface
is proportional to the difference of the surface specific humidity ``q_s`` and the saturation
specific humidity given ``T_{skin}`` and surface pressure ``p_s``.
This assumes that a very thin layer of air just above the ocean is saturated but over land this
assumption is less well justified as it should be a function of the soil moisture and how much
of that is available to evaporate given vegetation, hence we include a soil moisture availability.
We again make the simplification that ``q_s = q_N``, i.e. specific humidity of the surface is the
same as in the lowermost atmospheric layer above.

The surface humidity flux in ``kg/s/m^2``, positive up, is then

```math
F_q^\uparrow = \rho_s C V_0 (\alpha_{sw} q^* - q_s)
```

with ``q^*`` the saturation specific humidity calculated from the skin temperature ``T_{skin}``
and surface pressure ``p_s``. The availability of soil water over land is represented by
(over the ocean ``\alpha_{sw} = 1``)

```math
\alpha_{sw} = \frac{D_{top} W_{top} + f_{veg} D_{root} \max(0, W_{root} - W_{wil})}{
    D_{top}W_{cap} + D_{root}(W_{cap} - W_{wil})}
```

following the Fortran SPEEDY documentation [KucharskiMolteniBracco2006](@citep),
which follows [ViterboBeljaars1995](@citet). The variables (or spatially prescribed arrays) are water content in the top
soil layer ``W_{top}`` and the root layer below ``W_{root}``
(see [Land soil moisture](@ref) for how these are modelled) using the vegetation
fraction ``f_{veg} = veg_{high} + 0.8 veg_{low}`` composed of a (dimensionless)
high and low vegetation cover per grid cell ``veg_{high}, veg_{low}``.
The constants are depth of top soil layer ``D_{top} = 20~cm``, depth of root layer
``D_{root} = 2~m`` (defaults subject to change, check `model.land.geometry`),
soil wetness at field capacity (volume fraction)
``W_{cap} = 0.3``, and soil wetness at wilting point (volume fraction) ``W_{wil} = 0.17``. How this flux is integrated in time is described in [Time integration of surface fluxes](@ref surface_flux_time_integration).

## Land-sea mask

SpeedyWeather uses a _fractional_ land-sea mask, i.e. for every grid-point

- 1 indicates land
- 0 indicates ocean
- a value in between indicates a grid-cell partially covered by ocean and land

The land-sea mask determines solely how to weight the surface fluxes
coming from land or from the ocean. For the sensible heat fluxes this uses
land and sea surface temperatures and weights the respective fluxes
proportional to the fractional mask. Similar for humidity fluxes.
You can therefore define an ocean on top of a mountain, or vice versa.
Let ``F_L, F_S`` be the fluxes coming from land and sea, respectively.
Then the land-sea mask ``a \in [0,1]`` weights them as follows for the
total flux ``F``

```math
F = aF_L + (1-a)F_S
```

Note that this weighting is purely multiplicative: it does not protect
against `NaN` propagating into ``F`` should ``F_L`` or ``F_S`` be `NaN`.
Land and sea surface temperatures (and other quantities entering the flux
calculation, such as soil moisture) must therefore be finite everywhere,
even where the land-sea mask is 0 or 1, i.e. a `NaN` land or sea surface
temperature is not supported anymore. Use a fallback temperature for
grid cells without "real" land or ocean instead (e.g. `SlabOcean`'s
`land_temperature`, or `SeasonalLandTemperature`'s `ocean_temperature`).

Setting the land-sea mask to ocean therefore will disable any fluxes that
may come from land, and vice versa. However, with an ocean-everywhere land-sea mask
you must still also define sea surface temperatures everywhere, otherwise the fluxes
in those regions will be incorrect rather than simply zero.

For more details see [The land-sea mask](@ref) implementation section.
