# Surface fluxes

The surfaces fluxes in SpeedyWeather represent the exchange of momentum, heat,
and moisture between ocean and land as surface into the lowermost atmospheric
layer. Surface fluxes of momentum represent a drag that the boundary layer
wind experiences due to friction over more or less rough ground on land
or over sea. Surface fluxes of heat represent a sensible heat flux from 
a warmer or colder ocean or land into or out of the surface layer of the
atmosphere. Surface fluxes of moisture represent evaporation of sea water
into undersaturated surface air masses or, similarly, evaporation from
land with a given soil moisture and vegetation's evapotranspiration.

## Surface flux implementations

Currently implemented surface fluxes of momentum are

```@example surface_fluxes
using InteractiveUtils # hide
using SpeedyWeather
subtypes(SpeedyWeather.AbstractSurfaceWind)
```

!!! note "Interdependence of surface flux computations"
    `SurfaceWind` computes the surface fluxes of momentum but also the computation
    of the surface wind (which by default includes wind gusts) meaning that `NoSurfaceWind`
    will also effectively disable other surface fluxes unless custom surface fluxes
    have been implemented that do not rely on `column.surface_wind_speed`.

with more explanation below. The surface heat fluxes currently implemented are

```@example surface_fluxes
subtypes(SpeedyWeather.AbstractSurfaceHeatFlux)
```

and the surface moisture fluxes, i.e. evaporation (this does not include [Convection](@ref) or
[Large-scale condensation](@ref) which currently immediately removes humidity instead of
fluxing it out at the bottom) implemented are

```@example surface_fluxes
subtypes(SpeedyWeather.AbstractSurfaceEvaporation)
```

The calculation of thermodynamic quantities at the surface (air density, temperature, humidity)
are handled by

```@example surface_fluxes
subtypes(SpeedyWeather.AbstractSurfaceThermodynamics)
```

and the computation of drag coefficients (which is used by default for the surface fluxes above)
is handled through the `model.boundary_layer` where currently implemented are

```@example surface_fluxes
subtypes(SpeedyWeather.AbstractBoundaryLayer)
```

Note that `LinearDrag` is the linear drag following Held and Suarez (see [Held-Suarez forcing](@ref))
which does not compute a drag coefficient and therefore by default effectively disables other
surface fluxes (as the Held and Suarez forcing and drag is supposed to be used instead of
physical parameterizations).

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
with gravity ``g`` and layer-thickness ``\Delta p_k`` (see [Sigma coordinates](@ref)) so that
the right-hand side divides the absorbed flux by the mass of layer ``k`` (per unit area).
Tendencies for ``v, q`` equivalently with their respective absorbed fluxes.

The temperature tendency is calculated from the absorbed heat flux as

```math
\frac{\partial T_k}{\partial t} = \frac{g \Delta F_k}{c_p \Delta p_k}
```

with heat capacity of air at constant pressure ``c_p`` with a typical value of ``1004~J/kg/K``.
Because we define the heat flux as having units of ``W/m^2`` the conversion includes
the division by the heat capacity to convert to a rate of temperature change.

## Bulk Richardson-based drag coefficient

All surface fluxes depend on a dimensionless drag coefficient ``C`` which
we calculate as a function of the bulk Richardson number ``Ri`` following
Frierson, et al. 2006 [^Frierson2006] with some simplification as outlined below.
We use the same drag coefficient for momentum, heat and moisture fluxes.
The bulk Richardson number at the lowermost model layer ``k = N`` of height ``z_N`` is

```math
Ri_N = \frac{gz_N \left( \Theta_v(z_N) - \Theta_v(0) \right)}{|v(z_N)|^2 \Theta_v(0)}
```

with ``gz_N = \Phi_N`` the [Geopotential](@ref) at ``z = z_N``, ``\Theta =  c_pT_v + gz``
the virtual dry static energy and ``T_v`` the [Virtual temperature](@ref).
Then the drag coefficient ``C`` follows as

```math
C = \begin{cases}
        \kappa^2 \left( \log(\frac{z_N}{z_0}) \right)^{-2} \quad &\text{for} \quad Ri_N  \leq 0\\
        \kappa^2 \left( \log(\frac{z_N}{z_0}) \right)^{-2} \left(1 - \frac{Ri_N}{Ri_c}\right)^2 \quad &\text{for} \quad 0 < Ri_N < Ri_c \\
        0 \quad &\text{for} \quad Ri_N \geq Ri_c. \\
    \end{cases}
```

with ``\kappa = 0.4`` the von Kármán constant, ``z_0 = 3.21 \cdot 10^{-5}`` the 
roughness length. There is a maximum drag ``C`` for negative bulk Richardson numbers ``Ri_N``
but the drag becomes 0 for bulk Richardson numbers being larger than a critical
``Ri_c = 1`` with a smooth transition in between. The height of the ``N``-th
model layer is ``z_N = \tfrac{\Phi_N - \Phi_0}{g}`` with the geopotential
```math
\Phi_N = \Phi_{0} + T_N R_d ( \log p_{N+h} - \log p_N)
```
which depends on the temperature ``T_N`` of that layer. To simplify this calculation
and avoid the logarithm we use a constant ``Z \approx z_N`` from a reference temperature
``T_{ref}`` instead of ``T_N`` in the calculation of ``log(z_N/z_0)``.
While ``z_N`` depends on the vertical resolution (``\Delta p_N`` of the lowermost layer)
it is proportional to that layer's temperature which in Kelvin does not change much
in space and in time, especially with the logarithm applied.
For ``T_{ref} = 200K`` with specific gas constant ``R_d = 287 J/kg/K``
on sigma level ``\sigma_N = 0.95`` we have ``log(z_N/z_0) \approx 16.1`` for
``T_{ref} = 300K`` this increases to ``log(z_N/z_0) \approx 16.5`` a 2.5%
increase which we are happy to approximate. Note that we do
not apply this approximation within the bulk Richardson number ``Ri_N``.
So we calculate once a typical height of the lowermost layer
``Z = T_{ref}R_d \log(1/\sigma_N)g^{-1}`` for the given parameter choices
and then define a maximum drag constant

```math
C_{max} = \left(\frac{\kappa}{\log(\frac{Z}{z_0})} \right)^2
```

to simplify the drag coefficient calculation to

```math
C = C_{max} \left(1 - \frac{Ri_N^*}{Ri_c}\right)^2
```
with ``Ri_N^* = \max(0, \min(Ri_N, Ri_c))`` the clamped ``Ri_N`` which
is at least ``0`` and at most ``Ri_c``.

## Surface momentum fluxes

The surface momentum flux is calculated from the surface wind velocities
```math
u_s = f_w u_N, \quad v_s = f_w v_N
```

meaning it is scaled down by ``f_w = 0.95`` (Fortran SPEEDY default, [^SPEEDY])
from the lowermost layer wind velocities ``u_N, v_N``. A wind speed scale
accounting for gustiness with ``V_{gust} = 5~m/s`` (Fortran SPEEDY default, [^SPEEDY])
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
assuming adiabatic descent but that is currently not implemented.
In practice we use a flux limiter for numerical stability which limits the magnitude
(preserving the sign). Choosing ``F_{uv, max}^\uparrow = 0.5 Pa`` would mean
that the drag for winds faster than about ``33~m/s`` (typical ``C = 5\cdot 10^{-4}``)
does not further increase. This can help to prevent oscillations drag terms can produce
for sudden strong wind gusts.

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

and apply a similar flux limiter as for the momentum flux to prevent a sudden strong heating
or cooling. For ``F_{T, max}^\uparrow = 100~W/m^2`` 

## Surface evaporation

The surface moisture flux, i.e. evaporation of soil moisture over land and evaporation of
sea water over the ocean is proportional to the difference of the surface specific humidity
``q_s`` and the saturation specific humidity given ``T_{skin}`` and surface pressure ``p_s``.
This assumes that a very thin layer of air just above the ocean is saturated but over land
this assumption is less well justified as it should be a function of the soil moisture
and how much of that is available to evaporate given vegetation. We again make the
simplification that ``q_s = q_N``, i.e. specific humidity of the surface is the
same as in the lowermost atmospheric layer above. 

The surface evaporative flux is then (always positive)

```math
F_q^\uparrow = \rho_s C V_0 \max(0, \alpha_{sw} q^* - q_s)
```

with ``q^*`` the saturation specific humidity calculated from the skin temperature ``T_{skin}``
and surface pressure ``p_s``. The available of soil water over land is represented by
(over the ocean ``\alpha_{sw} = 1``)

```math
\alpha_{sw} = \frac{D_{top} W_{top} + f_{veg} D_{root} \max(0, W_{root} - W_{wil})}{
    D_{top}W_{cap} + D_{root}(W_{cap} - W_{wil})}
```

following the Fortran SPEEDY documentation[^SPEEDY] which follows Viterbo and Beljiars 1995
[^Viterbo95]. The variables (or spatially prescribed arrays) are water content in the top
soil layer ``W_{top}`` and the root layer below ``W_{root}`` using the vegetation
fraction ``f_{veg} = veg_{high} + 0.8 veg_{low}`` composed of a (dimensionless)
high and low vegetation cover per grid cell ``veg_{high}, veg_{low}``.
The constants are depth of top soil layer ``D_{top} = 7~cm``, depth of root layer
``D_{root} = 21~cm``, soil wetness at field capacity (volume fraction)
``W_{cap} = 0.3``, and soil wetness at wilting point (volume fraction) ``W_{wil} = 0.17``.

## Land-sea mask

SpeedyWeather uses a _fractional_ land-sea mask, i.e. for every grid-point

- 1 indicates land
- 0 indicates ocean
- a value in between indicates a grid-cell partially covered by ocean and land

The land-sea mask determines solely how to weight the surface fluxes
coming from land or from the ocean. For the sensible heat fluxes this uses
land and sea surface temperatures and weights the respective fluxes
proportional to the fractional mask. Similar for evaporation.
You can therefore define an ocean on top of a mountain, or a land without
heat fluxes when the land-surface temperature is not defined, i.e. `NaN`.
Let ``F_L, F_S`` be the fluxes coming from land and sea, respectively.
Then the land-sea mask ``a \in [0,1]`` weights them as follows for the
total flux ``F``

```math
F = aF_L + (1-a)F_S
```

but ``F=F_L`` if the sea flux is NaN (because the ocean temperature is not defined)
and ``F=F_S`` if the land flux is NaN (because the land temperature or soil moisture
is not defined, for sensible heat fluxes or evaporation), and ``F=0`` if both fluxes
are NaN.

Setting the land-sea mask to ocean therefore will disable any fluxes that
may come from land, and vice versa. However, with an ocean-everywhere land-sea mask
you must also define sea surface temperatures everywhere, otherwise the fluxes
in those regions will be zero.

For more details see [The land-sea mask](@ref) implementation section.

## References

[^Frierson2006]: Frierson, D. M. W., I. M. Held, and P. Zurita-Gotor, 2006: A Gray-Radiation Aquaplanet Moist GCM. Part I: Static Stability and Eddy Scale. J. Atmos. Sci., 63, 2548-2566. DOI: [10.1175/JAS3753.1](https://doi.org/10.1175/JAS3753.1). 

[^SPEEDY]: Franco Molteni and Fred Kucharski, 20??. Description of the ICTP AGCM (SPEEDY) Version 41. [https://users.ictp.it/~kucharsk/speedy_description/km_ver41_appendixA.pdf](https://users.ictp.it/~kucharsk/speedy_description/km_ver41_appendixA.pdf)

[^Viterbo95]: Viterbo, P., and A. C. M. Beljaars, 1995: An Improved Land Surface Parameterization Scheme in the ECMWF Model and Its Validation. J. Climate, 8, 2716-2748, DOI:[10.1175/1520-0442(1995)008<2716:AILSPS>2.0.CO;2](https://doi.org/10.1175/1520-0442(1995)008<2716:AILSPS>2.0.CO;2). 