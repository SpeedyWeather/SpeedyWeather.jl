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
The upward flux ``F^u`` out of layer ``k+1`` into layer ``k`` (vertical indexing ``k`` increases downwards)
is ``F^u_{k+h}`` (``h = \frac{1}{2}`` for half, as the flux sits on the cell face, the half-layer
in between ``k`` and ``k+1``) and similarly ``F^d_{k+h}`` is the downward flux. While for the sum of
all fluxes we have ``\sum_i F^d_{i, k+h} = - \sum_i F^u_{i, k+h}`` for clarity we may define
fluxes as either upward or downward depending on the process.
The absorbed flux in layer ``k`` is

```math
\Delta F_k = (F^u_{k+h} - F^u_{k-h}) + (F^d_{k-h} - F^d_{k+h})
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

The surface flux computation depends on a dimensionless drag coefficient ``C`` which
we calculate as a function of the bulk Richardson number ``Ri`` following
Frierson, et al. 2006 [^Frierson2006]. We use the same drag
coefficient for momentum, heat and moisture fluxes. The bulk Richardson number at
the lowermost model layer ``k = N`` of height ``z_N`` is

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

There's a maximum drag ``C`` for negative bulk Richardson numbers ``Ri_N``
but the drag becomes 0 for bulk Richardson numbers being larger than a critical
``Ri_c = 1`` with a smooth transition in between. To simplify this calculation
and avoid the logarithm we use a constant ``Z \approx z_N`` from a reference temperature
because while ``z_N`` depends on the vertical resolution (``\Delta p_N`` of the lowermost layer)
it is proportional to that layer's temperature which in Kelvin does not change much,
so we do ``Z = \tfrac{\Phi_n - \Phi_0}{g}`` with

```math
\Phi_n = \Phi_{0} + T_{ref}R_d ( \ln p_s - \ln p_n)
```

## Surface momentum fluxes



## Surface heat fluxes

## Surface evaporation

## References

[^Frierson2006]: Frierson, D. M. W., I. M. Held, and P. Zurita-Gotor, 2006: A Gray-Radiation Aquaplanet Moist GCM. Part I: Static Stability and Eddy Scale. J. Atmos. Sci., 63, 2548–2566. DOI: [10.1175/JAS3753.1](https://doi.org/10.1175/JAS3753.1). 