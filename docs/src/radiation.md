# Radiation

## Longwave radiation implementations

Currently implemented is

```@example radiation
using InteractiveUtils # hide
using SpeedyWeather
subtypes(SpeedyWeather.AbstractLongwave)
```

## Uniform cooling

Following Paulius and Garner[^PG06], the uniform cooling of the atmosphere
is defined as

```math
\frac{\partial T}{\partial t} = \begin{cases} - \tau^{-1}\quad&\text{for}\quad T > T_{min} \\
                                            \frac{T_{strat} - T}{\tau_{strat}} \quad &\text{else.} \end{cases}
```

with ``\tau = 16~h`` resulting in a cooling of -1.5K/day for most of the atmosphere,
except below temperatures of ``T_{min} = 207.5~K`` in the stratosphere where a
relaxation towards ``T_{strat} = 200~K`` with a time scale of ``\tau_{strat} = 5~days``
is present.

## Jeevanjee radiation

Jeevanjee and Zhou [^JZ22] (eq. 2) define a longwave radiative flux ``F`` for atmospheric cooling
as (following Seeley and Wordsworth [^SW23], eq. 1)

```math
\frac{dF}{dT} = α*(T_t - T)
```

The flux ``F`` (in ``W/m^2/K``) is a vertical upward flux between two layers (vertically adjacent)
of temperature difference ``dT``. The change of this flux across layers depends on the temperature
``T`` and is a relaxation term towards a prescribed stratospheric temperature ``T_t = 200~K`` with
a radiative forcing constant ``\alpha = 0.025 W/m^2/K^2``. Two layers of identical temperatures
``T_1 = T_2`` would have no net flux between them, but a layer below at higher temperature would
flux into colder layers above as long as its temperature ``T > T_t``. This flux is applied above
the lowermost layer and above, leaving the surface fluxes unchanged. The uppermost layer is tied
to ``T_t`` through a relaxation at time scale ``\tau = 6~h``

```math
\frac{\partial T}{\partial t} = \frac{T_t - T}{\tau}
```

The flux ``F`` is converted to temperature tendencies at layer ``k`` via

```math
\frac{\partial T_k}{\partial t} = (F_{k+1/2} - F_{k-1/2})\frac{g}{\Delta p c_p}
```

The term in parentheses is the absorbed flux in layer ``k`` of the upward
flux from below at interface ``k+1/2`` (``k`` increases downwards, see
[Vertical coordinates and resolution](@ref) and [Sigma coordinates](@ref)).
``\Delta p = p_{k+1/2} - p_{k-1/2}`` is the pressure thickness of layer ``k``,
gravity ``g`` and heat capacity ``c_p``.

## Shortwave radiation

Currently implemented schemes:

```@example radiation
using SpeedyWeather
subtypes(SpeedyWeather.AbstractShortwave)
```

### OneBandShortwave: Single-band shortwave radiation with diagnostic clouds

*Note: Stratocumulus cloud parameterization and ozone-radiation interactions are currently not implemented in the `OneBandShortwave` scheme.*

The `OneBandShortwave` scheme provides a single-band (broadband) shortwave radiation parameterization, including diagnostic cloud effects following [^KMB06].

**Cloud diagnosis:**
Cloud properties are diagnosed from the relative humidity and total precipitation in the atmospheric column. The cloud base is set at the interface between the lowest two model layers, and the cloud top is the highest layer where both

$$
\mathrm{RH}_k > \mathrm{RH}_{cl} \quad \text{and} \quad Q_k > Q_{cl}
$$

are satisfied. The cloud cover (CLC) in a layer is then given by

$$
\mathrm{CLC} = \min\left[1,\ w_{pcl} \sqrt{\min(p_{mcl}, P_{lsc} + P_{cnv})}+ \min\left(1, \left(\frac{\mathrm{RH}_k - \mathrm{RH}_{cl}}{\mathrm{RH}'_{cl} - \mathrm{RH}_{cl}}\right)^2\right)\right]
$$

where $w_{pcl}$ and $p_{mcl}$ are parameters, $P_{lsc}$ and $P_{cnv}$ are large-scale and convective precipitation, and $\mathrm{RH}_{cl}$ is a threshold.

**Stratocumulus clouds:**
Stratocumulus cloud cover over oceans is parameterized based on boundary layer static stability (GSEN):

$$
\mathrm{CLS} = F_{ST} \max(\mathrm{CLS}_{\max} - \mathrm{CLC}, 0)
$$

with

$$
F_{ST} = \max\left(0, \min\left(1, \frac{\mathrm{GSE}_N - \mathrm{GSES}_0}{\mathrm{GSES}_1 - \mathrm{GSES}_0}\right)\right)
$$

Over land, the stratocumulus cover $\mathrm{CLS}$ is further modified to be proportional to the surface relative humidity:

$$
\mathrm{CLS}_L = \mathrm{CLS} \cdot \mathrm{RH}_N
$$
where $\mathrm{RH}_N$ is the surface (lowest model layer) relative humidity.

**Radiative transfer:**
The incoming solar flux at the top of the atmosphere is computed from astronomical formulae. Ozone absorption in the lower and upper stratosphere is subtracted, yielding the downward flux into the first model layer:

$$
F_{h}^{\downarrow, SR} = F_{0}^{\downarrow, sol} - \Delta F_{ust}^{ozone} - \Delta F_{lst}^{ozone}
$$

Shortwave radiation is then propagated downward through each layer using a transmissivity $\tau_{k}^{SR}$, which depends on zenith angle, layer depth, humidity, and cloud properties:

$$
F_{k+h}^{\downarrow, SR} = F_{k-h}^{\downarrow, SR} \, \tau_k^{SR}
$$

In the cloud-top layer, cloud reflection is included:

$$
F_{k+h}^{\downarrow, SR} = F_{k-h}^{\downarrow, SR} (1 - A_{cl} \, \mathrm{CLC}) \, \tau_{k}^{SR}
$$

At the surface, stratocumulus reflection and surface albedo are applied:

$$
F_{N+h}^{\downarrow, SR} = F_{N-h}^{\downarrow, SR} (1 - A_{cls} \, \mathrm{CLS}) \, \tau_{N}^{SR}
$$

The upward flux at the surface is

$$
F_{s}^{\uparrow, SR} = F_{s}^{\downarrow, SR} \, A_s
$$

and is propagated upward as

$$
F_{k-h}^{\uparrow, SR} =  F_{k+h}^{\uparrow, SR} \, \tau_{k}^{SR}
$$

with cloud reflection added at the cloud-top layer. The upward part is only modeled for the visible band, as near-infrared is mostly absorbed downward.

**Key features:**

- Single-band (broadband) shortwave radiative transfer
- Cloud reflection and absorption parameterized using diagnosed cloud cover
- Surface albedo can vary between land and ocean
- Top-of-atmosphere insolation set by the solar constant and zenith angle
- Designed for use in idealized and moist atmospheric simulations

#### Usage

To use the new scheme, construct your model with:

```julia
model = PrimitiveWetModel(spectral_grid; shortwave_radiation=OneBandShortwave(spectral_grid))
```

and run as usual. The scheme will output surface and top-of-atmosphere shortwave fluxes, as well as cloud-modified radiative tendencies.

## References

[^PG06]: Paulius and Garner, 2006. JAS. DOI:[10.1175/JAS3705.1](https://doi.org/10.1175/JAS3705.1)

[^SW23]: Seeley, J. T. & Wordsworth, R. D. Moist Convection Is Most Vigorous at Intermediate Atmospheric Humidity. Planet. Sci. J. 4, 34 (2023). DOI:[10.3847/PSJ/acb0cb](https://doi.org/10.3847/PSJ/acb0cb)

[^JZ22]: Jeevanjee, N. & Zhou, L. On the Resolution‐Dependence of Anvil Cloud Fraction and Precipitation Efficiency in Radiative‐Convective Equilibrium. J Adv Model Earth Syst 14, e2021MS002759 (2022). DOI:[10.1029/2021MS002759](https://doi.org/10.1029/2021MS002759)

[^KMB06]: Kucharski, F., Molteni, F., & Bracco, A. SPEEDY: A simplified atmospheric general circulation model. ICTP, Trieste, Italy. Appendix A: Model Equations and Parameters (2006). [PDF](https://users.ictp.it/~kucharsk/speedy_description/km_ver41_appendixA.pdf)
