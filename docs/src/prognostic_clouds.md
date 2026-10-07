# Prognostic clouds

With [`ImplicitCondensation`](@ref) (see [Large-scale condensation](@ref)) humidity above a
threshold rains out in the same time step and column, and clouds are diagnosed for radiation
from relative humidity and precipitation. [`PrognosticCloudCondensation`](@ref) instead carries
the condensed water as a prognostic variable: it is advected, diffused, converted to rain and
snow over time, and handed to radiation as cloud water and ice in every layer. Cloud radiative
effects then follow from the model's condensate history, and clouds can be advected away from
where they formed.

The scheme follows the structure of the Zhao and Carr (1997) scheme that NCEP's spectral GFS
used from 2001 to 2019: one condensate for liquid and ice with the phase decided by
temperature, autoconversion after Sundqvist et al. (1989), precipitation diagnosed in one
top-down sweep, and a Xu and Randall (1996) cloud fraction for radiation. Its condensation
closure is replaced by the implicit relaxation of [`ImplicitCondensation`](@ref), so it is
not a published scheme as such.

!!! warning "Not tuned yet"
    The parameters are GFS and SpeedyWeather starting values. Cloud cover, water paths and the
    cloud radiative effects have not been tuned against observations.

## Usage

The condensation scheme replaces the default `large_scale_condensation`. To let radiation see
the clouds, use [`PrognosticClouds`](@ref) for the one-band shortwave radiation and
[`CloudyLongwaveTransmissivity`](@ref) for the one-band longwave radiation

```julia
using SpeedyWeather
spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)

large_scale_condensation = PrognosticCloudCondensation(spectral_grid)
shortwave = OneBandShortwave(spectral_grid; clouds = PrognosticClouds(spectral_grid))
longwave = OneBandLongwave(spectral_grid; transmissivity = CloudyLongwaveTransmissivity(spectral_grid))
radiation = Radiation(spectral_grid; shortwave, longwave)

model = PrimitiveWetModel(spectral_grid; large_scale_condensation, radiation)
add!(model, SpeedyWeather.CloudOutput()...)     # condensate, cloud fraction, water paths
simulation = initialize!(model)
run!(simulation, period = Day(10), output = true)
```

The condensate is `simulation.variables.prognostic.cloud_condensate` (spectral) and
`simulation.variables.grid.cloud_condensate` (grid), the cloud state is in
`simulation.variables.parameterizations` (`cloud_fraction`, `cloud_liquid_water`,
`cloud_ice_water`, `cloud_liquid_effective_radius`, `cloud_ice_effective_radius`,
`liquid_water_path`, `ice_water_path`).

## The condensate variable

The cloud condensate ``q_c`` [kg/kg] is a prognostic variable like humidity, declared by the
scheme and fused with the other atmospheric variables (see [Variable system](@ref)): its
spectral-to-grid transform and the grid-to-spectral transforms of its tendency and of the fluxes
``u q_c``, ``v q_c`` are part of the batched transforms that already exist, so on GPU it adds
work but no kernel launches. It is advected in flux form, vertically advected and
horizontally diffused exactly as humidity is.

Spectral transport makes an intermittent field like ``q_c`` slightly negative at cloud edges.
The grid copy is not clipped. Instead the scheme fills negative condensate from vapour in the
same cell with the corresponding latent heating (as NCEP's GFS does), which conserves water and
enthalpy.

## Processes in a column

The scheme reads temperature ``T``, humidity ``q`` and condensate ``q_c`` of the previous time
step (as all parameterizations do) and works from the top of the column down. All sinks of
``q_c`` are bounded with the time step ``\Delta t_p`` the state is advanced with, ``2\Delta t``
for leapfrog, so that the condensate stays non-negative.

**Phase.** The ice fraction of the condensate is a linear ramp from ``f_i = 0`` at freezing
``T_0`` to ``f_i = 1`` at `ice_temperature` (``-20˚C`` by default). The latent heat between
vapour and condensate is ``L = L_v + f_i L_f`` and the saturation humidity is a blend of
saturation over liquid and over ice (Clausius-Clapeyron with ``L_v + L_f``)

```math
q^\star = (1 - f_i) q^\star_l + f_i q^\star_i
```

**Condensation and evaporation.** Humidity above the threshold ``r q^\star`` condenses into
the condensate with the implicit relaxation of [`ImplicitCondensation`](@ref)

```math
C = \frac{\max(q - r q^\star, 0)}{\tau \left[ 1 + \frac{L r}{c_p} \frac{\partial q^\star}{\partial T} \right]}
```

with ``\tau`` = `time_scale` ``\times \Delta t``. Below the threshold, condensate evaporates
in the clear part ``1 - a`` of the cell (``a`` the cloud fraction) with its own, longer time
scale ``\tau_e`` = `evaporation_time_scale` ``\times \Delta t``, but not more than available

```math
E = \min \left( (1 - a) \frac{\max(r q^\star - q, 0)}{\tau_e \left[ 1 + \frac{L r}{c_p} \frac{\partial q^\star}{\partial T} \right]}, \frac{q_c}{\Delta t_p} \right)
```

**Cloud fraction.** After condensation the cloud fraction follows Xu and Randall (1996) in the
form and with the constants of NCEP's GFS

```math
a = \mathrm{RH}^{1/4} \left[ 1 - \exp \left( -\frac{\alpha q_c}{((1 - \mathrm{RH}) q^\star)^{1/4}} \right) \right]
```

with ``\alpha = 2000``, the denominator clamped to ``[10^{-4}, 1]``, no cloud for
``q_c < 10^{-6} p/1000~\mathrm{hPa}`` and cloud fractions below 0.001 set to zero. There is
no cloud without condensate.

**Autoconversion.** The liquid part ``(1-f_i) q_c`` above a minimum ``q_{min}`` converts to
rain with the Sundqvist rate

```math
k_l = C_0 F \left[ 1 - \exp \left( -\left( \frac{(q_l - q_{min}) F}{m_r a} \right)^2 \right) \right], \quad
F = (1 + c_1 \sqrt{P}) (1 + c_2 \sqrt{\min(\max(268~\mathrm{K} - T, 0), 20~\mathrm{K})})
```

enhanced by collection of the precipitation ``P`` [kg/m²/s] falling into the layer and by the
Bergeron-Findeisen process in supercooled cloud. The ice part converts to snow with
``k_i = c_i \exp(0.025 (T - T_0))`` (Lin et al. 1983). Both are applied exponentially in time,
``\Delta q = (q - q_{min}) (1 - \exp(-k \Delta t_p))``, so that they never remove more than
available, however large the rates.

**Precipitation.** Rain and snow are diagnosed: the fluxes from above melt (snow, energy-limited
above `melting_threshold`) and reevaporate (rain, proportional to the subsaturation, not beyond
saturation) as in [`ImplicitCondensation`](@ref), then this layer's autoconversion is added and
falls into the layer below. There is no storage and no fall speed.

**Conservation.** Per column, vapour and condensate lost equal the precipitation reaching the
surface, and the heating equals the latent heat of the vapour lost plus the latent heat of
fusion of the frozen water gained (ice in the condensate at its ice fraction, snow at the
surface). Both are tested to round-off. When the condensate changes its ice fraction because it
is advected to a different temperature, that phase change carries no latent heat, as in Zhao
and Carr (1997).

## Clouds in radiation

The scheme writes per layer the cloud fraction, the grid-mean cloud liquid and ice water and
their effective radii (10 µm over ocean, 5 to 10 µm over land increasing with the ice fraction,
50 µm for ice). Radiation schemes read this cloud state and nothing else from the cloud scheme;
without a prognostic cloud scheme it is zero and there are no clouds.

[`PrognosticClouds`](@ref) provides clouds for the one-band shortwave radiation: the column cloud
cover from the layer cloud fractions with maximum-random overlap, the cloud top as the highest
cloudy layer, and the cloud albedo from the in-cloud optical depth of the column

```math
\tau = \frac{3}{2} \left( \frac{\mathrm{LWP}}{\rho_w r_l} + \frac{\mathrm{IWP}}{\rho_i r_i} \right) / a_{col}, \quad
R = \frac{\beta \tau}{1 + \beta \tau}, \quad \beta = \frac{\sqrt{3}}{2} (1 - g)
```

the two-stream reflectivity of a non-absorbing layer with asymmetry factor ``g = 0.85``.

[`CloudyLongwaveTransmissivity`](@ref) multiplies the transmissivity of a clear-sky scheme
(by default `FriersonLongwaveTransmissivity`) in every layer by ``1 - a \varepsilon`` with
the emissivity of the cloudy part ``\varepsilon = 1 - \exp(-D(\kappa_l W_l + \kappa_i W_i))`` from the
in-cloud liquid and ice water paths ``W``, ``D = 1.66`` and mass absorption coefficients after
the NCAR CCM3 (Kiehl et al. 1998).

## Known limitations

- Not tuned: in a first 10-day run at T31 the global mean cloud cover is about 0.26 and the
  liquid water path about 6 g/m², far below observed values.
- No convective condensate: Betts-Miller convection detrains no cloud water, and
  [`PrognosticClouds`](@ref) has no separate convective or stratocumulus cloud.
- The shortwave reflects at the cloud top only; cloud absorption is that of the one-band scheme.
- The NumericalRadiation ecCKD scheme is clear-sky and does not use the cloud state yet.

## References

- Sundqvist, H., E. Berge, J. E. Kristjánsson (1989): Condensation and cloud parameterization
  studies with a mesoscale numerical weather prediction model. *Mon. Wea. Rev.*, 117, 1641–1657.
- Zhao, Q., F. H. Carr (1997): A prognostic cloud scheme for operational NWP models.
  *Mon. Wea. Rev.*, 125, 1931–1953.
- Xu, K.-M., D. A. Randall (1996): A semiempirical cloudiness parameterization for use in
  climate models. *J. Atmos. Sci.*, 53, 3084–3102.
- Lin, Y.-L., R. D. Farley, H. D. Orville (1983): Bulk parameterization of the snow field in a
  cloud model. *J. Climate Appl. Meteor.*, 22, 1065–1092.
- Kiehl, J. T., et al. (1998): The National Center for Atmospheric Research Community Climate
  Model: CCM3. *J. Climate*, 11, 1131–1149.
