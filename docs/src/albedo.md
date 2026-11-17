# Albedo

Albedo is the surface reflectivity to downward solar shortwave radiation.
A value of 1 indicates that all of the radiative flux is reflected at
the Earth's surface and sent back up through the atmospheric column.
In contrast, a value of 0 means no reflection and all of that radiative
flux is absorbed, typically heating ocean or land surface.

In SpeedyWeather the albedo is its own model component, `model.albedo`, that is computed
for every grid column on every time step before the [Radiation](@ref) so that the shortwave
radiation can use it. Albedo is therefore not a static boundary condition but recomputed each step
so that transient surface changes (e.g. snowfall, sea ice, wind roughening the ocean surface)
brighten or darken the surface without losing the underlying bare albedo.

## Ocean and land albedo

Albedo is computed separately for ocean and land as averaging those in grid cells which are
partially land, partially ocean would yield inaccurate results. Think 10% land having a lower
heat capacity than ocean but being treated with an albedo that comes from 90% ocean.
Not very realistic. The default albedo is therefore a container `OceanLandAlbedo` that
holds one albedo for the ocean and one for the land

```@example albedo
using SpeedyWeather
spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
albedo = OceanLandAlbedo(spectral_grid)
```

with `JinOceanAlbedo` used over the ocean and `LearnedLandAlbedo` used over land by default,
see below. They write into `vars.parameterizations.ocean.albedo` and `vars.parameterizations.land.albedo`,
respectively, and the shortwave radiation computes the reflected fluxes over ocean and land separately.
The shortwave radiation also stores the land-sea mask-weighted albedo in `vars.parameterizations.albedo`.

## Implemented albedos

The following albedos are currently implemented

```@example albedo
using InteractiveUtils # hide
subtypes(SpeedyWeather.AbstractAlbedo)
```

Roughly sorted by increasing complexity these are

| Albedo | Intended for | Depends on |
| :----- | :----------- | :--------- |
| `GlobalConstantAlbedo` | ocean, land or both | nothing, constant |
| `ManualAlbedo` | ocean, land or both | nothing, 2D field set by the user |
| `PrescribedAlbedo` | ocean, land or both | nothing, set externally (e.g. coupling) |
| `AlbedoClimatology` | ocean, land or both | nothing, 2D field read from file |
| `OceanSeaIceAlbedo` | ocean | sea ice concentration |
| `LandSnowAlbedo` | land | vegetation cover, snow depth |
| `JinOceanAlbedo` (default) | ocean | solar zenith angle, surface wind speed, sea ice concentration, cloud cover |
| `LearnedLandAlbedo` (default) | land | solar zenith angle, vegetation cover and leaf area, soil moisture and temperature, orography, snow, surface air temperature and wind speed, cloud cover |

and `OceanLandAlbedo` is the composite that combines one albedo for the ocean with one for the land.
Each albedo is briefly described in the following.

### GlobalConstantAlbedo

A globally constant albedo (default `albedo = 0.3`).
To be used for land and ocean or only one of them within an `OceanLandAlbedo`.

### ManualAlbedo

A manual albedo field, to be used with `set!`, that is copied into the diagnostic variables on every time step.
This way parameterizations can change the albedo at every time step (e.g. snow cover) without
losing the information of the original surface albedo. See [Customizing the albedo](@ref) for an example.

### PrescribedAlbedo

Like `ManualAlbedo` but the albedo is not held in the struct and not copied into the diagnostic variables
automatically. It's expected that the albedo variable (e.g `vars.parameterizations.land.albedo`)
is set externally at every time step, e.g. when coupling SpeedyWeather to another model.

### AlbedoClimatology

An albedo climatology loaded from a netCDF file (by default `albedo.nc` from
[SpeedyWeatherAssets](https://github.com/SpeedyWeather/SpeedyWeatherAssets)) and interpolated onto
the model grid at initialization. It then behaves like a `ManualAlbedo`, i.e. it is constant in time.

### OceanSeaIceAlbedo

Albedo that scales linearly between ocean and ice albedo depending on sea ice concentration ``ℵ``

```math
α = α_{ocean} + ℵ (α_{ice} - α_{ocean})
```

with ``α_{ocean} = 0.06`` (`albedo_ocean`) and ``α_{ice} = 0.6`` (`albedo_ice`) by default,
see also [Sea ice models](@ref).

### LandSnowAlbedo

Albedo over land based on bare soil, high and low vegetation and snow cover

```math
α = v_h α_h + v_l α_l + (1 - v_h - v_l) α_{land} + σ_s α_{snow}
```

with high and low vegetation cover ``v_h, v_l`` from the land model and their albedos
``α_h = 0.15`` (`albedo_high_vegetation`) and ``α_l = 0.2`` (`albedo_low_vegetation`),
bare land albedo ``α_{land} = 0.4`` (`albedo_land`) and an additive snow albedo
``α_{snow} = 0.4`` (`albedo_snow`). Without vegetation the albedo is that of bare land.
The snow cover ``σ_s`` is diagnosed from the snow depth ``S`` of the land model's snow bucket
(see [Snow model](@ref)) by either the `LinearSnowCover` ``σ_s = \min(S / S_0, 1)`` or the default
`SaturatingSnowCover` ``σ_s = S / (S + S_0)`` with `snow_depth_scale` ``S_0 = 0.05~m``.
The snow therefore brightens areas with low vegetation cover more than forests, which somewhat
reflects that in forests the snow cover is broken up and snow lies in between trees.

### JinOceanAlbedo

Ocean albedo following [Jin2011](@citet), the default ocean albedo. It depends on the
cosine of the solar zenith angle ``μ`` and the surface wind speed ``w``, which roughens the
ocean surface with a slope variance ``σ^2 = 0.003 + 0.00512 w`` [CoxMunk1954](@citep).
The albedo is the sum of the reflection at the air-water interface (surface scattering)
and the light scattered back from below the surface (volume scattering), with direct (black-sky)
and diffuse (white-sky) albedos computed separately, see [Black-sky, white-sky and blue-sky albedo](@ref).
It is composed of the following components

- `ocean_direct_reflection = OceanDirectReflection()`: the Fresnel reflectance ``r(n, μ)`` of a flat
  air-water interface with refractive index ``n = 1.34``, corrected for the roughness of a windy, wavy
  surface with a regression ``f(μ, σ)`` from [Jin2011](@citet), eq. (1),
  ``α_s^{dir} = r(n, μ) - \frac{r(n, μ)}{r(n_0, μ)} f(μ, σ)``.
- `ocean_diffuse_reflection = OceanDiffuseReflection()`: a linear regression in ``n`` and ``σ``,
  [Jin2011](@citet) eq. (5b), for the diffuse radiation under clouds, ``α_s^{dif} = c_0 + c_1 n + c_2 n σ``.
- `ocean_volume = ConstantOceanVolumeScattering()`: a constant broadband volume-scattering
  albedo ``α_w = 0.006`` for clear (case 1) waters, [Jin2011](@citet) eq. (15).
- `foam = WhitecapFoam()`: whitecaps cover a fraction ``f = \min(2.95 \cdot 10^{-6} w^{3.52}, 1)``
  of the ocean surface with an albedo of ``α_{foam} = 0.55`` [Koepke1984](@citep), such that
  ``α \to f α_{foam} + (1 - f) α``. Use `foam = NoFoam()` to switch the whitecaps off.

The black-sky albedo is then ``α_s^{dir} + α_w``, the white-sky albedo ``α_s^{dif} + α_w``,
both corrected for whitecaps. Sea ice is blended in linearly with its concentration
as in `OceanSeaIceAlbedo` with ``α_{ice} = 0.6`` (`albedo_ice`). Without a direct beam
(the sun below the horizon) the black-sky albedo is set to the white-sky albedo.
At low sun angles the Fresnel reflection makes the ocean much brighter than
the ``α = 0.06`` of `OceanSeaIceAlbedo`.

### LearnedLandAlbedo

Learned albedo over land based on surface properties and some ambient conditions,
the default land albedo. It uses the semi-empirical, kernel-driven bidirectional reflectance
distribution function (BRDF) model of [Lucht2000](@citet), which describes the reflectance of the
land surface as the sum of three scattering mechanisms with weights

- ``f_{iso}``, isotropic scattering: light at the surface is reflected equally in all directions,
  regardless of where the sun is,
- ``f_{vol}``, volumetric scattering: light interacting with a horizontally homogeneous vegetation
  canopy, it meets the canopy, bounces around a bit beneath the canopy and exits back towards the sky,
- ``f_{geo}``, geometric-optical scattering: shadow casting of 3D objects at the surface,
  like tree crowns or big rocks, and obscuration of shadows due to the viewing angle.

These weights are predicted for the visible and near-infrared separately, each by an expression that
was derived by distilling a trained neural network using symbolic regression (Munday et al., in prep.),
so six expressions in total, one for each `AbstractLandScattering` component of `LearnedLandAlbedo`.
The predictors are the high and low vegetation cover and their leaf area index, the soil moisture and
soil temperature in the top and bottom soil layers, the orography (as geopotential), the snow depth
and snow cover (via `LinearSnowCover` or `SaturatingSnowCover` as for `LandSnowAlbedo`), the surface
air temperature and the surface wind speed. Predictors that are not available, e.g. soil moisture
in a `PrimitiveDryModel`, are set to zero. Visible and near-infrared weights are combined into
broadband weights ``f = w_{vis} f^{vis} + (1 - w_{vis}) f^{nir}`` with ``w_{vis} = 0.5143`` (`vis_weight`).
From these the white-sky albedo and the black-sky albedo, depending on the solar zenith angle ``θ``,
are computed following [Lucht2000](@citet)

```math
\begin{aligned}
α_{ws} &= f_{iso} + 0.189184 f_{vol} - 1.377622 f_{geo} \\
α_{bs}(θ) &= f_{iso} + (-0.007574 - 0.070987 θ^2 + 0.307588 θ^3) f_{vol}
    + (-1.284909 - 0.166314 θ^2 + 0.04184 θ^3) f_{geo}
\end{aligned}
```

All coefficients of the learned expressions are parameters of the model, see `parameters(LearnedLandAlbedo(spectral_grid))`.

## Black-sky, white-sky and blue-sky albedo

The albedo of a surface depends on the direction of the incoming light.
The black-sky albedo ``α_{bs}`` is the albedo for direct radiation only, i.e. a sky without
any scattering ("black") and therefore depends on the solar zenith angle.
The white-sky albedo ``α_{ws}`` is the albedo for completely diffuse radiation, i.e. coming
from all directions of a uniformly bright ("white") sky, and is independent of the solar zenith angle.
The actual, "blue-sky", albedo is a blend of the two with the fraction of direct radiation ``f_{dir}``

```math
α = f_{dir} α_{bs} + (1 - f_{dir}) α_{ws}
```

`JinOceanAlbedo` and `LearnedLandAlbedo` compute both black-sky and white-sky albedo,
stored as `black_sky_albedo` and `white_sky_albedo` in the `ocean` or `land` namespace of
`vars.parameterizations`. The shortwave radiation then blends them as above with
``f_{dir} = 1 - \text{cloud cover}`` (`OneBandShortwave`, `TwoBandShortwave`) or
``f_{dir} = 1`` (`TransparentShortwave`). Shortwave radiation schemes that do not
know the fraction of direct radiation use the white-sky albedo. All other albedos
are independent of the direction of the incoming light and do not define
black-sky and white-sky albedos, nothing is blended in this case.

## Scattering and other components

`JinOceanAlbedo`, `LearnedLandAlbedo` and `LandSnowAlbedo` are composed of
components that compute parts of the albedo. These are not albedos themselves, i.e. they cannot be
used as `ocean` or `land` in `OceanLandAlbedo`, but they can be exchanged within their albedo,
see [Customizing the albedo](@ref). The scattering components for land

```@example albedo
subtypes(SpeedyWeather.AbstractLandScattering)
```

are the six learned expressions of `LearnedLandAlbedo`, see above.
The ocean has scattering components for the surface reflection

```@example albedo
subtypes(SpeedyWeather.AbstractOceanSurfaceScattering)
```

used as `ocean_direct_reflection` and `ocean_diffuse_reflection` in `JinOceanAlbedo`, and for
volume scattering, used as `ocean_volume`

```@example albedo
subtypes(SpeedyWeather.AbstractOceanVolumeScattering)
```

`ConstantOceanVolumeScattering` is described above. `SpectralOceanVolumeScattering` computes the
volume scattering per wavelength band from the backscattering and absorption of clear waters depending
on the chlorophyll concentration [Jin2011, MorelMaritorena2001](@citep). This is future work and cannot
be used yet as it requires a shortwave radiation scheme with a surface albedo per band,
`JinOceanAlbedo` throws an error if constructed with it.
The whitecaps of `JinOceanAlbedo` are its `foam` component

```@example albedo
subtypes(SpeedyWeather.AbstractFoam)
```

and the snow cover of `LandSnowAlbedo` and `LearnedLandAlbedo` is diagnosed from snow depth
with one of

```@example albedo
subtypes(SpeedyWeather.AbstractSnowCover)
```

## Customizing the albedo

All albedos and their components are constructed with the `spectral_grid` as first argument
and their parameters as keyword arguments. For example, to switch off the whitecaps
and to use a brighter sea ice in `JinOceanAlbedo` and to use the simpler `LandSnowAlbedo`
with a linear snow cover over land

```@example albedo
ocean = JinOceanAlbedo(spectral_grid, foam = NoFoam(), albedo_ice = 0.7)
land = LandSnowAlbedo(spectral_grid, snow_cover = LinearSnowCover(), albedo_snow = 0.5)
albedo = OceanLandAlbedo(spectral_grid; ocean, land)
```

The first argument of `OceanLandAlbedo` (without the `spectral_grid`) is used for ocean the second for land
but you can use keywords too as above. You can mix any albedos, they are internally
two independent albedos that are applied to fluxes separately, e.g.

```@example albedo
albedo = OceanLandAlbedo(GlobalConstantAlbedo(spectral_grid), AlbedoClimatology(spectral_grid))
```

constructs an albedo that is a global constant (default 0.3) for the ocean
but the `AlbedoClimatology` read from file used for the land.
Alternatively you can also drop the separation into ocean/land albedo
(e.g. idealised simulations or aqua planet, rocky planet) and just
use

```@example albedo
albedo = GlobalConstantAlbedo(spectral_grid)
```

and this definition of the albedo will be used for both ocean and land fluxes.
In all cases you can then pass on the albedo to the model constructor, e.g.

```@example albedo
albedo = OceanLandAlbedo(GlobalConstantAlbedo(spectral_grid, albedo = 0.1), ManualAlbedo(spectral_grid))
set!(albedo.land, (λ, φ) -> 0.2 + 0.3 * abs(φ) / 90)

model = PrimitiveWetModel(spectral_grid; albedo)
model.albedo
```

The albedo in the `model` is now the one defined just in the lines above,
using a globally constant albedo of 0.1 for the ocean but a higher albedo
over land which also increases to 0.5 towards the poles.

You can always output the land-sea mask weighted albedo with
`add!(model, OutputVariables.Albedo())` or inspect it as follows

```@example albedo
simulation = initialize!(model)
run!(simulation, steps = 1)   # run for a step to "diagnose" albedo = ocean/land weighted

using CairoMakie
(; albedo) = simulation.variables.parameterizations
heatmap(albedo, title = "Custom albedo, separately defined for ocean/land")
save("ocean_land_albedo.png", ans) # hide
nothing # hide
```
![Ocean-land albedo](ocean_land_albedo.png)

In comparison, the default albedo (`JinOceanAlbedo` over the ocean and `LearnedLandAlbedo` over land)
after one day of simulation

```@example albedo
model = PrimitiveWetModel(spectral_grid)
simulation = initialize!(model)
run!(simulation, period = Day(1))

(; albedo) = simulation.variables.parameterizations
heatmap(albedo, title = "Default albedo after 1 day")
save("default_albedo.png", ans) # hide
nothing # hide
```
![Default albedo](default_albedo.png)

## Defining a new albedo

A new albedo is defined as `MyAlbedo <: SpeedyWeather.AbstractAlbedo` and implements

- `initialize!(albedo::MyAlbedo, model)`
- `albedo!(ij, albedo_field, vars::Variables, albedo::MyAlbedo, model)`

`initialize!` is as for every model component but in contrast to other parameterizations albedos
should not implement `parameterization!` but rather `albedo!`, which is called from within `parameterization!`.
This is because every albedo needs to be used either as part of `OceanLandAlbedo` or the same albedo
is called over both land and ocean surfaces. If a surface namespace (`vars.parameterizations.land` or `.ocean`)
contains both `black_sky_albedo` and `white_sky_albedo`, a shortwave scheme that knows the direct-beam fraction
blends them into that namespace's albedo. Otherwise it leaves albedo alone.
See [Example: Albedo](@ref) for a step-by-step example.
