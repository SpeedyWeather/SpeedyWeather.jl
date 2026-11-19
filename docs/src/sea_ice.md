# Sea ice models

The following describes the currently implemented sea ice models, which are

```@example sea_ice
using InteractiveUtils # hide
using SpeedyWeather
subtypes(SpeedyWeather.AbstractSeaIce)
```

## No sea ice

You can set `sea_ice = nothing` which will not touch the initial
`variables.prognostic.ocean.sea_ice_concentration` (which when allocated is zero),
nor advance it. But you may use `set!(simulation, sea_ice_concentration=...)` to set the
sea ice concentration manually. No sea ice does not do anything on every time step
so your manual modifications will prevail after `initialize!`.

To be used like

```@example sea_ice
spectral_grid = SpectralGrid(truncation=32)
model = PrimitiveWetModel(spectral_grid; sea_ice=nothing)
model.sea_ice
```

Note that in SpeedyWeather the sea ice model and its [albedo](@ref "Albedo") are defined
independently, means you can have a sea ice model without affecting the albedo
and `sea_ice=nothing` but set the sea ice concentration manually and use its albedo
effect, this will be discussed in `ThermodynamicSeaIce` below. The sea ice concentration also
lowers the surface longwave emission of the ocean, with any sea ice model and also with
`sea_ice=nothing`, see [Radiation](@ref).

## Thermodynamic sea ice model

Thermodynamic sea ice models do not account for advection through wind and currents or rheological processes
but are only concerned with local freezing, ice growth and melt and the thermodynamic processes involved.
The `ThermodynamicSeaIce` model defined here has one prognostic variable, the sea ice volume per unit ocean area
``V`` [m], stored in `variables.prognostic.ocean.sea_ice_volume`. The sea ice concentration ``ℵ`` (the area fraction
of the ocean covered by ice, `sea_ice_concentration`) and the ice thickness ``h`` follow from it.

### Thickness and concentration

The ice thickness is tied to the concentration

```math
h(ℵ) = h_{thin} + (h_{pack} - h_{thin}) ℵ
```
with ``h_{thin} = 0.2~\text{m}`` (`thin_ice_thickness`) and ``h_{pack} = 2.5~\text{m}`` (`pack_ice_thickness`), such that
a low concentration represents new, thin ice and full cover thick pack ice. The volume is ``V = ℵ h(ℵ)``,
which is inverted exactly as

```math
ℵ(V) = \frac{2V}{h_{thin} + \sqrt{h_{thin}^2 + 4 (h_{pack} - h_{thin}) V}}
```
for ``0 \leq V \leq h_{pack}``. Beyond full cover, ``V > h_{pack}``, the concentration stays 1 and the ice thickens,
``h = h(1) + V - h_{pack} = V``.

### Ice surface temperature and conduction

Heat is conducted through the ice following a zero-layer model (Semtner, 1976): the conductive flux
from the ice base at the freezing temperature ``T_f`` to the ice surface at ``T_i`` is

```math
F_c = \frac{k}{h} (T_f - T_i)
```
with the conductivity ``k = 2~\text{W/m/K}`` (`conductivity`). The ice surface temperature follows from the ice
surface energy balance

```math
ϵσT_i^4 = (1 - α_i) R_{sd} + R_{ld} + F_c
```
with ``ϵ`` the ocean emissivity of the longwave scheme, ``α_i`` the albedo of sea ice (`albedo_ice` of the ocean albedo),
``R_{sd}`` and ``R_{ld}`` the downward shortwave and longwave radiation at the surface. It is solved with
four Newton steps in the longwave radiation and capped at the melting temperature ``T_m = 0˚C``
(`melting_temperature`). When capped, the surplus of the balance melts ice at its top.
`OneBandLongwave` and `JeevanjeeRadiation` use the downward longwave of the same time step, `AnalyticBandLongwave`
and `ClearSkyEcCKDRadiation` (which compute the upward emission first) that of the previous time step, and
`ClearSkyEcCKDRadiation` also the previous step's downward shortwave. The stored downward longwave is not part of a
restart file, so the first step after a restart (as the first step of a run) uses zero there; energy is conserved
either way, but a restart with these two schemes is not bitwise identical to a continued run. Without an albedo
scheme the ice albedo ``α_i`` stays zero.

Each ocean grid point is split into the ice-covered part ``ℵ`` and the open water ``1 - ℵ``, which emits at the
sea surface temperature. The ocean then emits at the flux-equivalent temperature

```math
T^* = \left( (1 - ℵ) SST^4 + ℵ T_i^4 \right)^{1/4}
```
such that ``ϵσT^{*4}`` is exactly the emission of the open water and the ice together. The turbulent heat and
humidity fluxes over ice remain insulated as before (see below), the albedo blends linearly with the concentration.

### Energy budget

With the radiation split this way, the net surface flux over the ocean part of a grid point equals what the ocean is
due, the open water's budget, the insulated turbulent fluxes and the conductive heat loss ``ℵF_c`` under the ice,
plus the top melt ``ℵ \times \text{surplus}``. The `SlabOcean` subtracts the top melt from its budget
and `ThermodynamicSeaIce` turns it into a loss of ice volume (with ``ρ_i = 917~\text{kg/m}^3`` and
``L_f = 3.34 \times 10^5~\text{J/kg}``, `density` and `latent_heat_fusion`)

```math
\frac{dV}{dt} = -\frac{ℵ \times \text{surplus}}{ρ_i L_f}
```

After every time step `ThermodynamicSeaIce` applies freezing and bottom melt as exact adjustments in `filter!`,
with ``C`` the heat capacity per area of the slab ocean's mixed layer

1. Top melt beyond the ice present (``V < 0``) warms the ocean instead: ``SST \leftarrow SST + ρ_i L_f \max(-V, 0)/C``, ``V \leftarrow \max(V, 0)``.
2. Freezing: the heat the ocean is below freezing becomes ice, ``V \leftarrow V + C \max(T_f - SST, 0)/(ρ_i L_f)``, and ``SST \leftarrow \max(SST, T_f)``.
3. Bottom melt: a basal flux ``γ \max(SST - T_f, 0)`` with ``γ`` = `basal_melt_coefficient` melts at most the ice present,
   ``ΔV = \min(γ \max(SST - T_f, 0) Δt/(ρ_i L_f), V)``, paid for by the ocean, ``SST \leftarrow SST - ρ_i L_f ΔV/C``.

Then ``ℵ = ℵ(V)``. Hence the ocean heat content minus the latent heat stored in sea ice, ``C \cdot SST - ρ_i L_f V``,
changes only by the surface fluxes, in every season. The default ``γ = 6.1~\text{W/m}^2\text{/K}`` matches the
previous melt rate of ``10^{-7}`` per second and Kelvin for thin ice. Over a prescribed (not interactive) ocean
the ice grows with a nominal mixed layer of 50 m (`prescribed_ocean_heat_capacity`) and the sea surface temperature is
not changed.

Air-sea fluxes ``F`` of heat and humidity are reduced following a simple sea ice concentration-dependent
scaling

```math
F \leftarrow \frac{F}{1 + \frac{ℵ}{c_s}} \\
```
With ``c_s = 0.01`` meaning that at a sea ice concentration of ``ℵ = 1`` only 1% of the
flux goes through the ice. The albedo ``a`` is chosen linearly with respect to sea ice concentration

```math
a = a_o + ℵ(a_i - a_o)
```
with ``a_o = 0.06`` the albedo of ocean, and ``a_i = 0.6`` the albedo of ice.

Setting the sea ice concentration directly, e.g. with `set!(simulation, sea_ice_concentration=...)`, is overwritten
by the concentration of the ice volume after the next time step: set the `sea_ice_volume` instead
(``V = ℵ h(ℵ)``). Initial concentrations are converted to volumes in `initialize!`.

### Filtering

SpeedyWeather calls `filter!(vars, model.sea_ice, model)`
(and similarly for `model.ocean`, `model.land`) once every time step, after all prognostic
variables have been updated, so that components can apply such corrections that don't fit
naturally into a tendency, here freezing and bottom melt. (See also [Ocean models](@ref) and
[Land surface model](@ref) for the respective `filter!` methods.)

### Usage

To be created like

```@example sea_ice
sea_ice = ThermodynamicSeaIce(spectral_grid)
```

and most often together with an albedo that scales linearly with sea ice concentration

```@example sea_ice
albedo = OceanLandAlbedo(ocean=OceanSeaIceAlbedo(spectral_grid), land=AlbedoClimatology(spectral_grid))
```

Using `ocean=GlobalConstantAlbedo(0.06)` instead would disable the effect of sea ice on
albedo (basically an ocean-coloured sea ice), then passed on to the model constructor

```@example sea_ice
model = PrimitiveWetModel(spectral_grid; sea_ice, albedo)
nothing # hide
```

Note that the insulating factor (``r`` above) of sea ice on air-sea fluxes is controlled
in `SlabOcean` which takes as keyword argument `sea_ice_insulation = (x) -> x` a function
of one argument (the sea ice concentration) yielding the reduction of air-sea fluxes,
see [Slab ocean](@ref). So you likely want to map 0 to 0 with that function to have
no insulation (i.e. no impact) in the case of no ice.

## Output

And output is added like

```@example sea_ice
add!(model, OutputVariables.SeaIceConcentration())
```

With `ThermodynamicSeaIce` also its sea ice volume per unit ocean area (the prognostic variable) and the
ice surface temperature

```@example sea_ice
add!(model, OutputVariables.SeaIceVolume(), OutputVariables.SeaIceSurfaceTemperature())
```

A restart file written before the sea ice volume became the prognostic variable has a concentration but no
volume. Set the volume from the concentration after such a restart, e.g. with
`V = SpeedyWeather.sea_ice_volume.(ℵ, Ref(model.sea_ice))`.

or as part of `OutputVariables.Ocean()`

```@example sea_ice
add!(model, OutputVariables.Ocean())
```

For more information see [Output variables](@ref).