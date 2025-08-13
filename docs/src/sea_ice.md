# Sea ice models

The following describes the currently implemented sea ice models, which are

```@example sea_ice
using InteractiveUtils # hide
using SpeedyWeather
subtypes(SpeedyWeather.AbstractSeaIce)
```

## No sea ice

`NoSeaIce` as the name says, initializes `prognostic_variables.sea_ice_concentration` 
to zero. But you may use `set!(simulation, sea_ice_concentration=...)` to set the
sea ice concentration manually. `NoSeaIce` does not do anything on every time step
so your manual modifications will prevail after `initialize!`.

To be used like

```@example sea_ice
spectral_grid = SpectralGrid(trunc=31)
sea_ice = NoSeaIce(spectral_grid)
```

which is then passed on to the model constructor like

```@example sea_ice
model = PrimitiveWetModel(spectral_grid; sea_ice)
model.sea_ice
```

Note that in SpeedyWeather the sea ice model and its albedo are defined
independently, means you can have a sea ice model without affecting the albedo
and `NoSeaIce` but set the sea ice concentration manually and use its albedo
effect, this will be discussed in `ThermodynamicSeaIce` below.

## Thermodynamic sea ice model

Thermodynamic sea ice models do not account for advection through wind and currents or rheological processes
but are only concerned with local freezing, ice growth and melt and the thermodynamic processes involved.
The `ThermodynamicSeaIce` model defined here uses only sea ice concentration ``c`` in units of
``[\text{m}^2/\text{m}^2]`` as a prognostic variable.

Melting and freezing is dependend on the sea surface temperature ``SST`` at the previous time step ``i-1``
and the current ``i`` and air-sea fluxes applied to get from one time step to the other.
Sea ice may modify the sea surface temperatures SST (particularly in the case of freezing) and so we denote the
uncorrected SST with ``SST^*`` meaning that the sea ice time step has not been applied yet.
The time step size is ``\Delta t``.

First determine an insulation factor ``r`` (``r=0`` equals full insulation, ``r=1`` no insulation) linearly
from sea ice concentration which reduces the air-sea flux ``F`` used as the tendency for the slab ocean SST
```math
\begin{aligned}
r &= 1 - c_{i-1} \\
SST_i^* &= SST_{i-1} + \Delta t~r~F \\
\end{aligned}
```
This happens inside the [Slab ocean](@ref) model.

Now determine a tendency in sea ice concentration ``\Delta c`` from the melting and freezing, both proportional to the
difference of the sea surface temperature to freezing temperature with freeze rate ``f`` in ``[\text{m}^2/\text{m}^2/K]``
and melt rate ``m`` in ``[\text{m}^2/\text{m}^2/s/K]``. The freeze rate misses a ``1/s`` in the units as it is applied with
``1/\Delta t`` below because ``\min(SST_i^* - T_f, 0)`` estimates ``\Delta t rF`` from above not purely the
(insulated) flux ``rF``

```math
\begin{aligned}
\Delta c &= -m \max(SST_i^* - T_f, 0) - \frac{f}{\Delta t} \min(SST_i^* - T_f, 0) \\
c_i^* &= c_{i-1} + \Delta t \Delta c
\end{aligned}
```

And we bound the uncorrected sea surface temperature ``SST_i^*`` by freezing at ``T_f = -1.8˚C``
and the uncorrected sea ice concentration ``c_i^*`` in ``[0, 1]`` by

```math
\begin{aligned}
SST_i &= max(SST_i^*, T_f) \\
c_i &= \max( \min( c_i^*, 1) , 0)
\end{aligned}
```

### Usage

To be created like

```@example sea_ice
sea_ice = ThermodynamicSeaIce(spectral_grid)
```

and most often together with an albedo that scales linearly with sea ice concentration

```@example sea_ice
albedo = Albedo(ocean=OceanSeaIceAlbedo(spectral_grid), land=AlbedoClimatology(spectral_grid))
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
add!(model, SpeedyWeather.SeaIceConcentrationOutput())
```

or as part of `SpeedyWeather.OceanOutput()` which however needs splatting `...`
to unpack the tuple

```@example sea_ice
add!(model, SpeedyWeather.OceanOutput()...)
```

