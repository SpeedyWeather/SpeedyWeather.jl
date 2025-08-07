# Ocean and sea ice models

The following describes the currently implement ocean and sea ice models.


## Slab ocean model

### Usage

To be used like

```julia
ocean = SlabOcean(spectral_grid, mixed_layer_depth=50)
model = PrimitiveWetModel(spectral_grid; ocean)
```

And output is added like

```julia
add!(model, SpeedyWeather.SeaSurfaceTemperatureOutput())
```


## Thermodynamic sea ice model

Thermodynamic sea ice models do not account for advection through wind and currents or rheological processes
but are only concern with local freezing, ice growth and melt and the thermodynamic processes involved.
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
This happens inside the [Slab ocean model](@ref).

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


To be used like

```julia
sea_ice = ThermodynamicSeaIce(spectral_grid)
albedo = Albedo(ocean=OceanSeaIceAlbedo(spectral_grid), land=AlbedoClimatology(spectral_grid))
model = PrimitiveWetModel(spectral_grid; sea_ice, albedo)
```

And output is added like

```julia
add!(model, SpeedyWeather.SeaIceConcentrationOutput())
```