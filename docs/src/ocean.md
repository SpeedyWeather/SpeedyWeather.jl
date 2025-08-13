# Ocean models

The following describes the currently implemented ocean models,
some prescribed sea surface temperature (not dependent on the state
of other variables) others are active (dependent on the atmospheric
state). All but SlabOcean force the atmosphere, as the sea surface temperatures are
used to calculate surface heat fluxes. SlabOcean interacts with both Atmosphere and SeaIce. All models can be used
with PrimitiveDry and PrimitiveWet models; for the former, only
surface heat fluxes are applied, albeit not the humidity fluxes.

```@example ocean
using InteractiveUtils # hide
using SpeedyWeather
subtypes(SpeedyWeather.AbstractOcean)
```

## Aqua planet

The `AquaPlanet` in SpeedyWeather is a prescribed sea surface temperature
that only depends on latitude, applying a cosine squared between the
Equator and the poles

```@example ocean
ocean = AquaPlanet(spectral_grid)
```

pole and equator temperatures can be modified and `mask` can mask the
sea surface temperatures according to `model.land_sea_mask`. Otherwise
sea surface temperatures are defined everywhere but the land-sea mask
will determine their proportional contribution to surface fluxes.

## Constant ocean climatology

`ConstantOceanClimatology` is like `SeasonalOceanClimatology` but constant in time.
At `initalize!(model)` the seasonal ocean climatology is read from file and interpolated
to the current time, but the sea surface temperature field is not further updated
thereafter. To be used like

```@example ocean
ocean = ConstantOceanClimatology(spectral_grid)
```

```@example ocean
model = PrimitiveWetModel(spectral_grid, ocean=ocean)
simulation = initialize!(model, time=DateTime(2000, 6, 1))
nothing # hide
```

which will use the sea surface temperature climatology from 1 June but
not change it thereafter. Note that because nothing happens in the ocean time step
you can use `set!(simulation, sea_surface_temperature=...)` to modify the 
sea surface temperatures further at any point.

## Seasonal ocean climatology

Created like

```@example ocean
spectral_grid = SpectralGrid(trunc=31)
ocean = SeasonalOceanClimatology(spectral_grid)
```

The `SeasonalOceanClimatology` reads monthly sea surface temperature
fields from file, and interpolates them in time on every time step
and writes them to the prognostic variables. Several options
exist to load another `file` from `path` etc. For a full list of
options type `?SeasonalOceanClimatology`. 
To be passed on to the model constructor like

```@example ocean
model = PrimitiveWetModel(spectral_grid, ocean=ocean)
nothing # hide
```

Or, equivalently, `; ocean` as Julia can match keyword arguments by name,
where `;` just means that every argument that follows is a keyword argument.

The time of the year is determined by the clock in `prognostic_variables.clock`
such that `initialize!(model, time=DateTime(2000, 1, 1))` would interpolate
the seasonal climatology onto the first of January.

## Slab ocean

The most complex ocean model implemented is a slab ocean model:
it can heat up and cool down given atmospheric fluxes and in turn
warm up or cool down the atmosphere from below and provide
sources or sink of humidity. The slab ocean model has
only sea surface temperature as prognostic variable so
the humidity flux from precipitation does not decrease
salinity as there is none. You can think of the `SlabOcean`
as a freshwater ocean. To be used like

```@example ocean
ocean = SlabOcean(spectral_grid)
```

`specific_heat_capacity`, `mixed_layer_depth` and `density` determine
multiplicatively the effective heat capacity of the mixed layer
``C_0``. Then the temporal evolution of sea surface temperature ``SST``
is given by


```math
C_0 \frac{d SST}{dt} = (1-r) \left( R_{sd} - R_{su} - R_{lu} + R_{ld} - L_v*E_v - S \right)
```

with ``r \in [0, 1]`` an insulating factor due to sea ice, see [Thermodynamic sea ice model](@ref),
by default chosen as ``r = c``, i.e proportional to sea ice concentration ``c``, insulating
the ocean from surface fluxes at full sea ice coverage (``c=1``) and no insulation for
ice-free ocean (``c=0``). The surface fluxes are ``R_{sd}`` for shortwave (``s``) radiative
flux downward (positive sign, heating the ocean); ``R_{su}`` reflected upward shortwave flux
given albedo (negative sign, cooling the ocean); ``R_{lu}`` a longwave upward flux through
which the ocean loses heat dependent on its sea surface temperature (negative sign);
``R_{ld}`` downward longwave flux where the atmosphere warms the ocean due to its own
tempreature (positive sign); ``E_v`` is the upward evaporative (or surface condensation) flux,
converted to ``W/m^2`` through multiplication with the latent heat of condensation
``L_v`` (negative sign as evaporation cools the ocean); and ``S`` the sensible heat flux
through turbulent transport of heat in the atmospheric boundary layer. All surface flux
terms have units of ``W/m^2``.

The slab ocean is then passed on to the model constructor as with all other
ocean models

```@example ocean
model = PrimitiveWetModel(spectral_grid; ocean)
nothing # hide
```

## Output

Sea surface temperature output is added like

```@example ocean
add!(model, SpeedyWeather.SeaSurfaceTemperatureOutput())
```

or collectively with sea ice concentration

```@example ocean
add!(model, SpeedyWeather.OceanOutput()...)
```

where the splatting operator `...` has to be applied to unpack all output variables in
the tuple, all output variable groups are defined as tuples.
