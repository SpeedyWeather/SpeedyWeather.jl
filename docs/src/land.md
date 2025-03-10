# Land surface model

The land surface in SpeedyWeather is represented through several
components each of which can be changed in a modular way.

- The actual [LandModel](@ref) or [DryLandModel](@ref) describing the equations for soil temperature and soil moisture, and vegetation or rivers.
- [The land-sea mask](@ref)
- The surface [Albedo](@ref)
- The [Orography](@ref)

As these components are largely independent of another,
it is possible to have a white (albedo high) ocean at the
top of Mt Everest, or to paint the Sahara black and moving
it below sea-level. The first one of these is discussed below
for the others see the respective sections.

# Dry vs wet land

The type hierarchy of theactual land surface model is defined as

```@example land
using SpeedyWeather
using InteractiveUtils # hide
subtypes(SpeedyWeather.AbstractLand)
```
So that a dry land does not have moisture whether the atmosphere
is a `PrimitiveWetModel` (with humidity so it can rain) or a
`PrimitiveDryModel` (without humidity). In contrast a wet land
does have some soil moisture (which can be zero) and in combination
with a `PrimitiveWetModel` will increase with precipitation for
example. You can combine dry and wet land and dry and wet atmosphere
freely.

# LandModel

The default `LandModel` in SpeedyWeather contains
(at the moment other than 2 soil layers are not supported or experimental)

```@example land
spectral_grid = SpectralGrid(trunc=31, nlayers=8, nlayers_soil=2)
land = LandModel(spectral_grid)
```

With `land.geometry` currently used to define the layer thickness

```@example land
land.geometry
```

And `land.thermodynamics` used to define thermodynamics-relevant
surface properties as conceptual "constants" used in the
other components of the land surface model `land`.

```@example land
land.thermodynamics
```

To change these you can either mutate the fields

```@example land
land.thermodynamics.heat_conductivity = 0.25
```

or create a new model component `thermodynamics` 
passed on to the land model constructor

```@example land
thermodynamics = LandThermodynamics(spectral_grid, heat_conductivity=0.25)
land = LandModel(spectral_grid; thermodynamics)
land.thermodynamics
```

(note that `; thermodynamics` is a Julia shortcut instead of `thermodynamics=thermodynamics`
we often use matching keyword arguments by their name).
Similarly you can change `land.geometry` which however is not mutable.

Generally, a `LandModel` or `DryLandModel` is similarly constructed to a `PrimitiveWetModel`
or `PrimitiveDryModel`. One starts by defining its non-default components,
always passes on the spectral grid as the first argument, and then constructs
the actual land model by passing on those components.
Finally, we use that land model to construct the atmospheric model

```@example land
model = PrimitiveWetModel(spectral_grid; land)
model.land
```

is now the land defined above used when integrating a SpeedyWeather `model`.

# DryLandModel

Alternatively, one can use the `DryLandModel` to explicitly disable any
functionality around soil moisture. By doing so, soil moisture will
always be zero, or treated like zero skipping unnecessary computations.
It is created in the same way (reusing some non-default `thermodynamics` from above)

```@example land
land = DryLandModel(spectral_grid; thermodynamics)
```

but it does not contain soil moisture, vegetation or rivers in contrast to the
`LandModel`.

# Land soil temperature

Currently implemented soil temperatures are

```@example land
subtypes(SpeedyWeather.AbstractLandTemperature)
```

You can use them by passing them on to a
`LandModel`/`DryLandModel` model constructor

```@example land
temperature = LandBucketTemperature(spectral_grid)
land = LandModel(spectral_grid; temperature)
```

such that

```@example land
land.temperature
```

is the `LandBucketTemperature` just defined. Similarly

```@example land
land = DryLandModel(spectral_grid; temperature)
```

if you do not want the land to hold any moisture, vegetation or rivers.

# LandBucketTemperature

`LandBucketTemperature` is a prognostic model of the soil temperature
in the land surface model, interacting two-way with the surface air
temperature. It can warm up through radiation and other surface heat fluxes,
retain thermal energy and release this back to the atmosphere either
in the form of longwave radiative fluxes or sensible heat fluxes
(latent heat fluxes depend on soil moisture, see [Surface fluxes](@ref)).
It is a bucket model such that interaction between neighbouring
grid cells ("buckets") of the land surface only interact through the 
atmosphere with another, there are no direct horizontal fluxes
between cells. In the sense of soil moisture, you can fill a bucket
from above with rainfall, it may leak/drain at the bottom but buckets
are laterally isolated from another. A similar concept applies to the
heat fluxes of the `LandBucketTemperature`.

The `LandBucketTemperature` here follows MITgcm's 2-layer model, as defined
[here](https://mitgcm.readthedocs.io/en/latest/phys_pkgs/land.html).
As this is a 2-layer model, `SpectralGrid(nlayers_soil=2)` is required.
The equations are

```math
\begin{aligned}
\Delta z_1 C_1 \frac{dT_1}{dt} &= F - \lambda\frac{T_1 - T_2}{(\Delta z_1 + \Delta z_2)/2} \\
\Delta z_2 C_2 \frac{dT_2}{dt} &= \lambda\frac{T_1 - T_2}{(\Delta z_1 + \Delta z_2)/2}
\end{aligned}
```

for two layers of thicknesses ``\Delta z_1 = 0.1~m`` (top) and ``\Delta z_2 = 4.0~m`` (layer below)
and respective temperatures ``T_1, T_2``. The total surface downward heat flux
`F` (in ``W/m^2``) forces the surface layer 1, and diffusion scales with the heat conductivity
``\lambda = 0.42 W/m/K`` in the opposite direction of the heat gradient between the layers.
For the parameter choices here it is typical that the surface layer is dominated by the daily cycle,
but the layer below by the seasonal cycle.

The heat capacities ``C_1, C_2`` are diagnosed from the heat capacity of water
``C_w = 4.2 \times 10^6 J/m^3/K`` and dry soil ``C_s = 1.13 \times 10^6 J/m^3/K``
given the soil moistures ``W_1, W_2`` (ratio of available water to field capacity)
of the respective layers.

```math
\begin{aligned}
C_1 &= C_w W_1 \gamma + C_s \\
C_2 &= C_w W_2 \gamma + C_s
\end{aligned}
```

with ``\gamma`` being the field capacity per meter soil.

# Land soil moisture

Currently implemented soil moistures are

```@example land
subtypes(SpeedyWeather.AbstractSoilMoisture)
```




# Albedo