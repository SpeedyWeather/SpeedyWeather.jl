# Land surface model

The land surface in SpeedyWeather is represented through several
components each of which can be changed in a modular way.

- The actual [`LandModel`](@ref) or [`DryLandModel`](@ref) describing the equations for soil temperature and soil moisture, and vegetation or rivers.
- [The land-sea mask](@ref)
- The surface [Albedo](@ref)
- The [Orography](@ref)

As these components are largely independent of another,
it is possible to have a white (albedo high) ocean at the
top of Mt Everest, or to paint the Sahara black and moving
it below sea-level.

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

# Land soil temperature

Currently implemented soil temperatures are

```@example land
subtypes(SpeedyWeather.AbstractLandTemperature)
```

You can use them by passing them on to a
`LandModel`/`LandDryModel` model constructor

```@example
temperature = LandBucketTemperature(spectral_grid)
land = LandModel(spectral_grid; temperature)
land.temperature
```

and similarly

```@example
land = LandDryModel(spectral_grid; temperature)
```

if you do not want the land to hold any moisture.

# Land soil moisture

Currently implemented soil moistures are

```@example land
subtypes(SpeedyWeather.AbstractSoilMoisture)
```



# Albedo