# Models

SpeedyWeather implements several models that define which equation are being solved

- `BarotropicModel` solves the 2D barotropic vorticity equations
- `ShallowWaterModel` solves the 2D shallow water equations
- `PrimitiveDryModel` solves the 3D primitive equations without humidity
- `PrimitiveWetModel` solves the 3D primitive equations with humidity

but each of these models has further modularity such that the equations being
solved also depend on the model components being defined. For example,
each of these models can have custom forcing and drag terms, or the
primitive equation models can have endless combinations of parameterizations
altering all sorts of components of the planetary climate system that is
simulated.

We briefly provide an overview of the components of each model. Constructing
a `model` follows the general syntax that is
`SomeModel(spectral_grid, component=my_component, ...)` which we call
the _model constructor_. The first argument is the
`spectral_grid` defining the resolution, then additional keyword arguments
can change model components for custom or non-default ones. The order
does not matter and you can provide as many as you like just the
`model.component` names have to fit the keywords of the arguments.

The model construction essentially gathers all default or non-default
model components (it will construct the default components) but doesn't to much else.
All the initialization of the model components and the allocation
of [Variables](@ref) happens in `initialize!(model)`.

## Model components

Many types and components in SpeedyWeather have
extended Julia's `show` function with some pretty printing to give you
a better overview of its contents, e.g. a `clock::Clock` is printed as

```@example models
using SpeedyWeather
clock = Clock()
```

illustrating the fields within a clock (one per row), their types (indicated by `::`)
and (unless they are an array) also their values. For structs (or NamedTuples) within
structs however, this information, is not be printed by default.
You can use Julia's autocomplete like `clock.<tab>` by hitting tab after the `.` to inspect
the fields of an instance and that way go down every branch of the simulation "tree".

## BarotropicModel

The `BarotropicModel` is the simplest model we have, which will not have many of
the model components that are needed to define the primitive equations for example.
If you create a model with non-default conponents they will show up here,
so your model may look different depending on what you have constructed!

```@example models
using SpeedyWeather
spectral_grid = SpectralGrid(nlayers=1)     # 2D models require nlayers=1
model = BarotropicModel(spectral_grid)
```

## ShallowWaterModel

The `ShallowWaterModel` is similar to the `BarotropicModel`, but it contains for example
orography, that the `BarotropicModel` doesn't have.

```@example models
spectral_grid = SpectralGrid(nlayers=1)     # 2D models require nlayers=1
model = ShallowWaterModel(spectral_grid)
```

## PrimitiveDryModel

The `PrimitiveDryModel` is a big jump in complexity compared to the 2D models, but
because it doesn't contain humidity, several model components like `surface_humidity_flux`
aren't needed.

```@example models
spectral_grid = SpectralGrid()
model = PrimitiveDryModel(spectral_grid)
```

## PrimitiveWetModel

The `PrimitiveWetModel` is the most complex model we currently have, hence its
field tree is the longest, defining many components for the physics parameterizations.

```@example models
model = PrimitiveWetModel(spectral_grid)
```