# Tree structure

At the top of SpeedyWeather's type tree sits the `Simulation`, containing
variables and model, which in itself contains model components with their
own fields and so on.
(Note that we are talking about the structure of structs within structs
not the type hierarchy as defined by subtyping abstract types.)
This can quickly get complicated with a lot of nested structs. The following
is to give users a better overview of how simulation, variables and model
are structured within SpeedyWeather. Many types in SpeedyWeather have
extended Julia's `show` function to give you an overview of its contents,
e.g. a `clock::Clock` is printed as

```@example structure
using SpeedyWeather
clock = Clock()
```

illustrating the fields within a clock, their types and (unless they are an array)
also their values. For structs within structs however, this information,
would not be printed by default. You could use Julia's autocomplete like
`clock.<tab>` by hitting tab after the `.` to inspect the fields of an instance
but that would require you to manually go down every branch of that tree.
To better visualise this, we have defined a `tree(S)` function for any instance
`S` defined in SpeedyWeather which will print every field, and also
its containing fields if they are also defined within SpeedyWeather.
The "if defined in SpeedyWeather" is important because otherwise
the tree would also show you the contents of a complex number or other types
defined in Julia Base itself that we aren't interested in here. 
But let's start at the top.

## Simulation

When creating a `Simulation`, its fields are
```@example structure
spectral_grid = SpectralGrid(nlayers = 1)
model = BarotropicModel(; spectral_grid)
simulation = initialize!(model)
```
the `prognostic_variables`, the `diagnostic_variables` and the `model` (that we just
initialized). We could now do `tree(simulation)` but that gets very lengthy and
so will split things into `tree(simulation.prognostic_variables)`,
`tree(simulation.diagnostic_variables)` and `tree(simulation.model)` for more
digestible chunks. You can also provide the `with_types=true` keyword to get
also the types of these fields printed, but we'll skip that here.

## Prognostic variables

The prognostic variables struct is parametric on the model type, `model_type(model)`
(which strips away its parameters), but this is only to dispatch over it.
The fields are for all models the same, just the barotropic model would not
use temperature for example (but you could use nevertheless). 

```@example structure
tree(simulation.prognostic_variables)
```

The prognostic variable struct can be mutated (e.g. to set new initial conditions) with the [`SpeedyWeather.set!`](@ref) function. 

## Diagnostic variables

Similar for the diagnostic variables, regardless the model type, they contain the
same fields but for the 2D models many will not be used for example.

```@example structure
tree(simulation.diagnostic_variables)
```

## BarotropicModel

The `BarotropicModel` is the simplest model we have, which will not have many of
the model components that are needed to define the primitive equations for example.
Note that forcing or drag aren't further branched which is because the default
`BarotropicModel` has `NoForcing` and `NoDrag` which don't have any fields. 
If you create a model with non-default conponents they will show up here. 
`tree` dynamicallt inspects the current contents of a (mutable) struct and
that tree may look different depending on what model you have constructed!

```@example structure
model = BarotropicModel(; spectral_grid)
tree(model)
```

## ShallowWaterModel

The `ShallowWaterModel` is similar to the `BarotropicModel`, but it contains for example
orography, that the `BarotropicModel` doesn't have.

```@example structure
model = ShallowWaterModel(; spectral_grid)
tree(model)
```

## PrimitiveDryModel

The `PrimitiveDryModel` is a big jump in complexity compared to the 2D models, but
because it doesn't contain humidity, several model components like evaporation
aren't needed.

```@example structure
spectral_grid = SpectralGrid()
model = PrimitiveDryModel(; spectral_grid)
tree(model)
```

## PrimitiveWetModel

The `PrimitiveWetModel` is the most complex model we currently have, hence its
field tree is the longest, defining many components for the physics parameterizations.

```@example structure
model = PrimitiveWetModel(; spectral_grid)
tree(model)
```

## Size of a Simulation

The `tree` function also allows for the `with_size::Bool` keyword (default `false`),
which will also print the size of the respective branches to give you an idea of
how much memory a SpeedyWeather simulation uses.

```@example structure
tree(simulation, max_level=1, with_size=true)
```

And with `max_level` you can truncate the tree to go down at most that many levels.
1MB is a typical size for a one-level T31 resolution simulation. In comparison,
a higher resolution `PrimitiveWetModel` would use

```@example structure
spectral_grid = SpectralGrid(trunc=127, nlayers=8)
model = PrimitiveWetModel(spectral_grid)
simulation = initialize!(model)
tree(simulation, max_level=1, with_size=true)
```