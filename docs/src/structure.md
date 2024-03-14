# Tree structure

At the top of SpeedyWeather's type tree sits the `Simulation`, containing
variables and model, which in itself contains model components and so on.
(Note that we are talking about the structure of structs within structs
not the type hierarchy as defined by subtyping abstract types.)
This can quickly get complicated with a lot of nested structs. The following
is to give users a better overview of how simulation, variables and model
are structured within SpeedyWeather. Many types in SpeedyWeather have
extended Julia's `show` function to give you an overview of its contents,
e.g. a `clock::Clock` is printed as

```@example structure
clock = Clock()
```

illustrating the fields within a clock, their types and (unless they are an array)
also their values. For structs within structs however, this information, however,
would not be printed by default. You could use Julia's autocomplete like
`clock.<tab>` by hitting tab after the `.` to inspect the fields of an instance
but that would require you to manually go down every branch of that tree.
To better visualise this, we have defined a `tree(S)` function for any instance
`S` that's defined in SpeedyWeather which will print every field, and also
its containing fields if they are also defined within SpeedyWeather.
But let's start at the top.

## Simulation

When creating a Simulation, it's fields are
```@example structure
model = BarotropicModel()
simulation = initialize!(model)
```
the `prognostic_variables`, the `diagnostic_variables` and the `model` (that we just
initialized.) We could now do `tree(simulation)` but that gets very lengthy and
so will will split things into `tree(simulation.prognostic_variables)`,
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

## Diagnostic variables

Similar for the diagnostic variables, regardless the model type, they contain the
same fields but for the 2D models many will not be used for example.

```@example structure
tree(simulation.diagnostic_variables)
```

## BarotropicModel

The `BarotropicModel` is the simplest model we have, which will not have many of
the model components that are needed to define the primitive equations for example.

```@example structure
model = BarotropicModel()
tree(model)
```

## ShallowWaterModel

The `ShallowWaterModel` is similar to the `BarotropicModel`, but it contains for example
orography, that the `BarotropicModel` doesn't have.

```@example structure
model = ShallowWaterModel()
tree(model)
```

## PrimitiveDryModel

The `PrimitiveDryModel` is a big jump in complexity compared to the 2D models, but
because it doesn't contain humidity, several model components like evaporation
aren't needed.

```@example structure
model = PrimitiveDryModel()
tree(model)
```

## PrimitiveWetModel

The `PrimitiveWetModel` is the most complex model we currently have, hence its
field tree is the longest, defining many components for the physics parameterizations.

```@example structure
model = PrimitiveWetModel()
tree(model)
```