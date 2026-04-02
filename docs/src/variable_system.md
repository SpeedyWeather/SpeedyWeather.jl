# Variable system

SpeedyWeather implements a dynamical variables system. This means there is no central list
of all variables being allocated by every model and every model component can declare
a set of variables as being required which will then be allocated when the
model is initialized.

In most cases you are advised to use and reuse as much as possible the variables
required, including scratch arrays you can always write to (but don't read from them
before you've overwritten their undefined state). However, there are also situation
where you may want to declare your own new variables or define new variable
"types" (some, typically mutable, object) which we call dimension. This is explained in
the following

## Declare variables

Say we define a new custom component

```@example variable_system
using SpeedyWeather
struct MyAlbedo <: SpeedyWeather.AbstractAlbedo end
```

and in order to compute it we need to have another 2D field available, called `my_albedo`.
Then we can extend the `variables` function as

```@example variable_system
function SpeedyWeather.variables(::MyAlbedo)
    return (
        ParameterizationVariable(:my_albedo, SpeedyWeather.Grid2D(), namespace=:radiation),
    )
end
```

This will now allocate `simulation.variables.parameterizations.radiation.my_albedo` as a 2D field
if the grid type, array type and number format as defined by `spectral_grid` (i.e. GPU-ready etc.).
The first argument, the `name::Symbol` and the 2nd argument the dimension are required arguments,
`namespace::Symbol`, `units::String`, `desc::String` are optional.

Let's see whether this actually works

```@example variable_system
spectral_grid = SpectralGrid()
my_albedo = MyAlbedo()
model = PrimitiveWetModel(spectral_grid, albedo=my_albedo)
simulation = initialize!(model)
simulation.variables.parameterizations.radiation.my_albedo
```

Here we go, the variable has been allocated at the required path and is initialized with zeros.

## Variable groups

In [Declare variables](@ref) we declared a `ParameterizationVariable` but depending on what
you do you may want to use another variable group. These are the hardcoded groups available

- `prognostic` via `PrognosticVariable`
- `grid` via `GridVariable`
- `tendencies` via `TendencyVariable`
- `dynamics` via `DynamicsVariable`
- `parameterizations` via `ParameterizationVariable`
- `particles` via `ParticleVariable`
- `scratch` via `ScratchVariable`

each of them require `name::Symbol` and `dim::AbstractVariabelDim` with `namespace::Symbol`,
`units::String`, `desc::String` optional. 

## Variable dimensions

Variable dimensions are being used to declare the type
(and explicitly or implicitly also the size and element type) of a variable required.
For example, `Grid2D` is a variable dimension declaring it to be of the grid as in
`spectral_grid` but only 2D with no additional e.g. vertical dimension.

Many variable dimensions have already been defined

```@example variable_system
using InteractiveUtils # hide
subtypes(SpeedyWeather.AbstractVariableDim)
```

So for example, following on from above, we could require a 10x10 matrix as

```@example variable_system
SpeedyWeather.variables(A::MyAlbedo) = (
    ParameterizationVariable(:albedo_matrix, SpeedyWeather.MatrixDim(10, 10), units="1", desc="What is it?", namespace=:radiation),
)
```

which would then be allocated at `simulation.variables.parameterizations.radiation.albedo_matrix` but the
size would be hardcoded to 10x10, you can however adapt this definition to use `A.m` or `A.n` in
case this information is in your albedo `A::MyAlbedo` or you can use any information from `model`
when using it as the 2nd argument

```julia
function SpeedyWeather.variables(A::MyAlbedo, model::AbstractModel) 
    n = model.spectral_grid.nlayers
    return (
        ParameterizationVariable(:albedo_matrix, SpeedyWeather.MatrixDim(n, n), units="1", desc="What is it?", namespace=:radiation),
    )
end
```

in which case the `albedo_matrix` is always allocated of size `nlayers`x`nlayers`
as determined in the spectral grid.

## Define new variable dimensions

When defining new custom model components you may need a new variable dimension
which you can define as follows. Here illustrated by introducing `MyDictDim` --
a dictionary dimension.

```@example variable_system
struct MyDictDim <: SpeedyWeather.AbstractVariableDim end
Base.zero(::SpeedyWeather.AbstractVariable{MyDictDim}, model::AbstractModel) = Dict()
```

and use like

```@example variable_system
SpeedyWeather.variables(::MyAlbedo) = (
    ScratchVariable(:my_dict, MyDictDim()),
)
```

So that when initializing a `model` we have

```@example variable_system
simulation = initialize!(model)
simulation.variables.scratch.my_dict
```

a dictionary as a variable!