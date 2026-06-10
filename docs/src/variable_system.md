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

## Prognostic variables

A variable is _prognostic_ if its contains information that carries on into the next time step.
_Diagnostic_ variables, in contrast, depend only on the prognostic variables and particularly
not on themselves.
Prognostic variables the solution to either a (continuous) differential equation or a discrete map
(in the sense of an evolution function).
But we also consider some time-varying boundary conditions like greenhouse gas concentrations
to be prognostic although their evolution might be prescribed rather than dynamically evolving.
While they could also be implemented as a time-dependent forcing, making them a prognostic
variable means there is (1) a user interface to change the state of this variable (think
climate scenarios like an abrupt 4xCO2 increase) and (2) this variable can be used more widely
across different components. Some examples for prognostic variables

- temperature as the solution to the partial differential equation for temperature
- the clock which "solves" the equation ``d/dt (clock) = 1``
- a random pattern following an autoregressive process
- CO2 concentration following a prescribed evolution

Note that there is an exception to our definition for the tendencies in multi-step methods.
These retain information from previous tendencies to form an average with the current tendency
for more accuracy. As such these tendencies qualify as prognostic variable as defined above,
which we, however, ignore in favour or grouping them into `variables.tendencies`.
It is arguably less confusing to have some memory in the tendencies than to move a
tendency to the prognostic variables for some time steppers but not for others.

## Time stepped variables

The prognostic variables _can but do not have to_ be subject to time stepping as defined in `model.time_stepping`.
In the examples above, one likely want temperature to be time stepped but the random pattern
should not in the same way as the time stepping is different (discrete vs continuous).
So a random pattern is prognostic but the random processes should determine its temporal evolution,
and not the general time stepper in `model.time_stepping` that would also be used for temperature.
Similarly, one may want to choose whether sea surface temperature is time stepped with the
general time stepper or whether the `model.ocean` component itself should be responsible for
evolving this prognostic variable in time.

To satisfy this flexibility a prognostic variable that is subject to time stepping is defined by
the existence of a tendency field (grid or spectral) in

- `simulation.tendencies`
- `simulation.tendencies.tracers`
- `simulation.tendencies.ocean`
- `simulation.tendencies.land`

Other namespace (e.g. `grid`) are ignored. This means that if you define a new prognostic variable `var` and you want
it to be time stepped then also define `simulation.tendencies.var`. Or choose namespaces `ocean` or `land` for both,
other namespaces are ignored.
A prognostic variable ignored by the general time stepper should not have a tendency defined.
Use a for example in-place updates instead or another namespace that is ignored. For example,

```julia
function SpeedyWeather.variables(::MyComponent)
    return (
        PrognosticVariable(:var1, SpeedyWeather.Grid2D()),                      # not time stepped
        PrognosticVariable(:var2, SpeedyWeather.Grid2D(), namespace=:ocean),    # time stepped
        PrognosticVariable(:var3, SpeedyWeather.Grid2D(), namespace=:other),    # not time stepped
        
        TendencyVariable(:var2, SpeedyWeather.Grid2D(), namespace=:ocean),
        TendencyVariable(:var3, SpeedyWeather.Grid2D(), namespace=:other),
    )
end
```

`var1` is _not_ time stepped as no tendency `var1` exists.
`var2` is time stepped as a tendency exists in the same name space `:ocean`.
`var3` is _not_ time stepped even though the tendency exists as the namespace `:other` is not `:ocean` or `:land`.

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
SpeedyWeather.allocate(::SpeedyWeather.AbstractVariable{MyDictDim}, model::AbstractModel) = Dict()
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

## Variable fusion 

Variables can be "fused" together to be allocated as one block of contiguous memory.
This is done primarly with GPU-optimiziation in mind, so that e.g. `transform!` calls can be batched together not just across levels,
but also across different variables. All variables that are defined with the same `fuse = :fuse_name` keyword argument are allocated together.
At the same time a `view` is defined that allows for the variable to be used as it would without the variable fusion.
As this is primarily a performance optimization, it is not required to run basic custom parameterizations or model components. 