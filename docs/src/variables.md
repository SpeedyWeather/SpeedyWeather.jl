# Variables

At the top of SpeedyWeather's type tree sits the `Simulation`, containing
`Variables` and model (e.g. `BarotropicModel`), which in itself contains
model components with their own fields and so on, see [Models](@ref).

The variables are split into

- groups: prognostic, grid, tendencies, dynamics, parameterizations, particles and scratch
- namespace, these are optional sub-categories, e.g. `prognostic.land` and `prognostic.ocean`

Such that their paths are fully defined by

```julia
simulation.variables.group.namespace.name  # if they are in a namespace
simulaiton.variables.group.name            # otherwise
```

Where to find the variables can quickly get complicated at that degre of nesting.
The following is to give users a better overview of how simulation, variables and model
are structured within SpeedyWeather.

## Simulation

Let's start at the top. When creating a `Simulation`, its fields are

```@example variables
using SpeedyWeather
spectral_grid = SpectralGrid(nlayers = 1)
model = BarotropicModel(spectral_grid)
simulation = initialize!(model)
```

the `variables` contain all arrays but also the clock or other scalars
that are supposed to be _variable_, so changing while a simulation is running.
Note that the `variables` depend on what variables are requested by the model and its
components, see [Variables](@ref) and [Variable system](@ref) for more information.

In contrast, we largely think of `model` as being constant after initialization.
This is not completely true, as `model` does contain mutuable structs
but mostly for output and feedback. Values that would influence the
variables are considered read-only after initialization but there is no
hard restriction on this, e.g. you can use [Intrusive callbacks](@ref) to change
the model during integration.

## The `Variables` struct

All simulation variables (prognostic and diagnostic) are stored in `simulation.variables`,
which is a `Variables` struct with 7 (hardcoded) groups `prognostic`, `grid`, `tendencies`,
`dynamics`, `parameterizations`, `particles` and `scratch`.
The variables are model-specific, each model only allocates the variables it needs.
The prognostic variables in `variables.prognostic` are generally in spectral coefficients,
`variables.grid` hold gridded variables, `variables.tendencies` the tendencies,
`variables.dynamics` work arrays that are computed by the dynamical core. `variables.parameterizations`
are those required by the parameterizations and `variables.particles` by the particle advection.
`variables.scratch` are scratch arrays: These can be used in any computation by should be considered
in an undefined state, so write to it before you read from it. Any other component can leave this
in any state. But you can use them to avoid allocations and hold intermediate results.

A full overview of all variables can be easily printed with:

```@example variables
simulation.variables
```

These are the default variables of the `BarotropicModel`, for the `ShallowWaterModel` we have

```@example variables
spectral_grid = SpectralGrid(nlayers = 1)
model = ShallowWaterModel(spectral_grid)
simulation = initialize!(model)
simulation.variables
```

The `PrimitiveDryModel` has the following default variables

```@example variables
spectral_grid = SpectralGrid(nlayers = 8)
model = PrimitiveDryModel(spectral_grid)
simulation = initialize!(model)
simulation.variables
```

And the most complex model, the `PrimitiveWetModel` allocates

```@example variables
spectral_grid = SpectralGrid(nlayers = 8)
model = PrimitiveWetModel(spectral_grid)
simulation = initialize!(model)
simulation.variables
```

## Setting variables

The prognostic variables can be mutated (e.g. to set new initial conditions) with the [`SpeedyWeather.set!`](@ref) function.
Other variables can be set too but they might be overwritten such that your changes may have a different
effect than you expect. You can specify `group` (default `=:prognostic`) and `namespace` (default `=nothing`)
in `set!` to set variables, e.g.

```@example variables
set!(simulation, sea_surface_temperature=300, namespace=:ocean)
```

For another example, see [Set tracers](@ref).

## Step dimension

As you'll notice when inspecting `simulation.variables` is that many prognostic and tendency variables
will have one more dimension than you may think they should have. This is the step dimension as
many time stepping schemes require either several steps in the prognostic variables (think Leapfrog)
or several tendency steps (think multi-step methods). In many cases this additional dimension
will just be a trailing singleton dimension and you can drop it by selecting it, e.g. `[:, 1]`.
As the step dimension essentially contains several "versions" of the same variable the error you
will make by selecting the wrong one is generally small. In Leapfrog, you would simply select
the previous time step for example. So prognostic and tendency variables have the following dimensions
in that order:

- horizontal dimension (unravalled longitude and latitude, or degree/order of spherical harmonics)
- vertical (optional, many 2D variables will drop this)
- step (all [Time stepped variables](@ref) and tendencies have this)

So a conceptually 2D horizontal-only variable may use an array of

- one dimesion (unravelled horizontal dimension),
- two dimensions (trailing singleton dimension for the vertical),
- two dimensions (trailing dimension for the step),
- or three dimesions (singleton dimension for the vertical, and dimension for step)

and similar with 3D variables (but the vertical dimension isn't singleton then).
Whether the step dimension is singleton or not depends on the time stepping scheme in use.

To ease the selection of the step dimension you can use the `get_step` function without any argument
which will automatically create a view onto the array selecting the last step as this
in many cases represents the "current" step and not any previous ones.
But this depends on your time stepping scheme (which the variable itself does not know about).
For example

```julia
size(simulation.variables.grid.u)
```

is a 3D array representing a 3D variable. The array has one more dimension (the last) for the step
but one less because horizontal dimensions are unravelled into one. The step dimension
here has two steps which represent with Leapfrog time stepping the previous and current step,
so to not select the previous time step, you should index `[ij, k, 2]` here.
This is what the `get_step` function will do for you

```julia
get_step(simulation.variables.grid.u)       # selects the last step index automatically
get_step(simulation.variables.grid.u, 1)    # select first step
get_steP(simulation.variables.grid.u, 2)    # select 2nd step
```

Then you can use `[ij, k]` indexing afterwards again. 