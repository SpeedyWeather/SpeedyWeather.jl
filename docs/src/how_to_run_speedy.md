# How to run SpeedyWeather.jl

Creating a SpeedyWeather.jl simulation and running it consists conceptually of 4 steps.
In contrast to many other models, these steps are bottom-up rather then top-down.
There is no monolithic interface to SpeedyWeather.jl, instead all options that a
user may want to adjust are chosen and live in their respective model components.

1. Create a [SpectralGrid](@ref) which defines the grid and spectral resolution.
2. [Create model components](@ref create_model_components) and combine to a [model](@ref create_model).
3. [Initialize the model to create a simulation](@ref initialize).
4. [Run that simulation](@ref run).

In the following we will describe these steps in more detail.
If you want to start with some examples first, see [Examples](@ref Examples)
which has easy and more advanced examples of how to set up
SpeedyWeather.jl to run the simulation you want.

## SpectralGrid

The life of every SpeedyWeather.jl simulation starts with a `SpectralGrid` object.
A `SpectralGrid` defines the physical domain of the simulation and its discretization.
This domain has to be a sphere because of the spherical harmonics, but it can have a different radius.
The discretization is for spectral, grid-point space and the vertical as this determines the size of many
arrays for preallocation, for which als the number format is essential to know.
That's why `SpectralGrid` is the beginning of every SpeedyWeather.jl simulation and that is why
it has to be passed on to (most) model components.

The default `SpectralGrid` is

```@example howto
using SpeedyWeather
spectral_grid = SpectralGrid()
```
You can also get the help prompt by typing `?SpectralGrid`.
Let's explain the details: The spectral resolution is T31, so the largest
wavenumber in spectral space is 31, and all the complex spherical harmonic
coefficients of a given 2D field (see [Spherical Harmonic Transform](@ref))
are stored in a [`LowerTriangularMatrix`](@ref lowertriangularmatrices)
in the number format Float32. The radius of the sphere is
6371km by default, which is the radius of the Earth.

This spectral resolution is combined with an
[octahedral Gaussian grid](@ref OctahedralGaussianGrid) of 3168 grid points.
This grid has 48 latitude rings, 20 longitude points around the poles
and up to 96 longitude points around the Equator. Data on that
grid is also stored in Float32. The resolution is therefore on average about 400km.
In the vertical 8 levels are used, using [Sigma coordinates](@ref).

The resolution of a SpeedyWeather.jl simulation is adjusted using the
`trunc` argument, this defines the spectral resolution and the grid
resolution is automatically adjusted to keep the aliasing between
spectral and grid-point space constant (see [Matching spectral and grid resolution](@ref)).
```@example howto
spectral_grid = SpectralGrid(trunc=85)
```
Typical values are 31, 42, 63, 85, 127, 170, ... although you can technically
use any integer, see [Available horizontal resolutions](@ref) for details.
Now with T85 (which is a common notation for `trunc=85`) the grid
is of higher resolution too. You may play with the `dealiasing` factor,
a larger factor increases the grid resolution that is matched with a given
spectral resolution. You don't choose the resolution of the grid directly,
but using the `Grid` argument you can change its type (see [Grids](@ref))
```@example howto
spectral_grid = SpectralGrid(trunc=85, dealiasing=3, Grid=HEALPixGrid)
```

## Vertical coordinates and resolution

The number of vertical layers or levels (we use both terms often interchangeably)
is determined through the `nlayers` argument. Especially for the
`BarotropicModel` and the `ShallowWaterModel` you want to set this to
```@example howto
spectral_grid = SpectralGrid(nlayers=1)
```
For a single vertical level the type of the vertical coordinates does not matter,
but in general you can change the spacing of the sigma coordinates
which have to be discretized in ``[0, 1]``
```@example howto
vertical_coordinates = SigmaCoordinates(0:0.2:1)
```
These are regularly spaced [Sigma coordinates](@ref), defined through their half levels.
The cell centers or called full levels are marked with an ×.
You have to provide this as an argument to `Geometry`,
i.e. `Geometry(spectral_grid, vertical_coordinates=σ)` and pass this on to the
model constructor if you want to use custom sigma coordinates. At the moment,
other vertical coordinates are not supported.

## [Creating model components](@id create_model_components)

SpeedyWeather.jl deliberately avoids the namelists that are the main
user interface to many old school models. Instead, we employ a modular approach
whereby every non-default model component is created with its respective
parameters to finally build the whole model from these components.

If you know which components you want to adjust you would create them right away,
however, a new user might first want to know which components there are,
so let's flip the logic around and assume you know you want to run a `ShallowWaterModel`.
You can create a default `ShallowWaterModel` like so and inspect its components
```@example howto
model = ShallowWaterModel(spectral_grid)
```

So by default the `ShallowWaterModel` uses a [Leapfrog `time_stepping`](@ref leapfrog),
which you can inspect like so
```@example howto
model.time_stepping
```

Model components often contain parameters from the `SpectralGrid` as they are needed
to determine the size of arrays and other internal reasons. You should, in most cases,
just ignore those. But the `Leapfrog` time stepper comes with `Δt_at_T31` which
is the parameter used to scale the time step automatically. This means at a spectral
resolution of T31 it would use 30min steps, at T63 it would be ~half that, 15min, etc.
Meaning that if you want to have a shorter or longer time step you can create a new
`Leapfrog` time stepper. All time inputs are supposed to be given with the help of 
`Dates` (e.g. `Minute()`, `Hour()`, ...). But remember that (almost) every model component
depends on a `SpectralGrid` as first argument.
```@example howto
spectral_grid = SpectralGrid(trunc=63, nlayers=1)
time_stepping = Leapfrog(spectral_grid, Δt_at_T31=Minute(15))
```
The actual time step at the given resolution (here T63) is then `Δt_sec`, there's
also `Δt` which is a scaled time step used internally, because SpeedyWeather.jl
[scales the equations](@ref scaled_swm) with the radius of the Earth,
but this is largely hidden (except here) from the user. With this new 
`Leapfrog` time stepper constructed we can create a model by passing
on the components (they are keyword arguments so either use `; time_stepping`
for which the naming must match, or `time_stepping=my_time_stepping` with
any name)
```@example howto
model = ShallowWaterModel(spectral_grid; time_stepping)
```
This logic continues for all model components, see [Examples](@ref Examples).
All model components are also subtype (i.e. `<:`) of some abstract supertype,
this feature can be made use of to redefine
not just constant parameters of existing model components but also
to define new ones. This is more elaborated in [Extending SpeedyWeather](@ref).

A model in SpeedyWeather.jl is understood as a collection of model components that
roughly belong into one of three groups. 

1. Components (numerics, dynamics, or physics) that have parameters, possibly contain some data (immutable, precomputed) and some functions associated with them. For example, a forcing term may be defined through some parameters, but also require precomputed arrays, or data to be loaded from file and a function that instructs how to calculate this forcing on every time step.
2. Components that are merely a collection of parameters that conceptually belong together. For example, `Earth` is an `AbstractPlanet` that contains all the orbital parameters important for weather and climate (rotation, gravity, etc).
3. Components for output purposes. For example, `OutputWriter` controls the NetCDF output, and `Feedback` controls the information printed to the REPL and to file.

## [Creating a model](@id create_model)

SpeedyWeather.jl currently includes 4 different models:

1. [`BarotropicModel`](@ref barotropic_vorticity_model) for the 2D barotropic vorticity equation,
2. [`ShallowWaterModel`](@ref shallow_water_model) for the 2D shallow water equations,
3. [`PrimitiveDryModel`](@ref primitive_equation_model) for the 3D primitive equations without humidity, and
4. [`PrimitiveWetModel`](@ref primitive_equation_model) for the 3D primitive equations with humidity.

Overall, all above-mentioned models are kept quite similar, but there are still fundamental
differences arising from the different equations. For example,
the barotropic and shallow water models do not have any physical
parameterizations. Conceptually you construct these different models with

```julia
spectral_grid = SpectralGrid(trunc=..., ...)
component1 = SomeComponent(spectral_grid, parameter1=..., ...)
component2 = SomeOtherComponent(spectral_grid, parameter2=..., ...)
model = BarotropicModel(spectral_grid; all_other_components..., ...)
```
or `model = ShallowWaterModel(spectral_grid; ...)`, etc.

## [Model initialization](@id initialize)

In the previous section the model was created, but this is conceptually
just gathering all its components together. However, many components
need to be initialized. This step is used to precompute arrays,
load necessary data from file or to communicate those between components.
Furthermore, prognostic and diagnostic variables are allocated.
It is (almost) all that needs to be done before the model can be run
(exception being the output initialization). Many model components
have a `initialize!` function associated with them that it executed here.
However, from a user perspective all that needs to be done here is
```@example howto
simulation = initialize!(model)
```
and we have initialized the `ShallowWaterModel` we have defined earlier.
As `initialize!(model)` also initializes the prognostic (and diagnostic)
variables, it also initializes the clock in
`simulation.prognostic_variables.clock`. To initialize with a specific
time, do
```@example howto
simulation = initialize!(model, time=DateTime(2020,5,1))
simulation.prognostic_variables.clock.time
```
to set the time to 1st May, 2020 (but you can also do that manually).
This time is used by components that depend on time, e.g. the solar
zenith angle calculation.

After this step you can continue to tweak your model setup but note that
some model components are immutable, or that your changes may not be
propagated to other model components that rely on it. But you can, for
example, change the output time step like so
```@example howto
simulation.model.output.output_dt = Second(3600)
```
Now, if there's output, it will be every hour. Furthermore the initial
conditions can be set with the `initial_conditions` model component
which are then set during `initialize!(::AbstractModel)`, but you can also
change them now, before the model runs 
```@example howto
simulation.prognostic_variables.vor[1][1, 1] = 0
```
So with this we have set the zero mode of vorticity of the first (and only)
layer in the shallow water to zero. Because the leapfrogging is a 2-step
time stepping scheme we set here the first. As it is often tricky
to set the initial conditions in spectral space, it is generally advised
to do so through the `initial_conditions` model component.

## [Run a simulation](@id run)

After creating a model, initializing it to a simulation, all that is left
is to run the simulation.
```@example howto
run!(simulation)
```
By default this runs for 10 days without output and returns a
[unicode plot](https://github.com/JuliaPlots/UnicodePlots.jl)
of surface relative vorticity (see [Shallow water with mountains](@ref) for more
on this). Now `period` and `output` are the only options to change, so with
```@example howto
model.output.id = "test" # hide
run!(simulation, period=Day(5), output=true)
```
You would continue this simulation (the previous `run!` call already integrated
10 days!) for another 5 days and storing default [NetCDF output](@ref).