# Tracer advection

A tracer is a property ``q`` of the fluid that is advected with the flow
``\mathbf{u} = (u, v, w)`` (for 3D, ``\mathbf{u} = (u, v)`` for 2D)
without changes along that trajectory.

```math
\frac{Dq}{Dt} = \frac{\partial q}{\partial t} + (\mathbf{u} \cdot \nabla)q = 0
```

If the tracer ``q`` does not impact the flow it is considered _passive_.
Humidity, for example, is an _active_ tracer as it changes the [Geopotential](@ref)
(and therefore the pressure gradient force) through the [Virtual temperature](@ref).
A tracer is conserved in the absence of sources or sinks (the zero on the right-hand side above).
Aerosols from wildfires might be considered to be a passive tracer, but the
source term should increase the aerosol concentration whereever and whenever there is 
a wildfire. And also a sink term should be added representing aerosols being washed out
in rainfall, or deposited on the ground. However, aersols should be considered
_active_ and not _passive_ if they influence the radiation and hence the temperature
which couples the tracer equation above two-way with the other equations.

# Eulerian advection

Numerically we solve ``Dq/Dt`` as 

```math
\frac{\partial q}{\partial t} = -\nabla\cdot(\mathbf{u}q) + q\mathcal{D} - W(q)
```

with ``\mathcal{D}`` being the horizontal divergence (see [Primitive equations](@ref primitive_equation_model)).
``\mathbf{u} = (u, v)`` is here the horizontal wind only because ``W(q)`` is the [Vertical advection](@ref)
operator (zero for 2D models). The products ``\mathbf{u}q, q\mathcal{D}`` are computed in grid space,
transformed to spectral space, where the divergence is taken for the former and added to the latter.
The time stepping is then performed in spectral space.

## Add/delete tracers

For every tracer in SpeedyWeather the tracer advection equation as outlined above is solved.
One can add a new tracer to the `model` _before_ it is initialized to a `simulation`

```@example tracers
using SpeedyWeather
spectral_grid = SpectralGrid(trunc=63, nlayers=1)
model = ShallowWaterModel(spectral_grid)

# add a tracer called :abc
add!(model, Tracer(:abc))
```

This returns `model.tracers`, a dictionary, which will always give you an overview of which tracers
are defined. Tracers are defined through a `key::Symbol` for which we use `Symbol`
(not strings, because Symbols are immutable). We just wrap the `key` here in
`Tracer` to define a tracer. You can add more tracers

```@example tracers
add!(model, Tracer(:co2), Tracer(:ch4))
```

or delete them again

```@example tracers
delete!(model, Tracer(:co2))
delete!(model, Tracer(:ch4))
```

Tracers are just defined through their `key`, e.g. `:co2`, so while you can do
`tracer1 = Tracer(:co2)` and `tracer2 = Tracer(:co2)`, they will be considered
the same tracer -- and no two tracers with the same key can exist inside `model`
(and `simulation`).

You can also add a tracer to a `simulation`, i.e. after the model is initialized.

```@example tracers
simulation = initialize!(model)
add!(simulation, Tracer(:xyz))
```

which will add the tracer to `model.tracers` as above but also add
it to the prognostic and diagnostic variables (otherwise done at `initialize!`).
What you should not do is add a tracer to the `model` _after_ it has been initialized.
Then you end up with an additional tracer in `model` without there
being variables for it, throwing an error. You can check that
the tracers exists in the variables with

```@example tracers
simulation.prognostic_variables
```

where both `:abc` and `:xyz` are listed. Tracers in SpeedyWeather are
based on dictionaries so the order of the tracers is arbitrary,
they are always defined by their `key` instead.

Note that a tracer can be added to or deleted from a simulation at _any time_.
So you can run a simulation, add a tracer, continute the simulation,
or delete a tracer and continue. You can also just activate them or
deactivate them, see below.

## (De)activate tracers

While a tracer is defined through its key, e.g.
```@example tracers
Tracer(:dust)
```
it also has a field `active` which can be changed any time.
An active tracer is advected, a deactivated tracer does not change in time
(=frozen) but continues to exist and all its variables remain in place.
You can (de)activate a tracer with

```@example tracers
activate!(model, Tracer(:abc))
deactivate!(model, Tracer(:abc))
```

which is equivalent to `model.tracers[:abc].active = true` (default, or `false`)
and also equivalent to (de)activating them in the `simulation` instead,
i.e. `activate!(simulation, Tracer(:abc))`.

## Set tracers

Tracers can be set to values by using the `set!` function, which
can take scalars, fields (spectral or grid) or functions as arguments,
e.g.

```@example tracers
set!(simulation, abc=1)

(; GridVariable3D, nlat_half, nlayers) = spectral_grid
set!(simulation, abc=randn(GridVariable3D, nlat_half, nlayers))
set!(simulation, abc=(λ, φ, σ) -> exp(-(λ-180)^2/10^2))
```

The first one sets `abc` to a global constant (not super exciting),
the second to some random values on a grid (transforms automatically!),
and the third sets the tracer to a Gaussian ridge that runs through
the Pacific (see [Tracer visualisation](@ref) below).

For more examples how to use `set!` see [Changing orography manually](@ref),
[Manual land-sea mask](@ref), and
[Rossby-Haurwitz wave in a BarotropicModel](@ref).
But note that because we are setting a (in general) 3D variable here
the vertical dimension must align: Hence `nlayers` for the grid,
and the anonymous function must take three arguments, including
the vertical coordinate `σ` even if it's independent of it.

## Tracer visualisation

Let us illustrate some tracer advection in practice

```@example tracers
using SpeedyWeather
spectral_grid = SpectralGrid(trunc=85, nlayers=1)
model = ShallowWaterModel(spectral_grid)
simulation = initialize!(model)

# add and set tracer and run a 0-day simulation
add!(simulation, Tracer(:abc))
set!(simulation, abc = (λ, φ, σ) -> exp(-(λ-180)^2/10^2))
run!(simulation, period=Day(0))

# visualise the initial conditions for this tracer
using CairoMakie
abc0 = simulation.diagnostic_variables.grid.tracers_grid[:abc][:, 1]

heatmap(abc0, title="Tracer abc, initial conditions")
save("tracer_abc.png", ans) # hide
nothing # hide
```
![Tracer abc](tracer_abc.png)

So we started with a north-south stripe of some tracer.
`[:, 1]` is used to pull out all values `:` on the one and only
layer `1`.
The `ShallowWaterModel` has by default a jet in the northern
hemisphere which will advect that tracer, after some days:

```@example tracers
run!(simulation, period=Day(3))

abc1 = simulation.diagnostic_variables.grid.tracers_grid[:abc][:, 1]
heatmap(abc1, title="Tracer abc, after 3 days")
save("tracer2.png", ans) # hide
nothing # hide
```
![Tracer after 3 days](tracer2.png)


## Output tracers

more to come...