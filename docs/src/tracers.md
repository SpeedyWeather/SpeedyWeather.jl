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
source term should increase the aerosol concentration wherever and whenever there is
a wildfire. And also a sink term should be added representing aerosols being washed out
in rainfall, or deposited on the ground. However, aerosols should be considered
_active_ and not _passive_ if they influence the radiation and hence the temperature
which couples the tracer equation above two-way with the other equations.
Active or passive tracers are treated equally in SpeedyWeather. By default a
tracer is passive, but if a [forcing](@ref "Custom forcing and drag") or [parameterization](@ref "Parameterizations") is defined that depends
on the tracer it becomes active, affecting the flow. A tracer definition
in itself therefore does not make distinction between active or passive
tracers but a forcing/parameterization definition can make existing tracers
active. (See also [Particle advection](@ref) for advecting individual
particles instead of continuous tracer fields.)

# Eulerian advection

Numerically we solve ``Dq/Dt`` as

```math
\frac{\partial q}{\partial t} = -\nabla\cdot(\mathbf{u}q) + q\mathcal{D} - W(q)
```

with ``\mathcal{D}`` being the horizontal divergence (see [Primitive equations](@ref primitive_equation_model)).
``\mathbf{u} = (u, v)`` is here the horizontal wind only because ``W(q)`` is the [Vertical advection](@ref)
operator (zero for 2D models). The products ``\mathbf{u}q, q\mathcal{D}`` are computed in grid space,
transformed to spectral space, where the divergence is taken for the former and added to the latter.
The [time stepping](@ref time_stepping) is then performed in spectral space.

## Define tracers

For every tracer in SpeedyWeather the tracer advection equation as outlined above is solved.
Tracers are defined when the `model` is constructed, through the `tracers` keyword,
as a `NamedTuple` of `Tracer`s

```@example tracers
using SpeedyWeather
spectral_grid = SpectralGrid(truncation=64, nlayers=1)

# define three tracers called abc, co2 and ch4
tracers = (; abc = Tracer(), co2 = Tracer(), ch4 = Tracer())
model = ShallowWaterModel(spectral_grid; tracers)
model.tracers
```

The key is the name of the tracer, e.g. `abc`, and there are no two tracers with the same
name (just like the keys of any `NamedTuple`). `model.tracers` will always give you an overview
of which tracers are defined. Without the `tracers` keyword a model has no tracers,
which is the same as `tracers = (;)` (an empty `NamedTuple`) or `tracers = nothing`.

The tracers of a model are fixed once it is constructed, because its type depends on them.
To use other tracers, construct another `model`, e.g. with only the tracer `abc`

```@example tracers
model = ShallowWaterModel(spectral_grid; tracers = (; abc = Tracer()))
```

The tracers then appear in the variables, which are allocated by `initialize!(model)`

```@example tracers
simulation = initialize!(model)
simulation.variables
```

Note that the tracer appear in several arrays, in the prognostic variables,
a tendency and so on. After initialization `simulation.variables`
is an immutable collection of (mutable) arrays (see [Variables](@ref)). If you do want to start or
pause a tracer advection at any time you have to activate or deactivate them,
see below.

## (De)activate tracers

A tracer has a field `active`
```@example tracers
Tracer()
```
which can be changed any time.
Do not confuse this with the discussion of active vs passive,
categorising whether a tracer impacts the flow or not.
Active versus deactivated here is solely used to describe
whether its time evolution is (temporarily) _deactivated_.
An activated tracer (the default) is advected,
a deactivated tracer does not change in time
(=frozen) but continues to exist and all its variables remain in place.
`Tracer(active = false)` defines a tracer that is deactivated from the start.
You can (de)activate a tracer by its name with

```@example tracers
activate!(model, :abc)
deactivate!(model, :abc)
```

which is equivalent to `model.tracers.abc.active = true` (default, or `false`)
and also equivalent to (de)activating them in the `simulation` instead,
i.e. `activate!(simulation, :abc)`.

## Set tracers

Tracers can be set to values by using the `set!` function, which
can take scalars, fields (spectral or grid) or functions as arguments.
But note that tracers are placed in the namespace `tracer` which
you have to pass on as a keyword argument (see [Setting variables](@ref))
e.g.

```@example tracers
set!(simulation, abc = 1, namespace = :tracers)

(; GridVariable3D, nlat_half, nlayers) = spectral_grid
set!(simulation, abc=randn(GridVariable3D, nlat_half, nlayers), namespace = :tracers)
set!(simulation, abc=(λ, φ, σ) -> exp(-(λ-180)^2/10^2), namespace = :tracers)
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
the [vertical coordinate](@ref vertical_coordinates_page) `σ` even if it's independent of it.

## Tracer visualisation

Let us illustrate some tracer advection in practice

```@example tracers
using SpeedyWeather
spectral_grid = SpectralGrid(truncation=86, nlayers = 1)
model = ShallowWaterModel(spectral_grid; tracers = (; abc = Tracer()))
simulation = initialize!(model)

# add and set tracer and run a 0-day simulation
set!(simulation, abc = (λ, φ, σ) -> exp(-(λ-180)^2/10^2), namespace =:tracers)
run!(simulation, period = Day(0))

# visualise the initial conditions for this tracer
using CairoMakie
abc0 = get_step(simulation.variables.grid.tracers.abc)[:, 1]

heatmap(abc0, title="Tracer abc, initial conditions")
save("tracer_abc.png", ans) # hide
nothing # hide
```
![Tracer abc](tracer_abc.png)

So we started with a north-south stripe of some tracer.
`[:, 1]` is used to pull out all values `:` on the one and only
layer `1`, see [Step dimension](@ref)
The [`ShallowWaterModel`](@ref) has by default a jet in the northern
hemisphere which will advect that tracer, after some days:

```@example tracers
run!(simulation, period=Day(3))

abc1 = get_step(simulation.variables.grid.tracers.abc)[:, 1]
heatmap(abc1, title="Tracer abc, after 3 days")
save("tracer2.png", ans) # hide
nothing # hide
```
![Tracer after 3 days](tracer2.png)


## Output tracers

This follows equivalent to other [Output variables](@ref) but the
tracer name as `Symbol` has to be provided (unlike the `Tracer()` of the model,
the output variable `OutputVariables.Tracer` takes a name)

```@example tracers
add!(model, OutputVariables.Tracer(:abc))
```

and other keyword arguments like `long_name::String` and `units::String`
can be passed on too. 
