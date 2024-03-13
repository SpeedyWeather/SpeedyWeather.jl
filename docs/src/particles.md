# Particle advection

All SpeedyWeather.jl models support particle advection. Particles are objects without mass or volume
at a location ``\mathbf{x} = (\lambda, \theta, \sigma)`` (longitude ``\lambda``, latitude ``\theta``,
vertical sigma coordinate ``\sigma``, see [Sigma coordinates](@ref)) that are moved with the wind
``\mathbf{u}(\mathbf{x})`` at that location. The location of the ``p``-th particle changes as follows

```math
\frac{d \mathbf{x}_p}{d t} = \mathbf{u}(\mathbf{x}_p)
```

This equation applies in 2D, i.e. ``\mathbf{x} = (\lambda, \theta)`` and ``\mathbf{u} = (u, v)`` or
in 3D, but at the moment only 2D advection is supported. In the [Primitive equation model](@ref primitive_equation_model)
the vertical layer on which the advection takes place has to be specified. It is therefore not
advected with the vertical velocity but maintains a constant pressure ratio compared to the
surface pressure (``\sigma`` is constant).

## Discretization of particle advection

The particle advection equation has to be discretized to be numerically solved. While the
particle location can generally be anywhere on the sphere, the velocity ``\mathbf{u}``
is only available on the discrete grid points of the simulation, such that
``\mathbf{u}(\mathbf{x}_p)`` requires an interpolation in order to obtain a velocity
at the particles location ``\mathbf{x}_p`` to move it around. Dropping the subscript
``p`` in favour a subscript ``i`` denoting the time step, with Euler forward 
the equation can be discretized as 

```math
\mathbf{x}_{i+1} = \mathbf{x}_i + \Delta t~\mathbf{u}_i (\mathbf{x}_i)
```

Meaning we have used the velocity field at both departure time ``i`` and departure
location ``\mathbf{x}_i`` to update a particle's location which makes this scheme
first order accurate. But only a single interpolation of the velocity field,
which, in fact, is one per dimension, is necessary. Note that the time step
``\Delta t`` here and the time step to solve the dynamics do not have to be
identical. We could use a larger time step for the particle advection then to
solve the dynamics inside the model, and because the stability criteria for these
equations are different, one is encouraged to do so. Also because the particles
are considered passive, meaning that their location does not influence the
other prognostic variables.

We can write down a more accurate scheme at the cost of a second interpolation step.
The Heun method, also called predictor-corrector is 2nd order accurate and uses an average of the
velocity at departure time ``i`` and location ``\mathbf{x}_i`` and at a
(_predicted_ meaning preliminary) arrival point ``x^\star_{i+1}`` and arrival time ``i+1``.

```math
\begin{aligned}
\mathbf{x}^\star_{i+1} &= \mathbf{x}_i + \Delta t~\mathbf{u}_i (\mathbf{x}_i) \\
\mathbf{x}_{i+1} &= \mathbf{x}_i + \frac{\Delta t}{2}~\left(
    \mathbf{u}_i (\mathbf{x}_i) + \mathbf{u}_{i+1} (\mathbf{x}^\star_{i+1})\right)
\end{aligned}
```

Because we don't have ``\mathbf{u}_{i+1}`` available at time ``i``, we perform this
integration retrospectively, i.e. if the other model dynamics have reached time ``i+1``
then we let the particle advection catch up by integrating them from ``i`` to ``i+1``.
This, however, requires some storage of the velocity ``\mathbf{u}_i`` at the previous
advection time step. Remember that this does not need to be the time step for the
momentum equations and could be much further in the past. We could either store
``\mathbf{u}_i`` as a grid-point field or only its interpolated values. In the case
of fewer particles than grid points the latter is more efficient and this is also
what we do in SpeedyWeather. Let square brackets ``[]`` denote
an interpolation then we perform the interpolation
``\mathbf{u}_i [\mathbf{x}_i]`` that's required to step from ``i`` to ``i+1`` already
on the time step that goes from ``i-1`` to ``i``. 

```math
\begin{aligned}
\mathbf{x}^\star_{i+1} &= \mathbf{x}_i + \Delta t~\mathbf{u}_i (\mathbf{x}_i) \\
\mathbf{x}_{i+1} &= \mathbf{x}_i + \frac{\Delta t}{2}~\left(
    \mathbf{u}_i (\mathbf{x}_i) + \mathbf{u}_{i+1} [\mathbf{x}^\star_{i+1}]\right) \\
\mathbf{u}_{i+1} (\mathbf{x}_{i+1}) &=  \mathbf{u}_{i+1} [\mathbf{x}_{i+1}]
\end{aligned}
```
Denoted here as the last line with the left-hand side becoming the last term of the
first line in the next time step ``i+1 \to i+2``. Now it becomes clearer that there
are two interpolations required on every time step.

We use for horizontal coordinates degrees, such that we need to scale the time
step ``\Delta`` with ``\frac{360˚}{2\pi R}`` (radius ``R``) for advection in
latitude and with ``\frac{360˚}{2\pi R \cos(\theta)}`` for advection in longitude
(because the distance between meridians decreases towards the poles). We move the
division by the radius conveniently into the time step as are also the momentum
equations scaled with the radius, see [Radius scaling](@ref scaling).

Technically, a particle moved with a given velocity follows a
[great circle](https://en.wikipedia.org/wiki/Great-circle_distance)
in spherical coordinates. This means that

```math
\theta_{i+1} \approx \theta_i + \frac{\Delta t}{R} \frac{360}{2\pi} v_i
```

becomes a bad approximation when the time step and or the velocity are large. However,
for simplicity and to avoid the calculation of the great circle we currently do
use this to move particles with a given velocity. We essentially assume a local
cartesian coordinate system instead of the geodesics in spherical coordinates.
However, for typical time steps of 1 hour and velocities not exceeding 100 m/s
the error is not catastrophic and can be reduced with a shorter time step.
We may switch to great circle calculations in future versions.

## Create a particle

So much about the theory

A [`Particle`](@ref) at location 10˚E and 30˚N (and ``\sigma = 0``) can be created as follows,

```@example particle
using SpeedyWeather
p = Particle(lon=10, lat=30, σ=0)
p = Particle(lon=10, lat=30)
p = Particle(10, 30, 0)
p = Particle(10, 30)
```

All of the above are equivalent. Unless a keyword argument is used, longitude is the first
argument, followed by latitude (necessary), followed by ``\sigma`` (can be omitted).
Longitudes can be -180˚E to 180˚E or 0 to 360˚E, latitudes have to be -90˚N to 90˚N.
You can create a particle with coordinates outside of these ranges (and no error or
warning is thrown) but during particle advection they will be wrapped into
[0, 360˚E] and [-90˚N, 90˚N], using the [`mod(::Particle)`](@ref) function, which
is similar to the modulo operator but with the second argument hardcoded to the
coordinate ranges from above, e.g.
```@example particle
mod(Particle(lon=-30, lat=0))
```
which also takes into account pole crossings which adds 180˚ in longitude
```@example particle
mod(Particle(lon=0, lat=100))
```
as if the particle has moved across the pole. That way all real values for longitude
and latitude are wrapped into the reference range [0, 360˚E] and [-90˚N, 90˚N].

!!! info "Particles are immutable"
    Particles are implemented as immutable `struct`, meaning you cannot change their
    position by `particle.lon = value`. You have to think of them as integers or floats
    instead. If you have a particle `p` and you want to change its position to the
    Equator for example you need to create a new one `new_particle = Particle(p.lon, 0, p.σ)`.

By default `Float32` is used, but providing coordinates in `Float64` will promote the
type accordingly. Also by default, particles are _active_
which is indicated by the 2nd parametric type of `Particle`, a boolean. Active particles
are moved following the equation above, but inactive particles are not. You can
[`activate`](@ref) or [`deactivate`](@ref) a particle like so
```@example particle
deactivate(p)
```
and so
```@example particle
activate(p)
```
or check its activity by [`active(::Particle)`](@ref) returning `true` or `false`.
The zero-element of the [`Particle`](@ref) type is
```@example particle
zero(Particle)
```
and you can also create a random particle which uses a
[raised cosine distribution](https://en.wikipedia.org/wiki/Raised_cosine_distribution)
in latitude for an equal area-weighted uniform distribution over the sphere
```@example particle
rand(Particle{Float32})         # specify number format
rand(Particle{Float32, true})   # and active/inactive
rand(Particle)                  # or not (defaults used instead)
```

## Advecting particles

The [`Particle`](@ref) type can be used inside vectors, e.g.
```@example particle
zeros(Particle{Float32}, 3)
rand(Particle{Float64}, 5)
```
which is how particles are represented inside a SpeedyWeather [`Simulation`](@ref).
Note that we have not specified whether the particles inside these vectors are
active (e.g. `Particle{Float32, true}`) or inactive (e.g. `Particle{Float64, false}`)
because that would generally force all particles in these vectors to be either
active or inactive as specified such that
```@example particle
v = zeros(Particle{Float32, false}, 3)
v[1] = Particle(lon = 134.0, lat = 23)      # conversion to inactive Particle{Float32, false}
v
```
would not just convert from `Float64` to `Float32` but also
from an active to an inactive particle. In SpeedyWeather all particles can be
activated or deactivated at any time.

First, you create a [`SpectralGrid`](@ref) with the `n_particles` keyword
```@example particle
spectral_grid = SpectralGrid(n_particles = 3)
```
Then the particles live as `Vector{Particle}` inside the prognostic variables
```@example particle
model = BarotropicModel(;spectral_grid)
simulation = initialize!(model)
simulation.prognostic_variables.particles
```
Which are placed in random locations (using `rand`) initially.
In order to change these (e.g. to set the initial conditions) you do
```@example particle
simulation.prognostic_variables.particles[1] = Particle(lon=-120, lat=45)
simulation.prognostic_variables.particles
```
which sets the first particle (you can think of the index as the particle identification)
to some specified location, or you could deactivate a particle with
```@example particle
first_particle = simulation.prognostic_variables.particles[1]
simulation.prognostic_variables.particles[1] = deactivate(first_particle)
simulation.prognostic_variables.particles
```

To actually advect these particles inside a SpeedyWeather simulation we have to create
a `ParticalAdvection2D` instance that lets you control the time step used
for particle advection and which vertical layer to use in the 3D models.

```@example particle
particle_advection = ParticleAdvection2D(spectral_grid, layer = 1)
```

we choose the first (=top-most) layer although this is the default anyway. Now we can
advect our three particles we have defined above

```@example particle
model = BarotropicModel(;spectral_grid, particle_advection)
simulation = initialize!(model)
simulation.prognostic_variables.particles
```

Which are the initial conditions for our three particles. After 10 days of simulation
they have changed

```@example particle
run!(simulation, period=Day(10))
simulation.prognostic_variables.particles
```

Woohoo! We just advected some particles. This is probably not as exciting as actually
tracking the particles over the globe and being able to visualise their trajectory
which we will do in the next section

## Tracking particles

A [`ParticleTracker`](@ref) is implemented as a callback, see [Callbacks](@ref), outputting
the particle locations via netCDF. We can create it like

```@example particle_tracker
using SpeedyWeather
spectral_grid = SpectralGrid(n_particles = 100)
particle_tracker = ParticleTracker(spectral_grid, schedule=Schedule(every=Hour(3)))
```

which would output every 3 hours (the default). This output frequency might be slightly adjusted
depending on the time step of the dynamics to output every `n` time steps (an `@info` is thrown if
that is the case), see [Schedules](@ref). Further options on compression are available
as keyword arguments `ParticleTracker(spectral_grid, keepbits=15)` for example.
The callback is then added after the model is created

```@example particle_tracker
particle_advection = ParticleAdvection2D(spectral_grid)
model = ShallowWaterModel(;spectral_grid, particle_advection)
add!(model.callbacks, particle_tracker)
```

which will give it a random key too in case you need to remove it again (more on this in 
[Callbacks](@ref)). If you now run the simulation the particle tracker is called on
`particle_tracker.every_n_timesteps` and it continuously writes into `particle_tracker.netcdf_file`
which is placed in the run folder similar to other [NetCDF output](@ref). For example,
the run id can be obtained after the simulation by `model.output.id`.

```@example particle_tracker
simulation = initialize!(model)
run!(simulation, period=Day(10))
model.output.id
```
so that you can read the netCDF file with

```@example particle_tracker
using NCDatasets
run_id = "run_$(model.output.id)"                    # create a run_???? string with output id
path = joinpath(run_id, particle_tracker.file_name)  # by default "run_????/particles.nc"
ds = NCDataset(path)
ds["lon"]
ds["lat"]
```
where the last two lines are lazy loading a matrix with each row a particle and each column a time step.
You may do `ds["lon"][:,:]` to obtain the full `Matrix`. We had specified
`spectral_grid.n_particles` above and we will have time steps in this file
depending on the `period` the simulation ran for and the `particle_tracker.Δt` output
frequency. We can visualise the particles' trajectories with

```@example particle_tracker
lon = ds["lon"][:,:]
lat = ds["lat"][:,:]

using PythonPlot
ioff() # hide
fig, ax = subplots(1, 1, figsize=(10, 6))
ax.plot(lon', lat')
ax.set_xlabel("longitude")
ax.set_ylabel("latitude")
ax.set_title("Particle advection")
tight_layout() # hide
savefig("particles.png", dpi=70) # hide
nothing # hide
```
![Particle trajectories](particles.png)

Instead of providing a polished example with a nice projection we decided to keep it simple here
because this is probably how you will first look at your data too. As you can see, some particles
in the Northern Hemisphere have been advected with a zonal jet and perform some wavy motions as
the jet does too. However, there are also some horizontal lines which are automatically plotted
when a particles travels across the prime meridian 0˚E = 360˚E. Ideally you would want to use
a more advanced projection and plot the particle trajectories as geodetics. 

With Makie.jl you can do

```julia
using GeoMakie, GLMakie

fig = Figure()
ga = GeoAxis(fig[1, 1]; dest = "+proj=ortho +lon_0=19 +lat_0=50")

lines!(ga, GeoMakie.coastlines())
ga.xticklabelsvisible[] = false
ga.yticklabelsvisible[] = false

n_particles = size(lon)[1]
[lines!(ga, lon[i,:], lat[i,:]) for i in 1:n_particles]
fig
```