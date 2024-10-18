# Initial conditions

The following showcases some examples of how to set the initial conditions
for the prognostic variables in SpeedyWeather.jl. In essence there
are three ways to do this

1. Change the arrays in `simulation.prognostic_variables`
2. Use the `set!` function
3. Set the `initial_conditions` component of a model

where 1 is a rather low-level, the `set!` function builds a convenient
interface around 1 (so you don't have to know about grid and spectral space details)
and 3 collects method 1 or 2 (or a combination of both) into a single struct
to "save" some initial conditions for one or several variables.
This lets you use predefined (inside SpeedyWeather or externally) initial conditions
as easy as `initial_conditions = RossbyHaurwitzWave()`.
Let us illustrate this with some examples where we will refer back those
methods simply as 1, 2, 3.

## Rossby-Haurwitz wave in a BarotropicModel

We define a `BarotropicModel` of some resolution but keep all its components
as default

```@example haurwitz
using SpeedyWeather
spectral_grid = SpectralGrid(trunc=63, nlayers=1)
model = BarotropicModel(; spectral_grid)
simulation = initialize!(model)
```

Now `simulation.prognostic_variables` contains already some
initial conditions as defined by `model.initial_conditions` (that's method 3).
Regardless of what those are, we can still mutate them
before starting a simulation. The Rossby-Haurwitz wave[^Williamson92] is
defined as an initial condition for vorticity ``\zeta`` (which is the sole
prognostic variable in the barotropic vorticity model) as

```math
ζ(λ, θ) = 2ω*\sin(θ) - K*\sin(θ)*\cos(θ)^R*(R^2 + 3R + 2)*\cos(R*λ)
```
with longitude ``\lambda`` and latitude ``\theta``. The parameters
are order ``R = 4``, frequencies ``\omega = 7.848e-6s^{-1}, K = 7.848e-6s^{-1}``.
Now setting these initial conditions is as simple as

```@example haurwitz
R = 4
ω = 7.848e-6
K = 7.848e-6

ζ(λ, θ, σ) = 2ω*sind(θ) - K*sind(θ)*cosd(θ)^R*(R^2 + 3R + 2)*cosd(R*λ)
set!(simulation, vor=ζ)
```

with only two difference from the mathematical notation. (1) SpeedyWeather's
coordinates are in degrees, so we replaced ``\sin, \cos`` with `sind` and `cosd`;
and (2) To generalise to vertical coordinates, the function `ζ(λ, θ, σ)` takes
exactly three arguments, with `σ` denoting the vertical [Sigma coordinates](@ref).
This is important so that we can use the same definition of initial conditions
for the 2D barotropic vorticity model also for the 3D primitive equations.

Some authors filter out low values of spectral vorticity with some cut-off amplitude
``c = 10^{-10}``, just to illustrate how you would do this (example for method 1)

```@example haurwitz
c = 1e-10       # cut-off amplitude

# 1 = first leapfrog timestep of spectral vorticity
vor = simulation.prognostic_variables.vor[1]
low_values = abs.(vor) .< c
vor[low_values] .= 0
nothing # hide
```
which is just treating `vor` as an array of something and tweaking the values within!

Let us illustrate these initial conditions. `set!` will set the initial conditions
in spectral space, taking care of the transform from the equation defined
in grid coordinates. So to show vorticity again in grid space we transform
back

```@example haurwitz
# [1] for first leapfrog time step, [:, 1] for all values on first layer
vor = simulation.prognostic_variables.vor[1][:, 1]
vor_grid = transform(vor)

using CairoMakie
heatmap(vor_grid, title="Relative vorticity [1/s] of Rossby-Haurwitz wave")

save("haurwitz.png", ans) # hide
nothing # hide
```
![Rossby-Haurwitz wave](haurwitz.png)

That's the Rossby-Haurwitz wave! This wave is supposed to travel
(without changing its shape) eastward around the globe, so let us run
a simulation for some days

```@example haurwitz
run!(simulation, period=Day(10))

# a running simulation always transforms spectral variables
# so we don't have to do the transform manually but just pull 
# layer 1 (there's only 1) from the diagnostic variables
vor = simulation.diagnostic_variables.grid.vor_grid[:, 1]

heatmap(vor, title="Relative vorticity [1/s], Rossby Haurwitz wave after 10 days")
save("haurwitz_day10.png", ans) # hide
nothing # hide
```
![Rossby-Haurwitz wave after day 10](haurwitz_day10.png)


## Rossby-Haurwitz wave in primitive equations

We could apply the same to set the Rossby-Haurwitz for a primitive equation
model, but we have also already defined `RossbyHaurwitzWave` as
`<: AbstractInitialConditions` so you can use that directly, regardless
the model. Note that this definition currently only includes vorticity
not the initial conditions for other variables. Williamson et al. 1992
define also initial conditions for height/geopotential to be used
in the shallow water model (eq. 146-149) -- those are currently not included,
so the wave may not be as stable as its supposed to be.

The following shows that you can set the same `RossbyHaurwitzWave` initial
conditions also in a `PrimitiveDryModel` (or `Wet`) but you probably
also want to set initial conditions for temperature and pressure
to not start at zero Kelvin and zero pressure. Also no orography,
and let's switch off all physics parameterizations with `physics=false`.

```@example haurwitz
spectral_grid = SpectralGrid(trunc=42, nlayers=8)
initial_conditions = InitialConditions(
                        vordiv=RossbyHaurwitzWave(),
                        temp=JablonowskiTemperature(),
                        pres=PressureOnOrography())

orography = NoOrography(spectral_grid)
model = PrimitiveDryModel(; spectral_grid, initial_conditions, orography, physics=false)
simulation = initialize!(model)
run!(simulation, period=Day(10))
nothing # hide
```

Note that we chose a lower resolution here (T42) as we are simulating
8 vertical layers now too. Let us visualise the surface vorticity
(`[:, 8]` is on layer )

```@example haurwitz
vor = simulation.diagnostic_variables.grid.vor_grid[:, 8]
heatmap(vor, title="Relative vorticity [1/s]")

save("haurwitz_primitive.png", ans) # hide
nothing # hide
```
![Rossby-Haurwitz wave in primitive equations](haurwitz_primitive.png)

As you can see the actual Rossby-Haurwitz wave is not as stable anymore
(because those initial conditions are not a stable solution of the primitive equations)
and so the 10-day integration looks very different from the barotropic model!

## References

[^Williamson92]: DL Williamson, JB Drake, JJ Hack, R Jakob, PN Swarztrauber, 1992. *A standard test set for numerical approximations to the shallow water equations in spherical geometry*, **Journal of Computational Physics**, 102, 1, DOI: [10.1016/S0021-9991(05)80016-6](https://doi.org/10.1016/S0021-9991(05)80016-6)