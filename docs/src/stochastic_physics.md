# Stochastic physics

Stochastic physics introduces stochasticity into the parameterizations of physical process.
There is conceptually several classes of how this can be done

    - Stochastic perturbations of the tendencies
    - Stochastic parameter perturbations
    - Stochastic perturbations to the inputs of parameterizations

all of these use random numbers created from some random processes
(white noise or with some autocorrelation in space or time or sampled from some other distribution)
and are applied additive or multiplicative.


## Stochastic physics implementations

Currently implemented is

```@example radiation
using InteractiveUtils # hide
using SpeedyWeather
subtypes(SpeedyWeather.AbstractStochasticPhysics)
```

## Stochastically perturbed parameterization tendencies (SPPT)

SPPT is based on the idea that the [dynamics](@ref primitive_equation_model) ``D`` have a higher certainty but that
the [Parameterization](@ref) ``P`` are more uncertain and hence any stochasticity should scale
with the size of the tendencies coming from the parameterizations. 
Conceptually an atmospheric state ``\mathbf{x}`` is integrated in time ``t`` as

```math
\frac{\partial \mathbf{x}}{\partial t} = D(\mathbf{x}) + P(\mathbf{x}, t; p)
```

with the parameterizations being a function of the atmospheric state (column only for the single column 
parameterizations in SpeedyWeather), time ``t`` and parameters ``p``. Now SPPT changes this to

```math
\frac{\partial \mathbf{x}}{\partial t} = D(\mathbf{x}) + (1+r)P(\mathbf{x}, t; p)
```

with ``r \in [-1, 1]`` being a random value changing with time and space but not in the vertical
(but you could define some tapering). The idea to constrain ``r \in [-1, 1]`` is such that

- the sign up the tendencies ``P`` is never reversed
- the tendencies ``P`` are at most 2x larger after perturbation

So the stochasticity is multiplicative, in regions/times when the parameterizations are
small, so is the perturbation, but with larger tendencies form parameterizations the
perturbation is also stronger; satisfying the conceptual idea of uncertainty from the
parameterizations scaling with the size of those.

You can use SPPT as follows. Start by defining a `random_process`, this controls how
the random values ``r`` are created

```@example SPPT
using SpeedyWeather

spectral_grid = SpectralGrid()
random_process = SpectralAR1Process(spectral_grid)
```

and you could change the length scale via `wavenumber` the `time_scale` or set
a `seed` for reproducibility. `seed=0` (default) means that a random seed is taken
from Julia's global random number generator. Now we define SPPT as


```@example SPPT
stochastic_physics = StochasticallyPerturbedParameterizationTendencies(spectral_grid)
```

`tapering` can be used to change vertically the amplitude of `r`, e.g
`tapering = σ < 0.8 ? 1 : 1 - (σ - 0.8)/0.2` (in [Sigma coordinates](@ref)) would reduce
the SPPT perturbation towards the surface. A tapering ``\tau(\sigma)`` is applied like
``(1 + \tau r)``, where ``r = r(\lambda, \varphi, t)`` is a function of horizontal
coordinates longitude ``\lambda``, latitude ``\varphi`` and time ``t`` only.

Now we pass these on to the model constructor and run a simulation

```@example SPPT
model = PrimitiveWetModel(spectral_grid; random_process, stochastic_physics)
simulation = initialize!(model)
run!(simulation, period=Day(10))

# surface humidity, kg/kg -> g/kg
humid = simulation.diagnostic_variables.grid.humid_grid[:, 8]*1000

heatmap(model.orography.orography, title="Surface humidity [g/kg], with SPPT", colormap=:oslo)
save("humid_sppt.png", ans) # hide
nothing # hide
```
![Surface humidity with SPPT](humid_sppt.png)

if we now run the same simulation again (but with a different seed for the random process,
that's randomly drawn at `initialize!`)

```@example SPPT
simulation = initialize!(model)
run!(simulation, period=Day(10))

# surface humidity, kg/kg -> g/kg
humid = simulation.diagnostic_variables.grid.humid_grid[:, 8]*1000

heatmap(model.orography.orography, title="Surface humidity [g/kg], with SPPT, other seed", colormap=:oslo)
save("humid_sppt2.png", ans) # hide
nothing # hide
```
![Surface humidity with SPPT](humid_sppt2.png)

