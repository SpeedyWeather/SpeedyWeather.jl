# Stochastic physics

Stochastic physics introduces stochasticity into the parameterizations of physical process.
There is conceptually several classes of how this can be done

- Stochastic perturbations of the tendencies
- Stochastic parameter perturbations
- Stochastic perturbations to the inputs of parameterizations

all of these use random numbers created from some random processes
(white noise or with some autocorrelation in space or time or sampled from some other distribution)
and are applied additive or multiplicative. 

Separate from this is the [`StochasticLeapfrog`](@ref stochastic_leapfrog) that integrates the complete model as a SDE with a Wiener noise process added to the tendencies.

## [Random processes](@id random_processes)

The random numbers for stochastic physics come from a random process, a model component
`model.random_process` (`nothing` by default, i.e. no random numbers are drawn).
Once per time step the random process draws a new random pattern after `update_random!` is called. The stochastic
physics then read it from `simulation.variables.grid.random_pattern`. Currently implemented are

```@example random_process
using InteractiveUtils # hide
using SpeedyWeather
subtypes(SpeedyWeather.AbstractRandomProcess)
```

All of them have the following options

- `standard_deviation` ``\sigma`` of the pattern at every grid point, by default ``\sigma = 1/3``
  so that values beyond ``\pm 1`` (``3\sigma``) are rare,
- `seed` for reproducibility, `seed = 0` (default) takes a new random seed at every `initialize!`,
- `n_patterns` the number of independent patterns (default 1),
- `draw_on_cpu = true` to draw the random numbers on the CPU also if the model runs on GPU
  (slower, but the same random numbers on every architecture).

### First-order autoregressive process (AR1)

`SpectralAR1Process` is a random process with memory. Every time step ``\Delta t`` it updates the
spherical harmonic coefficients ``\xi_{lm}`` of the pattern as

```math
\xi_{lm}(t + \Delta t) = a \, \xi_{lm}(t) + F_l \, \eta_{lm}, \qquad a = \exp\left(-\frac{\Delta t}{\tau}\right)
```

with new (complex) standard normal random numbers ``\eta_{lm}`` every time step. Its  correlation decays exponentially with the `time_scale` ``\tau`` (default 6 hours). The noise
factors ``F_l`` decrease with the degree ``l`` as

```math
F_l \propto \exp\left(-\frac{l(l+1)}{2k(k+1)}\right)
```

so that the pattern is also correlated in space: scales larger than the `wavenumber` ``k``
(default 12) keep their variance, smaller scales are damped. Features have a size of about
``R/k`` with ``R`` the radius of the planet, e.g. about 500 km for ``k = 12`` on Earth.
The amplitude of ``F_l`` is chosen (following ECMWF Technical Memorandum 598) such that the pattern
has the standard deviation ``\sigma`` at every grid point once it is statistically stationary
(starting from zero that takes a few ``\tau``). `wavenumber = 0` stands for ``k \to \infty``,
the same ``F_l`` for all degrees ``l`` (spectrally flat).
Neither the global mean (``l = 0``) nor the last degree (the truncation) get noise.
The ``m = 0`` coefficients are real, with their real part scaled by ``\sqrt{2}``, so that the
pattern is statistically homogeneous on the sphere. After the transform to grid space the pattern is clamped into `clamp = (-1, 1)`
(default), as required for [SPPT](@ref "Stochastically perturbed parameterization tendencies (SPPT)").

### White noise in spectral space

`SpectralWhiteNoiseProcess` is the AR1 process above with `time_scale = 0`, so that
``a = 0``: it has no memory, the pattern is drawn anew every time step
independent of the previous one (white in time)

```math
\xi_{lm}(t) = F_l \, \eta_{lm}
```

With `wavenumber = 0` (default) ``F_l`` is the same for all degrees ``l``, the pattern is
then also white in space (spectrally flat). With a `wavenumber` ``k > 0`` it uses the same ``F_l``
as the AR1 process above, so that the pattern is white in time but correlated in space.
As for the AR1 process, neither the global mean (``l = 0``) nor the last degree get noise
and the ``m = 0`` coefficients are real. There is no clamping by default (`clamp = (-Inf, Inf)`) to keep the pattern Gaussian,
use e.g. `clamp = (-1, 1)` together with SPPT.

### White noise in grid space

`GridWhiteNoiseProcess` draws an independent normal random number with standard deviation
``\sigma`` at every grid point and every time step. It is drawn directly in grid space,
no spectral transform is involved. As grid cells generally differ in area, the pattern is the
same at every grid point but not statistically homogeneous on the sphere.

## Stochastic physics implementations

Currently implemented is

```@example radiation
using InteractiveUtils # hide
using SpeedyWeather
subtypes(SpeedyWeather.AbstractStochasticPhysics)
```

## Stochastically perturbed parameterization tendencies (SPPT)

SPPT is based on the idea that the [dynamics](@ref primitive_equation_model) ``D`` have a higher certainty but that
the [Parameterizations](@ref) ``P`` are more uncertain and hence any stochasticity should scale
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
the random values ``r`` are created (see [Random processes](@ref random_processes))

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

`tapering` (by default an
[anonymous function](https://docs.julialang.org/en/v1/manual/functions/#man-anonymous-functions)
hence it looks like `var"#269#273" = #269` as the compiler assigns a "name")
can be used to change vertically the amplitude of `r`, e.g
`tapering = σ -> σ < 0.8 ? 1 : 1 - (σ - 0.8)/0.2` (in [Sigma coordinates](@ref sigma_coordinates_physics)) could
be passed on a (keyword) arugment to reduce the SPPT perturbation towards the surface.
A tapering ``\tau(\sigma)`` is applied like ``(1 + \tau r)``, where ``r = r(\lambda, \varphi, t)``
is a function of horizontal coordinates longitude ``\lambda``, latitude ``\varphi`` and time ``t`` only.

Now we pass these on to the model constructor and run a simulation

```@example SPPT
model = PrimitiveWetModel(spectral_grid; random_process, stochastic_physics)
simulation = initialize!(model)
run!(simulation, period=Day(10))

# surface humidity, kg/kg -> g/kg
k = spectral_grid.nlayers   # surface layer
humid = get_step(simulation.variables.grid.humidity)[:, k]*1000

using CairoMakie
heatmap(humid, title="Surface humidity [g/kg], with SPPT", colormap=:oslo)
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
humid = get_step(simulation.variables.grid.humidity)[:, k]*1000

heatmap(humid, title="Surface humidity [g/kg], with SPPT, other seed", colormap=:oslo)
save("humid_sppt2.png", ans) # hide
nothing # hide
```
![Surface humidity with SPPT](humid_sppt2.png)

## Implementing your stochastic parameterizations with their own random processes

If you want to implement a new stochastic parameterization with an independent noise realization, you can just make a random process one of the fields of your parameterization. In this case however, you have to manually update/draw new numbers with `update_random!` within your parameterization,
e.g. in its `parameterization!(vars, ...)` method that is called once per time step before the
column parameterizations. A minimal example that heats or cools the surface layer randomly
(see [Define your own parameterizations](@ref) for the general interface)

```@example own_random_process
using SpeedyWeather

struct StochasticHeating{NF, RP} <: SpeedyWeather.AbstractStochasticPhysics
    amplitude::NF           # heating amplitude [K/s]
    random_process::RP      # its own random process
end

StochasticHeating(SG::SpectralGrid; amplitude = 1.0e-4) =
    StochasticHeating(SG.NF(amplitude), GridWhiteNoiseProcess(SG))

# its own pattern to draw into, separate from the `random_pattern` of `model.random_process`
SpeedyWeather.variables(::StochasticHeating) =
    (ParameterizationVariable(:heating_noise, SpeedyWeather.Grid2D(), desc = "Random heating pattern", units = "1"),)

# initialize (and seed) the random process with the parameterization
SpeedyWeather.initialize!(heating::StochasticHeating, model::PrimitiveEquation) =
    SpeedyWeather.initialize!(heating.random_process, model)

# draw new random numbers once per time step
SpeedyWeather.parameterization!(vars::Variables, heating::StochasticHeating, model::PrimitiveEquation) =
    SpeedyWeather.update_random!(vars.parameterizations.heating_noise, heating.random_process)

# then use them in every column
Base.@propagate_inbounds function SpeedyWeather.parameterization!(ij, vars::Variables, heating::StochasticHeating, model)
    surface = size(vars.tendencies.grid.temperature, 2)
    vars.tendencies.grid.temperature[ij, surface] += heating.amplitude * vars.parameterizations.heating_noise[ij]
end

spectral_grid = SpectralGrid(truncation = 32, nlayers = 8)
stochastic_physics = StochasticHeating(spectral_grid)
model = PrimitiveWetModel(spectral_grid; stochastic_physics)
simulation = initialize!(model)
run!(simulation, period = Day(1))
```

Passed on as `stochastic_physics` it is initialized and executed (after all other parameterizations)
by default. It does not need a `model.random_process`, and with one, both draw independent noise
(unless given the same `seed`).


