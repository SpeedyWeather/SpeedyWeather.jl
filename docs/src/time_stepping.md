# [Time stepping](@id time_stepping)

This page describes the mathematics and internal implementation of SpeedyWeather.jl's
time integration schemes. For how to create and configure a time stepper and pass it to
the model constructor, see [Time integration](@ref) in the Usage section.

SpeedyWeather.jl supports several time integration schemes, selected by passing a
`time_stepping` component to the model constructor:

- [`Leapfrog`](@ref leapfrog), a 2-step leapfrog scheme with a Robert-Asselin and Williams filter
  (the default for `ShallowWaterModel`, `PrimitiveDryModel` and `PrimitiveWetModel`), described below.
- [`StochasticLeapfrog`](@ref stochastic_leapfrog), the leapfrog scheme with additive Wiener noise.
- [`NCycleLorenz`](@ref ncycle), a family of semi-implicit Lorenz N-cycle schemes
  ([Hotta2016](@citep); the default for `BarotropicModel`).

All schemes share a common framework (see [Time steppers and variable steps](@ref steps)) in which
the time stepper decides, for every model component, which stored *step* of each variable to read
or write. This decouples the dynamical core and parameterizations from the time-stepping
bookkeeping: e.g. leapfrog stores two steps of the prognostic variables, while the Lorenz N-cycle
stores only one but keeps a second tendency for its weighted accumulation.

## [Leapfrog](@id leapfrog)

SpeedyWeather.jl's default time integration is the
[Leapfrog time integration](https://en.wikipedia.org/wiki/Leapfrog_integration),
which, for relative vorticity ``\zeta``, is
in its simplest form
```math
\frac{\zeta_{i+1} - \zeta_{i-1}}{2\Delta t} = RHS(\zeta_i),
```
meaning we step from the previous time step ``i-1``, leapfrogging over the current time step``i``
to the next time step ``i+1`` by evaluating the tendencies on the right-hand side ``RHS``
at the current time step ``i``. The time stepping is done in spectral space
(see [LowerTriangularArrays](@ref lowertriangularmatrices) for how spectral coefficients are stored).
Once the right-hand side ``RHS`` is evaluated, leapfrogging is a linear operation, meaning
that its simply applied to every spectral coefficient ``\zeta_{lm}`` as one would evaluate
it on every grid point in grid-point models.

For the Leapfrog time integration two time steps of the prognostic variables have to be stored,
``i-1`` and ``i``. Time step ``i`` is used to evaluate the tendencies which are then added
to ``i-1`` in a step that also swaps the indices for the next time step ``i \to i-1`` and ``i+1 \to i``,
so that no additional memory than two time steps have to be stored at the same time.

## Leapfrog initialisation

The Leapfrog time integration has to be initialized with an Euler forward step in order
to have a second time step ``i+1`` available when starting from ``i`` to actually leapfrog over.
SpeedyWeather.jl therefore does two initial time steps that are different from
the leapfrog time steps that follow and that have been described above.

1) an Euler forward step with ``\Delta t/2``, then
2) one leapfrog time step with ``\Delta t``, then
3) leapfrog with ``2 \Delta t`` till the end

This is particularly done in a way that after 2. we have ``t=0`` at ``i-1`` and ``t=\Delta t`` at ``i``
available so that 3. can start the leapfrogging without any offset from the intuitive spacing
``0, \Delta  t, 2\Delta t, 3\Delta t, ...``. The following schematic can be useful

|                    | time at step ``i-1`` | time at step ``i`` | time step at ``i+1`` |
| ------------------ | -------------------- | ------------------ | -------------------- |
| Initial conditions | ``t = 0``            |                    |                      |
| 1: Euler           | (T) ``\quad t = 0``  |  ``t=\Delta t/2``  |                      |
| 2: Leapfrog with ``\Delta t``|``t = 0``|(T) ``\quad t = \Delta t/2``| ``t = \Delta t``|
| 3 to ``n``: Leapfrog with ``2\Delta t``|``t-\Delta t``|(T) ``\qquad \quad \quad t``| ``t+\Delta t`` |

The time step that is used to evaluate the tendencies is denoted with (T).
It is always the time step furthest in time that is available.

Before the time integration starts, the initial conditions -- which only occupy step ``i-1`` --
are also copied into step ``i``. This is purely a bookkeeping convenience: it lets
`which_prognostic_step` always return step ``i`` (2) as the one to evaluate the right-hand side on,
without a special case for the very first (Euler) step where step ``i`` would otherwise not yet
exist.

The initial Euler step is not filtered, see [Robert-Asselin and Williams filter](@ref) below.
On the first two steps (Euler, and the first leapfrog step) `Leapfrog`'s `update_prognostic!`
therefore disables the filter weights (both are 0), and enables them for every step after that.

## Robert-Asselin and Williams filter

The standard leapfrog time integration is often combined with a Robert-Asselin filter
[Robert1966, Asselin1972](@citep)
to dampen a computational mode. The idea is to start with a standard leapfrog step to obtain
the next time step ``i+1`` but then to correct the current time step ``i`` by applying a filter
which dampens the computational mode. The filter looks like a discrete Laplacian in time
with a ``(1, -2, 1)`` stencil, and so, maybe unsurprisingly, is efficient to filter out
a "grid-scale oscillation" in time, aka the computational mode. Let ``v`` be the unfiltered
variable and ``u`` be the filtered variable, ``F`` the right-hand side tendency,
then the standard leapfrog step is
```math
v_{i+1} = u_{i-1} + 2\Delta tF(v_i)
```
Meaning we start with a filtered variable ``u`` at the previous time step ``i-1``, evaluate
the tendency ``F(v_i)`` based on the current time step ``i`` to obtain an
unfiltered next time step ``v_{i+1}``. We then filter the current time step ``i``
(which will become ``i-1`` on the next iteration)
```math
u_i = v_i + \frac{\nu}{2}(v_{i+1} - 2v_i + u_{i-1})
```
by adding a discrete Laplacian with coefficient ``\tfrac{\nu}{2}`` to it, evaluated
from the available filtered and unfiltered time steps centred around ``i``:
``v_{i-1}`` is not available anymore because it was overwritten by the filtering
at the previous iteration, ``v_i, v_{i+1}`` are not filtered yet when applying
the Laplacian. The filter parameter ``\nu`` is typically chosen between 0.01-0.2,
with stronger filtering for higher values.

[Williams2009](@citet) then proposed an additional filter step to regain accuracy
that is otherwise lost with a strong Robert-Asselin filter [Amezcua2011, Williams2011](@citep).
Now let ``w`` be unfiltered, ``v`` be once filtered, and ``u`` twice filtered, then
```math
\begin{aligned}
w_{i+1} &= u_{i-1} + 2\Delta tF(v_i) \\
u_i &= v_i + \frac{\nu\alpha}{2}(w_{i+1} - 2v_i + u_{i-1}) \\
v_{i+1} &= w_{i+1} - \frac{\nu(1-\alpha)}{2}(w_{i+1} - 2v_i + u_{i-1})
\end{aligned}
```
with the Williams filter parameter ``\alpha \in [0.5, 1]``. For ``\alpha=1``
we're back with the Robert-Asselin filter (the first two lines).

The Laplacian in the parentheses is often called a *displacement*,
meaning that the filtered value is displaced (or corrected) in the direction
of the two surrounding time steps. The Williams filter now also applies
the same displacement, but in the opposite direction to the next time
step ``i+1`` as a correction step (line 3 above) for a once-filtered
value ``v_{i+1}`` which will then be twice-filtered by the Robert-Asselin
filter on the next iteration. For more details see the referenced publications.

The initial Euler step (see [Leapfrog initialisation](@ref)) is not filtered.
Both the the Robert-Asselin and Williams filter are then switched on for all
following leapfrog time steps.

Implementation-wise, `update_prognostic!` for `Leapfrog` folds both filters into two
weights ``w_1 = \tfrac{\nu\alpha}{2}`` and ``w_2 = \tfrac{\nu(1-\alpha)}{2}`` (both zero
on the unfiltered initial steps) and applies the leapfrog step, the Robert-Asselin and the
Williams filter in a single kernel per spectral coefficient (and vertical layer).

## [Stochastic leapfrog](@id stochastic_leapfrog)

`StochasticLeapfrog` is the leapfrog scheme above with additive Wiener noise, it integrates
the stochastic differential equation
```math
dx = F(x)dt + a dW
```
for every prognostic variable ``x`` with the deterministic tendencies ``F`` of the model and
white noise of a fixed amplitude ``a`` (`AdditiveWienerNoise`: independent in time, between
variables and between layers, spectrally flat by default). Every time step a new pattern
``\xi_i`` of unit variance is drawn and ``a\xi_i / \sqrt{\Delta t}`` is added to the tendencies,
so that the state receives increments ``\propto a\sqrt{\Delta t}\xi_i``, whose variance grows
linearly in time independent of the time step.

The scheme is *drift-implicit*: the noise is treated explicitly and is added to the complete
explicit tendencies right before the implicit terms of the time stepping. Splitting the
tendencies ``F = N + L`` into the terms ``N`` that are treated explicitly and the linear terms
``L`` that are treated implicitly (the [horizontal diffusion](@ref diffusion) with ``\alpha = 1``
and the gravity waves of the [semi-implicit time stepping](@ref implicit_primitive) with the
centering ``\alpha = \tfrac{1}{2}`` by default), a time step is
```math
x_{i+1} = x_{i-1} + 2\Delta t (1 - 2\Delta t \alpha L)^{-1}
\left[ N(x_i) + L x_{i-1} + \frac{a}{\sqrt{\Delta t}} \xi_i \right]
```
followed by the [Robert-Asselin and Williams filter](@ref "Robert-Asselin and Williams filter").
The noise is therefore on the right-hand side of the implicit solves like every explicit
tendency and the implicit terms damp it within the time step, as they do in the exact solution
of the linear problem, in which the noise enters as ``\int e^{L(t-s)} a dW(s)``.
Implementation-wise, `time_step!` calls `random_process!(vars, time_stepping, model)` between the
dynamical core and `diffusion_and_implicit!`, which for `StochasticLeapfrog` adds the noise to
the spectral tendencies and does nothing for the other time steppers.

This is how stiff linear terms are treated in the standard schemes for stochastic
differential equations

- The implicit Euler schemes of [KloedenPlaten1992](@citet) (chapter 12) and the stochastic
  theta method [Higham2000](@citep) evaluate the drift (partly) at the new time step and add
  the Wiener increment ``\Delta W_i`` explicitly,
  ``x_{i+1} = x_i + \Delta t [\theta F(x_{i+1}) + (1-\theta) F(x_i)] + a \Delta W_i``.
- The stochastic semi-implicit leapfrog scheme analysed for geophysical models by
  [EwaldTemam2005](@citet), in its weak version [EwaldWeakSchemes](@citep)
  ``x_{i+1} = x_{i-1} + 2\Delta t N(x_i) + 2\Delta t L x_{i+1} + a(\Delta W_{i-1} + \Delta W_i)``.

Two properties of the scheme here are worth knowing

- One increment ``2a\sqrt{\Delta t}\xi_i`` is added per leapfrog step, not the sum of two
  increments as in [EwaldTemam2005](@citet). The physical mode (the mean of the two interleaved
  leapfrog sequences) receives the same ``a\sqrt{\Delta t}(\xi_{i-1} + \xi_i)``, but the noise
  also excites the computational mode, which the time filter damps. The filter is hence
  part of the stochastic scheme.
- The alternative of adding the noise after the implicit terms (a split-step scheme,
  [HighamMaoStuart2002](@citealp)) converges to the same equation for ``\Delta t \to 0`` and
  injects exactly the specified noise covariance. However, it leaves the noise undamped in the
  stiff modes, which are implicit because they are stiff: the smallest scales
  (``2\Delta t \nu \approx 0.3`` at the truncation wavenumber with the default resolution and diffusion) and the
  fast gravity waves. The two variants agree where ``2\Delta t L`` is small, i.e. on the large scales.

The diffusion the noise injects is the nominal ``a^2\sigma^2/2`` (``\sigma`` is the standard
deviation of the pattern) only for modes for which neither the implicit terms nor the time
filter matter, it should be measured for the scales of interest. Ocean and land variables have
no implicit terms and are stepped with `EulerForward` (see [Ocean, sea ice and land](@ref)),
for them the scheme is the Euler-Maruyama method with grid-point white noise.

````@docs; canonical=false
StochasticLeapfrog
AdditiveWienerNoise
````

## [Time steppers and variable steps](@id steps)

Different time integration schemes need to store a different number of past states of the
prognostic variables and/or tendencies. Leapfrog needs the two spectral steps ``i-1`` and ``i``;
the Lorenz N-cycle needs only one prognostic step but a second tendency to accumulate weighted
tendencies. SpeedyWeather handles this generically: prognostic variables and tendencies carry an
extra [Step dimension](@ref) (the last dimension of their underlying array), and the time stepper
decides how many steps to allocate and which step each model component should read or write.

The number of steps is requested by the time stepper through `prognostic_steps` and
`tendency_steps` (with `prognostic_grid_steps`, `prognostic_spectral_steps`,
`tendency_grid_steps`, `tendency_spectral_steps` to distinguish grid/spectral and dispatch over
the model). For example leapfrog requests two spectral prognostic steps but, in the
primitive-equation models, also two grid steps (because the parameterizations are evaluated at the
previous grid state):

```julia
prognostic_steps(::AbstractLeapfrog) = 2
prognostic_grid_steps(::AbstractLeapfrog, ::PrimitiveEquation) = 2
tendency_steps(::AbstractLeapfrog) = 1
```

The Lorenz N-cycle, by contrast, only needs one prognostic step (the state is updated in place)
but two spectral tendency steps -- one for the explicit tendency ``F``, one for the weighted
accumulation ``G`` -- and only one grid tendency step (``G`` is only ever needed in spectral
space):

```julia
prognostic_steps(::NCycleLorenz) = 1
tendency_steps(::NCycleLorenz) = 2          # also for ocean and land variables, stepped directly on the grid
tendency_grid_steps(::NCycleLorenz) = 1     # atmospheric grid tendencies only hold F
```

Throughout the dynamical core and parameterizations a variable is then accessed with
`get_prognostic_step` / `get_tendency_step`, which return a view of the appropriate step:

```julia
# in a model component, e.g. the spectral→grid transform or a tendency term
vor  = get_prognostic_step(vars.prognostic.vorticity, time_stepping, component)
ζtend = get_tendency_step(vars.tendencies.vorticity,  time_stepping, component)
```

Which step is returned is decided by the time stepper via `which_prognostic_step` /
`which_tendency_step`, dispatched on the variable, the time stepper, the *component* and
(optionally) the model — so a scheme can choose a different step per process. The fallback is step
1, and leapfrog for instance overrides it to read the current (2nd) step for transforms and the
nonlinear dynamical core, but the previous (1st) step for the linear terms and horizontal
diffusion:

```julia
which_prognostic_step(var, ::AbstractLeapfrog, ::AbstractSpectralTransform)       = 2  # current
which_prognostic_step(var, ::AbstractLeapfrog, ::LinearDynamicalCore)             = 1  # previous
which_prognostic_step(var, ::AbstractLeapfrog, ::AbstractHorizontalDiffusion)     = 1  # previous
```

For the Lorenz N-cycle, only the explicit tendency ``F`` is reset to zero every substep for
accumulation into ``G``; the fallback step 1 covers this and every other component that only
ever needs one tendency step:

```julia
which_tendency_step(var, ::AbstractNCycleLorenz, ::ResetTendencies) = 1
```

`get_step(var, i)` is the low-level accessor used by all of the above. Which array dimension is
the step dimension is taken from the variable's dimension tag (see
[Array dimensions](@ref array_dimensions)): `get_step` selects the dimension tagged as time `T`,
so for a variable tagged `XYT`, `XYZT`, `LMT` or `LMZT` it returns a view of step `i`. A variable
*without* a time tag (`XY`, `XYZ`, `LM`, `LMZ`) has no step dimension to select from — the full
variable is returned as a view and `i` is ignored. Model components can therefore call it
unconditionally, whether or not the variable they are handed has steps. `get_steps(var)` returns
all steps as a tuple of views (a 1-tuple for a variable without a step dimension), and
`get_steps(var, Val(N))` does the same with a compile-time length, which is what Enzyme needs.

Dispatching on whether a variable has a step dimension at all is done with the
`AbstractFieldWithTime` / `LowerTriangularArrayWithTime` aliases, as `nsteps` does:

```julia
nsteps(var::Union{AbstractFieldWithTime, LowerTriangularArrayWithTime}) = size(var, ndims(var))
nsteps(::Union{AbstractField, LowerTriangularArray}) = 1
```

````@docs; canonical=false
get_step
SpeedyWeather.get_steps
````

When writing a new time stepper you implement the `*_steps` methods (how many steps to store), the
`which_*_step` methods (which step each component reads/writes) and an `update_prognostic!` method
(how a tendency advances the state); the rest of the model is agnostic to the scheme.

Time steppers can also reorder how diffusion and the semi-implicit correction are applied via
`diffusion_and_implicit!`: Leapfrog applies the implicit correction first, then horizontal
diffusion, while the Lorenz N-cycle applies diffusion first, then the implicit correction
(consistent with the ordering in [Hotta2016](@citet)). Right before, a stochastic time stepper
adds its noise to the complete explicit tendencies with `random_process!`, see
[Stochastic leapfrog](@ref stochastic_leapfrog).

### Ocean, sea ice and land

Dynamic ocean, sea ice and land components (e.g. `SlabOcean`, `ThermodynamicSeaIce`,
`LandBucketTemperature`, `LandBucketMoisture`, `SnowModel`) write tendencies into `vars.tendencies.ocean`
and `vars.tendencies.land`. The variables in these namespaces are time stepped with their own
time stepper within the time stepping setup `model.time_stepping`,
`SpeedyWeather.time_stepper(model.time_stepping, :ocean)` and `(model.time_stepping, :land)`
(sea ice lives in the `:ocean` namespace). By default this is the time stepper of the atmosphere,
but [`Leapfrog`](@ref) has `ocean` and `land` fields which default to
[`EulerForward`](@ref), which steps the variables with **Euler forward over ``\Delta t``**:

```julia
time_stepping = Leapfrog(spectral_grid)                             # Euler forward for ocean and land
time_stepping = Leapfrog(spectral_grid, land = nothing)             # leapfrog land as the atmosphere
model = PrimitiveWetModel(spectral_grid; time_stepping)
```

[`NCycleLorenz`](@ref) has the same `ocean` and `land` fields, defaulting to `nothing`, i.e. the
same N-cycle as the atmosphere. Any time stepper can be used for ocean and land, e.g.
`Leapfrog(spectral_grid, land = NCycleLorenz(spectral_grid))`.

The variables are allocated with the number of steps of their time stepper
(`prognostic_steps` and `tendency_steps`, see `SpeedyWeather.get_nsteps(time_stepping, :land)`), so with
`EulerForward` ocean and land variables have only one step, with `NCycleLorenz` one prognostic
step and two tendency steps. The time step of `ocean` and `land` is always set
to the ``\Delta t`` of the leapfrog (also with `set!(model, Δt=...)`), and on the first two
leapfrog steps, which each only advance the clock by ``\Delta t/2``, also Euler steps with ``\Delta t/2``.

Why not leapfrog? The tendencies are computed from surface fluxes (sensible and latent heat,
longwave radiation) that the parameterizations evaluate at the previous (1st) step, and they
contain stiff relaxation terms, e.g. the sensible heat flux feedback on the thin top soil layer has
timescales of about an hour in strong winds. Leapfrogging those over ``2\Delta t`` from
``t - \Delta t`` is an Euler step of ``2\Delta t``, which halves the stability limit compared to an
Euler step of ``\Delta t`` and made long integrations unstable. Terms that are meant to act
"within one time step" (e.g. restoring sea surface temperatures to freezing when sea ice forms)
use `SpeedyWeather.time_step(model.time_stepping, :ocean, clock)` for the same reason.

````@docs; canonical=false
EulerForward
````

## [Lorenz N-cycle](@id ncycle)

The Lorenz N-cycle [`NCycleLorenz`](@ref) is a semi-implicit time integration following
[Hotta2016](@citet). Over a cycle of ``N`` substeps it advances the state ``x``,
per substep, as

```math
\begin{aligned}
G &= w\,F_E(x) + (1-w)\,G \\
dx &= (I - \alpha\,\Delta t\,L_I)^{-1} (G + L_I x) \\
x &= x + \Delta t\,dx
\end{aligned}
```

where ``F_E`` are the explicit tendencies, ``L_I`` the linear (gravity-wave) terms treated
implicitly, and ``w`` a substep-dependent weight coefficient. Only one prognostic step is stored
(the state is updated in place), but two tendencies are kept: ``F`` (the current tendency) and
``G`` (the weighted accumulation that carries memory across substeps), hence
`tendency_spectral_steps(::NCycleLorenz) = 2`. In SpeedyWeather's implementation the implicit
solve ``(I - \alpha\,\Delta t\,L_I)^{-1}(G + L_I x)`` is folded into the same `implicit_correction!`
step used by Leapfrog (see [Semi-implicit time stepping](@ref implicit_primitive) for the
primitive equations and [Semi-implicit time integration](@ref implicit_swm) for the shallow water
model -- both implicit solves are shared between time steppers), so `update_prognostic!` for
`NCycleLorenz` itself only performs the weighted accumulation and the explicit state update

```math
\begin{aligned}
G &= w\,F + (1-w)\,G \\
x &= x + \Delta t\,G
\end{aligned}
```

for every spectral coefficient (and vertical layer), with the implicit correction applied to
``G`` beforehand by `diffusion_and_implicit!`.

Four weight variants are available, following the naming in [Hotta2016](@citet), selected via the
`variant` option:

| Variant | Weight ``w_k`` for substep ``k = 1, \dots, N-1`` (``w_0 = 1``) | Steps per period | Notes |
| ------- | --------------------------------------------------------------- | ----------------- | ----- |
| `NCycleLorenzA` (default) | ``w_k = N / (N - k)`` | ``N`` | |
| `NCycleLorenzB` | ``w_k = N / k`` | ``N`` | |
| `NCycleLorenzAB` | alternates a full cycle of `A`, then a full cycle of `B` | ``2N`` | 2 cycles per period |
| `NCycleLorenzABBA` | sequence `A`, `B`, `B`, `A` | ``4N`` | 4th-order accurate for ``N=4``; other ``N`` run but without that accuracy guarantee |

"Steps per period" is the number of substeps after which the weight sequence repeats. These
substeps can be interpreted in two ways: as sub-stages of a single step, the way a Runge-Kutta
method would use them, where the full time step ``\Delta t`` is only completed once all substeps
of a period have run; or as a multi-step method with changing weights that restarts every period,
where every substep already advances the state by a full ``\Delta t`` (as in the update
``x = x + \Delta t\,G`` above). SpeedyWeather.jl uses the latter, so the computational cost of the
Lorenz N-cycle is independent of ``N`` but, like leapfrog, inversely proportional to ``\Delta t``.

The cycle length ``N`` is set with `steps` (3 or 4 recommended, 4 is more stable). The current
substep ``k`` is `current_substep(L, clock) = mod(clock.step_counter, L.steps)`, and the weight for
that substep is computed by `weight_coefficient`, dispatched on the `variant`.

````@docs; canonical=false
NCycleLorenz
NCycleLorenzA
NCycleLorenzB
NCycleLorenzAB
NCycleLorenzABBA
````
## [Periodic schedules](@id periodic_schedules)

Many components don't need to do something on every time step but only periodically,
e.g. the semi-implicit solver of the primitive equations, which recomputes its reference
temperature profile (and with it the implicit operators) once a day by default.
SpeedyWeather.jl has three schedules for this, all subtypes of `AbstractSchedule`
and all used with `isscheduled(schedule, clock)`, which returns `true` on the time
steps on which to execute

| Schedule | Mutable | Periodic | Events | Initialized with | Primarily used for |
| -------- | ------- | -------- | ------ | ---------------- | ------------------ |
| [`Schedule`](@ref) | yes | `every::Second` | `times::Vector{DateTime}` | `initialize!(schedule, clock)` at every `run!` | callbacks, see [Schedules](@ref) |
| [`PeriodicSchedule`](@ref) | yes | `every::Second` | – | `initialize!(schedule, time_stepping)` in `initialize!(model)` | model components |
| [`StepSchedule`](@ref) | no | `every_n_steps::Int` | – | nothing to initialize | model components |

The `Schedule` precomputes a `BitVector` of all time steps of a simulation, so that it can
combine periodic execution with events at given times. But for that it needs the start
time and the length of the simulation, which are only known in `run!`.
`PeriodicSchedule` and `StepSchedule` are purely periodic and only need the time step, so
they can be part of model components that are initialized with `initialize!(model)`
before `run!` is called. A model component would typically use them as

```julia
@kwdef struct MyComponent <: SpeedyWeather.AbstractModelComponent
    schedule::PeriodicSchedule = PeriodicSchedule(every = Day(1))
    # other fields
end

function SpeedyWeather.initialize!(component::MyComponent, model::AbstractModel)
    initialize!(component.schedule, model.time_stepping)    # converts every to every_n_steps
    # initialize other things here
end

function SpeedyWeather.timestep!(vars::Variables, component::MyComponent, model::AbstractModel)
    isscheduled(component.schedule, vars.prognostic.clock) || return nothing
    # do something every day only
end
```

### PeriodicSchedule vs StepSchedule

A `PeriodicSchedule` is defined by a period of time, `every`, which is converted into a
number of time steps `every_n_steps` with the model's time step on `initialize!`.
`every` is rounded to a multiple of the time step, and is at least one time step;
`every = Second(0)` means never. It also counts its executions in `counter`.

```@example periodic_schedules
using SpeedyWeather
spectral_grid = SpectralGrid(truncation = 32, nlayers = 1)
time_stepping = Leapfrog(spectral_grid)
periodic_schedule = PeriodicSchedule(every = Hour(6))
initialize!(periodic_schedule, time_stepping)
```

so with a time step of `time_stepping.Δt_millisec` at T32 executing every 6 hours means every
`periodic_schedule.every_n_steps` time steps. Changing the resolution therefore changes
`every_n_steps` but not (up to rounding) `every`. A `StepSchedule` on the other hand is
immutable and directly defined by the number of time steps `every_n_steps`, independent of the
time step; `every_n_steps = 0` means never

```@example periodic_schedules
step_schedule = StepSchedule(every_n_steps = 4)
```

which is useful if something should execute every n-th time step regardless of resolution.
`StepSchedule(1)` executes on every time step (but the first), e.g. for the semi-implicit solver
to always recompute its reference temperature profile

```julia
reference_temperature = GlobalMeanReferenceTemperature()
implicit = ImplicitPrimitiveEquation(spectral_grid; schedule = StepSchedule(1), reference_temperature)
model = PrimitiveWetModel(spectral_grid; implicit)
```

The semi-implicit corrections are centred in time by default, `centering = 0.5` (Crank-Nicolson).
`centering = 1` is backward Euler, which damps the gravity waves. With (near-)uniformly spaced
σ levels, i.e. `SigmaCoordinates(spectral_grid)`, Crank-Nicolson can be unstable over steep
orography (a warning is raised), use the default `FriersonSigmaCoordinates` or `centering = 1` there.
The schedule is only used with a reference temperature profile that is recomputed from the
global-mean temperature. The semi-implicit solver linearises about a reference temperature
profile ``T_k``, chosen with `reference_temperature`

| Reference temperature | Profile ``T_k`` | Recomputed |
|---|---|---|
| `IsothermalReferenceTemperature(temperature = 300)` (default) | constant `temperature` | only when the time step changes |
| `GlobalMeanReferenceTemperature(offset = 0)` | global-mean temperature + `offset` | on the `schedule` |
| `MaximumReferenceTemperature(margin = 30)` | constant, maximum temperature + `margin` | on the `schedule` |

A warm isothermal reference profile is the most stable choice (Simmons, Hoskins and Burridge, 1978),
as also used in other models, e.g. 300 K in Isca and the IFS. The warmer the reference temperature,
the more the gravity waves treated implicitly are slowed down. `MaximumReferenceTemperature`
therefore follows the model's climate, e.g. for a snowball Earth or a hot climate, with an
isothermal reference temperature just above the warmest temperature in the model.

### Scheduling and the leapfrog spin-up

All schedules count the time steps with `clock.time_step_counter` which excludes the
initial conditions and the spin-up of the [Leapfrog initialisation](@ref): the first
Euler forward step with ``\Delta t/2`` does not count. With every 4th time step, i.e.
`Schedule(every = 4Δt)`, `PeriodicSchedule(every = 4Δt)` or `StepSchedule(4)`, all three
schedules execute on the same steps, after 4``\Delta t``, 8``\Delta t``, ... of simulated time.
For the first 10 leapfrog time steps (clock shown after each step)

| Time stepper step | Leapfrog step | `clock.step_counter` | `clock.time_step_counter` | Clock time | `isscheduled` |
|---|---|---|---|---|---|
| (initial) | - | 0 | 0 | 0 | |
| 1 | Euler, ``\Delta t/2`` (spin-up) | 1 | 0 | ``\tfrac{1}{2}\Delta t`` | |
| 2 | Leapfrog, ``\Delta t`` (spin-up end) | 2 | 1 | ``\Delta t`` | |
| 3 | Leapfrog, ``2\Delta t`` | 3 | 2 | ``2\Delta t`` | |
| 4 | Leapfrog, ``2\Delta t`` | 4 | 3 | ``3\Delta t`` | |
| 5 | Leapfrog, ``2\Delta t`` | 5 | 4 | ``4\Delta t`` | ✅ |
| 6 | Leapfrog, ``2\Delta t`` | 6 | 5 | ``5\Delta t`` | |
| 7 | Leapfrog, ``2\Delta t`` | 7 | 6 | ``6\Delta t`` | |
| 8 | Leapfrog, ``2\Delta t`` | 8 | 7 | ``7\Delta t`` | |
| 9 | Leapfrog, ``2\Delta t`` | 9 | 8 | ``8\Delta t`` | ✅ |
| 10 | Leapfrog, ``2\Delta t`` | 10 | 9 | ``9\Delta t`` | |

Callbacks are evaluated after the clock is stepped forward, so they see the clock as in the table.
Model components inside the time step, like the semi-implicit solver, are evaluated before the
clock is stepped forward. They see the same clock one time stepper step later, e.g. the semi-implicit
solver recomputes its reference temperature profile at the beginning of time stepper step 6,
computing the tendencies at time ``4\Delta t``. With events, `Schedule(times = ...)`, a `Schedule`
executes on the first time step on or after the given time, see [Scheduling logic](@ref).

````@docs; canonical=false
PeriodicSchedule
StepSchedule
````
