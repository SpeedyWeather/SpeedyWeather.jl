# [Time integration](@id leapfrog)

SpeedyWeather.jl supports several time integration schemes, selected by passing a
`time_stepping` component to the model constructor:

- [`Leapfrog`](@ref leapfrog), a 2-step leapfrog scheme with a Robert-Asselin and Williams filter
  (the default for `ShallowWaterModel`, `PrimitiveDryModel` and `PrimitiveWetModel`), described below.
- [`NCycleLorenz`](@ref ncycle), a family of semi-implicit Lorenz N-cycle schemes
  (Hotta et al. 2016[^Hotta2016]; the default for `BarotropicModel`).

All schemes share a common framework (see [Time steppers and variable steps](@ref steps)) in which
the time stepper decides, for every model component, which stored *step* of each variable to read
or write. This decouples the dynamical core and parameterizations from the time-stepping
bookkeeping: e.g. leapfrog stores two steps of the prognostic variables, while the Lorenz N-cycle
stores only one but keeps a second tendency for its weighted accumulation.

## Leapfrog

SpeedyWeather.jl's default time integration is the
[Leapfrog time integration](https://en.wikipedia.org/wiki/Leapfrog_integration),
which, for relative vorticity ``\zeta``, is
in its simplest form
```math
\frac{\zeta_{i+1} - \zeta_{i-1}}{2\Delta t} = RHS(\zeta_i),
```
meaning we step from the previous time step ``i-1``, leapfrogging over the current time step``i``
to the next time step ``i+1`` by evaluating the tendencies on the right-hand side ``RHS``
at the current time step ``i``. The time stepping is done in spectral space.
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


## Robert-Asselin and Williams filter

The standard leapfrog time integration is often combined with a Robert-Asselin filter[^Robert66][^Asselin72]
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

Williams[^Williams2009] then proposed an additional filter step to regain accuracy
that is otherwise lost with a strong Robert-Asselin filter[^Amezcua2011][^Williams2011].
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

The initial Euler step (see [Time integration](@ref leapfrog), Table) is not filtered.
Both the the Robert-Asselin and Williams filter are then switched on for all
following leapfrog time steps.

## Leapfrog options

The leapfrog time integration is controlled by creating a custom `Leapfrog`
component and passing it on to the model constructor

```@example leapfrog
using SpeedyWeather
spectral_grid = SpectralGrid()
time_stepping = Leapfrog(spectral_grid)
```

and with `?Leapfrog` you see a summary of the fields, only manually change those marked `[OPTION]`.
We will discuss some options in the following.

````@docs; canonical=false
Leapfrog
````

## Change the time step

SpeedyWeather chooses the time step automatically based on the resolution.
A default time step of

```@example leapfrog
time_stepping.Δt_at_T31
```

is used at T31 (`trunc=31`) spectral resolution (see [Available horizontal resolutions](@ref))
which is then (almost) linearly scaled to higher (or lower) resolution. Creating a simulation
at twice the resolution (T63) will approximately half the time step (20min if T31 runs at 40min).
This is such that in most cases the user does need to know what time step is stable. But if
you want a shorter time step the easiest is to choose `Δt_at_T31` (write `\Delta` then hit tab,
works in the Julia REPL and other interfaces) relative to its default. If you half that time step
you'll half the time step for all resolutions. The other "time steps" in `time_stepping` are
explained in the docstring (`?Leapfrog`). Note that internally SpeedyWeather uses a time step
scaled by the radius, so `time_stepping.Δt` will be in units of second divided by meter.

You can also choose the time step manually with

```@example leapfrog
set!(time_stepping, Δt=Minute(10))
```

which you can do before `initialize!(model)` or after -- it will change the other time step information
consistently, as shown here. You can provide any `Second`, `Minute`, `Hour`, but note that there is a
stability limit above which your simulation quickly blows up.

## Adjust with output

By default the time step is (slightly) adjusted to match the [Output frequency](@ref).
See that section for more information.

## Restart with Leapfrog

As a 2-step scheme, leapfrog time stepping has to be initialised with an Euler forward step to have
information for the 2nd step, see [Leapfrog initialisation](@ref). This Euler spin-up step
(`spin_up_steps(::AbstractLeapfrog) == 1`) is taken at the start of every integration and does not
count towards the clock or the output frequency.

SpeedyWeather also allows the user to issue several `run!(simulation)` calls one after another to
continue a simulation, possibly after some modification by the user. Each `run!` re-initialises the
clock (its step counter is reset to 0), so every `run!` -- including continued ones -- begins again
with the Euler spin-up step; the integration is restarted from the currently available state rather
than relying on a stored 2nd step.

For time steppers that need no such initialisation `spin_up_steps` is `0`, e.g. the
[Lorenz N-cycle](@ref ncycle).

## Passing `time_stepping` to the model constructor

Here we just called our time stepping scheme `time_stepping` but this needs to be passed on to the model constructor,
e.g. for the `PrimitiveDryModel`

```@example leapfrog
model = PrimitiveDryModel(spectral_grid; time_stepping)
nothing # hide
```

where `;` matches the `time_stepping` keyword argument by name. If you name `leapfrog = Leapfrog(spectral_grid)` then you
would need to change this to `time_stepping=leapfrog` in the function call arguments. The same
keyword takes any time stepper, e.g. `time_stepping = NCycleLorenz(spectral_grid)` for the
[Lorenz N-cycle](@ref ncycle).

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
prognostic_spectral_steps(::AbstractLeapfrog) = 2
prognostic_grid_steps(::AbstractLeapfrog, ::PrimitiveEquation) = 2
tendency_steps(::AbstractLeapfrog) = 1
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

`get_step(var, i)` is the low-level accessor used by all of the above; for a variable with a step
dimension it returns a view of step `i`.

````@docs; canonical=false
get_step
````

When writing a new time stepper you implement the `*_steps` methods (how many steps to store), the
`which_*_step` methods (which step each component reads/writes) and an `update_prognostic!` method
(how a tendency advances the state); the rest of the model is agnostic to the scheme.

## [Lorenz N-cycle](@id ncycle)

The Lorenz N-cycle [`NCycleLorenz`](@ref) is a semi-implicit time integration following
Hotta et al. (2016)[^Hotta2016]. Over a cycle of ``N`` substeps it advances the state ``x``,
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
`tendency_spectral_steps(::NCycleLorenz) = 2`.

Four weight variants are available, following the naming in Hotta et al. (2016):
`NCycleLorenzA` (default), `NCycleLorenzB`, `NCycleLorenzAB` (alternating A and B), and
`NCycleLorenzABBA` (an A-B-B-A sequence that is 4th-order accurate for ``N=4``). The cycle length
``N`` is set with `steps` (3 or 4 recommended, 4 is more stable).

```@example ncycle
using SpeedyWeather
spectral_grid = SpectralGrid()
time_stepping = NCycleLorenz(spectral_grid, steps=4, variant=NCycleLorenzAB())
model = PrimitiveWetModel(spectral_grid; time_stepping)
nothing # hide
```

````@docs; canonical=false
NCycleLorenz
````

## References

[^Robert66]: Robert, André. "The Integration of a Low Order Spectral Form of the Primitive Meteorological Equations." Journal of the Meteorological Society of Japan 44 (1966): 237-245.
[^Asselin72]: ASSELIN, R., 1972: Frequency Filter for Time Integrations. Mon. Wea. Rev., 100, 487-490, doi:[10.1175/1520-0493(1972)100<0487:FFFTI>2.3.CO;2](https://doi.org/10.1175/1520-0493(1972)100<0487:FFFTI>2.3.CO;2.)
[^Williams2009]: Williams, P. D., 2009: A Proposed Modification to the Robert-Asselin Time Filter. Mon. Wea. Rev., 137, 2538-2546, [10.1175/2009MWR2724.1](https://doi.org/10.1175/2009MWR2724.1).
[^Amezcua2011]: Amezcua, J., E. Kalnay, and P. D. Williams, 2011: The Effects of the RAW Filter on the Climatology and Forecast Skill of the SPEEDY Model. Mon. Wea. Rev., 139, 608-619, doi:[10.1175/2010MWR3530.1](https://doi.org/10.1175/2010MWR3530.1).
[^Williams2011]: Williams, P. D., 2011: The RAW Filter: An Improvement to the Robert-Asselin Filter in Semi-Implicit Integrations. Mon. Wea. Rev., 139, 1996-2007, doi:[10.1175/2010MWR3601.1](https://doi.org/10.1175/2010MWR3601.1).
[^Hotta2016]: Hotta, D., E. Kalnay, and P. Ullrich, 2016: A Semi-Implicit Modification to the Lorenz N-Cycle Scheme and Its Application for Integration of Meteorological Equations. Mon. Wea. Rev., 144, 2215-2233, doi:[10.1175/MWR-D-15-0330.1](https://doi.org/10.1175/MWR-D-15-0330.1).