# Time integration

SpeedyWeather.jl supports several time integration schemes, selected by passing a
`time_stepping` component to the model constructor:

- [`Leapfrog`](@ref leapfrog), a 2-step leapfrog scheme with a Robert-Asselin and Williams filter
  (the default for `ShallowWaterModel`, `PrimitiveDryModel` and `PrimitiveWetModel`).
- [`NCycleLorenz`](@ref ncycle), a family of semi-implicit Lorenz N-cycle schemes
  (the default for `BarotropicModel`).

This page covers how to create and configure a time stepper and pass it to the model
constructor. For the mathematics behind each scheme and how the time stepper decides which
stored step of a variable to read or write, see [Time stepping](@ref time_stepping) in the
Numerics section.

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
time_stepping.Δt_at_T32
```

is used at T31 (`truncation=32`) spectral resolution (see [Available horizontal resolutions](@ref))
which is then (almost) linearly scaled to higher (or lower) resolution. Creating a simulation
at twice the resolution (T63) will approximately half the time step (20min if T31 runs at 40min).
This is such that in most cases the user does need to know what time step is stable. But if
you want a shorter time step the easiest is to choose `Δt_at_T32` (write `\Delta` then hit tab,
works in the Julia REPL and other interfaces) relative to its default. If you half that time step
you'll half the time step for all resolutions. The other "time steps" in `time_stepping` are
explained in the docstring (`?Leapfrog`).

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

## Lorenz N-cycle options

Like `Leapfrog`, the Lorenz N-cycle is controlled by creating a custom `NCycleLorenz`
component and passing it on to the model constructor. The cycle length `steps` (3 or 4
recommended) and the weight `variant` ([`NCycleLorenzA`](@ref) (default), [`NCycleLorenzB`](@ref),
[`NCycleLorenzAB`](@ref) or [`NCycleLorenzABBA`](@ref), see [Lorenz N-cycle](@ref ncycle) for what
distinguishes them) are the two options you are most likely to change, the time step options
(`Δt_at_T32`, `adjust_with_output`) behave the same as for `Leapfrog` above.

```@example ncycle
using SpeedyWeather
spectral_grid = SpectralGrid()
time_stepping = NCycleLorenz(spectral_grid, steps=4, variant=NCycleLorenzAB())
```

````@docs; canonical=false
NCycleLorenz
````

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
model = PrimitiveWetModel(spectral_grid; time_stepping)
nothing # hide
```
