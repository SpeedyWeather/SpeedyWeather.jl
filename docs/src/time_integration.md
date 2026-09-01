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
keyword takes any time stepper, e.g. the `NCycleLorenz` created above

```@example ncycle
model = PrimitiveWetModel(spectral_grid; time_stepping)
nothing # hide
```
