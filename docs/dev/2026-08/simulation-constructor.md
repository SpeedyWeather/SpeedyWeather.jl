# `Simulation(model)` as an alternative to `initialize!(model)`

> Status: **completed**. `Simulation(model; kwargs...)` is added as an alias for
> `initialize!(model; kwargs...)`, documented in the "Model initialization" section of
> `docs/src/how_to_run_speedy.md`, and covered by a unit test in
> `SpeedyWeather/test/dynamics/simulation_constructor.jl`.

Date of initial draft: 2026-08-18

Base revision: `b65590e7` (`mk/simulation`)

## Originating prompt

> This branch allows the user to write `simulation = Simulation(model)` as an alternative to
> `simulation = initialize!(model)`. The latter highlights that model is actually mutated in
> this process, the former highlights that a Simulation object is returned. While we widely
> use the initialize! formulation and do not want to change that we want to mirror the user
> interface of Oceananigans.jl too that uses Simulation as interface. Can you write this into
> the documentation were it best fits, maybe in the how to run section?

## Revision log

- **2026-08-18, initial draft.** The constructor itself already existed on the branch
  (uncommitted change in `SpeedyWeather/src/models/simulation.jl`); this plan covers the
  documentation asked for, plus the unit test that the testing convention in `CLAUDE.md`
  requires for every new method. The user handles the `CHANGELOG.md` entry themselves.

## Problem description

SpeedyWeather.jl initializes a model into a runnable simulation with

```julia
simulation = initialize!(model)
```

The `!` correctly signals that `model` is mutated during this step: components precompute
arrays, load data from file, and communicate information between each other. What the name
does *not* convey is that the return value is a `Simulation` object, distinct from the
`model` that was passed in.

Oceananigans.jl, the other major Julia geophysical fluid dynamics model, spells the same step
as a constructor call, `Simulation(...)`. Users moving between the two packages, or coming to
SpeedyWeather from Oceananigans, benefit from the SpeedyWeather interface accepting the same
spelling.

## Background

`Simulation{V, M}` is defined in `SpeedyWeather/src/models/simulation.jl` and is already
exported. It is a plain two-field struct holding `variables` and `model`, so the default
constructor `Simulation(variables, model)` exists but is not the user-facing entry point —
`initialize!(::AbstractModel)` is, and it does the actual work of initializing components and
allocating variables.

Adding a single-argument method `Simulation(model::AbstractModel; kwargs...)` therefore does
not conflict with the default constructor (different arity) and can simply forward to
`initialize!`.

## Summary of changes

1. **`SpeedyWeather/src/models/simulation.jl`** — a one-line method

   ```julia
   Simulation(model::AbstractModel; kwargs...) = initialize!(model; kwargs...)
   ```

   with a docstring stating it is the same as `initialize!(model)`. Because it forwards
   `kwargs...`, every keyword `initialize!` accepts (notably `time`) works unchanged.

2. **`docs/src/how_to_run_speedy.md`** — a `!!! note` admonition in the
   [Model initialization](@id initialize) section, placed directly after the paragraph
   introducing `initialize!(model, time=...)`, i.e. after both forms of `initialize!` a user
   would meet have been shown. The note states that the two are identical including keyword
   arguments, explains the difference in emphasis (mutation vs. returned object), names
   Oceananigans.jl as the motivation, and says explicitly that the documentation continues to
   use `initialize!`. The example inside the admonition is an `@example howto` block, so the
   new method is exercised at documentation build time.

3. **`SpeedyWeather/test/dynamics/simulation_constructor.jl`** — new test file, picked up
   automatically by `find_tests`. It asserts that `Simulation(model)` returns a `Simulation`
   wrapping the same model object, that keyword arguments are forwarded (checking the clock
   time for `time=DateTime(2020,5,1)`), and that `initialize!` and `Simulation` produce
   equivalent initial state for a `ShallowWaterModel`.

No version bump beyond what the branch already carries: `SpeedyWeather/Project.toml` is
already at `0.22.0+DEV`, and this is an additive, non-breaking change, so the accumulating
`+DEV` tag already covers it.

## Testing and verification

```bash
julia --project=SpeedyWeather --check-bounds=yes -e \
  'using Test, SpeedyWeather, Dates; include("SpeedyWeather/test/dynamics/simulation_constructor.jl")'
```

passes (5/5). The documentation change is additionally verified by the docs build, since the
admonition contains a live `@example howto` block.

## Documentation changes

`docs/src/how_to_run_speedy.md`, section "Model initialization", as described above. No other
documentation page is changed: the deliberate decision is that `initialize!` remains the form
used throughout the documentation, examples, and docstrings, with `Simulation(model)`
mentioned once where a user first learns about initialization.

## Known limitations

- The alias is documented in exactly one place. A user who only ever reads, say,
  `docs/src/examples_2D.md` will not learn that `Simulation(model)` exists. This is
  intentional — showing both forms everywhere would suggest a distinction that does not exist
  — but it does mean discoverability rests on the how-to-run page and the docstring.
- The similarity to Oceananigans.jl is in spelling only. Oceananigans' `Simulation`
  constructor takes a model plus run-control keywords such as `Δt` and `stop_time`, whereas
  SpeedyWeather's takes only `initialize!` keywords and controls run length through
  `run!(simulation, period=...)`. Users porting code across cannot expect keywords to carry
  over.

## Future work

- If the `Simulation(model)` spelling turns out to be what users reach for, consider whether
  run-control keywords (`period`, `output`) should be accepted at construction time to close
  the remaining gap with Oceananigans' interface. That would be a genuine interface change
  rather than an alias, and would need its own plan.
