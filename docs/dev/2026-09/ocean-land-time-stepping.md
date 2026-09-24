# Separate time steppers for ocean and land

> Status: **completed** (pending CI). `Leapfrog` gets `ocean` and `land` fields holding the time steppers
> for these namespaces, defaulting to a new `EulerForward`. Ocean and land variables are
> allocated with the step counts of their own time stepper.

Date of initial draft: 2026-09-24

Base revision: `c14d88f3` (`mc/landstepping-euler`, PR #1264, targeting #1183)

## Originating prompt

> Can you look at PR1264? The target PR1183 has the aim to replace the hardcoded euler forward of
> land prognostic variables with an interface that would use a centralised time stepper. I
> understand that this is unstable with Leapfrog and so obviously it should be possible to choose
> a different time stepper for land and ocean. Before implementing anything can you present me a
> plan that would allow a different time stepping for atmosphere/ocean/land? Maybe there could be
> a field inside every time stepper .ocean .land to have an Euler fallback or chosen differently,
> or we could have a CompoundTimeStepper which distinguishes between atmosphere/ocean/land at a
> higher level? Try to not make this change too intrusive

## Revision log

- **Fields on the time stepper, not a `CompoundTimeStepper`.** `model.time_stepping` is used
  throughout as the atmospheric time stepper (`which_prognostic_step` methods, step counts, clock,
  `set!`, `Δt`, Adapt). A wrapper would have to forward all of that.
- **Distinguish ocean and land**, no "surface" umbrella, as that does not correspond to a
  namespace. Sea ice lives in the `:ocean` namespace and follows the ocean time stepper.
- **Allocate variables with the namespace's time stepper.**
- **New `EulerForward` type** rather than reusing `NCycleLorenz(steps = 1)`, which is Euler
  forward mathematically but needs 2 tendency steps and is confusing to read.
- **`namespace_time_stepping` instead of `time_stepping`** as function name, to avoid clashes with
  local variables of that name.
- **Namespace allocation with the child's own step counts.** `NCycleLorenz` for land failed with
  a `BoundsError` as land tendencies were allocated with the grid tendency steps (1, only F),
  but NCycleLorenz needs F and G. Ocean and land variables are stepped directly, so they are
  now allocated with `get_namespace_nsteps(model, namespace)`, i.e. `prognostic_steps` and
  `tendency_steps` of their time stepper, without the grid/spectral distinction of the
  atmosphere (Leapfrog: `prognostic_steps = 2`, NCycleLorenz: `tendency_steps = 2`).
- **`NCycleLorenz` also gets `ocean` and `land` fields** (default `nothing`, i.e. itself) so that
  the interface is not only the fallback. The Δt sync moved into the generic `calculate_Δt!`
  and `set!` so that it works for every time stepper with these fields.
- With Leapfrog for the atmosphere, an NCycleLorenz for ocean/land starts its first cycle at
  step counter 0 (w = 1, Euler), the cycle phase is then shifted by leapfrog's start-up step,
  which is fine.

## Problem description

#1183 moved soil temperature and moisture into the general `update_prognostic!`, so with
`Leapfrog` they were stepped over 2Δt from t−Δt, which is unstable for the stiff surface flux
relaxation terms. #1264 fixed this with an `update_prognostic_surface!` hook that hardcodes Euler
forward for Leapfrog, but keeps 2 (identical) prognostic steps allocated and offers no way to
choose the time stepping of ocean or land.

## Background

- Allocation is per component: each ocean, sea ice and land `variables(component, model)` calls
  `get_nsteps(model.time_stepping, model)`.
- `update_prognostic_namespaces!` is a generated loop over the `:ocean` and `:land` namespaces,
  the namespace is a compile-time literal there.
- Leapfrog starts with two steps that advance the clock by Δt/2 each.

## Summary of changes

- New `EulerForward <: AbstractTimeStepper` (`time_stepping/steppers/euler_forward.jl`) with the
  usual `Δt_at_T32`, `adjust_with_output`, `Δt_millisec`, `Δt` fields so that `set!` and
  `calculate_Δt!` work generically. 1 prognostic step, 1 tendency step,
  `update_prognostic!` does `var += Δt/scale * tendency`.
- `namespace_time_stepping(time_stepping, Val(namespace))` returns the time stepper of a namespace,
  defaulting to the time stepper itself; `namespace_time_stepping(model, namespace)` as
  convenience. (Not called `time_stepping` to avoid clashes with local variables of that name.)
- `Leapfrog` gets `ocean` and `land` fields (default `EulerForward(spectral_grid)`, `nothing`
  means leapfrog like the atmosphere). `initialize!(::Leapfrog, model)` and `set!(::Leapfrog, …)`
  set the Δt of the children to the Leapfrog's Δt.
- `time_step_scale(parent, child, clock)`: Leapfrog passes `scale = 2` to non-Leapfrog children
  on its first two steps (clock advances by Δt/2). `time_step(model, namespace, clock)` gives the
  resulting time step [s], used by `ThermodynamicSeaIce`.
- `update_prognostic_namespaces!` calls `update_prognostic!` with the namespace's time stepper;
  `update_prognostic_surface!` and `surface_time_step` from #1264 are removed.
- Ocean, sea ice and land `variables` use `get_nsteps(namespace_time_stepping(model, namespace), model)`,
  and every `get_prognostic_step`/`get_tendency_step` on ocean/land variables (also in surface
  fluxes and longwave radiation) passes the namespace's time stepper.
- Leapfrog's `copy_step_forward!` is a no-op for variables with a single step.

## Testing and verification

- `test/parameterizations/surface_time_stepping.jl` rewritten: `EulerForward` step and `scale`,
  Δt synced from Leapfrog via `initialize!` and `set!`, `time_step(model, namespace, clock)` on
  the start-up steps, 1-step allocation with the default and 2 steps with `land = nothing`, a short
  run with `land = nothing, ocean = nothing`.
- Only small tests locally; full test suite and long integrations on CI.

## Documentation changes

- `docs/src/time_integration.md`: the "Ocean, sea ice and land" subsection describes the
  `ocean`/`land` fields and `EulerForward`.

## Known limitations

- Ocean and sea ice share the `:ocean` namespace and cannot be stepped differently.
- Child time steppers are not nested (a child's own `ocean`/`land` fields are ignored).
- `Leapfrog` as ocean/land time stepper under a non-leapfrog atmosphere is not supported: its
  steps are only initialised (copied 1 → 2) when the atmosphere uses Leapfrog.
- No sub- or super-cycling: the children's Δt is always the Leapfrog's Δt.
- `land = nothing` / `ocean = nothing` with Leapfrog keeps the instability #1264 fixed; it is an
  opt-in for experiments.

## Future work

- `EulerForward` as a full atmospheric time stepper (clock `time_step!`, `which_prognostic_step`).
- Sub-/super-cycling of ocean and land.
