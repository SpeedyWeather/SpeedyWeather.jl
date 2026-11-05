# Compute the geopotential once, in grid space

> Status: **in progress**. Grid geopotential wired into the dynamical core, stability being checked.

Date of initial draft: 2026-11-05

Base revision: e82a6571

## Originating prompt

> We currently calculate the geopotential twice, once in spectral space for the dynamics and
> once in grid space for the parameterizations as we need to know the height for some of them.
> Ideally we only compute it once, in grid space, and also use that formulation in the
> dynamical core. Can you investigate whether this is possible? There's some issues around the
> time step to use, as the dynamics may use a different step as the parameterizations but in
> that case, maybe use the dynamics step as that's important for the stability and the implicit
> solver, whereas for the parameterizations it might be fine to use another time step, or retain
> the previous step so both dynamics/parameterizations can choose what do use.

> Yes please go ahead, and if you're able to get this stable (I would start in a setup without
> parameterizations) do open a PR, that can be a draft if there's still open questions but so I
> can better see what you're working on

## Revision log

- 2026-11-05: initial draft.
- 2026-11-05: comparison against main showed bitwise-close results for 2 steps, then
  ~1 m/s differences from step 3. Cause: during the Leapfrog start (`lf == 1`, Euler step and
  first leapfrog step) the spectral previous step stays at x₀ while the grid previous step was
  moved forward to x½. Fixed: `move_prognostic_grid_variables_back!` now leaves the previous
  grid step unchanged when the leapfrog update did not advance the previous spectral step.
  With that, branch and main agree to round-off over the first steps, with or without the
  Robert/Williams filter, i.e. the filtered-vs-unfiltered previous step is a negligible effect.

## Problem description

Primitive equation models compute the geopotential Φ twice per time step:

1. `geopotential!(vars, model)` in `transform!` (grid space) for the parameterizations,
   with the full virtual temperature `Tᵥ = T(1+μq)`, supporting σ and σ-pressure coordinates.
2. `geopotential!(vars, G, orography)` in `dynamics_tendencies!` (spectral space) for
   `∇²Φ` in the divergence tendency, with a linearised `Tᵥ = T + Tₖμq` and σ coordinates only.

## Background

- Both read step 1 under Leapfrog (`AbstractGeopotential` and `LinearDynamicalCore`).
  But spectral step 1 holds the Robert/Williams-filtered state, whereas grid step 1 is the
  copy of the previous grid state (`move_prognostic_grid_variables_back!`), i.e. unfiltered.
  The two Φ therefore differ by the filter increment ~ν/2·(xₙ₋₁ − 2xₙ + xₙ₊₁).
- The implicit solver's geopotential operator `R` acts on temperature tendencies only and is
  built from the `Geopotential` integration constants; it does not read the spectral Φ.
- ShallowWater already adds the grid Φ to the grid Bernoulli potential before the batched
  grid→spectral transform; the same can be done for primitive equation models at no extra
  transform.

## Summary of changes

- `bernoulli_grid_potential!` for primitive equation models writes `½(u²+v²) + Φ` with the
  grid geopotential, transformed in the existing batched transform.
- The linear `R·Tₖ·lnpₛ` pressure-gradient term stays spectral (it is the term the implicit
  solver treats at the `LinearDynamicalCore` step) and is added to the spectral Bernoulli
  potential after the transform in `linear_pressure_gradient!`.
- Remove the spectral geopotential (`geopotential!(vars, G, orography)`, its kernel,
  `spectral_geopotential`), `linear_virtual_temperature!` and the spectral `virtual_temperature`
  variable.
- The geopotential step remains chosen by `which_prognostic_step(var, TS, ::AbstractGeopotential)`
  (= 1 for Leapfrog, same as `LinearDynamicalCore`), so the dynamics and parameterizations share
  one Φ. A time stepper wanting a different step for parameterizations could later give the grid
  geopotential a step dimension.

## Testing and verification

- Compare Φ (grid vs. transformed old spectral) and the divergence tendency after one step
  against main.
- Long runs without parameterizations (`PrimitiveDryModel` with physics off, Held-Suarez,
  with orography), at default and larger Δt, check for stability and closeness to main.
- PrimitiveWet default run.
- Update `test/dynamics/geopotential.jl`.

## Documentation changes

None beyond docstrings.

## Known limitations

- Dynamics now sees the unfiltered previous state in Φ (see Background).
- Φ is nonlinear in grid space (`T(1+μq)`), its transform is subject to the usual aliasing.

## Future work

- With Φ in grid space the dynamical core no longer assumes `p ∝ σ` for the geopotential,
  a step towards hybrid coordinates (the implicit operator still does).
