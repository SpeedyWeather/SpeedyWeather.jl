# Semi-Lagrangian advection for the barotropic model

> Status: **in progress**. Draft implementation of a two-time-level semi-Lagrangian time stepper
> for `BarotropicModel`, plus an exponential (ETD1) treatment of hyperdiffusion that it depends on.
> Not yet tested — opened as a draft PR to review how the changes sit in the code.

Date of initial draft: 2026-09-20

Base revision: 35d764b6a9029b32399bb536906ce17732de118d

## Originating prompt

> I have written a 2D semi Lagrangian advection scheme for passive tracers before. Can you explain
> to me what changes when you have an additional forcing term, what would the algorithm look like
> given some time stepping method for the forcing term

> Could you have a look how feasible it would be to implement a semi Lagrangian advection scheme
> just for the 2D barotropic model to start with. Don't write any code yet, just check the dynamical
> core and modularity of the code structure and wrt to the time steppers whether this could work
> without too many code changes. RingGrids has 2D interpolation already, that should be used.

## Revision log

- **2026-09-20** — Initial feasibility survey. Concluded the `AbstractTimeStepper` interface is the
  right seam and that `ParticleAdvection2D` already contains a reusable spherical trajectory solver.
- **2026-09-20** — Asked whether the "SL as an effective tendency" variant could be written as an
  analytical exponential, as proposed for relaxation terms. Answer: not for advection (the advection
  operator has purely imaginary eigenvalues, so its exponential is a phase shift with no timescale;
  and the SL shift already *is* that exponential). But hyperdiffusion is diagonal in spectral space
  with a genuine per-degree timescale, and the existing backward-Euler treatment can be made exact
  by replacing one line. Added as a prerequisite to this plan.
- **2026-09-20** — Asked whether SL is compatible with `NCycleLorenz`. Answer: no, not for the
  transport term (see "Why not NCycleLorenz" below). Decision: implement a dedicated two-time-level
  `SemiLagrangian` stepper.
- **2026-09-21** — Review of PR #1266. Changes made in response:
  - Departure point coordinates and winds there are `Grid2D` `DynamicsVariable`s rather than bare
    `VectorDim` `ScratchVariable`s — there is exactly one departure point per grid point, so entry
    `ij` is the departure point of the trajectory *arriving* at cell `ij`.
  - Arrival coordinates are no longer stored on the stepper, they are `model.geometry.londs/latds`.
  - `SemiLagrangian` is `@kwdef` with the option defaults on the type.
  - `ζ+f` is formed in place in `vars.grid.vorticity` instead of in a dedicated array (nothing reads
    it again before the closing `transform!` overwrites it), removing one grid allocation.
  - `get_step(..., 1)` for the spectral vorticity replaced by `get_prognostic_step(..., DynamicalCore())`.
  - `STEP_COMPONENT` in the `which_prognostic_step` signatures spelled out as `AbstractModelComponent`
    and `AbstractSpectralTransform`.
  - The first-step special case in `extrapolate_winds!` is gone. The barotropic `transform!` now calls
    the existing `move_prognostic_grid_variables_back!` hook (no-op by default, and an explicit no-op
    for `Leapfrog` on 2D models which keeps a single grid step), so `uⁿ⁻¹ = uⁿ` going into the first
    step. Previously step 1 held *zeros*, not a copy of step 2, so the guard was load-bearing; this
    is the cleaner fix.
  - Verified `displace` does not allocate: `Particle` is `isbits`, `move`/`mod` are `@inline`, and
    the round trip infers `Tuple{NF, NF}` with 0 allocations over a 10⁴-point loop.

  Not changed, with reasoning in the PR thread: making `exponential` the *default* for
  `HyperDiffusion` (agreed it is generally better, but it changes results for every existing
  simulation so it wants its own PR), and turning `f` into a `Field` instead of a latitude vector
  (agreed, but it touches `Coriolis` everywhere — see Future work).
- **2026-09-21** (second review round):
  - **Trajectories moved to 3D Cartesian great circles.** The lat/lon form needs
    `dlon = u*Δt/(radius*cos(lat))`, singular at the poles. On T128 the outermost ring is at
    89.28˚N where `cos(lat) = 0.0125`, so a 50 m/s wind over a 45 min step gave a **97˚ longitude
    jump in one step** — the departure point near the poles was meaningless. The Cartesian form has
    no coordinate singularity and is exact for solid-body rotation.
  - **Diffusion fully decoupled from the time stepper.** `SemiLagrangian` no longer has an
    `exponential_diffusion` option and no longer mutates `model.horizontal_diffusion`. ETD1 is a
    property of the diffusion component, not of the scheme: the Lie–Trotter splitting error between
    transport and diffusion is `O(Δt²)`, the same order as the difference between `φ₁(z)` and
    `1/(1-z)`, so making the damping factor exact improves a constant, not the order. The earlier
    framing ("ETD1 needs `exponential=true`") overstated it.
  - `@kwdef` reverted to a plain mutable struct plus generator, matching the other time steppers.

  Measured at T128, 30 days, default `KolmogorovFlow` forcing (statistical steady state, so these
  reflect the forcing/dissipation balance rather than pure decay, but the forcing is identical):

  | configuration | enstrophy | KE |
  |---|---|---|
  | Eulerian `NCycleLorenz` 30min | 4.055e-5 | 2.851e7 |
  | SL 60min, default diffusion | 3.582e-5 | 2.354e7 |
  | SL 60min, `time_scale=Hour(100)` | 3.742e-5 | 2.494e7 |
  | SL 60min, `n_iterations=4` | 3.854e-5 | 2.493e7 |
  | SL 180min, default diffusion | 2.764e-5 | 6.506e7 |

  Two things to read off: weakening the diffusion 25× recovers only ~4% of enstrophy, confirming
  the damping is the interpolation and not the diffusion; and `n_iterations=4` recovers ~7%, so
  trajectory error contributes more than expected. The 180min KE is 2.3× the Eulerian value, which
  looks like the SETTLS extrapolation going unstable at long steps — to be checked with
  `extrapolate_winds=false`.
- **2026-09-21** (third round): added `CubicInterpolator` to RingGrids and made the interpolation
  machinery interpolator-flexible. Measured in isolation on an octahedral Gaussian grid, shifting a
  smooth field by a third of a grid cell 100 times (what semi-Lagrangian does every step):
  bilinear-class `AnvilInterpolator` retains **88.7%** of the amplitude, `CubicInterpolator`
  **99.87%**. Off-grid RMS error on a smooth field is 114x lower. `SemiLagrangian` now takes an
  `Interpolator` keyword and defaults to cubic.

  Flexibility fixes needed to get there: `find_grid_indices!` claimed `::AbstractLocator` but
  destructured Anvil-specific fields, so it now dispatches on the concrete locator; the
  `Interpolator(grid, npoints)` constructor is generic over `AbstractInterpolator` via `Locator(I)`
  and `nonparametric_type`; and the `AbstractLocator` contract (`npoints_output`, `js`, `Δys` filled
  generically by `find_rings!`, everything else the locator's own business) is now documented.

- **2026-09-21** (fourth round): cubic interpolation worked but cost ~2.3x Eulerian. Profiling at
  T128 showed the cost was **not** the interpolation but the locator update: `find_grid_indices!`
  was 6.79 ms against `interpolate!` at 1.03 ms and `find_rings!` at 0.56 ms. Three fixes:

  1. `find_cubic_indices_kernel!` rewritten — the four latitude weights are computed in one pass
     instead of four calls to a per-node helper doing nested `ifelse` selection (40 `ifelse` and 16
     divisions per point, down to ~12 and 4), the ring loop is unrolled, and the longitude wrap uses
     `ifelse` rather than `mod` (16 integer divisions per point removed). 6.79 → 4.08 ms.
  2. The first trajectory iteration was redundant: the first guess for the departure point *is* the
     arrival point, where the wind is already known, so the locator update and two interpolations
     were recomputing `u_trajectory` (then named `u_star`). Seed the departure winds directly instead. The winds at the *final*
     departure points are never read either, so that refresh is skipped too. `n_iterations = 2` now
     costs one locator update and two wind interpolations rather than two of each.
  3. Split the interpolators: `trajectory_interpolator` (default `AnvilInterpolator`) places the
     departure point, `interpolator` (default `CubicInterpolator`) shifts the transported field. Only
     the latter's damping accumulates — the winds are rebuilt from the spectral state every step.

  Measured at T128 (same machine, relative to Eulerian `NCycleLorenz` at 30min):

  | configuration | SYPD | vs Eulerian |
  |---|---|---|
  | Eulerian `NCycleLorenz` 30min | 213.3 | 1.00x |
  | SL 120min, cubic field / anvil trajectory | 232.4 | **1.09x** |
  | SL 120min, cubic everywhere | 188.0 | 0.88x |
  | SL 120min, anvil everywhere | 277.1 | 1.30x |

  Accuracy is unaffected by the cheap trajectory interpolator: at T128/30d, enstrophy 5.288e-5
  (cubic field / anvil trajectory) against 5.073e-5 (cubic everywhere) at 60min, and 5.093e-5
  against 5.159e-5 at 120min — within the scatter of a forced statistical steady state, and not
  systematically worse.


- **2026-09-21 (sixth review round).** milankl: *"can you double check that the great circle
  calculation is correct? Because I still see velocities blow up towards the poles, negative
  -500m/s near the north pole, 500m/s near the south pole, quickly growing there from 0 over 30
  days"*, plus comments on the cubic locator storage, the interpolator constructor constraint,
  where the interpolator should live, and field docstring length.

  The great-circle calculation itself is correct — `displace` agrees with the independent
  direct-geodesic formula to 1e-12 at every latitude except exactly ±90˚, where the local
  east/north frame is undefined and the two conventions differ in longitude only. But the review
  uncovered **three genuine bugs**, two of which were cancelling each other:

  1. **Trajectory winds were averaged across mismatched frames.** `u, v` are components in the
     *local* east/north frame. `_departure_point_kernel!` averaged the arrival and departure winds
     componentwise, which is only valid where the two frames coincide. The mismatch is the
     meridian convergence, O(Δλ sinφ) — negligible in midlatitudes, O(1) next to the poles, where
     a trajectory spans a large Δλ. This made the fixed-point iteration *diverge*: at T31/30d the
     polar ring-mean `u` was +218/-222 m/s at `n_iterations = 2` and -8800/+8766 at
     `n_iterations = 4`. Both winds are now mapped to a common 3D Cartesian frame
     (`cartesian_velocity`) before averaging, and `displace` projects the mean back onto the
     tangent plane. This is the reported symptom.

  2. **The trajectory was half the correct length.** The trapezoidal rule is
     `x_d = x_a - Δt (V(x_a) + V(x_d))/2`. The kernel formed the mean `(V_a + V_d)/2` and then
     *also* halved the time step, applying `Δt/4 (V_a + V_d)`. Measured against the exact
     trajectory for a uniform wind, the departure points sat at a ratio of 0.5000 at every
     latitude; the docstring already stated the correct formula.

  3. **SETTLS was a silent no-op.** `move_prognostic_grid_variables_back!` was called at the *end*
     of the barotropic `transform!`, after the current grid step had already been overwritten, so
     step 1 was copied from the *new* state and `u_old == u_new`. `3/2 uⁿ - 1/2 uⁿ⁻¹` therefore
     reduced to `uⁿ` exactly. It is now called before the transforms, and only afterwards on
     `initialize`, where there is no previous step yet.

  Bug 3 is why the scheme was unstable, and bug 2 was partly masking it by shortening every
  trajectory. The mechanism: the trajectory carries the Rossby term through `f(x_d) - f(x_a)`, so
  **the time level of the trajectory wind is the time discretisation of the Rossby wave**. With
  the wind at `tⁿ` it is forward Euler, unconditionally unstable on the imaginary axis. Measured
  amplification per step for a small-amplitude m=2, n=3 mode (free decay, no forcing, no drag,
  negligible diffusion):

  | Δt | measured, wind at tⁿ | `1 + (ωΔt)²/2` | measured, SETTLS |
  |---|---|---|---|
  | 60 min | 1.0038185 | 1.003827 | 1.0000468 |
  | 30 min | 1.0009562 | 1.000957 | 1.0000049 |
  | 15 min | 1.0002394 | 1.000239 | 1.0000007 |

  The middle column uses ω = 2Ωm/(n(n+1)) = Ω/3 for this mode; forward Euler is confirmed to
  within 3e-6. With SETTLS working, growth drops by 82x at 60 min and the scaling changes from
  Δt² to ≈Δt⁴, the AB2-like behaviour the extrapolation is supposed to give. `extrapolate_winds`
  is therefore not an optional refinement and is documented as such; `extrapolate_winds = false`
  is a diagnostic that blows up within days.

  End-to-end at T31, RMS wind over 60 days (`KolmogorovFlow` forcing + `LinearVorticityDrag`,
  and an unforced free-decay control):

  | | day 10 | day 20 | day 30 | day 40 | day 50 | day 60 |
  |---|---|---|---|---|---|---|
  | Leapfrog, forced | 17.4 | 21.9 | 24.8 | 27.4 | 29.3 | 30.6 |
  | SL, forced (before) | 677 | 685 | 642 | 587 | 569 | 547 |
  | SL, forced (after) | 14.3 | 18.5 | 21.4 | 22.3 | 22.6 | 24.7 |
  | Leapfrog, free decay | 13.1 | 12.9 | 12.6 | 12.7 | 12.6 | 12.5 |
  | SL, free decay (before) | 747 | 719 | 724 | 709 | 701 | 696 |
  | SL, free decay (after) | 11.8 | 11.5 | 11.5 | 12.3 | 13.9 | 16.7 |

  Polar ring-mean `u` at T31/30d is now -15.2/+12.3 m/s against Leapfrog's +2.1/-1.7, and global
  max|u| 44.9 against Leapfrog's 54.1 — the semi-Lagrangian run is now quieter than the Eulerian
  one rather than 8x louder. This also retracts the earlier hypothesis in this log that the
  120-180min KE excess was SETTLS instability: SETTLS was never running, and the excess was the
  forward-Euler Rossby term.

  The four other review points:

  - `CubicLocator.ijs`/`weights` are now `NSTENCIL_CUBIC x npoints_output` matrices rather than
    flat vectors with a stride, so the kernel indexes `[s, k]` and each point's stencil is one
    contiguous column. `RingGrids.extrema_in` was widened from `AbstractVector` to `AbstractArray`
    for the bounds check.
  - `CubicInterpolator(geometry, locator)` now requires a `CubicLocator`; `AnvilInterpolator` was
    equally unconstrained and now requires an `AnvilLocator`.
  - The locators moved from the time stepper into `Variables`, following the particle-advection
    split: `geometry` is a property of the grid and constant, so it stays in the stepper, while
    the locators are rewritten every step and are now `DynamicsVariable`s under the
    `:semi_lagrangian` namespace. `LocatorDim` was generalised from hardcoding `AnvilLocator` and
    `model.particle_advection.nparticles` to a type parameter and its own `n`.
  - Field docstrings in `SemiLagrangian` reduced to one line each.

## Problem description

Transport in SpeedyWeather is Eulerian and spectral. For `BarotropicModel`,
`vorticity_flux_grid_tendencies!` forms `u_tend = Fᵤ + v(ζ+f)`, `v_tend = Fᵥ - u(ζ+f)` on the grid
and `curl!` turns these into `∂ζ/∂t`. This is accurate and conservative, but the advective CFL
condition ties the time step to the grid spacing.

A semi-Lagrangian scheme removes that restriction: along a trajectory the barotropic equation is

```
D(ζ+f)/Dt = F - cζ + 𝓓ζ
```

so absolute vorticity is materially conserved up to sources. Transport becomes "interpolate the old
field at the departure point" rather than "evaluate a flux divergence", and the stability limit
becomes a trajectory-accuracy limit rather than a CFL limit.

The goal here is a first, deliberately narrow implementation: **2D, barotropic only**, reusing
`RingGrids`' existing interpolation, and leaving shallow water / primitive equations for later.

## Background

### Why the time stepper is the right seam

`time_step!(vars, ::AbstractTimeStepper, ::Union{Barotropic, ShallowWater})` already dispatches on
the stepper type. `NCycleLorenz` demonstrates that a fundamentally different scheme can be added
purely additively: it declares its own step counts, overrides `which_tendency_step`,
`diffusion_and_implicit!` and `update_prognostic!`, and nothing else in the core changed.

Critically, spectral and grid step counts are **already decoupled and dispatched on the model**
(`leapfrog.jl`):

```julia
prognostic_spectral_steps(::AbstractLeapfrog) = 2
prognostic_grid_steps(::AbstractLeapfrog, ::Union{<:Barotropic, <:ShallowWater}) = 1
prognostic_grid_steps(::AbstractLeapfrog, ::PrimitiveEquation) = 2
```

A two-time-level SL scheme wants the mirror image — one spectral state, two grid time levels of
`u, v` for the wind extrapolation — which needs no new mechanism:

```julia
prognostic_spectral_steps(::AbstractSemiLagrangian) = 1
prognostic_grid_steps(::AbstractSemiLagrangian, ::Barotropic) = 2
```

`variables(::Type{<:Barotropic}, nsteps)` already threads `ps` into spectral vorticity and `pg` into
the grid `u`/`v`/`vorticity` declarations, so the allocation follows automatically.

### Why not NCycleLorenz

One cycle of `NCycleLorenz` with N=3, variant A (weights 1, 3/2, 3) has amplification polynomial

```
P(z) = 1 + 3z + (9/2)z² + (9/2)z³
```

matching `exp(3z)` through third order. Feeding it an SL effective tendency `Δt·F = μζ` with
`μ = exp(-iθ) - 1` (θ = **u**·**k**Δt the per-step phase shift) gives `P(μ)`, whereas three exact SL
steps give `(1+μ)³ = exp(-3iθ)`, of modulus exactly 1 for any θ. The two agree only to first order
in μ:

```
P(μ) - (1+μ)³ = (3/2)μ² + (7/2)μ³
```

Since `|μ| = 2|sin(θ/2)|` grows to 2, this is not a small correction. At the grid scale
(θ = πC, C the Courant number):

| C | θ | \|P(μ)\| | exact |
|---|---|---|---|
| 0.25 | π/4 | 0.53 | 1 |
| 0.5 | π/2 | 7.6 | 1 |
| 1.0 | π | 23 | 1 |

Unstable between C = 0.25 and 0.5, and heavily damped where stable — a *worse* stability limit than
the Eulerian core it would replace. The general reason: any tendency-based multistep/multistage
method approximates `exp(z)` by a polynomial, and polynomials have bounded stability regions. SL's
value is an amplification factor on the unit circle for arbitrarily large θ; wrapping it in a
polynomial approximant re-imposes a bounded region. The group-property rescue that works for
leapfrog (`𝒜(Δt)² = 𝒜(2Δt)`, i.e. re-trace the trajectory over 2Δt) has no analogue here, because
`prognostic_steps(::NCycleLorenz) = 1` leaves no earlier state to re-trace from.

`NCycleLorenz` and SL are competing answers to the same question — how to take larger stable steps —
not complementary ones.

### Exponential hyperdiffusion (prerequisite)

`HyperDiffusion` is diagonal in spectral space with a genuine inverse timescale per degree
(`horizontal_diffusion.jl`):

```julia
∇²ⁿ[l+1, k] = -eigenvalue_norm^power / time_scale      # = -1/τ_l
∇²ⁿ_implicit[l+1, k] = 1 / (1 - Δt * ∇²ⁿ[l+1, k])      # = 1/(1 + Δt/τ_l), backward Euler
```

`1/(1-z)` is the [0/1] Padé approximant of `exp(z)`. The exact answer is available at the same cost.
Writing `z = Δt∇²ⁿ`, the exact ETD1 step is

```
ζⁿ⁺¹ = exp(z)ζⁿ + Δt φ₁(z) T,    φ₁(z) = (exp(z) - 1)/z
```

and substituting into the existing kernel form `tendency = (tendency + expl*var) * impl` followed by
`var += Δt*tendency` shows that **`expl` is unchanged** and only

```
impl = φ₁(z)    instead of    1/(1 - z)
```

Check with `tendency = 0`: `ζⁿ⁺¹ = ζⁿ(1 + zφ₁(z)) = ζⁿexp(z)`, exact. Same kernel, same arrays, same
shapes — only the precomputed values differ.

This is a modest accuracy gain for existing steppers (at T32 with the 4-hour default `time_scale`
and Δt = 30 min, `Δt/τ ≈ 0.125` at the truncation limit: 0.8889 vs 0.8825, ~0.7% per step at the
smallest resolved scale only) and buys **accuracy, not stability** — both forms are unconditionally
stable and monotone. It matters more in the SL context, where larger Δt means larger `z`, and it is
what makes the SL update exact for the diffusion operator.

## Summary of changes

### 1. `φ₁` hyperdiffusion (`SpeedyWeather/src/dynamics/horizontal_diffusion.jl`)

New option `HyperDiffusion.exponential::Bool` (default `false`, so existing behaviour is bit-identical
unless opted in). When `true`, the `impl` arrays are filled with `φ₁(z)` instead of `1/(1-z)`. A
numerically safe `φ₁` is used (series expansion for `|z| < eps^(1/3)` to avoid `0/0`).

`SemiLagrangian` sets it to `true` by default via `initialize!`.

### 2. Advection term made dispatchable (`SpeedyWeather/src/dynamics/tendencies.jl`)

`_vorticity_flux_kernel!` gains an `advection` multiplier:

```julia
u_tend_grid[ij, k] = (u_tend_grid[ij, k] + advection * v[ij, k] * ω) * coslat⁻¹j
```

With `advection = 1` (the Eulerian default, dispatched via `advection_factor(::AbstractTimeStepper)`)
behaviour is unchanged; with `advection = 0` (`::AbstractSemiLagrangian`) the flux term is switched
off and the remaining pipeline — forcing, drag, `coslat⁻¹` scaling, transform, `curl!` — produces
exactly the source term `S = ∇×(Fᵤ, Fᵥ) - cζ` that the SL update needs. This means
`dynamics_tendencies!` is reused unchanged.

### 3. Departure points (`SpeedyWeather/src/dynamics/semi_lagrangian.jl`, new)

`SemiLagrangianTrajectory` holds the `GridGeometry`, an `AnvilLocator` sized to `npoints(grid)`, and
the derived trajectory time step. Departure points are found by fixed-point iteration of the
trapezoidal ("iterated backward") trajectory,

```
x_d ← x_a - Δt/2 (u*(x_a) + u*(x_d))
```

with `u*` the time-extrapolated wind `3/2 uⁿ - 1/2 uⁿ⁻¹` (SETTLS-style; falls back to `uⁿ` on the
first step and when `extrapolate_winds = false`). Two iterations by default.

The spherical displacement helper mirrors `advect_2D` in `particle_advection.jl` but works on plain
`(lon, lat)` vectors rather than `Particle`s, and converts explicitly via the planetary radius:

```
dlon = u Δt / (R cos φ) · 180/π,   dlat = v Δt / R · 180/π
```

### 4. Transport (`semi_lagrangian.jl`)

Absolute vorticity `ζ + f·scale` is formed on the grid, interpolated to the departure points with
`RingGrids.interpolate!` through the existing locator, and `f·scale` subtracted at the arrival point.
`f` is scaled on the fly exactly as the Eulerian flux kernel does, so this is consistent with the
radius-scaled prognostic state.

### 5. The stepper (`SpeedyWeather/src/time_stepping/steppers/semi_lagrangian.jl`, new)

```julia
function time_step!(vars, TS::AbstractSemiLagrangian, model::Barotropic)
    reset_tendencies!(vars, TS)
    dynamics_tendencies!(vars, model)       # sources only — advection switched off by dispatch
    departure_points!(vars, TS, model)      # backward trajectories
    semi_lagrangian_transport!(vars, TS, model)   # shift ζ+f, write ζ* into the spectral state
    horizontal_diffusion!(vars, model)      # tendency = (S + ∇²ⁿζ*) φ₁
    update_prognostic!(vars, model)         # ζⁿ⁺¹ = ζ* + Δt·tendency  ⇒  exp(z)ζ* + Δt φ₁(z) S
    transform!(vars, model)
    particle_advection!(vars, model)
end
```

The ordering matters: the transport writes `ζ*` into `vars.prognostic.vorticity` *before*
`horizontal_diffusion!` reads it, so the `(tendency + expl*var)*impl` form picks up the
already-transported state and the final Euler update is exactly ETD1 applied to the shifted field.

`update_prognostic!` for this stepper is a plain forward Euler on the spectral state — correct
because transport is already baked into `ζ*` and diffusion into `φ₁`.

### 6. Grid time-level bookkeeping

`move_prognostic_grid_variables_back!` already exists as a generic hook with a no-op default
(`transform.jl`). The SL stepper defines its own method (copying the `:uv_grid` and `:grid` fuse
parents from step 2 to step 1) and calls it from its own `time_step!`. The barotropic `transform!`
is **not** modified, so `Leapfrog`'s single-grid-step behaviour on 2D models is untouched.

## Testing and verification

Not run yet — this is a draft. Planned:

- Solid-body rotation of a localised vorticity blob: after one full revolution the field should
  return to its initial position. Compare phase error against the Eulerian core.
- A stability sweep over Courant number, confirming SL remains stable well past C = 1 where the
  Eulerian core does not.
- `φ₁` unit test: with zero tendency, one diffusion step must reproduce `exp(-Δt/τ_l)` per degree to
  machine precision; and `φ₁(z) → 1` as `z → 0` without cancellation.
- Cross-check against a three-time-level SL variant built on the existing `Leapfrog`
  (trace over 2Δt from spectral step 1), which is exactly consistent by the group property. If the
  two-time-level and three-time-level results disagree, the trajectory or interpolation code is wrong.
- Conservation diagnostics (mean vorticity, enstrophy) to quantify the interpolation damping.

## Documentation changes

None yet. If the scheme graduates from draft, `docs/src/barotropic.md` and the time stepping docs
need a section, and `HyperDiffusion`'s new `exponential` option needs documenting.

## Known limitations

- **Rossby waves are only weakly damped, not neutral.** With SETTLS the trajectory wind makes the
  Rossby term AB2-like, which is still weakly unstable on the imaginary axis: amplification
  1 + 4.7e-5 per step at Δt = 60 min against forward Euler's 1 + 3.8e-3. Hyperdiffusion absorbs it
  at the resolutions tested, but free decay does creep up over 60 days (11.8 -> 16.7 m/s RMS at
  T31 against Leapfrog's 13.1 -> 12.5). A genuinely centred trajectory wind would need a
  predictor-corrector pass — recompute `u` from the provisional `ζⁿ⁺¹`, average with `uⁿ`, redo
  the trajectory — at roughly double the cost per step. Not attempted here.
- **`extrapolate_winds = false` is unstable**, not a cheaper first-order option. It is kept as a
  diagnostic and documented as such.
- **Not conservative.** SL is not conservative, and the grid↔spectral round trip adds to the drift.
- **Tracers remain Eulerian.** `tracer_advection!` is untouched, so tracers still use the spectral
  flux form and keep their CFL limit. Inconsistent with the vorticity transport; deliberate for a
  first pass.
- **Cost.** `update_locator!` over `npoints(grid)` points every iteration of every step, plus one
  extra grid→spectral transform, are real costs the Eulerian core does not pay. Measured at T128:
  1.09x Eulerian throughput at a 4x longer time step, see the fifth revision log entry.
- **Stencil memory.** `CubicLocator` stores 16 indices and 16 weights per output point, ~7.7 MB at
  T128. Leaner layouts are possible (the four ring offsets and one Δ per ring would do) at the
  cost of recomputing the Lagrange weights inside the interpolation kernel.
- **Barotropic only.** `ShallowWater` and `PrimitiveEquation` fall back to their existing steppers;
  constructing them with `SemiLagrangian` is not supported.
- **No GPU verification.** Kernels are written with `KernelAbstractions` and should be device-
  agnostic, but this has not been exercised.

## Future work

- Higher-order (cubic / quasi-monotone) interpolation on reduced `RingGrids`.
- Extend to `ShallowWater` (divergent trajectories), then to 3D with vertical trajectories.
- Semi-Lagrangian tracer transport sharing the same departure points — nearly free once the
  trajectories exist, and removes the tracer/vorticity inconsistency.
- A conservative (SLICE/CSLAM-style) cell-integrated variant if conservation turns out to matter.
- Make `Coriolis.f` a `Field` rather than a latitude vector, so `ζ+f` is a plain broadcast
  `vor .+ scale*f` instead of needing a `whichring` lookup in a kernel. Duplicates `f` across
  longitudes but is unlikely to matter for performance, and would simplify several call sites
  beyond this one.
- Flip `HyperDiffusion.exponential` to `true` by default once it has been validated against the
  existing test suite — it is strictly more accurate at identical cost, but changes results.
