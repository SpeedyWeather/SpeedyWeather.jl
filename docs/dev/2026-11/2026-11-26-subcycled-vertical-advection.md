# SubcycledVerticalAdvection: Courant-adaptive vertical advection for many layers

> Status: **in progress**. Revised after review and benchmarked; GPU efficiency of the adaptive mode is the main open point.

Date of initial draft: 2026-11-26

Base revision: fc7a46906a0b6dc6730b89b8f672a37f831d2c3c (first version), 7626f058bf27fb0a43472830fb8fcb955780de68 (revision)

## Originating prompt

> [From the developers:] You open two PRs separately for the cloud fix and the new vertical advection scheme.

> OK, got you! Please go ahead with addressing #1436 when suitable.

## Revision log

- First version (#1436): work arrays allocated in the component, WENO5 hard-coded as the stage
  operator, the integration interval stored in a `Ref` field set by the per-variable entry point,
  adaptive substeps only.
- Review: work arrays moved to the scratch variables `vars.scratch.grid.a`, `.b`, `.c` declared by the
  scheme; the stage operator is the option `scheme` (any `DiffusiveVerticalAdvection`, default WENO);
  `courant_target = nothing` gives a fixed number of `substeps`; the `Ref` was replaced by
  `advection_interval(scheme, vars, model)`, computed once in `vertical_advection!(vars, model)` and
  passed to `_vertical_advection!` together with `vars.scratch` (also on the tracer path).
  `advection_interval` returns `nothing` for all other schemes, so they are unchanged.
- `max_substeps` renamed `substeps`, as it is also the fixed number of substeps.
- Tests now go through the public constructor and `_vertical_advection!`.

## Problem description

With 100 layers (Frierson spacing) the lowest layer is ~19 m deep and Δσ ≈ 0.001–0.01 near the
surface. With Δt = 300 s the vertical Courant number |σ̇|Δt/Δσ reaches 1 at only ~0.3 m/s of ascent,
which ordinary deep convection exceeds. All existing schemes evaluate the advection once per time
step and blew up within 1–2 days at T85 L100 with a 2Δz checkerboard in ascending columns: centred,
upwind-3 and WENO5, also with the time step halved to 150 s.

## Background

- Diffusive schemes (upwind, WENO) advect ξ at t−Δt over the leapfrog interval 2Δt
  (`which_prognostic_step(..., ::DiffusiveVerticalAdvection) = 1`), so their effective Courant number
  is twice that of the time step Δt.
- SSP-RK3 (Shu and Osher, 1988) is linearly stable up to C ≈ 1.43 with the linear WENO5 (fifth-order
  upwind) operator and up to √3 ≈ 1.73 with centred operators (purely imaginary eigenvalues). Both
  values were checked numerically from the amplification factor `1 + z + z²/2 + z³/6`.
- An Oceananigans-style split (explicit up to a Courant limit, implicit first-order upwind beyond it)
  is designed as a guard against rare Courant spikes. At L100 violations near the surface are routine,
  so the implicit first-order part would be active in many columns every step. Subcycling keeps the
  advection explicit and high-order everywhere.

## Summary of changes

- `SpeedyWeather/src/dynamics/subcycled_vertical_advection.jl` (new):
  `SubcycledVerticalAdvection(spectral_grid; scheme = WENOVerticalAdvection(spectral_grid), courant_target = 1, substeps = 16)`,
  a `DiffusiveVerticalAdvection{NF, B}` with `B` from the inner scheme. One GPU thread per column:
  the number of substeps `n` comes from the column's maximum Courant number over the interval
  (clamped to 1 ≤ n ≤ `substeps`, non-finite values clamped before `ceil` so NaNs reach the NaN check
  instead of throwing in the kernel), or is fixed. Each substep is SSP-RK3 of the inner scheme's flux-form
  operator. The tendency `(ξ_end − ξ₀)/interval` is added to the tendencies.
- `SpeedyWeather/src/dynamics/vertical_advection.jl`: `advection_interval` (default `nothing`) and two
  trailing arguments `interval, scratch` for `_vertical_advection!`, ignored by the existing schemes.

## Testing and verification

- `test/dynamics/vertical_advection.jl`: the scheme is in the "runs" and "reference temperature
  profile" testsets. New testset: constructors and options, scratch variable declared, at a tiny
  interval it reproduces the WENO tendency, at Courant number 9.5 the result stays finite and bounded
  within 1 % and the bump moves by the right distance (WENO and upwind-3 stage operators), fixed
  10 substeps equal adaptive ceil(9.5) = 10 substeps, interval taken from the time stepper.
- Used in all our L100 runs from T31 to T512 (A100), including 90 days at T512 L100, all stable.
- Cost: `SpeedyWeather/benchmark/subcycled_vertical_advection.jl`. `PrimitiveWetModel` with default
  physics, 12-hour spin-up with the subcycled scheme, then each scheme timed on the same state
  (ms per `vertical_advection!` call for u, v, T, q). "step" is a full time step with the adaptive scheme.

  CPU (Apple M2, 4 threads):

  | T | L | WENO | adaptive | n = 1 | n = 4 | columns n = 1 | max n | max C | step | adaptive / step |
  |---|---|---|---|---|---|---|---|---|---|---|
  | 31 | 8 | 0.58 | 2.28 | 2.29 | 8.6 | 100 % | 1 | 0.14 | 14.2 | 16 % |
  | 31 | 100 | 8.5 | 28.1 | 24.9 | 97.0 | 99.4 % | 2 | 1.73 | 102 | 28 % |
  | 127 | 8 | 6.6 | 22.7 | 22.7 | 86.3 | 100 % | 1 | 0.14 | 145 | 16 % |
  | 127 | 100 | 95.0 | 301.5 | 283.3 | 1056 | 98.4 % | 7 | 6.08 | 1498 | 20 % |

  GPU (A100-SXM4-40GB):

  | T | L | WENO | adaptive | n = 1 | n = 4 | step | adaptive / step |
  |---|---|---|---|---|---|---|---|
  | 31 | 8 | 0.28 | 0.28 | 0.26 | 0.93 | 1.4 | 19 % |
  | 31 | 100 | 0.09 | 4.44 | 2.57 | 9.70 | 23.4 | 19 % |
  | 127 | 8 | 0.13 | 0.37 | 0.36 | 1.29 | 7.4 | 5 % |
  | 127 | 100 | 1.25 | 17.8 | 5.26 | 18.2 | 42.0 | 42 % |

  On CPU the cost is 3-4× WENO, as expected from 3 stages, and the adaptive mode costs about the same
  as n = 1 because 98-100 % of the columns take one substep. On GPU the adaptive mode at T127 L100
  costs as much as 4 substeps everywhere: the few columns with up to 7 substeps hold back their warps.
  At T31 L100 (~4.6k columns) the column-parallel kernel under-occupies the GPU, while WENO is
  parallel over columns × layers.

## Documentation changes

- Paragraph in `docs/src/primitiveequation.md`: options, stability limit, fixed mode, operator splitting.

## Known limitations

- Parallel over columns only: at low resolution a GPU is under-occupied, and columns with different
  `n` in one warp diverge (one convective column sets the time of its warp).
- The cost depends on the flow unless `courant_target = nothing`.
- Not yet tested with Reactant or Enzyme; the data-dependent loop count of the adaptive mode may need
  the fixed mode there.

## Future work

- GPU efficiency: one workgroup per column with the layers as work items and the column in local
  memory (`@localmem`, `@synchronize` between stages). Divergence would then only be between
  workgroups, and occupancy would match WENO. Alternatively sort or bin columns by `n`.

- Output the number of substeps (maximum over columns) as a diagnostic.
- With the Lorenz N-cycle time stepping for `PrimitiveWet`, revisit whether subcycling is still needed.
