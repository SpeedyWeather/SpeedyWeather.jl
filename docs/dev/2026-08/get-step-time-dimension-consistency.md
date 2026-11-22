# Making `get_step`/`get_steps` consistently select the time dimension

> Status: **completed**. `get_step` now selects the step dimension by dispatching on the time
> dimension `T` in the variable's `ArrayDimensions` (returning a full view when there is none),
> `get_steps` always splats that dimension into a tuple, and a new `nsteps` helper backs both.
> Docstrings were updated to match and the previously untested step machinery gained a test file
> with 122 tests. One pre-existing, unrelated failure in `dynamics/set.jl` was identified and
> left alone; `ColumnField` (`ZXY`/`ZXYT`) remains outside the step machinery, as before.

Date of initial draft: 2026-08-21

Base revision: `5b1e832c` (`mk/landstepping`), the HEAD the source change was drafted against.

Landed as `0373a8d9` "Claude made get_step/get_steps more consistent" and `9186dd39`
"unpack functions for local tests too".

## Originating prompt

> Made some changes to this file, can you update the docstrings?

(referring to `SpeedyWeather/src/time_stepping/steps.jl`, which had gained a runtime `hastime`
branch in the `{T, 2}` methods of `get_step` for both `Field` and `LowerTriangularArray`)

## Revision log

- **2026-08-21, initial request — docstrings.** Updated the two docstrings whose behaviour had
  changed, and dropped a stale "3D field" opening line on a method that takes a 2D field.
- **2026-08-21, "check whether the step functionality is actually tested".** It was not: `get_step`
  appeared in four test files only as a *helper* for reaching into prognostic variables, never as
  the subject of a test, so the new `hastime` branch was entirely uncovered. Added
  `SpeedyWeather/test/variables/steps.jl`. Writing it surfaced two inconsistencies, reported back
  rather than fixed unilaterally:
  1. `get_step(var)` (no step argument) threw a `BoundsError` on any 2D field/LTA, because
     `size(var, ndims(var))` reads the *horizontal* extent when there is no step dimension.
  2. `get_steps` never consulted `hastime`, so a vertical-only `XYZ` field reported its 3 layers
     as 3 "steps" — contradicting `get_step`, which after the `hastime` change returned the full
     array for exactly that type.
- **2026-08-21, "implement the two flagged points ... get_step should always return a view
  selecting the step dimension indicated by T in the ArrayDimensions returning a full view if no
  time exists, and get_steps should always splat the step dimension into a tuple".** Implemented
  as described below. Mid-implementation the dispatch strategy had to change: see *Known
  limitations*, the `WithTime` unions do not dispatch.
- **2026-08-21, "don't bother running the full test suite, that takes too long".** Switched to
  running targeted test files only; recorded as a standing preference.
- **2026-08-21, feedback on module access in tests.** A bare `using SpeedyWeather` in a test file
  loads a *second* copy of the package when the module is already locally loaded as
  `.SpeedyWeather`. The test file now resolves the loaded module at runtime instead.

## Problem description

`get_step(var, step)` selects one step (e.g. one leapfrog time level) from a variable as a view,
and `get_steps(var)` returns all of them as a tuple. Which array dimension is "the step
dimension" was decided inconsistently across the methods:

- For a matrix-like variable the 2nd dimension is ambiguous: it is the step dimension for a 2D
  variable (`XYT`, `LMT` — horizontal + time) but the *vertical* for a 3D one (`XYZ`, `LMZ`).
  The `{T, 2}` methods had just gained a runtime `ArrayDimensions.hastime(var)` branch to tell
  these apart.
- The `{T, 1}` methods (`XY`, `LM` — no step dimension at all) did `view(var, :, step)`, relying
  on a trailing singleton dimension, which only works for `step == 1`.
- `get_steps` ignored the dimension tags entirely and counted the last array dimension, so it
  disagreed with `get_step` for vertical-only variables.
- `get_step(var)` without an explicit step took `size(var, ndims(var))` as "the last step", which
  is the number of grid points or harmonics for a variable with no step dimension — an
  out-of-bounds index rather than a sensible default.

The unifying rule the codebase already had available, but was not applying uniformly, is that the
step dimension is exactly the time dimension `T` in the variable's `ArrayDimensions` tag.

## Background

Both `Field` and `LowerTriangularArray` carry an `ArrayDimensions` tag as their 5th type
parameter, naming each array dimension: `XY`, `XYZ`, `XYT`, `XYZT` for fields and `LM`, `LMZ`,
`LMT`, `LMZT` for spectral coefficients (plus `ZXY`/`ZXYT` for column-major `ColumnField`s).
`ArrayDimensions.hastime` is a compile-time trait over these types, so the distinction the
`hastime` branch made at runtime is available to dispatch.

Note the two array types number their `N` parameter differently, which matters for writing the
method signatures: for `Field`, `N` counts *physical* dimensions with the horizontal unravelled
(`XY` → 1, `XYT`/`XYZ` → 2, `XYZT` → 3), whereas for `LowerTriangularArray` `N` is the ndims of
the backing data array (`LM` → 1, `LMT`/`LMZ` → 2, `LMZT` → 3).

## Summary of changes

All changes are in `SpeedyWeather/src/time_stepping/steps.jl` plus a new test file.

### `get_step` — always selects the time dimension

The runtime `hastime` branch is gone. Every method now names its `Dims` parameter explicitly, so
the step dimension is chosen at compile time from the dimension tag:

| Dimensions | Behaviour |
|------------|-----------|
| `XYT`, `XYZT`, `LMT`, `LMZT` (with time) | view selecting the last (step) dimension |
| `XY`, `XYZ`, `LM`, `LMZ` (without time)  | full view, `step` ignored (should be 1) |

This makes the no-time cases uniform: `XY`/`LM` previously did `view(var, :, step)` and now return
a plain full view (`field_view(var, :)` / `lta_view(var, :)`) like the vertical-only cases already
did.

Plain `AbstractArray`s keep their previous behaviour — they carry no dimension information, so
`step` indexes the last dimension. This path exists because `Adapt` unwraps a `Field` to its bare
device array inside GPU kernels, where a `MethodError` would abort kernel compilation.

### `get_steps` — always splats the step dimension

The three ndims-based methods for `Field`/`LowerTriangularArray` collapse into one:

```julia
get_steps(var::Union{AbstractField, LowerTriangularArray}) = ntuple(step -> get_step(var, step), nsteps(var))
```

A variable with no time dimension yields a 1-tuple holding the full variable as a view, consistent
with `get_step`. `get_steps(field_xyz)` on a 3-layer vertical-only field now returns 1 element
instead of mistaking the 3 layers for 3 steps.

The `Val(N)` compile-time-length method is untouched — it is load-bearing for Enzyme (a runtime
`ntuple` returns a small union of tuple types that breaks Enzyme's type analysis on Julia ≥ 1.11,
[EnzymeAD/Enzyme.jl#3275](https://github.com/EnzymeAD/Enzyme.jl/issues/3275)).

### New `nsteps(var)` helper

Length of the step dimension, or 1 when the variable has no time dimension. It backs `get_steps`
and also fixes the no-argument `get_step(var)`, which is now `get_step(var, nsteps(var))` and
returns the full view for a variable without steps instead of throwing.

### Docstrings

All `get_step`/`get_steps` docstrings now state the rule in terms of the time dimension `T` in the
`ArrayDimensions` rather than "the last dimension", and say explicitly what happens when there is
no such dimension.

## Testing and verification

New file `SpeedyWeather/test/variables/steps.jl`, auto-discovered by `find_tests`, 122 tests,
~7 s. Six testsets:

- **`get_step` on fields** and **on LowerTriangularArrays** — all four dimension tags each, the
  with-time vs without-time split asserted against `ArrayDimensions.hastime`, correct step
  selection checked against explicit indexing, and confirmation that the returned views alias the
  parent rather than copying.
- **`get_step` on plain arrays** — the path taken inside GPU kernels after `Adapt` unwrapping.
- **`get_steps`** — runtime-length and `Val(N)` forms; that no-time variables give a 1-tuple.
- **Number of steps per time stepper** — `DummyTimeStepper` fallbacks, `DEFAULT_NSTEPS`, and
  `get_nsteps` for `Leapfrog` (2 spectral prognostic steps; 2 grid steps for primitive equations,
  1 for 2D models) and `NCycleLorenz` (1 prognostic, 2 spectral tendency steps for the F and G
  terms of Hotta et al. 2016). Note `BarotropicModel` defaults to `NCycleLorenz`, not `Leapfrog`,
  so both steppers are constructed explicitly.
- **`which_step` dispatch** — the fallback of 1, the model-dispatch → no-model fallback chain, and
  `Leapfrog`'s actual choices (step 2 for the spectral transform, step 1 for a parameterization),
  plus `get_prognostic_step`/`get_tendency_step` resolving through them.

Regression coverage for the dispatch change, all passing: `dynamics/time_stepping`,
`dynamics/geopotential`, `dynamics/horizontal_diffusion`, `dynamics/tracer_advection`,
`dynamics/vertical_advection`, `dynamics/copy_variables`, `variables/variables`, `variables/fuse`.

`dynamics/set.jl` errors with a `DimensionMismatch` inside `interpolate!` (an `XYZT` output field
against an `XY` input at `set.jl:3`). This was verified to fail identically at the base revision
with the source change reverted — **pre-existing and unrelated**. `dynamics/forcing_drag` was not
run (a batch invocation timed out before reaching it); the full suite was not run, by request.

### Test-file conventions applied

The test file deliberately contains no `using SpeedyWeather`. `runtests.jl` already does that once
and injects it into each worker via `init_code`, so exported names are in scope; and a bare
`using SpeedyWeather` in a file loads a *second* copy of the package when the module is already
locally loaded as `.SpeedyWeather` (Revise / `include`d into `Main`), whose types then do not
match. Non-exported names are reached by resolving the loaded module at runtime:

```julia
const SW = parentmodule(SpectralGrid)   # SpectralGrid is exported, hence already in scope
using .SW: get_step, get_steps, nsteps, which_step, ...
```

`using .SW: ...` does work against a `const` module binding — `using` resolves the leading `.`
name as a binding in the current module. Verified under both `using SpeedyWeather` and a
locally-bound module.

## Documentation changes

Docstrings in `steps.jl` only; no changes to the user-facing docs, as `get_step`/`get_steps` are
internal to the time-stepping machinery (`get_step` is exported but not documented in the manual).

## Known limitations

**The `...WithTime` type aliases do not work for dispatch.** `LowerTriangularArrayWithTime` and
`AbstractFieldWithTime` are defined as `T{...,Dims} where {..., Dims <: DimensionsWithTime}`. A
concrete `LMT` array *is* a subtype of that alias, but as a method signature the bound does not
make it more specific than the unbounded `LowerTriangularArray`, so the two are ambiguous and
Julia resolves them by definition order. The first implementation of `nsteps` used the aliases and
silently returned 1 for `LMT` variables. All step methods therefore dispatch on the **concrete**
dimension types (`ArrayDimensions.LMT`, `.XYZT`, …), with a comment in the source recording why.
The aliases remain fine for `isa` tests. This is worth checking wherever else in the codebase the
unions are used for dispatch rather than as predicates.

**`ColumnField` is not covered by the step machinery.** Column-major fields carry `ZXY`/`ZXYT`
(vertical first) and now hit a `MethodError` in `get_step`, where they previously matched the
generic `{T, 2}` method. That old path did `field_view(var, :, step)`, which indexed the
*horizontal* dimension — wrong for a column field regardless — and `field_view` cannot construct a
`ColumnField` at all (it builds a `Field`, and `ZXY[:, :]` is not defined). No caller in the
codebase invokes `get_step` on one, so speculative semantics were not invented for it. The failure
mode changed from silently-wrong to a `MethodError`, which is an improvement, but it is a change.

## Future work

- Add `ZXY`/`ZXYT` support to `get_step`/`nsteps` if column fields ever need step selection. This
  requires `field_view` to preserve the `ColumnField` wrapper and `ArrayDimensions` to define the
  corresponding `getindex` reductions for those tags.
- Audit other uses of `DimensionsWithTime`/`DimensionsWithVertical`-bounded aliases for the
  dispatch pitfall described above.
- Investigate the pre-existing `dynamics/set.jl` failure (`interpolate!` with an `XYZT` output
  against an `XY` input), which is independent of this work.
