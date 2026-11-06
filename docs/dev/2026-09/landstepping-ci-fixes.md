# CI fixes for the land time stepping PR (#1183)

> Status: **completed** (pending CI). Fixes the `Variables set!` test, the docs build and the
> Metal GPU spectral transform test failing on #1183.

Date of initial draft: 2026-09-25

Base revision: `fbed747d` (`mk/landstepping`, PR #1183)

## Originating prompt

> Great, can you look at the other fails too? Don't run any code, just summarize what's going on
> and propose a fix if easy

> Yes can you create this PR please

## Revision log

- Diagnosis from the CI logs of #1183, then the proposed fixes as described below.

## Problem description

1. `test/dynamics/set.jl` hard-coded `step = 2` also for ocean and land variables, which have only
   one step with the new default `EulerForward` (#1264): `BoundsError`.
2. Docs build: the merge of `main` (#1222 moved "Time steppers and variable steps" and "Lorenz
   N-cycle" from `time_integration.md` to `time_stepping.md`, #1239 replaced footnotes by
   citations) kept #1183's copies of these sections in `time_integration.md`, so `@id steps`
   and `@id ncycle` were defined twice and the `[^Hotta2016]` footnote undefined.
3. Metal GPU "single-layer step view" test: the test's `LowerTriangularArray` had the default `LM`
   tag with 2D data, so the dims-aware `get_step` of #1183 fell back to a plain `SubArray`.
   Additionally, `on_architecture` dropped the dimension tag of a `LowerTriangularArray`.

## Background

Since #1183 `get_step` dispatches on the dimension tag (`LM`, `LMZ`, `LMT`, `LMZT`), so the tag
has to be correct and preserved when moving arrays between architectures.

## Summary of changes

- `test/dynamics/set.jl`: ocean and land use their last (=current) step, `nsteps` of the variable.
- `docs/src/time_integration.md`: removed the duplicated sections from "Time steppers and variable
  steps" onwards. `docs/src/time_stepping.md`: added the "Ocean, sea ice and land" subsection,
  #1183's `get_step`/`nsteps` paragraph and `get_steps` docstring, and updated the Leapfrog and
  NCycleLorenz step snippets (`prognostic_steps(::AbstractLeapfrog) = 2`,
  `tendency_steps(::NCycleLorenz) = 2`).
- `LowerTriangularArrays`: `on_architecture(arch, L)` keeps `L.dims` (version already `+DEV`).
- `test/GPU/spectral_transform.jl`: the step-view test tags its array `LMT`.

## Testing and verification

- `test/dynamics/set.jl` and `LowerTriangularArrays/test/lta_trait_dispatch.jl` (new testset:
  `on_architecture` keeps `LM`, `LMZ`, `LMT`, `LMZT` for CPU and JLArrays) pass locally.
- Docs build and GPU tests on CI.

## Documentation changes

See above.

## Known limitations

- Buildkite (CUDA) was not inspected, assumed to be the same GPU test.

## Future work

None.
