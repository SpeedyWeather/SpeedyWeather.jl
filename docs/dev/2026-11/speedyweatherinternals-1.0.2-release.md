# Release SpeedyWeatherInternals 1.0.2 to fix unsatisfiable DocStringExtensions compat

> Status: **completed**. SpeedyWeatherInternals bumped to 1.0.2, dependants require it as lower bound.

Date of initial draft: 2026-11-10

Base revision: 03800ef2

## Originating prompt

> With PR1394 we now have the following problem: Unsatisfiable requirements detected for package
> SpeedyWeatherInternals [...] restricted by compatibility requirements with DocStringExtensions
> to versions: uninstalled — no versions left. Can you have a look what's going on and how this could be solved?
> — Yes please [open the PR]

## Revision log

- Initial draft.
- The other packages are not released yet (unreleased v0.24.0), so no `+DEV` bumps for them; CHANGELOG line goes into the 0.24.0 section.

## Problem description

#1310 raised `DocStringExtensions` compat to `0.10` in all five packages. The v0.24.0 release (#1394)
registered LowerTriangularArrays/RingGrids/SpeedyTransforms 1.0.1 and SpeedyWeather 0.24.0 with that
compat, but SpeedyWeatherInternals stayed at `1.0.1+DEV` and was not released. Every registered
SpeedyWeatherInternals version requires `DocStringExtensions = "0.9"`, so outside the monorepo
(no `[sources]`) the resolver finds no solution.

## Background

Inside the monorepo `[sources]` path entries mask the problem, so CI did not catch it.

## Summary of changes

- `SpeedyWeatherInternals` version `1.0.1+DEV` → `1.0.2` (to be registered).
- `SpeedyWeatherInternals = "1.0.2"` compat lower bound in LowerTriangularArrays, RingGrids,
  SpeedyTransforms and SpeedyWeather, so they can never resolve to an old SpeedyWeatherInternals.
  These versions (1.0.1, 0.24.0) are not registered yet, so they are released with this bound, no `+DEV`.

## Testing and verification

Compat-only change; CI. Effective fix once SpeedyWeatherInternals 1.0.2 is registered.

## Documentation changes

None.

## Known limitations

None.

## Future work

Release PRs should include every package whose `[compat]` changed since its last release.
