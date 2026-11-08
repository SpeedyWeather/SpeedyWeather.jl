# FriersonSigmaCoordinates as default, vertical diffusion off by default

> Status: **in progress**. Defaults switched, PR opened so that CI tests the new defaults across the suite.

Date of initial draft: 2026-10-23

Base revision: 0e9e1275533f0702b789644d807e524acdbfa3a9

## Originating prompt

> Can you open a PR now that makes the FriersonSigmaCoordinates the default so that CI can test whether this causes any other problems? In PR1283 you flagged that the vertical diffusion doesn't actually do anything for now, so instead of dealing with stability issues from there in this new PR, I'd set the vertical diffusion to nothing for both models, and we'll overhaul the vertical diffusion scheme in a follow-up instead.

## Revision log

## Problem description

The default vertical coordinates are equi-spaced sigma levels, which resolve the boundary
layer poorly. `FriersonSigmaCoordinates` refine near the surface and in the stratosphere.
Their thin lowermost layer made explicit surface fluxes unstable (#1335), fixed by the
damped surface flux time integration (#1336), so they can now become the default.

## Background

`BulkRichardsonDiffusion` currently has no meaningful effect (flagged in #1283) and thin
layers are a potential source of instability for an explicit vertical diffusion. Rather than
tuning it here, it is switched off and overhauled in a follow-up.

## Summary of changes

- `Geometry(spectral_grid)` defaults to `vertical_coordinates = FriersonSigmaCoordinates(spectral_grid)`
  (affects all models, for 1 layer this is identical to before).
- `PrimitiveDryModel` and `PrimitiveWetModel` default to `vertical_diffusion = nothing`.

## Testing and verification

- Frierson test in `surface_flux_time_integration.jl` now runs both primitive models with
  the default geometry and checks the default equals `FriersonSigmaCoordinates`.
- Full CI to catch tests that implicitly relied on equi-spaced levels or vertical diffusion.

## Documentation changes

- `vertical_coordinates.md`: Frierson spacing noted as the default.
- `vertical_diffusion.md`: vertical diffusion off by default, how to enable it.

## Known limitations

Changes results of all multi-layer simulations with default settings.

## Future work

Overhaul of the vertical diffusion scheme.
