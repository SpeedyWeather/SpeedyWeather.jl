# Precipitable water output variable

> Status: **in progress**. Implemented with unit tests; awaiting review.

Date of initial draft: 2026-11-01

Base revision: c147266e5ddbc234d1fcd3a1917b220386e6c593

## Originating prompt

> No that's fine but I'm wondering whether we can have an output variable for precipitable water that
> computes the vertical integral on the fly before the interpolation, Can you add that?

## Revision log

## Problem description

Precipitable water (column-integrated water vapor) is a standard diagnostic for a model's moisture. In
#1279 it was estimated offline from model-grid humidity with a fixed 1000 hPa surface pressure. Computing
it from interpolated 3D humidity output would be inaccurate and requires writing the 3D field.

## Background

PW = (1/g) ∑ₖ qₖ Δpₖ, with Δpₖ from `pressure_thickness(k, pₛ, coordinate)`, which supports both sigma
and hybrid sigma-pressure coordinates. Mean sea-level pressure, surface temperature and 10 m winds are
already computed on the model grid within `output!` before the interpolation.

## Summary of changes

- `OutputVariables.PrecipitableWater` (`pw`, kg/m², 2D + time). Its `output!` integrates over the
  vertical in a kernel over grid points into a scratch array, then interpolates and writes it, following
  the 10 m wind output.
- Included in `OutputVariables.Dynamics()` and therefore `All()`, but not in the default output of
  `PrimitiveWet`.
- Silently skipped for models without humidity (`path` points to `grid.humidity`).

## Testing and verification

- Unit test of the kernel: constant q gives q·pₛ/g, including a column with a different pₛ.
- NetCDF output test: `pw` is written with units kg/m² and has a plausible range.

## Documentation changes

The surface variables section in `docs/src/output.md` covers the formula.

## Known limitations

- Water vapor only, no condensate (the model has no prognostic cloud water).

## Future work

- Possibly add it to the default output of `PrimitiveWet`.
