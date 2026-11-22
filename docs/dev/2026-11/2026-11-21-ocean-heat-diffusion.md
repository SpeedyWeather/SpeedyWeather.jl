# Horizontal heat diffusion for the slab ocean

> Status: **in progress**. Implemented with unit tests; awaiting review.

Date of initial draft: 2026-11-21

Base revision: 32fbfa49d5c9e244340109cdb403e8642b3f90f2

## Originating prompt

> Build stage 0 of a dynamic slab ocean: a conservative horizontal diffusion of mixed-layer heat on ring grids,
> off by default so that existing behaviour is unchanged.

## Revision log

## Problem description

The `SlabOcean` has no horizontal heat transport, so tropics are too warm and poles too cold unless a
q-flux is prescribed.

## Background

Codron (2012, Clim. Dyn., doi:10.1007/s00382-011-1031-3) pairs diffusion (eddies) with Ekman transport (tropics).
This is the diffusion part only, which is useful on its own (Donnadieu et al. 2006 use the same idea).

## Summary of changes

- `SlabOcean` has a `heat_transport` component: `NoOceanHeatTransport` (default, a no-op that compiles away) or
  `OceanHeatDiffusion`.
- `OceanHeatDiffusion` precomputes at `initialize!` a padded sparse stencil (columns, weights) on the ring grid:
  zonal faces between neighbours, meridional faces between rings with both rings linearly interpolated in longitude
  to the midpoints of the face segments, every face weighted by the smaller ocean fraction of its two cells.
  The stencil is applied in a kernel that adds to the SST tendency.
- Diffusivity `25000 * 50 / (H + 150)` m²/s by default, a stability number is checked at `initialize!`.

## Testing and verification

Heat conservation, uniform field, no flux into land, accuracy of the Laplacian of a spherical harmonic, default
output unchanged, stability number; see `test/parameterizations/ocean_heat_transport.jl`.

## Documentation changes

A section "Ocean heat transport" in `docs/src/ocean.md`.

## Known limitations

- Explicit in time (stability number reported), mixed-layer depth is uniform.
- The land-sea mask has to be set before `initialize!`.
- Rings at the poles (e.g. a Clenshaw grid with poles) are not supported.

## Future work

Ekman transport with an upwelling deep layer, energy-based sea ice melt.
