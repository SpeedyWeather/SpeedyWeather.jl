# Horizontal heat diffusion for the slab ocean

> Status: **in progress**. Implemented with unit tests; awaiting review.

Date of initial draft: 2026-11-21

Base revision: 7626f058bf27fb0a43472830fb8fcb955780de68

## Originating prompt

> Build stage 0 of a dynamic slab ocean: a conservative horizontal diffusion of mixed-layer heat on ring grids,
> off by default so that existing behaviour is unchanged.

## Revision log

- Diffusivity parametrized by the depth-integrated D H, Codron's rule that D varies inversely with the total depth
  (D = 25000·50/(Hs+Hd) in his 2 layers), applied to the single layer; automatic substepping instead of a warning.
- After review: work arrays are scratch variables, `set!(model; land_sea_mask, Δt)` rebuilds the stencil, substeps
  use the ocean's time step, a bounds check catches a missing `initialize!`, land cells don't enter the interpolation,
  documented that sea ice breaks heat conservation of the coupled system, tests also on the HEALPix grid.

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
- The transport is set by the depth-integrated diffusivity `diffusivity_depth` = D H = 1.25e6 m³/s (Codron's
  25000 m²/s at 50 m), flux = ρ c D H ∇T, so it does not depend on the mixed-layer depth; `diffusivity` [m²/s] overrides.
  A later 2-layer scheme can split D H between layers.
- Explicit time stepping with a fixed number of substeps chosen at `initialize!` so that the stability number is below 0.5.

## Testing and verification

Heat conservation, uniform field, no flux into land, accuracy of the Laplacian of a spherical harmonic, default
output unchanged, stability number; see `test/parameterizations/ocean_heat_transport.jl`.

## Documentation changes

A section "Ocean heat transport" in `docs/src/ocean.md`.

## Known limitations

- Explicit in time (substeps fixed at `initialize!`, rebuilt by `set!(model; Δt)`), mixed-layer depth is uniform.
- Changing `model.land_sea_mask` other than through `set!` after `initialize!` needs a new `initialize!`.
- `ThermodynamicSeaIce` has no energy budget, so under sea ice heat is not conserved by the coupled system.
- Rings at the poles (e.g. a Clenshaw grid with poles) are not supported.

## Future work

Ekman transport with an upwelling deep layer, energy-based sea ice melt.
