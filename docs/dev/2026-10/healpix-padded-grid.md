# HEALPixPaddedGrid: HEALPix with at least 20 longitudes per ring, nested with HEALPix

> Status: **completed**. Grid, nested `copyto!` both directions and tests implemented in RingGrids.

Date of initial draft: 2026-10-17

Base revision: e410eed05c967a4e117b91c91591d17da155a774

## Originating prompt

> Can you identify the PR with the healpix weight optimisation and read up on the suggested padded
> healpix grid. Can you open a new PR and implement that grid. It should have 20, 20, 20, 20, 20,
> 24, 28, … longitude points per ring towards the poles and those points should nest with the
> original healpix. The latter is important so that we can identify a subset of grid points that
> correspond exactly to those on Healpix. Write a copyto! Method that can convert between these two
> grids but subsampling one way and doing a linear interpolation the other way.

## Revision log

- 2026-10-17: Initial draft. The requested 20 points on rings 2–4 cannot nest (see below);
  implemented the closest nesting alternative, the smallest *multiple* of HEALPix's `nlon` ≥ 20.

## Problem description

PR [#1240](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/1240) found that HEALPix
instability is caused by the polar cap rings having only 4, 8, 12, … longitudes, too few to
represent the zonal wavenumbers fields carry there (Courtier and Naughton 1994: ~20 points needed
on the polar ring). Its experimental `HEALPixPaddedGrid` used `4j + 16` longitudes, which does not
nest with HEALPix: HEALPix points are generally not grid points of the padded grid, so the
two can only be related by Fourier resampling or interpolation.

## Background: why 20 points on rings 2–4 cannot nest

All HEALPix points of ring `j` are equispaced with spacing `360˚/n`, `n = 4j` on the caps.
An equispaced ring with `n'` points contains them only if its spacing divides `360˚/n`, i.e.
`n' = k·n` for an integer `k`. 20 is a multiple of 4 (ring 1) and 20 (ring 5) but not of
8, 12, 16 (rings 2–4). So 20, 20, 20, 20, 20, 24, … is not nestable.

## Summary of changes

New reduced grid `HEALPixPaddedGrid` in `RingGrids/src/grids/healpix_padded.jl`:

- Same latitudes as `HEALPixGrid`.
- Ring `j` has `k_j · n_j` longitudes, where `n_j` is HEALPix's and `k_j = cld(20, n_j)` the
  smallest integer refinement with at least `npoints_min = 20` points. Polar rings are therefore
  20, 24, 24, 32, 20, 24, 28, … — beyond ring 5 identical to HEALPix.
- Longitudes follow HEALPix's convention (`dlon/2` offset or not). A refined ring has an offset iff
  the HEALPix ring has one and `k` is odd; with that choice all HEALPix points are padded points,
  HEALPix point `i` being padded point `k(i-1) + first`, `first = (k + 1)÷2` (odd `k`, offset),
  `k÷2 + 1` (even `k`, offset) or `1` (no offset). `reverse!` keeps working as offsets are
  0 or `dlon/2`.
- Quadrature: HEALPix's equal-area ring weights (each ring covers the same band of colatitude),
  the per-point solid angle divides that by the padded ring's `nlon` (generic `get_solid_angles`).
- `copyto!(::HEALPixField, ::HEALPixPaddedField)` subsamples (exact, nested points);
  `copyto!(::HEALPixPaddedField, ::HEALPixField)` linearly interpolates in longitude along each
  ring (latitudes coincide), exact on the nested points. Both are KernelAbstractions kernels, one
  thread per grid point and layer, for CPU and GPU, with matching `nlat_half` asserted.
- Export from RingGrids and SpeedyWeather; default dealiasing as for HEALPix (3).

## Testing and verification

`RingGrids/test/healpix_padded.jl`: nlons per ring, nesting of coordinates, total npoints,
sum of solid angles = 4π, subsample∘interpolate = identity, interpolation exact for zonally
constant fields and for wavenumber-1 linear profiles on nested points, multi-layer fields.
Added to the generic grid tests list.

## Documentation changes

Docstring of the new grid; CHANGELOG.

## Known limitations

- Not equal-area: refined cap pixels are smaller.
- The non-monotonic ring lengths (20, 24, 24, 32, 20) are the price of exact nesting.
- Linear interpolation HEALPix → padded is only first-order accurate between nested points.

## Future work

- Test model stability with this grid (the purpose of #1240's experiment 2).
- Optionally a Fourier (spectral) resampling as higher-order alternative to linear interpolation.
