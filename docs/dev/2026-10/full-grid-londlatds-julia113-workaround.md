# Work around a Julia v1.13.2 miscompile in `get_londlatds` for full grids

> Status: **completed**. Full-grid coordinates are built with explicit loops; the reproducer no longer crashes on Julia v1.13.2.

Date of initial draft: 2026-10-24

Base revision: e83b6e1f

## Originating prompt

> I'm getting segfaults when trying to initialize small primitive wet models, like T16, T8, T6 etc
> nlayers=2, 3, 4. And some seem to originate from interpolation. Can you test those and if you see
> some issues, maybe you can nail them down, and if you find the culprit post an issue and a PR to fix?

> I also saw this in the get_londlatds functions. On Julia v1.13.2 can you test these functions for no segfault?

> Can you do 2 only but post the issue in the SpeedyWeather repo with a MWE. I don't want to escalate it higher up yet

## Revision log

- Scope limited to a workaround in RingGrids (no upstream Julia issue yet), see #1347.
- Root cause identified (via @vchuravy): JuliaLang/julia#60924, GC-pointer null stores sunk past a safepoint, fixed on master but not yet backported to 1.13 (JuliaLang/julia#63831). Stress-tested other candidates with `--heap-size-hint=64M`; no further crashes, so no further workarounds.

## Problem description

On Julia v1.13.2, initializing small models segfaults or aborts with `GC error (probable corruption)`,
usually while allocating in `get_londlatds` / `GridGeometry` / `interpolator`. Repeatedly calling
`get_londlatds(FullClenshawGrid, n)` for `n in 1:48` reproduces it with RingGrids alone.

## Background

The full-grid method built the coordinates as `reshape(ones(nlat)' .* lond, :)` and
`reshape(latd' .* ones(nlon), :)`. Reduced to Base Julia (MWE in #1347), this crashes on every run on
Julia v1.13.2 but passes with `-O0`, with `--check-bounds=yes` and on Julia v1.12.7. That points to a
miscompile in Julia rather than a bug in RingGrids.

## Summary of changes

- `get_londlatds(::Type{<:AbstractFullGrid}, nlat_half)` preallocates `londs`, `latds` and fills
  them ring by ring with loops, mirroring the reduced-grid method. Same values and order as before.

## Testing and verification

- New test in `RingGrids/test/grids.jl`: full-grid `get_londlatds` matches the ring-by-ring
  coordinates for several `nlat_half`.
- Reproducer (looping `get_londlatds` over all grid types, `nlat_half` 1:48, 20 times, then
  `GC.gc(true)`) passes on Julia v1.13.2 without bounds checks.

- Stress tests on Julia v1.13.2 with `--heap-size-hint=64M` (makes the original MWE crash in 3 of 3
  runs), 2 runs each, all clean: RingGrids coordinates (`get_londlatds`, `get_lonlats`, `get_latd`,
  `get_lond`), grid construction, `GridGeometry`, `interpolator`, `interpolate`, `average_on_poles`,
  `grid_cell_average!` over all grid types and `nlat_half` 1:32; `SpectralTransform`,
  `initialize!(BarotropicModel)` and `initialize!(PrimitiveWetModel)` over 7 grid types,
  T5/8/16/31, 1/2/4 layers; and 3-day `PrimitiveWetModel` runs at small resolutions.
- Static search for the known triggers (Union-typed ternaries over tuples with GC pointers,
  reshaped broadcasts): the only ternary-tuple hits are bits-only (scalars, indices, ranges),
  which this bug cannot affect.

## Documentation changes

None.

## Known limitations

Other code paths may hit the same miscompile, so this may not remove every crash on Julia v1.13.2.

## Future work

Revert to the vectorised form once the backport of JuliaLang/julia#60924 (JuliaLang/julia#63831) is in a 1.13 release.
