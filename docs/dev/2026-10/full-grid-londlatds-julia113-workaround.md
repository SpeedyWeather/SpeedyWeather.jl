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

## Documentation changes

None.

## Known limitations

Other code paths may hit the same miscompile, so this may not remove every crash on Julia v1.13.2.

## Future work

Report the MWE to JuliaLang/julia if it persists in later 1.13 patch releases.
