# Fix flaky "Interpolate meridionally-constant field" test in RingGrids

> Status: **completed**. Test switched from relative to absolute tolerance; interpolation itself is correct.

Date of initial draft: 2026-10-30

Base revision: 28069552f493ae7122c90588304b6b7ea427c743

## Originating prompt

> Have a look at https://github.com/SpeedyWeather/SpeedyWeather.jl/actions/runs/41041281830/job/124014420572?pr=1283 I've seen this a few times that the "Interpolate meridionally-constant field" test in RingGrids is flaky, can you investigate why, improve the test and/or the interpolation itself if necessary and push this as a new PR?

## Revision log

> None.

## Problem description

The test occasionally fails for Float32, e.g. `0.00018556092f0 ≈ 0.00018310547f0 (rtol=0.0025)`.
The second half of the test interpolates a field equal to longitude in [-180˚, 180˚E] at random
longitudes in [-120˚, 120˚E], compared with `rtol = 2.5e-3`. When a random longitude lands very
close to 0˚E, the relative tolerance shrinks to nothing while the (absolute) Float32 rounding
error of order `eps(Float32) * 360 ≈ 4e-5` stays the same, so the test fails.

## Background

The field is linear in longitude on every ring, so bilinear interpolation reproduces it exactly
up to rounding. Measured over 10^5 points per grid: max abs error ≤ 6.1e-5 (Float32, ≤ 2 ulp of
360) and ≤ 6e-14 (Float64) for all grids. The interpolation is therefore correct.

## Summary of changes

- `RingGrids/test/interpolation.jl`: compare with `atol = 360 * 4 * eps(NF)` in both halves of
  the test instead of `rtol`. This is both robust near 0˚E and much stricter than before
  (previous rtol allowed errors up to ~0.3˚).

## Testing and verification

Ran the test file repeatedly locally, passing; error statistics above.

## Documentation changes

None.

## Known limitations

None. No source change so no version bump of RingGrids.

## Future work

None.
