# Year angle from the actual length of the year, and dilated rotation/orbit time

> Status: **completed**. `year_angle` now uses a 0-based day of year and the actual length of
> that year (365 or 366 days), removing a discontinuity of ~0.5 % of a full seasonal cycle at
> every New Year. Two further bugs found while investigating were fixed: the equinox phase in
> `SinSolarDeclination` carried the same off-by-one, and `rotation_time`/`orbit_time` ran a
> permanent Δt/2 ahead of the model time. Tests added and passing; user-facing documentation
> added to `docs/src/radiation.md`.

Date of initial draft: 2026-08-17

Base revision: `a8c31181` (`hw/solartime`)

## Originating prompt

> In this branch we changed the SolarDeclination, year_angle and solar_hour_angle functions to
> take fixed LENGTH_OF_DAY, LENGTH_OF_YEAR in as the functionality from Dates is relative to the
> Earth calendar. Instead we obtain flexibility wrt length of day/year by having an orbit_time
> and a rotation_time that can run faster, slower or even backwards. However, in the functions
> mentioned above, the length of year should be the length of the year of that respective time,
> as some years have 366 instead of 365 days. Can you adapt that and write some tests in
> test/parameterizations/zenith.jl to check that the different times run correctly and that the
> year/hour angles correctly go from 0 to 2pi without jumps

## Revision log

- **2026-08-17, initial draft.** Two decisions were put to the user before implementation:
  - *The Δt/2 offset.* Found during investigation, outside the original request. The user chose
    to fix it in this PR so the clock test could assert exact synchronization rather than
    tolerate a half-timestep of slack.
  - *The equinox phase.* Deriving it from the fixed `year_angle` would have been the simplest
    code, but makes the *year* of `planet.equinox` matter by 0.8 days (leap vs common), against
    its documented "year irrelevant". The user chose to keep the phase on the fixed 365.25-day
    mean year, so `LENGTH_OF_YEAR` was retained rather than deleted.
- **2026-08-17, during implementation.** Two test-authoring corrections, both mine not the
  code's: `Day(365) ÷ 2` truncates to 182 days and so is half a day short of mid-year (switched
  to `Hour(24 * daysinyear) ÷ 2`), and Jun-21 vs Dec-21 are both at midnight so they do not
  exercise the daily cycle at all (switched to `june + Hour(6)` for the rotation/orbit
  independence check).
- **2026-08-17, "add the user-facing controls and concept to the documentation".** Added a
  *Solar zenith angle, length of day and year* section to `docs/src/radiation.md`.

## Problem description

`year_angle` divided by the fixed 365.25-day `LENGTH_OF_YEAR` while its numerator counted actual
calendar days, and used the 1-based `Dates.dayofyear`. Measured on the base revision, as `g/2π`:

| time | before | after |
|---|---|---|
| 2000-01-01T00:00 | 0.00274 | 0.0 |
| 2000-12-31T23:59:59 (leap) | 1.00479 | 0.99999997 |
| 2001-01-01T00:00 | 0.00274 | 0.0 |
| 2001-12-31T23:59:59 (common) | 1.00205 | 0.99999997 |

Three consequences: the year angle was non-zero on Jan-01; it overshot 2π by a year-dependent
amount; and it jumped discontinuously at each New Year — by ~0.998 of a full cycle after a leap
year and ~0.0007 after a common year. The seasonal cycle therefore lurched once a year.

`solar_hour_angle` was checked and found **already correct** — a clean sawtooth from −π at
midnight through 0 at noon to +π just before the next midnight, continuous mod 2π across the day
boundary. It was left semantically unchanged.

Two further defects surfaced during investigation:

1. **Equinox phase.** `SinSolarDeclination` duplicated the same broken formula
   (`LENGTH_OF_DAY * dayofyear(equinox) / LENGTH_OF_YEAR`), so the declination peak sat about a
   day off the equinox and would have become inconsistent with a corrected `year_angle`.
2. **Δt/2 offset in the dilated times.** `leapfrog.jl` rewinds `clock.time -= Δt ÷ 2` to undo the
   Euler spin-up step, but `rotation_time` and `orbit_time` were never rewound. They therefore
   ran a constant 1200000 ms (at default Δt) ahead of the model time — measured to be constant,
   not accumulating, and present for every simulation on this branch.

## Background

`Dates` is inherently tied to the Earth calendar, so this branch deliberately fixes the zenith
angle calculations to Earth's day and year and expresses non-Earth day/year lengths through the
clock instead: `rotation_time` and `orbit_time` advance each step by `Δt` scaled by
`rotation_dilation`/`orbit_dilation`, derived from `planet.length_of_day`/`length_of_year`. That
split is sound; the bug was that the *denominator* of the year angle stayed fixed while the
*numerator* followed the real calendar.

## Summary of changes

`SpeedyWeather/src/parameterizations/radiation/zenith.jl`

- `year_angle` uses a 0-based day of year and divides by `Dates.daysinyear(time) * LENGTH_OF_DAY`,
  i.e. the actual length of that year. The division is done in `Float64` before converting to
  `T`, which also avoids the Float32 precision loss of the old `T(2π)/LENGTH_OF_YEAR`.
- `SinSolarDeclination` fixes the off-by-one but keeps its phase on the mean 365.25-day year, so
  the year of `planet.equinox` stays irrelevant — up to the one day by which a Feb-29 shifts
  `dayofyear`, which is inherent to any date-based phase and is documented at the call site.
- Both angle functions use `secondofday` instead of `Dates.second(Dates.Time(time).instant)`.
  Not cosmetic: the old form only worked via the type-piracy method `Dates.second(::Nanosecond)`
  in `clock.jl`, and `secondofday` is the form already overloaded for Reactant.

`SpeedyWeather/src/time_stepping/clock.jl` and `steppers/leapfrog.jl`

- New `dilate(Δt, dilation)` helper replaces the rounding expression that was duplicated in
  `time_step!`, and the leapfrog spin-up branch now rewinds `rotation_time`/`orbit_time` by their
  dilated Δt/2 alongside `clock.time`.

## Testing and verification

`SpeedyWeather/test/parameterizations/zenith.jl`, 72 assertions, ~1.5 min:

- **`year_angle`** — 0 at Jan-01T00:00, 2π at Dec-31T23:59:59, monotonic across the year, π at
  mid-year, and no jump at the New Year (`mod(g_after - g_before, 2π) < 1e-4`), for Float32 and
  Float64 and for both a leap (2000) and a common (2001) year.
- **`solar_hour_angle`** — −π at midnight, 0 at noon, π at 23:59:59, monotonic through the day,
  no jump across midnight, and the longitude offset adding directly.
- **`cos_zenith!`** — the pre-existing testset was stale and would have errored (it passed 4
  arguments to the now 5-argument signature); fixed, and extended with a case where rotation and
  orbit time differ, to cover that the two are used independently.
- **Clock** — dilations of 1×, 2×/3.6525× and −1× (backwards rotation); all three times
  synchronized after `initialize!`; `set!` re-synchronizing a deliberately desynced clock; and,
  now that the Δt/2 bug is fixed, the exact invariant
  `rotation_time - start ≈ dilation * (time - start)` at `rtol = 1e-5`.

Each new `year_angle` assertion was checked against the old formula and **all fail on it**, so
the tests pin the fix rather than merely passing. The Δt/2 rewind was verified to undo the
forward step exactly for dilations 1.0, 2.0, 3.6525, 0.5 and −1.0 across several Δt including
odd millisecond counts.

No regressions: `all_parametrizations` 42/42 and `shortwave_radiation` pass, consistent with no
hard-coded declination or `cos_zenith` values existing anywhere in the test suite or docs.

```bash
julia --project=SpeedyWeather --check-bounds=yes \
  -e 'using Test, SpeedyWeather, Dates; include("SpeedyWeather/test/parameterizations/zenith.jl")'
```

## Documentation changes

New section *Solar zenith angle, length of day and year* in `docs/src/radiation.md` — previously
the planet, the solar zenith angle and the day/year lengths were undocumented anywhere in
`docs/src`. It covers the `length_of_day`/`length_of_year` controls, the concept of dilated
`rotation_time`/`orbit_time` (with a worked 10-day-long-day example), backwards rotation and
tidal locking, the orthogonal `daily_cycle`/`seasonal_cycle` switches, and `axial_tilt`,
`equinox` and `solar_constant`.

## Known limitations

- The equinox phase still shifts by up to one day depending on the year of `planet.equinox`,
  because `dayofyear(Mar-20)` is 79 in a common year and 80 in a leap year. This is inherent to a
  date-based phase and is roughly 4× smaller than the pre-fix spread.
- `year_angle` still assumes an Earth calendar by construction. That is the deliberate design of
  this branch — non-Earth year lengths come from `orbit_dilation`, not from this function.
- `Dates.second(::Nanosecond)` and friends in `clock.jl` remain type piracy. Nothing in
  `zenith.jl` needs them after this change, but the Reactant extension still does.

## Future work

- `SolarZenith`/`SolarZenithSeason` are exported but have no `@docs` entry in `docs/src/api.md`,
  so they cannot be `@ref`-linked (the new docs section names them in plain code font instead).
  Adding them there would also satisfy `checkdocs = :exports`.
- The dilated times are advanced incrementally, rounding to whole milliseconds each step. This is
  exact for the dilations tested, and even a pathological one drifts only ~8 s over 100k steps
  (measured; dominated by the `Float32` dilation itself rather than the rounding), so it is
  negligible in practice. Computing them as `start + dilation * (time - start)` would remove the
  drift entirely, at the cost of no longer supporting a dilation that changes mid-run.
