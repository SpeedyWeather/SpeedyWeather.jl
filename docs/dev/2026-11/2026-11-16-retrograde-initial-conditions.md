# Initial conditions for retrograde rotation

> Status: **completed**. The Jablonowski–Williamson and balanced zonal initial conditions use a jet in the sense of rotation, so a retrograde planet starts from the mirror image of the prograde state and no longer blows up.

Date of initial draft: 2026-11-16

Base revision: b06af911

## Originating prompt

(From Paleo Planet Sim, a worldbuilding app that runs SpeedyWeather for user-designed planets,
where a retrograde planet NaN'd in the warm-up even when started from rest.)

> Can you scope out the retrograde fix? Start with understanding where the failure mode is,
> scanning the GitHub issues and PRs to see if anyone has addressed it, and if not scope out how
> you would fix it as a PR.

> Great, let's build the branch and draft the PR text.

## Revision log

- 2026-11-16: initial plan and implementation.

## Problem description

A planet with `rotation < 0` (retrograde) and the default initial conditions starts with its
poles warmer than the tropics. With `StartFromRest` the wet model blows up within 10 days
(T32, 8 layers):

| rotation | initial conditions | surface temperature at start | 10 days |
|---|---|---|---|
| +Ω | `BalancedZonalState` (default) | 233–300 K | fine |
| −Ω | `BalancedZonalState` (default) | 258–323 K, warm poles | survives |
| +Ω | `StartFromRest` | 225–309 K | fine |
| −Ω | `StartFromRest` | 266–346 K, warm poles | **NaN** |
| −Ω | `StartFromRest` with `JablonowskiTemperature(u₀ = 0)` | 287 K | fine |
| −Ω | `BalancedZonalState(u₀ = -56)` | 233–300 K, as prograde | fine |

So the dynamics and physics handle a negative rotation. The initial temperature does not.

No issue or PR mentions retrograde or negative rotation.

## Background

`JablonowskiTemperature`, `BalancedZonalState` and the `ZonalRidge` orography are in balance with a
zonal jet of speed `u₀ > 0`. The temperature anomaly has two terms, the curvature term `2A(φ)U`
(quadratic in the wind) and the Coriolis term `B(φ)RΩ` (linear in wind × rotation). With Ω < 0 and
`u₀` unchanged only the Coriolis term flips, which reverses the equator-to-pole gradient: near the
surface the poles start about 60 K too warm instead of about 60 K too cold. With
`ConstantRelativeHumidity` that makes a very moist 346 K polar column.

(Ω, u) → (−Ω, −u) is a symmetry of the equations (the mirror image λ → −λ), so the balanced
state for a retrograde planet is the prograde one with the jet reversed: the same temperature,
easterlies instead of westerlies.

`StartFromRest` uses `JablonowskiTemperature` without its jet, so it is not in balance on any planet.
On a 12-hour day the Coriolis term doubles and the poles start too cold. That is a separate question
(see Future work).

## Summary of changes

- `jet_speed(u₀, rotation) = copysign(u₀, rotation)`: westerly for prograde, easterly for retrograde.
- Used for `u₀` in `ZonalWind`, `JablonowskiTemperature`, `BalancedZonalState` (vorticity and
  temperature, now passed to the temperature kernel as an argument) and `ZonalRidge`.
- Prograde planets are bit-identical (`copysign(u₀, Ω) = u₀` for Ω > 0).
- The Jablonowski–Williamson perturbation is unchanged, so it is not mirrored.

## Testing and verification

- `test/dynamics/initial_conditions.jl`, "Retrograde rotation":
  - with the perturbation off, `BalancedZonalState` (wet) and `ZonalWind` + `JablonowskiTemperature`
    (dry) give the same temperature and the negated vorticity for −Ω as for +Ω;
  - `StartFromRest` gives the same temperature for −Ω;
  - `ZonalRidge` gives the same orography;
  - a retrograde `PrimitiveWetModel` from rest runs 10 days without NaNs (before: NaN).
- Offline mirror check (aquaplanet, no daily cycle, T32, 120 days, prograde vs retrograde with the
  mirrored jet): zonal-mean winds mirrored and zonal-mean temperatures within 1.4 K RMS for a single
  snapshot, as expected from the chaotic eddies.

## Documentation changes

- `initial_conditions.md`: a paragraph on retrograde rotation in the `BalancedZonalState` section.

## Known limitations

- The solar hour angle doesn't depend on the sign of `rotation`, so on a retrograde planet the sun
  still rises in the east. That only affects the phase of the daily cycle.

## Future work

- `StartFromRest` without the jet's thermal wind terms (`JablonowskiTemperature(u₀ = 0)`), so that
  a resting start is balanced on any planet. This changes the default rest start, so it's left
  for discussion. Possibly related to the polar ringing in #1079.
- The sun rising in the west on retrograde planets.
