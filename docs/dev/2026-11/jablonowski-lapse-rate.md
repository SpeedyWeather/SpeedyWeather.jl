# Jablonowski-Williamson temperature with a stable lapse rate

> Status: **completed**. `JablonowskiTemperature` uses Γ = 5 K/km (new option `lapse_rate`), and the baroclinic wave no longer breaks up.

Date of initial draft: 2026-11-06

Base revision: c96ec480

## Originating prompt

> Back to the instability of the baroclinic wave I just opened issue 1388, v0.23 was stable with centering = 0.5 and 1 but since then some PR introduced something that made this example fail. Can you go systematically through the merged PRs since v0.23 and identify which one introduced this?

> Okay I came to the same conclusion that PR1286 is the culprit. There has been some related changes since in the PRs 1332, 1369 and 1346 that we don't want to revert entirely (e.g. FriersonSigmaCoordinates were generally a good move). Can you investigate what's going on? You should be able to test by simply checking that vorticity is within plusminus 2e-5 on the southern hemisphere

## Revision log

> Yes please, would you recommend to still use the dry adiabatic lapse rate of 9.8K/km in model.atmosphere or should this be changed?

Kept `atmosphere.dry_lapse_rate` unchanged; the fix is local to `JablonowskiTemperature`.

## Problem description

Since #1286 the Jablonowski-Williamson baroclinic wave docs example (`PrimitiveDryModel`, T32, 8 layers)
breaks up within about 3 days in both hemispheres, independent of the semi-implicit centering
(0.5 or 1), the time step, or the vertical advection scheme. See #1388.

## Background

A run at v0.23.0 and at each of the 29 merged commits since then that change package source
shows that the example is bit-identical to v0.23.0 up to #1286 and breaks up from then on.

`JablonowskiTemperature` built its vertical profile T₀σ^(R_dΓ/g) with Γ = `lapse_rate(model.atmosphere)`,
which for `PrimitiveDryModel` is the dry-adiabatic `dry_lapse_rate = 9.8 K/km`. So the
troposphere is neutrally stratified (dlnT/dlnσ = κ), whereas Jablonowski and Williamson (2006)
use Γ = 5 K/km. Before #1286 the vertical advection acted only on the anomaly from the global-mean
profile, so the model never saw the stratification of the mean profile
(effectively stability κT/σ, much too stable). #1286 correctly advects the full temperature,
the model then sees the neutral stratification: the slowest internal gravity modes of the
linearised operator `R*L + U*W'` have c² < 0, and everything grows about 10× per day.

Diagnostics (max |ζ| on the southern hemisphere, surface layer):

| run on main | day 9 |
|---|---|
| as is, centering 0.5 or 1, Δt halved, upwind advection | ~2e-4 (by day 3) |
| no mean-profile vertical advection (v0.23 behaviour) | 4.8e-6 |
| Γ = 5 K/km, centering 0.5 | 5.4e-6 |
| Γ = 5 K/km, centering 1 | 5.9e-6 |

## Summary of changes

- `JablonowskiTemperature` gets an option `lapse_rate = 0.005` [K/m] used for its profile instead of `lapse_rate(model.atmosphere)`.
- Docs example back to the default `centering`, and the gravity-wave paragraph from #1386 replaced.

## Testing and verification

- Test in `test/dynamics/initial_conditions.jl` that the dry model's JW profile is stably stratified (dlnT/dlnσ < 0.8κ) and that the `lapse_rate` option is used.
- Docs example: southern-hemisphere |ζ| < 2e-5 through day 9 for centering 0.5 and 1.

## Documentation changes

Jablonowski-Williamson example in `docs/src/examples_3D.md`.

## Known limitations

`BalancedZonalState` and `PressureOnOrography` also use `lapse_rate(model.atmosphere)`, so for
`PrimitiveDryModel` their reference profile is neutrally stratified, too. Not changed here.

## Future work

Decide whether initial conditions in general should use a stable reference lapse rate (e.g. 6.5 K/km)
instead of the dry-adiabatic one.
