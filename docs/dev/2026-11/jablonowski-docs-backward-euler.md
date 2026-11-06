# Use backward Euler in the Jablonowski-Williamson docs example

> Status: **completed**. Docs example now uses `centering = 1`, with an explanation.

Date of initial draft: 2026-11-06

Base revision: 2894ecfb8554d8daa031af5f045c76bc6bd63067

## Originating prompt

> Find out why the baroclinic wave docs example breaks up quickly on main (suspect #1369).
> [...] can you open a PR with 2 please (= keep Crank-Nicolson default, set `centering=1`
> in the docs example and explain why)

## Revision log

- None.

## Problem description

Since #1369 changed the default `centering` of `ImplicitPrimitiveEquation` from 1
(backward Euler) to 0.5 (Crank-Nicolson), the southern-hemisphere jet of the
Jablonowski-Williamson example becomes unstable by day 3, making the 9-day plot noisy.

## Background

Diagnostic: max zonal-eddy surface vorticity in the SH (should stay ≈ 0) and max surface divergence.

| centering | SH eddy day 2 | SH eddy day 3 | max div day 2 |
|---|---|---|---|
| 0.5 | 8.0e-6 | 1.6e-4 (saturated) | 8.7e-5 |
| 0.51 | 7.1e-6 | 1.6e-4 | 8.0e-5 |
| 0.6 | 2.8e-6 | 1.4e-4 | 4.1e-5 |
| 0.75 | 7.9e-7 | 2.0e-5 | 1.6e-5 |
| 1.0 | 1.6e-7 | 2.9e-6 | 4.3e-6 |

The SH seed at day 1 (~3e-8) is identical across centering, in Float32 and Float64,
and with or without orography: it is gravity waves radiated by the unbalanced NH
perturbation. Crank-Nicolson does not damp them; the baroclinically unstable SH jet
then amplifies them by ~10-20x/day. With `centering = 1` the SH still breaks up by
day 4 (as before #1369), so the example was always marginal.

## Summary of changes

- `docs/src/examples_3D.md`: pass `ImplicitPrimitiveEquation(spectral_grid, centering = 1)`,
  explain why; fix the text saying `ZeroInitially` (code uses `ConstantPressure`).
- No change to model defaults.

## Testing and verification

Ran the example with `centering = 1` on main (table above).

## Documentation changes

As above.

## Known limitations

The SH jet still becomes unstable around day 4.

## Future work

Check whether other setups (e.g. Held-Suarez) get noisier with the Crank-Nicolson default.
