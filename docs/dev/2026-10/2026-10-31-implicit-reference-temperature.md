# Fix the implicit temperature operator and make its reference temperature selectable

> Status: **completed**. Operator fixed and tested, reference temperature selectable, isothermal 300 K default.

Date of initial draft: 2026-10-31

Base revision: `c147266e`

## Originating prompt

> Yes please separate PR for those. Or could you get reference pressure from mode 0 of spectral surface
> pressure? As it's done similarly for the temperature profile? If it helps they could be scaled from those?

> Can you externalize into one PR all changes you're proposing to implicit, pressure gradient and
> geopotential that can be made independently of the hybrid coordinates, e.g. with a flag that either uses
> an isothermal atmosphere for the vertical profile in implicit (in which case it doesn't need to be
> reinitialized I guess) or the given global profile? Or maybe also the possibility to have global profile
> plus x Kelvin?

## Revision log

- 2026-10-31: initial draft. The reference pressure from mode 0 of ln(pₛ) only changes the operators for
  hybrid coordinates (σ-coordinate weights do not depend on pₛ), so it stays in #1236 together with the
  pressure gradient force and geopotential changes, which are hybrid-only too.

## Problem description

Since #1286 the explicit vertical advection advects the full temperature T = T' + Tₖ and relies on the
semi-implicit operator L to contain the linearised vertical advection of the reference profile Tₖ.
The vertical advection part of L (`L1`), ported from SPEEDY, does not match the explicit tendency:

1. the ΔT factor multiplies only the first of the two terms of the linearised vertical velocity,
2. full-level σ is used where the vertical velocity is defined on half levels (σₖ₊½),
3. the sign is opposite to the explicit `-1/(2Δσₖ)(σ̇ₖ₊½ΔT_below + σ̇ₖ₋½ΔT_above)`.

Cf. SPEEDY `implicit.f90` (`xc = ya + xa*xb`, `xb` the full half-level σ̇ per unit divergence).
For a resting atmosphere T = Tₖ the explicit temperature tendency from a divergence D differs from L·D
by 30-90% on `main`, it agrees to round-off with the fix.

With the fixed operator, linearising about the (non-isothermal, cold-topped) global-mean temperature
profile is unstable for Held-Suarez at T63 (NaN on day 1), consistent with Simmons, Hoskins & Burridge
(1978). A warm isothermal reference (300 K, as in Isca/GFDL; IFS uses ~350 K) is stable at T31 and T63.

## Background

- Simmons, Hoskins & Burridge (1978), Simmons & Burridge (1981): a warm isothermal reference state
  makes the semi-implicit scheme stable for the fast gravity modes.
- Isca `spectral_dynamics.F90`: `ref_temperature_implicit = 300`, full T advected explicitly.
- `docs/dev/2026-10/2026-10-19-implicit-reference-update.md`: periodic recomputation of the global-mean
  reference profile (#1332).

## Summary of changes

- Fix `L1` in `initialize!(::ImplicitPrimitiveEquation, ...)` (ΔT factor, half levels, sign).
- New option `reference_temperature` of `ImplicitPrimitiveEquation`:
  - `IsothermalReferenceTemperature(temperature = 300)`: constant profile, only reinitialized when
    the time step changes (no scheduled recomputation needed),
  - `GlobalMeanReferenceTemperature(offset = 0)`: the global-mean profile plus an offset in Kelvin,
    recomputed on the existing `schedule`.
- Default `IsothermalReferenceTemperature(300)` from the stability tests below.

Stability with the fixed operator, 20 days at nlayers = 8:

| Reference | Held-Suarez T31 | Held-Suarez T63 | PrimitiveWet T31 | PrimitiveWet T63 |
|---|---|---|---|---|
| isothermal 300 K | stable | stable | stable | stable |
| isothermal 250 K | stable | stable | stable | stable |
| global mean | stable | NaN day 1 | stable | stable |
| global mean + 30 K | stable | NaN day 2 | stable | stable |

## Testing and verification

- Unit test: explicit temperature tendency of a resting atmosphere equals L·D (isothermal and
  non-isothermal Tₖ), fails on `main`.
- Unit tests for both reference temperature options (profile set, isothermal not rescheduled).
- Stability: Held-Suarez and default PrimitiveWet at T31/T63, 20 days, for each option.

## Documentation changes

- `docs/src/time_stepping.md`: reference temperature options.

## Known limitations

- Only the centred 2nd-order vertical advection is linearised exactly by L1; for other schemes L1 is
  the centred linearisation.

## Future work

- Hybrid coordinate linearisation (reference surface pressure), in #1236.
