# Fix the implicit temperature operator and make its reference temperature selectable

> Status: **completed**. Operator fixed and tested, reference temperature selectable, isothermal 300 K and Crank-Nicolson default, warning for uniform σ spacing.

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

- 2026-11-02: `centering = 0.5` (Crank-Nicolson) as default again, stable with the default
  `FriersonSigmaCoordinates`; a warning for centering < 1 with (near-)uniformly spaced σ levels, the CI
  test with uniform levels uses `centering = 1`.

> Update the plan doc and PR, revert to 0.5 if that's possible with Frierson spacing but maybe raise a
> warning for uniform spacing if combined with 0.5 centering?

> I don't necessarily mind retiring the uniform spacing in favour of Frierson spacing if that's always
> stable, but I'm worried the Frierson spacing is masking another problem

- 2026-11-02: reverted to `centering = 1` (backward Euler). Crank-Nicolson blows up on day 1 for Held-Suarez
  with Earth orography at T32 with 16 uniformly spaced σ levels (CI, `vertical_coordinates.jl`), for both a
  300 K and a 350 K reference; backward Euler is stable there.

- 2026-11-01: default `centering = 0.5` (Crank-Nicolson) instead of 1 (backward Euler), stable in all tests.

> In that PR could you also test whether we can set the centering now finally to 0.5 to have Crank Nicolson
> instead of backward Euler?

- 2026-11-01: `MaximumReferenceTemperature(margin)`, isothermal at the model's maximum temperature plus a
  margin, recomputed on the schedule, to follow the model's climate (e.g. snowball Earth). Corrected the
  IFS reference temperature: IFS documentation CY49R1 Part III, Ch. 2: "Tref = 300 K and prefs = 800 hPa".

> Could you search the ECMWF documentation and tech memorandums to find a reference why IFS uses 350K? Is the
> idea here to use isothermal but slightly warmer because … ? Ideally I would like to tie this to the actual
> temperature in the model, e.g. to do snowball simulations too

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
(1978). A warm isothermal reference (300 K, as in Isca/GFDL and the IFS) is stable at T31 and T63.

## Background

- Simmons, Hoskins & Burridge (1978), Simmons & Burridge (1981): a warm isothermal reference state
  makes the semi-implicit scheme stable for the fast gravity modes.
- Isca `spectral_dynamics.F90`: `ref_temperature_implicit = 300`, full T advected explicitly.
- IFS documentation CY49R1, Part III, Ch. 2: "The reference values adopted for the semi-implicit scheme
  are Tref = 300 K and prefs = 800 hPa" (no rationale given).
- Bénard (2003, MWR), summarising Simmons, Hoskins & Burridge (1978): for an isothermal atmosphere at T̄
  and isothermal reference T*, stability requires 0 ≤ T̄ ≤ 2T*; "a large static-stability for the
  reference-state is necessary"; a warm isothermal reference "was then widely adopted for NWP".
- `docs/dev/2026-10/2026-10-19-implicit-reference-update.md`: periodic recomputation of the global-mean
  reference profile (#1332).

## Summary of changes

- Fix `L1` in `initialize!(::ImplicitPrimitiveEquation, ...)` (ΔT factor, half levels, sign).
- New option `reference_temperature` of `ImplicitPrimitiveEquation`:
  - `IsothermalReferenceTemperature(temperature = 300)`: constant profile, only reinitialized when
    the time step changes (no scheduled recomputation needed),
  - `GlobalMeanReferenceTemperature(offset = 0)`: the global-mean profile plus an offset in Kelvin,
    recomputed on the existing `schedule`,
  - `MaximumReferenceTemperature(margin = 30)`: isothermal at the model's maximum temperature plus
    `margin`, recomputed on the `schedule`.
- Default `IsothermalReferenceTemperature(300)` from the stability tests below.

Stability with the fixed operator, 20 days at nlayers = 8:

| Reference | Held-Suarez T31 | Held-Suarez T63 | PrimitiveWet T31 | PrimitiveWet T63 |
|---|---|---|---|---|
| isothermal 300 K | stable | stable | stable | stable |
| isothermal 250 K | stable | stable | stable | stable |
| global mean | stable | NaN day 1 | stable | stable |
| global mean + 30 K | stable | NaN day 2 | stable | stable |
| maximum + 30 K (332-338 K) | stable | stable | stable | stable |
| maximum + 10 K (312-317 K) | stable | stable | stable | stable |

Crank-Nicolson (`centering = 0.5`) vs backward Euler (`centering = 1`) with the isothermal
300 K reference, 30 days, nlayers = 8: stable for Held-Suarez and PrimitiveWet at T31, T63, T127 with
both. Held-Suarez T127 over 200 days: stable with both, same temperature statistics (0.1% quantiles per
layer within ~2 K). The occasional ~160 K minimum occurs with both centerings, in the top layer over
the Andes (~22°S, 292°E), i.e. it is caused by the orography, not by the centering. However, with 16
uniformly spaced σ levels (Held-Suarez, Earth orography, T32) Crank-Nicolson is unstable (NaN on day 1),
so backward Euler remains the default.

### Crank-Nicolson with uniform σ spacing

All checks with one uninterrupted `run!`: repeated short `run!` calls restart the leapfrog with damping
Euler steps and masked the instability in earlier tests.

Crank-Nicolson, Held-Suarez with Earth orography at T32, 10 days: Frierson spacing stable at 8, 16 and
24 levels; uniform spacing stable at 8, NaN at 16 and 24 levels. Defaults (Frierson, 8 levels) stable
for Held-Suarez and PrimitiveWet at T31, T63, T127 over 30 days. Smallest stable centering with uniform
spacing: 0.55 at 16 levels, 0.9 at 24 levels.

Isolation, 24 uniform levels, centering 0.7, 5 days (NaN unless stated):

| Change | Result |
|---|---|
| none / no forcing and drag / Robert-Williams 0.2, 1.0 | NaN |
| reference temperature 250 K, 300 K, 350 K, global mean | NaN |
| advect only T − Tₖ (before #1286) | NaN |
| WENO vertical advection / no vertical advection of u, v | NaN |
| no vertical advection of T | stable |
| 1st-order upwind vertical advection | stable |
| 4× stronger horizontal diffusion | stable |
| flat orography / half the time step | stable |

The implicit operators R, U, W and L all equal the linearised explicit tendencies to round-off for
uniform and Frierson spacing, so the implicit solver is consistent. The instability grows from the
explicit centred vertical advection of temperature over steep orography (New Guinea first, mid-levels,
explosive after ~15 h), it does not depend on the reference temperature. Backward Euler damps it, and so
does Frierson spacing, which may mask the underlying cause (see Future work).

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
- Root cause of the Crank-Nicolson instability with uniform σ spacing: the explicit centred vertical
  advection of temperature over steep orography. Check at which layers it starts and whether 4th-order
  centred advection or T's vertical advection evaluated at the previous time step avoid it.
