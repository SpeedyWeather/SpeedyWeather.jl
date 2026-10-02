# ImplicitCondensation: melting heat, melting units, Clausius-Clapeyron derivative

> Status: **completed**. All three issues from #1276 confirmed by failing tests and fixed.

Date of initial draft: 2026-09-25

Base revision: fd1ef0d9

## Originating prompt

> Now, relatedly can you look at issue 1276 where Urs raises some concerns on the condensation
> scheme? If you agree with his points, open a new PR to fix them and then respond?

## Revision log

- 2026-09-25: initial draft, executed as planned.

## Problem description

@ursho2552 reports in #1276 three issues in `large_scale_condensation!`:

1. The melting cooling `δT = -Lᵢ/cₚ δq_melt` is overwritten by `δT = -Lᵥ/cₚ δq` before it is
   added to `temp_tend`. Freezing heat is released when snow forms, but it is never taken up
   when that snow melts lower down, so the column gains energy.
2. The melting `δT` is a temperature change over the time step [K], but it is used as a rate:
   `T = temp + Δt * δT`, which can give nonsensical (even negative) `T` in `dqsat_dT ∝ 1/T²`.
3. `dqsat_dT = qsat * RH * Lᵥ/cₚ / (Rᵥ T²)` has an extra `1/cₚ`. Clausius-Clapeyron gives
   `dqsat/dT = qsat Lᵥ / (Rᵥ T²)` [1/K]; the implicit factor (Frierson et al. 2006, eq. 21) is
   `1 + Lᵥ/cₚ dqsat/dT`. With the extra `1/cₚ` the implicit correction is ~1000x too weak.

## Background

All three were confirmed by reading the code. Point 3 is a behaviour change: the implicit
correction now actually damps condensation in a single time step (factor ~2-3 in the tropics),
as the scheme was designed to.

## Summary of changes

- Melting cooling as a tendency: `δT = -Lᵢ/cₚ * δq_melt / Δt` [K/s].
- `T = temp + Δt * δT` is now the temperature after melting [K].
- Condensation heating is added to (not overwriting) the melting cooling: `δT += -Lᵥ/cₚ * δq`.
- `dqsat_dT` uses `Lᵥ` instead of `Lᵥ/cₚ`.

## Testing and verification

New tests in `test/parameterizations/large_scale_condensation.jl` on a prescribed column
(cold, supersaturated layers aloft producing snow, falling into warm layers below that melt it):

- Column enthalpy budget: `cₚ/g Σ dT/dt Δp = Lᵥ/g Σ (-dq/dt) Δp + Lᵢ ρ snow_rate_surface`
  (fails before the fix, points 1 and 2).
- Implicit correction: in a layer without incoming precipitation the humidity tendency equals
  `δq_cond / ((1 + Lᵥ/cₚ dqsat/dT) τ Δt)` with `dqsat/dT = RH qsat Lᵥ/(Rᵥ T²)` (point 3).

## Documentation changes

`docs/src/large_scale_condensation.md` already had the correct equations; added a sentence that
the melting cooling is added to the temperature tendency.

## Known limitations

The budget ignores the small heat content of precipitation itself, as the model does elsewhere.

## Future work

None.
