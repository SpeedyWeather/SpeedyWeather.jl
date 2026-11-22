# Reevaporation per unit mass in ImplicitCondensation

> Status: **completed**. Reevaporation scales with layer thickness Δσ, making it independent of the number of layers.

Date of initial draft: 2026-11-22

Base revision: 32fbfa49d5c9e244340109cdb403e8642b3f90f2

## Originating prompt

> Can you have a look at issue 1426 and see what you think? Then respond? I believe suggestion 2 is good, if you agree, go open a PR with the fixes.

## Revision log

- Test uses `relative_humidity_threshold = 1` so that the subsaturated layers below do not condense.

## Problem description

`ImplicitCondensation` evaporates a fraction `min(c * max(q* - q, 0), 1)` of the falling rain in
every layer (#1426). The total reevaporation through a column therefore grows with the number of
layers: at L100 rain passes ~12x more layers than at L8 and almost all of it evaporates.

## Background

Suggestion 2 of #1426: scale the evaporated fraction with the layer's mass, `1 - exp(-c dq Δp/Δp_ref)`.

## Summary of changes

- Evaporated fraction is `1 - exp(-c * max(q* - q, 0) * Δσ_k / Δσ_ref)` with `Δσ_k = Δp_k / pₛ`.
  For constant `q* - q` the product over layers depends only on the total Δσ, not on nlayers.
- New parameter `reevaporation_layer_thickness = 1/8` (Δσ_ref), so that L8 keeps roughly its previous
  behaviour (for small `c dq` the two formulations agree, `1 - exp(-r) ≈ r`; for larger `r` the
  exponential evaporates slightly less than the old linear cap, e.g. 0.26 vs 0.3 at r = 0.3).
- The exponential replaces the `min(…, 1)` cap: the fraction is always in [0, 1).

## Testing and verification

New test: rain formed in the top layer falls through layers all `dq` below saturation; the surface rain
rate matches `exp(-c dq Σ_{k≥2} Δσ_k / Δσ_ref)` times the rain formed, for both L8 and L64.

## Documentation changes

`docs/src/large_scale_condensation.md` formula and text updated.

## Known limitations

Results change slightly at L8 too (see above). Other layer-count dependencies from #1426 (entrainment,
parcel start, cloud cover) are not addressed here.

## Future work

Suggestions 1, 3, 4, 5 of #1426.
