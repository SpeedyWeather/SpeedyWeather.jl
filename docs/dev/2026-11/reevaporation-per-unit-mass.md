# Reevaporation per unit mass and implicit in ImplicitCondensation

> Status: **completed**. Reevaporation scales with layer thickness Δσ (independent of nlayers), is implicit separately from condensation and conserves water.

Date of initial draft: 2026-11-22

Base revision: 32fbfa49d5c9e244340109cdb403e8642b3f90f2

## Originating prompt

> Can you have a look at issue 1426 and see what you think? Then respond? I believe suggestion 2 is good, if you agree, go open a PR with the fixes.

## Revision log

> Make a separate PR for convection and clarify this in the issue. For the large scale condensation PR can you also include the bug raised by issue 1420? The idea is to treat also reevaporation implicitly but separate from the implicit condensation and then add them at the end. At the moment we mix the two in the implicit formulation which isn't consistent. Verify this inconsistency for yourself first and see whether you agree.

- Added #1420: reevaporation was removed from the rain flux in full but entered humidity only after division
  by the implicit factor `(1 + Lᵥ/cₚ dqsat/dT) * time_scale ≥ 3`, losing ≥ 2/3 of it (water and latent cooling).
  Also, with condensation and evaporation in the same layer, evaporation reduced the new rain a second time.
  Now: condensation implicit as before (relative to threshold); reevaporation implicit separately,
  `δq_e = δq_e* dq / (dq + γ δq_e*)` with `γ = 1 + Lᵥ/cₚ dqsat/dT` (at 100%), so evaporation cannot oversaturate,
  no `time_scale` division; exactly `δq_e` is removed from the rain flux; tendencies added at the end.

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

- Column water conservation (#1420): humidity change + surface precipitation = 0 to 1e-10 relative,
  and no layer evaporated beyond saturation, at L8 and L64.

New test: rain formed in the top layer falls through layers all `dq` below saturation; the surface rain
rate matches `exp(-c dq Σ_{k≥2} Δσ_k / Δσ_ref)` times the rain formed, for both L8 and L64
(little rain so the implicit limiter is negligible, rtol 1e-3).

## Documentation changes

`docs/src/large_scale_condensation.md` formula and text updated.

## Known limitations

Results change slightly at L8 too (see above). Other layer-count dependencies from #1426 (entrainment,
parcel start, cloud cover) are not addressed here.

## Future work

Suggestions 1, 3, 4, 5 of #1426.
