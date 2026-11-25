# Reevaporation in ImplicitCondensation: close the water budget and scale with layer thickness

> Status: **completed**. Reevaporation conserves water, cannot oversaturate a layer, and scales with layer thickness Δσ so it is independent of the number of layers.

Date of initial draft: 2026-11-23

Base revision: babf64cb33d1fd80c300777c822c5d97d1d7202c

## Originating prompt

> Can you take a look at this fix and compare it to mine? Which is better? (PR #1427)

> Can you please give me a step by step of how to implement the layer thickness update into my PR?

> Have I done everything correctly? Please update the docs for me and include a plan doc.

## Revision log

- Started as a fix for the water budget only (#1420). Then added the layer-thickness scaling of
  reevaporation (#1426, suggestion 2), following the approach of #1427.
- Kept a hard cap on reevaporation (evaporate at most up to saturation, accounting for evaporative
  cooling) instead of the smooth backward-Euler factor `δq* dq / (dq + γ δq*)` in #1427.
  Away from saturation the cap leaves the tuned `reevaporation` rate unchanged.
- Saturation humidity is recomputed after snow melting and used for both condensation and
  reevaporation. #1427 uses only the post-melt temperature, in `dqsat_dT`.
- Copied the unit tests from #1427 (layer-count independence, column water conservation, no
  supersaturation from evaporation).
- The hard cap was later replaced by the smooth limiter `E = E_trial C / (E_trial + C)`.
- Review (#1431): docs now describe the smooth limiter, clarify units [m/s] of the rain flux,
  the meaning of `c` with the layer-thickness scaling, and why reevaporation is limited but not
  relaxed over `Δt_c`; the water-budget test moved into its own testset, now also with snow melting.

## Problem description

1. **Water not conserved (#1420).** The evaporated rain was removed from the rain flux in full, but
   was added to the humidity tendency only after division by the implicit condensation factor
   `(1 + Lᵥ/cₚ ∂q*/∂T) * time_scale ≥ 3`. So at least 2/3 of it was lost, both the water and its
   latent cooling. Also, `δq_rain = min(0, δq)` used the sum of condensation and evaporation,
   so in a layer with both, evaporation reduced the new rain a second time.
2. **Reevaporation depends on the number of layers (#1426).** In every layer, a fraction
   `min(c max(q* - q, 0), 1)` of the falling rain evaporates. More layers means more evaporation
   steps, so at high vertical resolution almost all rain evaporates before reaching the ground.

## Background

`ImplicitCondensation` relaxes humidity towards `relative_humidity_threshold * q*` with an implicit
correction for the latent heat feedback (Frierson et al. 2006, eq. 21). Reevaporation is not a
relaxation with a time scale. It is limited by the rain falling into the layer, so it should not be
divided by that implicit factor.

## Summary of changes

`SpeedyWeather/src/parameterizations/large_scale_condensation.jl`:

- New parameter `reevaporation_layer_thickness = 1/8` (Δσ_ref).
- Evaporated fraction is `1 - exp(-c max(q* - q, 0) Δσ_k / Δσ_ref)` with `Δσ_k = Δp_k / pₛ`.
  For constant `q* - q` the fraction reaching the surface is `exp(-c (q* - q) ΣΔσ / Δσ_ref)`, which
  is independent of the number of layers. The exponential stays below 1, so `min(…, 1)` is no
  longer needed.
- Reevaporation is capped at `max(q* - q, 0) / (1 + Lᵥ/cₚ ∂q*/∂T)`, converted to a rain flux, so a
  layer cannot be evaporated into beyond saturation.
- Exactly the evaporated rain is removed from the flux and added to humidity. The implicit
  correction applies only to condensation.
- Precipitation in a layer comes from condensation only.
- Saturation humidity is recomputed after the latent cooling from snow melt.

## Testing and verification

New testset "Reevaporation independent of number of layers" in
`SpeedyWeather/test/parameterizations/large_scale_condensation.jl`, at L8 and L64:

- Rain from a supersaturated top layer falls through layers all `dq` below saturation. The surface
  rain matches the analytic exponential (rtol 1e-3).
- Column water: humidity change + surface precipitation = 0 to 1e-10 relative.
- No layer ends up supersaturated by evaporation.

## Documentation changes

`docs/src/large_scale_condensation.md`, section Re-evaporation: new formula with Δσ scaling, the
evaporation cap, water conservation, and the post-melt saturation humidity.

## Known limitations

- Results change at L8 too: for large `r` the exponential evaporates slightly less than the old
  linear fraction (0.26 vs 0.30 at r = 0.3). At L16 and above, reevaporation is much weaker than
  before.
- With `relative_humidity_threshold < 1`, a layer between the threshold and 100% humidity can
  condense and reevaporate in the same time step. Water is conserved, but the two partly cancel.
- The `min` in the cap has a kink, so it is only piecewise differentiable.

## Future work

- Other layer-count dependencies from #1426 (convective entrainment, parcel start, cloud cover).
- Possibly a consistent treatment of the band between the condensation threshold and saturation.
