# Bulk Richardson number with the surface temperature

> Status: **completed**. The surface bulk Richardson number now uses the surface (skin) temperature in both `BulkRichardsonDrag` and `BulkRichardsonDiffusion`.

Date of initial draft: 2026-09-29

Base revision: c7d631c658feefc851eb7af5b93aa5754430faf5

## Originating prompt

> Can you have a look at the humidity initial conditions we use by default. [...] Overall
> problem is that there's heavy precipitation within the first 1-3days particularly in the
> tropics [...] Can you figure out why and find some initial conditions for humidity that do not
> cause this initial heavy rain dump?

> [...] can we also define a wet profile for the temperature initial conditions for the
> PrimitiveWetModel that would be closer to the equilibrium profile?

> Fix the Richardson bug first and open a PR with that alone, then stack another PR with the
> BalancedZonalState, I like your plan

## Revision log

- 2026-09-29: found while analysing the precipitation spin-up for new initial conditions (see
  `humidity-initial-conditions.md`, a follow-up PR). With a near-equilibrium initial state,
  evaporation over the tropical ocean stayed at about 0.2 mm/day for days, although
  qsat(SST) = 23.5 g/kg was far above q_air = 8.8 g/kg. The boundary layer drag was at its
  minimum, 1e-5.
- 2026-09-29: also fixed the surface Richardson number in `BulkRichardsonDiffusion`. While
  testing, found that the vertical diffusion has no effect at all (neighbour indices
  `k₋ = max(k, 1)`, `k₊ = min(k, nlayers)` give zero gradients since 80926431, March 2024). Left
  for a separate issue (#1282) and PR, because switching it on also needs the conversion of K from m²/s
  to σ-space.

## Problem description

`bulk_richardson_surface` computed

    Θ₀ = cₚ T_v(T_N)          # "virtual dry static energy at surface", but with the air temperature T_N
    Θ₁ = Θ₀ + gz              # at the lowermost layer N, height z
    Ri = gz (Θ₁ − Θ₀) / (Θ₀ V²) = (gz)² / (cₚ T_v V²)

so Ri > 0 (stable) always, whatever the surface temperature. For the lowermost layer at about
500 m (8 layers) Ri exceeds Ri_c = 10 for any wind below about 3 m/s. The drag then drops to
`drag_min = 1e-5`, roughly 30–100 times too small. This shuts down the surface sensible heat,
humidity and momentum fluxes in calm conditions, even over a tropical ocean 10 K warmer than the
air. The surface Ri in `BulkRichardsonDiffusion` was built the same way.

## Background

Frierson et al. (2006), eq. (15): Ri_a = g z_a (θ_v(z_a) − θ_v(0)) / (θ_v(0) |v(z_a)|²). The
surface value θ_v(0) uses the surface temperature. In dry static energy form this is
Θ₀ = cₚ T_v(Tₛ) at z = 0 and Θ₁ = cₚ T_v(T_N) + gz. Tₛ comes from the ocean (SST) or land
(uppermost soil layer), as already used by the surface heat fluxes.

## Summary of changes

- New `surface_skin_temperature(ij, vars, land_sea_mask, time_stepping, component, T_air)`:
  SST over ocean and the uppermost soil layer temperature over land, weighted by the land
  fraction. Where a surface type is undefined (NaN) the other one is used. If neither is
  available it falls back to the air temperature, which recovers the previous behaviour.
- `bulk_richardson_surface` (drag): Θ₀ = cₚ T_v(Tₛ), Θ₁ = cₚ T_v(T_N) + gz. `land_sea_mask` is
  passed through.
- `bulk_richardson!` (diffusion), surface layer: Θ₀ = cₚ T_v(Tₛ) + Φₛ, Θ₁ = cₚ T_v(T_N) + Φ_N,
  Ri_N = (Φ_N − Φₛ)(Θ₁ − Θ₀)/(cₚ T_v(Tₛ) V²). The height is now measured above the surface
  rather than as absolute geopotential. The layers above are unchanged.

## Testing and verification

New `test/parameterizations/boundary_layer.jl`:
- `surface_skin_temperature` returns SST on ocean points, soil temperature on land points and
  the weighted value on coastal points.
- For calm wind at an ocean point, a surface 10 K warmer than the air gives Ri < 0 and the
  neutral maximum drag (κ/ln(z/z₀))². A surface 10 K colder gives Ri > 0 and `drag_min`.
- The diffusion finds a boundary layer over the warm surface and none over the cold one.

Effect on the default `PrimitiveWetModel` (start 2000-01-01, zonal and time mean over days
60–120, main → this PR):

| | global T, lowest layer [K] | T 0–10°, lowest layer [K] | T 0–10°, σ ≈ 0.45 [K] | q 0–10°, lowest layer [g/kg] | max zonal-mean u [m/s] |
|---|---|---|---|---|---|
| T31L8 | 277.4 → 281.0 | 287.0 → 291.3 | 250.9 → 259.0 | 7.2 → 10.0 | 50.4 → 57.0 |
| T63L8 | 277.9 → 280.7 | 287.9 → 291.2 | 252.4 → 258.5 | 7.6 → 9.8 | 45.3 → 51.6 |
| T31L16 | 282.2 → 283.6 | 292.1 → 294.1 | 257.1 → 261.3 | 9.9 → 11.6 | 51.5 → 52.8 |

Global mean precipitation over days 10–15: T31L8 2.24 → 2.58, T63L8 1.95 → 2.49 and T31L16
2.13 → 2.31 mm/day. The model's tropics warm, which reduces the known cold bias (#856). The
effect is smaller at 16 layers, where the lowermost layer is closer to the surface. There z is
smaller, so the old Ri was less often above Ri_c.

## Documentation changes

None. The docstrings of `bulk_richardson_surface` and `surface_skin_temperature` describe the
formulation.

## Known limitations

- Tₛ is the SST or the top soil layer temperature, not a proper skin temperature. The surface
  heat fluxes make the same assumption.
- Sea ice is not considered separately. The SST is used under sea ice, as in the heat flux.

## Future work

- The vertical diffusion is a no-op (neighbour indices), and its K would need the conversion to
  σ-space before it can be switched on, #1282.
- The land humidity flux ρC_DV(α·qsat(T_soil) − q_air) is negative over dry soil whenever
  q_air > α·qsat, so dry soil takes up moisture from unsaturated air.
