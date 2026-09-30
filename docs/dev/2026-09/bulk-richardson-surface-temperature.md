# Bulk Richardson number with the surface temperature and a working vertical diffusion

> Status: **in progress**. Richardson fix, working implicit vertical diffusion and Ri_c = 1 implemented; stable up to T85L16, but T31L24 blows up over Tibet through an explicit land–atmosphere coupling instability (limiter proposed, awaiting decision).

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

- 2026-09-30: review of the drag behaviour:
  > Can we talk about PR1283, you replaced the temperature used in the bulk Richardson from
  > lowermost layer to surface ocean/land temperatures but not the surface air temperatures. Now
  > the boundary layer drag is always maxed out [...] maybe variables.parameterization.surface_temperature
  > should be used instead because that's actually at z=0 but is still air [...]

  Kept the ocean/land temperature. The surface air temperature is T_N σ^(−κ), the lowermost
  layer's potential temperature, so Θ₁ − Θ₀ ≈ 0 and Ri ≈ 0 (always neutral) by construction.
  The drag is near its maximum at ~90% of points because the surface is on average 3.3 K
  (ocean) and 5.5 K (land) warmer than the lowermost layer's potential temperature.
- 2026-09-30:
  > Can you turn PR1283 into fixing both drag and vertical diffusion as suggested?

  PR extended to fix #1282 (see Summary of changes, items 4–8). Explicit diffusion with K in
  σ-space blew up at T85L16, so it's now implicit. Solving over Δt with Leapfrog applying the
  tendency over 2Δt overshot and blew up at T31L16, so it's now solved over 2Δt. The stability
  factor used Ri at the boundary-layer top (≈ Ri_c by construction), which made K ≈ 1 m²/s, so
  it now uses the surface Ri_N (Frierson eq. 20). T31L24 still blows up; diagnosed as explicit
  land–atmosphere coupling (see Known limitations).

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

Surface bulk Richardson number (#1281):

1. New `surface_skin_temperature(ij, vars, land_sea_mask, time_stepping, component, T_air)`:
   SST over ocean and the uppermost soil layer temperature over land, weighted by the land
   fraction. Where a surface type is undefined (NaN) the other one is used. If neither is
   available it falls back to the air temperature, which recovers the previous behaviour.
2. `bulk_richardson_surface` (drag): Θ₀ = cₚ T_v(Tₛ), Θ₁ = cₚ T_v(T_N) + gz.
3. `bulk_richardson!` (diffusion), surface layer: Θ₀ = cₚ T_v(Tₛ) + Φₛ, Θ₁ = cₚ T_v(T_N) + Φ_N,
   Ri_N = (Φ_N − Φₛ)(Θ₁ − Θ₀)/(cₚ T_v(Tₛ) V²).

Vertical diffusion (#1282):

4. Neighbour indices in `_vertical_diffusion!`: `k₋ = max(k - 1, 1)`, `k₊ = min(k + 1, nlayers)`
   (were both `k`, so all gradients were zero).
5. K [m²/s] converted to σ coordinates: K̃ = K (gσ/(RT))² [1/s], from ∂z = −(gσ/(RT)) ∂σ.
6. Implicit (backward Euler) solve per column with the Thomas algorithm, over
   `implicit_vertical_diffusion_time_step` = 2Δt for Leapfrog (tendencies are applied over 2Δt
   from the previous step), Δt otherwise. No flux through the boundary-layer top. Dry static
   energy is diffused and its tendency converted with 1/cₚ, rather than dividing K by cₚ. Two
   new scratch arrays `vertical_diffusion_c`, `vertical_diffusion_d`.
7. Stability factor, Frierson eq. (20), uses the surface Ri_N instead of Ri at the
   boundary-layer top.

Both:

8. Critical Richardson number 10 → 1 in `BulkRichardsonDrag` and `BulkRichardsonDiffusion`
   (Frierson 2006, still to be checked against the paper). 10 was probably tuned to live with
   the always-stable bug.
9. `BulkRichardsonDrag(SG::SpectralGrid; kwargs...)`, previously `, kwargs...`, which rejected
   keyword arguments.

## Testing and verification

`test/parameterizations/boundary_layer.jl`:
- `surface_skin_temperature` on ocean, land and coastal points.
- Calm wind at an ocean point: a surface 10 K warmer than the air gives Ri < 0 and the neutral
  maximum drag (κ/ln(z/z₀))². A surface 10 K colder gives `drag_min`. The diffusion finds a
  boundary layer over the warm surface and none over the cold one.
- The vertical diffusion changes u, conserves the column integral Σ tendency·Δσ, and reduces
  the vertical variance within the boundary layer.

Diffusion before and after (T31L8, 3 days): identical results with and without
`vertical_diffusion` on main. With this PR, K ≈ 50–100 m²/s near the surface in a typical
column, and tendencies of a few K/day.

Drag with this PR (T31L8, day 20): ≥ 99% of the maximum at 91% of ocean and 83% of land points,
at `drag_min` at 1% (ocean) and 8% (land). Surface minus lowermost-layer potential temperature
is +3.5 K (ocean) and +5.1 K (land) on average; the diffusion doesn't change this.

Stability, default `PrimitiveWetModel`, 20 days:

| | main | this PR |
|---|---|---|
| T31L4, T31L8, T63L8, T31L16, T85L16 | stable | stable |
| T31L24 | stable | NaN after ~12 h (over Tibet) |
| T31L24, 20-min time step | – | stable |
| T31L24, no orography | – | stable |
| T31L32 | NaN | – |

The climate impact numbers from the first version (Richardson fix only) are outdated and will
be redone once the stability question is settled.

## Documentation changes

None. The docstrings of `bulk_richardson_surface` and `surface_skin_temperature` describe the
formulation.

## Known limitations

- **T31L24 blows up over Tibet.** The soil top layer swings between 218 and 342 K from hour to
  hour and the drag flips between minimum and maximum. Momentum diffusion brings jet-level winds
  (~40 m/s) to the 2.5 km plateau surface. The thin lowermost layer (z ≈ 190 m) raises the
  maximum drag to ~4.7e-3. Together ρcₚC_DV ≈ 150 W/m²/K, an exchange time scale of about
  20–25 min for both the top soil layer and the lowermost air layer. That is shorter than the
  2 × 40 min Leapfrog step, so the explicit surface heat flux overshoots. Proposed: limit
  ρC_DV so that the exchange time scale can't fall below the time step. Awaiting decision.
- Tₛ is the SST or the top soil layer temperature, not a proper skin temperature. The surface
  heat fluxes make the same assumption.
- Sea ice is not considered separately.
- The critical Richardson number of 1 is from memory of Frierson (2006) and not yet checked.

## Future work

- Implicit coupling of the surface fluxes with land and the lowermost layer, as an
  alternative to a limiter.
- The land humidity flux ρC_DV(α·qsat(T_soil) − q_air) is negative over dry soil whenever
  q_air > α·qsat, so dry soil takes up moisture from unsaturated air.
- Refit `BalancedZonalState` (#1284) to the new equilibrium.
