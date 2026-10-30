# Bulk Richardson number relative to the land and ocean surface

> Status: **completed**. `BulkRichardsonDrag` computes separate land and ocean Richardson numbers from the
> soil and sea surface temperatures; vertical diffusion and the critical Richardson number are unchanged.

Date of initial draft: 2026-09-29 (as the plan for #1283, adapted to #1363 on 2026-10-28)

Base revision: e292baf012d70813f9eab23081661bf1e168c3fb

## Originating prompt

> Can you have a look at the humidity initial conditions we use by default. [...] Overall
> problem is that there's heavy precipitation within the first 1-3days particularly in the
> tropics [...] Can you figure out why and find some initial conditions for humidity that do not
> cause this initial heavy rain dump?

> Fix the Richardson bug first and open a PR with that alone [...]

(#1283, which led to #1281.) #1363 was opened independently:

> We were hitting the `drag_min` quite a lot after I fixed a bug in the surface roughness in #1349
> and was trying to work out why - this seems to fix it, and improve the surface fluxes as a result.

## Revision log

- 2026-09-29: found while analysing the precipitation spin-up for new initial conditions. With a
  near-equilibrium initial state, evaporation over the tropical ocean stayed at about 0.2 mm/day
  for days, although qsat(SST) = 23.5 g/kg was far above q_air = 8.8 g/kg. The boundary layer drag
  was at its minimum, 1e-5.
- 2026-09-30: considered the surface *air* temperature for Θ₀ and rejected it. It is T_N σ^(−κ),
  the lowermost layer's potential temperature, so Θ₁ − Θ₀ ≈ 0 and Ri ≈ 0 everywhere by construction.
- 2026-09-30: #1283 was extended to also fix the vertical diffusion (#1282), with the critical
  Richardson number lowered to 1 (Frierson 2006). the 24-layer configuration then blew up over Tibet.
- 2026-10-28: #1363 fixes the same bug independently, with separate land and ocean Richardson
  numbers and a saturated surface humidity. It replaces the Richardson part of #1283. Decisions in
  review:
  - Ri_c stays at 10. `(1 − Ri/Ri_c)²` only sets the transition; with Ri_c = 1 the drag is already
    at 25% for Ri = 0.5, where shear and buoyancy are comparable. Ri_c = 10 shuts the drag off only
    in strongly stable boundary layers.
  - Vertical diffusion is left alone (it is `nothing` by default).
  - No NaN guard: land and ocean must be initialised without NaNs; mixing them is the job of the
    land-sea mask.

## Problem description

`bulk_richardson_surface` computed

    Θ₀ = cₚ T_v(T_N)          # "virtual dry static energy at surface", but with the air temperature T_N
    Θ₁ = Θ₀ + gz              # at the lowermost layer N, height z
    Ri = gz (Θ₁ − Θ₀) / (Θ₀ V²) = (gz)² / (cₚ T_v V²)

so Ri > 0 (stable) always, whatever the surface temperature. For the lowermost layer at about
500 m (8 layers) Ri exceeds Ri_c = 10 for any wind below about 3 m/s. The drag then drops to
`drag_min = 1e-5`, roughly 30–100 times too small. This shuts down the surface sensible heat,
humidity and momentum fluxes in calm conditions, even over a tropical ocean 10 K warmer than the air.

## Background

Frierson et al. (2006), eq. (15): Ri_a = g z_a (θ_v(z_a) − θ_v(0)) / (θ_v(0) |v(z_a)|²). The
surface value θ_v(0) uses the surface temperature. In dry static energy form this is
Θ₀ = cₚ T_v(Tₛ, qₛ) at z = 0 and Θ₁ = cₚ T_v(T_N, q_N) + gz. Since land and ocean drag coefficients
are already computed separately (with separate roughness lengths), Ri is computed per surface too.

## Summary of changes

1. `bulk_richardson_surface` returns `(Ri_land, Ri_ocean)`:
   - Θ₁ = cₚ T_v(T_N, q_N) + gz from the lowermost layer.
   - Ocean: Tₛ = SST, qₛ = q⋆(SST, pₛ), as in `SurfaceOceanHumidityFlux`.
   - Land: Tₛ = uppermost soil layer temperature, qₛ = α q⋆(Tₛ, pₛ) with the soil moisture
     availability α, as in `SurfaceLandHumidityFlux`.
   - Dry models: qₛ = q_N = 0.
2. `boundary_layer_drag!` applies the stability factor `(1 − clamp(Ri, 0, Ri_c)/Ri_c)²` per surface
   to the land and ocean momentum, heat and humidity drag coefficients.
3. Documentation in `docs/src/surface_fluxes.md` describes the surface Θ₀ over land and ocean.

Not changed: Ri_c = 10, `BulkRichardsonDiffusion` (still uses the old surface Θ₀), version
(`SpeedyWeather` is already at `0.24.0-DEV`).

## Testing and verification

`test/parameterizations/drag.jl`, for `PrimitiveWetModel` and `PrimitiveDryModel`:
- Land and ocean drag coefficients are finite and ≥ `drag_min` after a few steps.
- Surface 20 K warmer than the air (Ri < 0): neutral transfer coefficient.
- Surface 10 K colder at 10 m/s (0 < Ri < Ri_c): between `drag_min` and neutral.
- Surface 20 K colder at 0.5 m/s (Ri ≥ Ri_c): `drag_min`.

Default `PrimitiveWetModel` with `truncation = 31`, 10 days: 8 and 24 layers stable; heat drag at `drag_min` at ~0.03%
(L8) and 0% (L24) of grid points.

## Documentation changes

`docs/src/surface_fluxes.md`, section on the drag coefficient.

## Known limitations

- Tₛ is the SST or the top soil layer temperature, not a proper skin temperature. The surface heat
  fluxes make the same assumption.
- Sea ice is not considered separately.
- `BulkRichardsonDiffusion` is inconsistent with the drag: its surface Ri still uses the air
  temperature, and its neighbour indices give zero gradients (#1282).

## Future work

- Fix the vertical diffusion (#1282): neighbour indices, K in σ coordinates, implicit solve, and
  the surface Ri from `bulk_richardson_surface`. A first attempt is on the #1283 branch.
- With a working diffusion, explicit land–atmosphere coupling can become unstable for thin
  lowermost layers (seen with 24 layers over Tibet in #1283).
