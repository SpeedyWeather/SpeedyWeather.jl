# BalancedZonalState: moist initial conditions without an initial precipitation burst

> Status: **completed**. `BalancedZonalState` and a height-dependent `ConstantRelativeHumidity` are the new default initial conditions of the `PrimitiveWetModel`. Destacked from the bulk Richardson fix (#1283), this PR contains only the initial conditions.

Date of initial draft: 2026-09-29

Base revision: c7d631c658feefc851eb7af5b93aa5754430faf5

## Originating prompt

> Can you have a look at the humidity initial conditions we use by default. PR1279 already
> implemented a relative humidity depending on the vertical which I would like to implement
> independently of that PR. Overall problem is that there's heavy precipitation within the first
> 1-3days particularly in the tropics due to some adjustment of temperature, pressure, ... I'm not
> sure about. Can you figure out why and find some initial conditions for humidity that do not
> cause this initial heavy rain dump?

## Revision log

- 2026-09-29: initial draft with diagnosis and spin-up experiments.
- 2026-09-29: follow-up prompt:
  > I like the idea to make the humidity profile both temperature (as now) and height dependent
  > (new) but can we also define a wet profile for the temperature initial conditions for the
  > PrimitiveWetModel that would be closer to the equilibrium profile? Technically this depends
  > on the radiative heating but maybe there's a good default profile you could implement?

  Plan changed from "lower the default RH to 0.5" to a new balanced temperature and wind
  initial condition fitted to the model's equilibrium, plus RH that depends on height. While
  analysing the slow spin-up that remained, found a bug in the bulk Richardson number (surface
  temperature not used), proposed as a separate PR.
- 2026-09-29: signed off:
  > Fix the Richardson bug first and open a PR with that alone, then stack another PR with the
  > BalancedZonalState, I like your plan

  The Richardson fix is #1283 (issue #1281). While fixing it, found that the vertical diffusion
  has no effect at all (#1282, not addressed here). Refitted the parameters to the equilibrium
  with the Richardson fix: T₀ 281 → 283 K, Γ 6 → 5.7 K/km, T_min 222 → 226 K, u₀ 52 → 56 m/s,
  H = 2 unchanged. `PressureOnOrography` and `StartFromRest` stay unchanged (see Known
  limitations).
- 2026-10-19: follow-up prompt:
  > Look at PR1284, can you destack this from PR1283? In that PR I just want the
  > BalancedZonalState and the new relative humidity initial conditions

  Rebased onto main without the #1283 commits. The parameters stay fitted to the equilibrium
  with #1283, so without it the initial state is slightly off the model's equilibrium until
  #1283 is merged.
- 2026-10-19: follow-up prompt:
  > Can you investigate the initial conditions of pressure too? JW starts with p=p_ref
  > everywhere but they also have no orography, so we currently just correct surface pressure
  > to be reduced to orographic height but there's still a wave bouncing around in the first
  > days, any suggestions?

  Diagnosed (see "Surface pressure on orography" below) and, signed off with "Just do 1 here
  please!", `BalancedZonalState` now sets the surface pressure itself and evaluates the
  temperature at p/p₀ = σpₛ/p₀. Digital filter initialization for the remaining wave is future work.
- 2026-10-19: follow-up prompt:
  > When outputting MeanSeaLevelPressureOutput() the initial conditions are always Inf,
  > presumably because the surface_temperature is still 0, wondering whether we could either fix
  > this in the output variable itself or whether it's actually a good idea to run the
  > parameterizations also at initialize! simulation but setting the tendencies back to 0

  > Can you just fix it in this PR please given it's also directly related to initial conditions?

  Fixed in the output variable (see Summary of changes, 5) rather than running the
  parameterizations at initialization, which also update accumulators, ocean/land fluxes and
  diagnostics that resetting the tendencies would not undo.
- 2026-10-20: follow-up prompt:
  > Can you look at PR1284 and add some better initial conditions for snow depth? Starting with 0
  > everywhere just doesn't cut it because you're on the wrong side of the snow-albedo feedback
  > [...] So could you have a good estimated guess at how much snow depth there should be with
  > upper-layer soil temperature as predictor and then use that as initial condition?

  Added as Summary of changes, 6. Soil temperature is now initialized before snow.
- 2026-10-23: follow-up prompt:
  > Can you expand the documentation on how they have been derived? [...] as long as it's just
  > about a few parameters that's fine, let's just state them and highlight that in order to
  > start with a balanced state for a climatology they would need to be adapted.

  Docs now derive the temperature from gradient-wind balance, state that the parameters were
  read off the zonal/time-mean equilibrium (no formal optimisation) with a table of how each
  was chosen, explain how to refit them, and contrast with DCMIP-2016 (Ullrich et al. 2016).
- 2026-10-23: follow-up prompt:
  > can you use the reference temperature please as this is the climate we ultimately would
  > like to have (roughly) but do feel free to use your lapse rate and write this into
  > EarthAtmosphere as the new default

  Removed the `T₀` and `lapse_rate` fields of `BalancedZonalState`, which now uses
  `atmosphere.reference_temperature` (288 K) and `lapse_rate(atmosphere)`. Default
  `moist_lapse_rate` of `EarthAtmosphere` 5 → 5.7 K/km (also affects JW and `PressureOnOrography`
  in wet models). Notation: z = -ln σ, radius R.
- 2026-10-28: follow-up prompt:
  > At the moment the initial conditions for the zonal velocity peak at sigma=0, can this be
  > lowered? [...] the jet should be around sigma=0.25 or so. This likely requires a temperature
  > distribution where the equator-pole gradient changes sign [...] As long as there's a knob to
  > get the balanced as is right now back

  U = u₀ tanh(z/H) f(z) with f = exp(-((z - z_jet)/D)²) above z_jet = -ln σ_jet, new options
  `σ_jet = 0.3`, `jet_decay = D = 1.5`; `σ_jet = 0` recovers the previous profile. Balance
  stays exact (the temperature formula holds for any U(z)); above the jet dU/dz < 0 reverses the
  meridional temperature gradient. Scanned (σ_jet, D) ∈ {(0.25, 1), (0.25, 1.5), (0.25, 2),
  (0.3, 1.5), (0.3, 2), (0.35, 2.5)} at T31L16: lowering the jet and keeping the reversal small
  compete; D = 1 gives a 54 K reversal, D = 2-2.5 peaks the jet at σ ≈ 0.15 with 27-35 K.
  Chose (0.3, 1.5): jet ≈ 35 m/s at σ ≈ 0.2, reversal up to 40 K at σ ≈ 0.09. In a 20-day
  T31L8 run the jet at 45° stays at σ ≈ 0.19 (33 m/s at day 20).

## Problem description

The default `PrimitiveWetModel` rains about 10 times its equilibrium rate in the first
12–18 hours. At T31L8 the global mean is 21 mm/day for the first 6 h and 18 mm/day for the next
6 h. In the tropics (|φ| < 20°) it is 39 and 34 mm/day. After one day the global mean settles at
about 2 mm/day. Almost all of this rain comes from convection. Large-scale condensation gives
about 5 mm/day, and only in the first 6 h.

## Background

The default initial conditions are `ZonalWind`, `PressureOnOrography`, `JablonowskiTemperature`
and `ConstantRelativeHumidity(relhumid_ref = 0.7)`. The temperature is the dry
Jablonowski–Williamson (2006) baroclinic-wave profile. It was designed for dry dynamical cores,
not moist physics. In the tropics it is too warm near the surface and too unstable:

| σ | T, 0–10°N/S, initial [K] | T, 0–10°N/S, model after 30 days [K] |
|---|---|---|
| 0.938 | 307.2 | 290.0 |
| 0.812 | 298.9 | 282.5 |
| 0.688 | 287.8 | 275.5 |
| 0.562 | 274.7 | 267.2 |
| 0.312 | 244.1 | 241.7 |

The lapse rate between σ = 0.94 and 0.56 is about 7–8 K/km. The moist adiabat is about 4 K/km
there, so the column is conditionally very unstable. With 70% relative humidity at 307 K, the
lowest layer holds q ≈ 26 g/kg. Tropical precipitable water starts at 80 mm, while the model's
own equilibrium is about 32 mm. The surface parcel therefore has large CAPE. Betts–Miller
convection removes it within its 4-hour time scale. Tropical precipitable water drops from 80 to
64 mm on day 1. The lowest layer cools by 10 K and the upper troposphere warms by up to 12 K.

So the burst is not an adjustment of pressure or winds. It is the convection scheme removing
the instability that JW temperature plus 70% RH creates at the start. The humidity in the
tropical boundary layer is the main control. Humidity in the upper troposphere matters little.

## Experiments

Global mean precipitation in mm/day in consecutive 12-hour intervals for 15 days. The scripts
are not committed. Equilibrium is about 2 mm/day.

| Initial humidity | T31L8 first 2 days | T63L8 first 2 days | T31L16 first 2 days |
|---|---|---|---|
| const RH 0.7 (current default) | **17.8** 3.3 1.9 1.8 | **18.3** 3.2 1.9 1.7 | 3.1 **7.6** 4.5 2.6 |
| const RH 0.6 | **8.9** 2.7 1.7 1.5 | **9.0** 2.6 1.7 1.5 | 0.5 1.2 2.8 3.2 |
| const RH 0.55 | 5.1 2.7 1.7 1.4 | 4.8 2.6 1.7 1.4 | 0.3 0.6 0.7 1.9 |
| **const RH 0.5** | 1.1 1.5 1.9 1.4 | 1.4 1.3 1.6 1.3 | 0.2 0.4 0.3 0.8 |
| RH 0.7, linear in σ to 0 at σ = 0.02 (#1279) | **6.1** 2.3 1.6 1.4 | (run crashed) | 0.5 2.5 2.3 1.8 |
| RH 0.6, linear in σ | 1.6 0.9 1.1 1.0 | 1.8 0.7 1.0 0.9 | 0.2 0.4 0.4 1.0 |
| RH 0.5, linear in σ | 0.3 0.2 0.1 0.1 (dry for 3 days) | 0.6 0.1 0.1 0.2 | 0.1 0.2 0.1 0.2 |
| DCMIP-2016 q, q₀ = 18 g/kg | 0.3 1.5 2.1 1.6 | 0.6 1.3 1.9 1.5 | 0.1 0.3 0.4 0.8 |
| DCMIP-2016 q, q₀ = 14 g/kg | 0.1 0.1 0.2 0.3 (dry for 3 days) | 0.2 0.1 0.2 0.3 | 0.1 0.1 0.1 0.2 |
| DCMIP-2016 q, q₀ = 22 g/kg | **8.7** 2.8 1.9 1.6 | **8.8** 2.7 1.9 1.5 | 0.6 1.1 1.7 2.8 |

The DCMIP-2016 moist baroclinic wave (Ullrich et al. 2016) prescribes specific humidity
q = q₀ exp(−(φ/φ_w)⁴) exp(−((σ−1)p₀/p_w)²), with φ_w = 40°, p_w = 340 hPa, and q = 10⁻¹²
above σ = 0.1.

All runs reach the same equilibrium within about 5 days.

- The burst appears once tropical boundary-layer RH is above about 0.55.
- The #1279 profile only reduces RH aloft. Its boundary-layer RH is 0.66, so the burst drops
  from 18 to 6 mm/day but does not go away.
- With JW temperature: const RH 0.5 and DCMIP q₀ = 18 g/kg both remove the burst at all three resolutions. They
  rain at 1–2 mm/day from the first 12 h at L8, and within about 2 days at L16.
- DCMIP dries everything poleward of about 50° (RH 0.01–0.4). The model settles at RH ≈ 0.7
  there. Const RH 0.5 has no such problem, and it can't supersaturate under any temperature
  initial condition.

## A near-equilibrium, balanced temperature and wind initial condition

### Equilibrium climate of the model

Zonal and time means over days 60–120, starting 2000-01-01, at T31L8, T63L8 and T31L16. The
three resolutions agree within about 2 K:

- **Global mean temperature**: follows T₀σ^(RΓ/g) with T₀ ≈ 280–282 K and Γ ≈ 6 K/km. The
  tropopause is near σ ≈ 0.2 at about 223 K, and the stratosphere is roughly isothermal at
  226–230 K.
- **Tropical lowest layer**: 287–292 K, with a lapse rate of about 5.5–6 K/km.
- **Meridional contrast**: equator to pole is about 45 K near the surface, about 40 K up to
  σ ≈ 0.45, and it shrinks above.
- **Wind**: at 45° it grows roughly linearly with −ln σ, from about 5 m/s near the surface to
  about 40 m/s at σ = 0.2.
- **Relative humidity**: about 0.6–0.75 throughout the troposphere, 0.25 at σ ≈ 0.16 and
  about 0.05 above σ ≈ 0.1.

JW differs mainly in its vertical wind profile. u₀cos^{3/2}(ηᵥ) puts most of the shear into the
lower troposphere, which gives a 78 K surface contrast and 307 K in the tropics.

### Generalising the Jablonowski–Williamson balance

In JW the temperature follows analytically from the wind through thermal-wind (gradient-wind)
balance. For any vertical wind profile U(σ), with x = −ln σ,

    u(φ, σ) = U(σ) sin²(2φ)
    T(φ, σ) = T̄(σ) + (1/R) dU/dx [2 A(φ) U + B(φ) aΩ]
    A(φ) = −2 sin⁶φ (cos²φ + 1/3) + 10/63
    B(φ) = 8/5 cos³φ (sin²φ + 2/3) − π/4

This reproduces JW's eq. (6) for U = u₀cos^{3/2}(ηᵥ). A and B have zero global mean, so T̄(σ) is
the global-mean profile and can be chosen freely without breaking balance. The vorticity is JW's
eq. (3) with U(σ) in place of u₀cos^{3/2}(ηᵥ). Divergence and the JW perturbation stay as they
are.

Prototype values, fitted to the equilibrium before #1283 and refitted afterwards (see Summary of changes):

- U(σ) = u₀ tanh(−ln σ / H), with u₀ = 52 m/s and H = 2. U is 0 at the surface, so uniform
  surface pressure stays balanced. It grows about linearly in −ln σ through the troposphere and
  levels off in the stratosphere.
- T̄(σ) = max(T₀ σ^(RΓ/g), T_min), with T₀ = 281 K, Γ = 6 K/km and T_min = 222 K.

The initial state (T31L8, 0–10° and 65–90° bands) then matches the equilibrium closely:

| | tropical T at σ = 0.94 | polar T at σ = 0.94 | tropical T at σ = 0.56 | tropical precipitable water |
|---|---|---|---|---|
| JW (current) | 307.2 K | 229 K | 274.7 K | 80 mm |
| new | 289.7 K | 247.5 K | 265.3 K | 32 mm |
| model equilibrium | 287.0 K | 240 K | 262.9 K | 32 mm |

### Humidity: RH depends on both temperature and height

RH = relhumid_ref below σ_moist, decreasing linearly to 0 at σ_dry. With relhumid_ref = 0.7,
σ_moist = 0.25 and σ_dry = 0.1 this matches the equilibrium: 0.26 at σ = 0.156 (equilibrium
0.25) and 0.56 at σ = 0.219 (equilibrium 0.55). Manabe and Wetherald's profile (#1279) is the
special case σ_moist = 1, σ_dry = 0.02. It is too dry in the mid-troposphere: RH 0.7 gives 2–3
days of almost no rain.

### Spin-up with the new initial conditions

Global mean precipitation in mm/day for each 12-hour interval of the first 4 days.

| | T31L8 | T63L8 | T31L16 |
|---|---|---|---|
| JW + RH 0.7 (current) | **17.8** 3.3 1.9 1.8 1.9 1.8 1.7 1.8 | **18.3** 3.2 1.9 1.7 1.8 1.6 1.6 1.8 | 3.1 **7.6** 4.5 2.6 1.7 1.9 1.6 2.1 |
| new T/wind + RH 0.7 tapered | 1.4 0.8 0.9 0.8 0.8 0.9 0.9 1.0 | 1.5 0.7 0.8 0.7 0.8 0.8 0.8 1.0 | 0.5 0.3 0.5 0.8 1.2 1.4 1.5 1.7 |
| same, with the Ri fix (below) | 1.3 1.2 1.4 1.5 1.6 1.8 1.8 2.1 | – | 0.5 0.7 1.3 1.4 1.4 1.8 1.8 2.1 |

Equilibrium is about 2.0–2.3 mm/day, or about 2.6 with the Ri fix.

The burst is gone at every resolution. Instead precipitation ramps up over about 4–5 days.
Random temperature noise of 0.1–1 K to seed eddies makes no difference.

### Why the ramp-up is slow: bulk Richardson number ignores the surface temperature

The ramp-up is limited by surface evaporation. It is about 0.2 mm/day over the tropical ocean
at first, although qsat(SST) = 23.5 g/kg is far above q_air = 8.8 g/kg. The reason is
`bulk_richardson_surface` (`surface_fluxes/boundary_layer.jl`). It sets Θ₀ = cₚT_v from the
**lowest-layer air temperature** and Θ₁ = Θ₀ + gz, so Θ₁ − Θ₀ = gz always. Ri = (gz)²/(cₚT_v V²)
is therefore always positive (stable) and does not depend on SST or land temperature. With the
lowest level at about 500 m (L8), Ri > Ri_c = 10 for any wind below about 3 m/s. The drag then
drops to `drag_min = 1e-5`, about 30 times too small, even over a tropical ocean 10 K warmer
than the air. Frierson (2006) eq. 15 uses the surface virtual potential temperature.

A scratch override that uses the skin temperature (SST and top soil layer, weighted by land
fraction) raises tropical ocean drag to 9×10⁻⁴ and evaporation to 2.3–2.9 mm/day from the
first hours. It also changes the equilibrium climate substantially. At T31L8 the tropical
column warms by 4–8 K and holds up to 40% more water, and the global mean warms by 3–6 K. At
T31L16 the warming is 2–5 K. This partly addresses the "tropics too cold" bias (#856). The fix
belongs in its own PR. Because it shifts the equilibrium, the defaults T₀, Γ, u₀ and H should
be refitted once it lands. At T31L8 with the fix, T₀ rises by about 4 K.

A second observation that is not addressed here: the land humidity flux
ρC_DV(α·qsat(T_soil) − q_air) turns negative whenever q_air > α·qsat. Dry soil then takes
moisture out of unsaturated air (−1.8 mm/day over land in the first day, −3.5 with the Ri fix).
SPEEDY had max(…, 0) there. It was removed to allow dew, but that also allows deposition onto
dry soil.

### Surface pressure on orography

The default initial state launches a gravity wave over orography: the RMS of the hourly
surface pressure change is 3 hPa/h in the first hours vs about 0.1 hPa/h in equilibrium.
Two inconsistencies cause it. `PressureOnOrography` uses the atmosphere's reference
temperature and lapse rate (288 K, 5 K/km), not T₀ and Γ of `BalancedZonalState`. And the
temperature is a function of σ, so air on a 3 km mountain starts at the sea-level 283 K
instead of about 266 K, out of hydrostatic balance with the reduced surface pressure.

RMS of the hourly surface pressure change [hPa/h], T31L8 `PrimitiveWetModel`:

| Variant | hours 1–6 | day 2 | day 5 |
|---|---|---|---|
| A: before | 3.0 2.9 2.0 1.2 0.9 0.9 | 0.16 | 0.12 |
| B: pₛ with T₀, Γ of `BalancedZonalState` | 2.7 2.6 1.8 1.1 0.8 0.9 | 0.16 | 0.12 |
| C: B and temperature at p = σpₛ (**implemented**) | 1.9 1.8 1.2 0.7 0.5 0.5 | 0.13 | 0.10 |
| G: C and the wind at p = σpₛ too | 1.9 1.8 1.3 0.7 0.4 0.4 | 0.11 | 0.09 |
| D: at rest, T̄(p) only | 0.07 0.06 … | 0.08 | 0.06 |
| F: at rest, T̄(σ) only | 1.5 1.3 1.0 0.8 0.7 0.7 | 0.14 | 0.09 |
| E: no orography | 0.07 0.09 … | 0.10 | 0.09 |

At rest the consistent state (D) has no wave at all. The remainder in C and G is the wind
blowing into orography with zero initial divergence, which no analytic profile removes.

## Summary of changes

1. New `BalancedZonalState` initial condition (`dynamics/initial_conditions.jl`). It sets
   vorticity, divergence and temperature from one set of parameters, so wind and temperature
   can't be made inconsistent:
   - `u₀ = 56` m/s and `H = 2` for U(σ) = u₀ tanh(−ln σ / H).
   - `T₀ = 283` K, `lapse_rate = 5.7e-3` K/m and `Tmin = 226` K for T̄(σ).
   - The JW perturbation (`perturb_lat`, `perturb_lon`, `perturb_uₚ`, `perturb_radius`).

   The vorticity is a `BalancedZonalVorticity` functor used with `set!`, and the divergence
   reuses `JablonowskiDivergence`. The temperature comes from a kernel with per-layer T̄, U and
   dU/dx/R precomputed on the CPU, following `JablonowskiTemperature`. The JW vorticity
   perturbation is factored out into `jablonowski_vorticity_perturbation` and shared by
   `JablonowskiVorticity` and `BalancedZonalVorticity`.
2. `ConstantRelativeHumidity` gets `σ_moist = 0.25` and `σ_dry = 0.1`. RH is `relhumid_ref`
   below σ_moist and decreases linearly to 0 at σ_dry. `σ_moist = σ_dry = 0` gives vertically
   constant RH, and `σ_moist = 1, σ_dry = 0.02` gives Manabe–Wetherald (#1279).
3. `InitialConditions(spectral_grid, PrimitiveWet)` becomes
   `(; zonal = BalancedZonalState, humid = ConstantRelativeHumidity)`.
   PrimitiveDry keeps JW.
4. `BalancedZonalState` also sets the surface pressure, hydrostatic with T₀ - Γz on orography,
   and evaluates T̄, U and dU/dx in its temperature at x = -ln(σpₛ/p₀) instead of -ln(σ).
   The wind stays a function of σ (zero at the surface), as G gave no improvement.
5. `MeanSeaLevelPressureOutput` computes the surface air temperature on the fly from the lowest
   layer (dry adiabatic, as `SurfaceCondition`) instead of reading
   `vars.parameterizations.surface_air_temperature`, which is 0 for the initial conditions
   and gave Inf. It also no longer overwrites that variable when `dynamics_only`. Tested in
   `test/output/netcdf_output.jl` on the first output time step.
6. `SnowModel` initial snow depth (`parameterizations/land/snow.jl`) estimated from the top-layer
   soil temperature T₁ instead of 0: `clamp(initial_snow_depth_rate * (melting_threshold - T₁), 0,
   initial_snow_depth_max)`, zero over ocean, with 5 mm/K and a 0.2 m cap (water equivalent).
   With the 0.05 m albedo `snow_depth_scale`, land 10 K below the 275 K melting threshold starts
   fully snow covered, and ~250 K land (Siberia, Canada, Greenland in winter) starts with ~0.12 m,
   close to observed winter maxima of 0.1–0.2 m. The land model now initializes soil temperature
   before snow.

## Testing and verification

`test/dynamics/initial_conditions.jl`:
- `BalancedZonalState` is the default. Its area-weighted global mean temperature matches
  T̄(σ) (rtol 1e-3) on every layer. Its vorticity matches −4U/a sinφ cosφ (2 − 5sin²φ) at 30°S
  (rtol 2%). The tropics are more than 30 K warmer than the poles near the surface.
- `ConstantRelativeHumidity` profile: mean RH is 0.7 below σ_moist, humidity is 0 above σ_dry,
  and RH is linear in between.
- No initial burst: the default PrimitiveWet at T31L8 rains less than 3 mm in the first day,
  area-weighted global mean. JW gives about 10 mm.

Also passing locally: `output/netcdf_output.jl`, `dynamics/{dispatch,time_stepping,set,simulation_constructor,copy_variables}.jl`,
`output/boundary_layer_output.jl`, `parameterizations/{all_parametrizations,boundary_layer}.jl`.

Initial state at T31L8, tropics 0–10°, compared with the equilibrium (days 60–120) with #1283:

| σ | T initial | T equilibrium | q initial [g/kg] | q equilibrium [g/kg] |
|---|---|---|---|---|
| 0.938 | 292.8 | 291.3 | 10.8 | 10.0 |
| 0.438 | 257.9 | 259.0 | 1.9 | 2.2 |

Spin-up, global mean precipitation in mm/day per 12 hours, JW → BalancedZonalState, both with
#1283:

| | first 3 days, JW | first 3 days, BalancedZonalState |
|---|---|---|
| T31L8 | **17.3** 3.3 1.9 1.8 1.7 1.9 | 1.8 1.4 1.4 1.6 1.5 1.8 |
| T63L8 | **17.8** 3.2 1.9 1.7 1.7 1.8 | 1.9 1.3 1.4 1.4 1.5 1.7 |
| T31L16 | 4.0 **8.2** 4.2 2.4 1.7 2.1 | 0.6 0.8 1.1 1.6 1.3 1.8 |

Equilibrium is about 2.3–2.6 mm/day and is reached after about 4 days.

Balance, with `dynamics_only = true`, no orography, no perturbation and constant surface
pressure, T31L8, after 1 day (5 days in brackets):

| | rms divergence | max \|v\| |
|---|---|---|
| JW | 2.0e-7 (1.7e-8) 1/s | 0.98 (0.04) m/s |
| BalancedZonalState | 3.4e-8 (1.6e-9) 1/s | 0.14 (0.006) m/s |

## Documentation changes

New section "Default initial conditions for the PrimitiveWetModel" in
`docs/src/initial_conditions.md`. It covers the derivation, the humidity profile, how to change
parameters and how to go back to JW.

## Known limitations

- Fitted to the equilibrium after #1283 at T31L8, T63L8 and T31L16. The fit depends on
  radiation, the surface fluxes and the season, and it is hemispherically symmetric while the
  default start date is 1 January. The parameters will need refitting when the physics changes
  substantially, e.g. once #1282 switches on the vertical diffusion.
- The jet sits at 45° and the zonal wind is westerly everywhere. There are no trades or surface
  easterlies. At the surface U = 0, which keeps the constant surface pressure balanced.
- The polar stratosphere is colder than in equilibrium, about 205 K vs 218 K at σ ≈ 0.2–0.3.
- Combining `BalancedZonalState` with `PressureOnOrography` overwrites the consistent surface
  pressure. `PressureOnOrography` itself still uses the atmosphere's reference values.
- The wind over orography still launches a gravity wave, about 60% of the original amplitude,
  mostly decayed after a day.
- `StartFromRest` still uses JW temperature, with the new humidity profile.

## Future work

- #1282: vertical diffusion, then refit.
- Digital filter initialization (Lynch and Huang, 1992) to remove the remaining gravity wave,
  for any initial condition.
- Surface winds with trades and midlatitude westerlies, balanced by a surface pressure pattern.
- Land humidity flux into dry soil (see the #1283 plan).
