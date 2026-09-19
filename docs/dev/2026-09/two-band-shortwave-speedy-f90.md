# Two-band shortwave radiation, seasonal ozone and cloud fixes following speedy.f90

> Status: **in progress**. Implemented in PR #1262, awaiting review; `TwoBandShortwave`
> is the new default for `PrimitiveWetModel`.

Date of initial draft: 2026-09-19

Base revision: `35d764b6` (`main`)

## Originating prompt

> Could you clone Sam Hatfield's speedy.f90 repository so we can make some comparisons
> between their shortwave radiation scheme and our shortwave radiation scheme? Can you create
> a new PR that lists the algorithmic differences, I believe we treat the ozone differently but
> I think speedy.f90 uses a UV band that bypasses the cloud reflection. Can you also check how
> they implement clouds. The PR should make our shortwave radiation scheme at least as
> sophisticated as theirs to start with.

## Revision log

- **2026-09-19, initial draft.** Written retroactively after the implementation was pushed
  as PR #1262 (the plan should have preceded it). Records the comparison with speedy.f90, the
  bugs found in `DiagnosticClouds`, and the global-mean energy budget used for verification.
- **2026-09-19.** Correction to the originating prompt: speedy.f90 has no UV band. The band
  that bypasses cloud reflection is a near-infrared band (5% of incoming flux) absorbed only by
  water vapor. Ozone absorption is subtracted from the visible band in the stratosphere.
- **2026-09-19.** Cloud albedo reset from 0.6 to speedy.f90's 0.43. Once the cloud cover
  bugs were fixed, 0.6 gave a planetary albedo of 0.48. 0.6 had been compensating for the
  too-low cloud cover.
- **2026-09-19.** Ozone seasonal cycle synchronized with the solar zenith angle
  (prompt: "The ozone seasonal cycle seems to hardcode Earth. The zenith angle calculation we
  have is now completely flexible, can the seasonal cycles be synced somehow?"). The SPEEDY
  factor `max(0, cos α)`, with α the year angle from the northern winter solstice (hardcoded as
  Jan 1 + 10 days on a 365-day year), is replaced by `max(0, -δ/|δₘₐₓ|)` with the solar
  declination δ from `model.solar_zenith` and δₘₐₓ the planet's `axial_tilt`. This is identical
  for Earth's sinusoidal declination, but it now follows `length_of_year`, `equinox`,
  `axial_tilt` and `seasonal_cycle`. Removed the `solstice_offset` field. New helpers
  `year_angle(NF, zenith, orbit_time)` and `solar_declination(NF, zenith, orbit_time)` are
  shared with the zenith calculation.
- **2026-09-19.** Cloud scheme compared with speedy.f90 by zonal means of every cloud term
  (prompt: "Can you also compare the cloud scheme similarly?"). This found two more
  algorithmic differences, now fixed: 5. the q > q_min threshold was also applied to layer N−1
  (speedy.f90 exempts it), so there were no humidity clouds in cold polar air; 6. snowfall was
  missing from the precipitation term. The remaining differences come from the model climate,
  see "Cloud comparison" below.
- **2026-09-19.** Generalized from speedy.f90's fixed 8 layers to any number of layers
  (prompt: "Speedy.f90 is hardcoded to 8 layers, the lowermost is the planetary boundary layer,
  the top two are stratosphere. We want to generalise this to any amount of layers. [...] I think
  it is reasonable to assume certain processes never go to the surface [...] Have a check that we
  indeed generalise in a reasonable way"). The remaining layer-index assumptions (`nlayers - 1`
  as cloud base, `nlayers` as boundary layer) are replaced by two σ-boundaries in
  `DiagnosticClouds`: `σ_tropopause = 0.14` and `σ_boundary_layer = 0.9`, speedy.f90's half
  levels between its stratosphere and troposphere and between its free troposphere and boundary
  layer. `free_troposphere_layers` returns `layer_top:cloud_base`, and `clouds!` now also returns
  `cloud_base`. Clouds (humidity search, cloud top from convection/condensation, cloud absorption
  in both transmissivities) stay within the free troposphere and never reach the surface layer.
  Stratocumulus stability is the Δs/ΔΦ between the cloud base and the surface layer, and
  stratocumulus reflects at the top of the boundary layer. Nothing changes for 8 equally spaced
  layers, except that layer 1 (σ = 0.06) is no longer searched for clouds (it never had any).
- **2026-09-19.** CI failure (prompt: "Some CI on github failed, can you have a look?"):
  `long_integrations/default_primitive_wet` (1 year, Julia 1.13) hit NaNs at step 7651 (31 July).
  This reproduces locally, while `main` passes the same test. The top layer over Greenland
  overheated (312 K) before the blow-up. Cause 1: ozone was distributed by σ, so over high
  orography the thinner top σ-layer took the full ozone absorption (≈45% higher heating rate at
  pₛ = 700 hPa). speedy.f90 avoids this by scaling with pₛ/p₀. Fix: ozone is now distributed by
  pressure overlap with [0, 50 hPa] and [50, 140 hPa] (`pressure_upper`, `pressure_lower` replace
  `σ_upper`, `σ_lower`). Heating rates no longer depend on pₛ, and the ozone column is retained.
  After the fix the CI test passes locally, but the top layer still reaches 310 K in August
  (265 K on `main`). Cause 2: the default `FriersonLongwaveTransmissivity` gives the top layer an
  optical depth of ≈0.02 at the poles, so it can hardly emit. speedy.f90 has a CO₂ band and a
  stratospheric correction term in its longwave scheme for this.
- **2026-09-19.** Version of `SpeedyWeather` bumped to `0.23.0-DEV` (new public types,
  changed default).

## Problem description

Our `OneBandShortwave` is a single broadband scheme with blended absorptivities. It uses a
fixed vertical ozone profile with no latitude or season dependence. speedy.f90 (Fortran SPEEDY,
Molteni 2003) instead has two bands and a latitude- and season-dependent ozone absorption.
Comparing the diagnostic clouds also showed that ours contained several bugs, so cloud cover
was about half of speedy.f90's.

## Background

speedy.f90 (`source/shortwave_radiation.f90`, `mod_radcon.f90`, `physics.f90`):

- **Bands**: visible (95%) and near-infrared (5%, `fband2`). The near-IR band is absorbed only
  by water vapor (`abswv2 = 15` per g/kg per 10⁵ Pa). It is not reflected by clouds or the
  surface and is fully absorbed at the surface.
- **Ozone**: `epssw = 0.02` of the TOA flux. `0.5ε` is absorbed in the upper stratosphere
  (layer 1). `0.4ε(1 + max(0, cos α) sin φ + 1.8 P₂(sin φ))` is absorbed in the lower
  stratosphere (layer 2), where α is the year angle from the northern winter solstice. Both are
  multiplied by the zenith correction factor. The seasonal phase is hardcoded for Earth.
- **Clouds**: cloud cover is `wpcl √min(10, P[mm/day]) + min(1, (RHmax − 0.3)/0.7)²`. RHmax is
  the maximum over tropospheric layers with q > 0.2 g/kg, and the cloud top is at that level
  (or the precipitation top if higher). Cloud albedo is 0.43. Clouds absorb in the visible
  band from the cloud top down to layer N−1.
- **Stratocumulus**: stability is `GSE = Δs/ΔΦ` between the lowest two layers. Cover is
  `clamp((GSE − 0.25)/0.15, 0, 1) · max(0.6 − 1.2 CLC, 0)`, and over land
  `max(CLS, 0.15) · RH_N`. The albedo is 0.5, applied at the top of the surface layer.

Bugs found in our `DiagnosticClouds`:

1. Rain rate was converted from m/s with `86400 rain / 1000` instead of `rain · 86400 · 1000`,
   so the precipitation term was ~0.
2. Stratocumulus stability used a raw `s_N − s_{N−1}` in J/kg (negative when stable) against
   dimensionless thresholds, so stratocumulus almost never formed.
3. The RH term used the lowest qualifying layer, while the cloud top was the highest one.
   speedy.f90 uses the maximum RH, with the cloud top at the same level.
4. The land stratocumulus minimum (`clsminl = 0.15`) was missing.
5. The specific humidity threshold was also applied to the layer directly above the surface
   layer. speedy.f90 always considers that layer, so ours had no humidity clouds wherever
   q < 0.2 g/kg (13% of columns, mostly polar). Arctic cloud cover was 0.37 vs 0.63.
6. Precipitation only included rain. speedy.f90 uses all precipitation, so snow is now included.

## Summary of changes

- `TwoBandShortwave <: AbstractShortwave` with components `clouds`, `transmissivity`, `ozone`
  and `radiative_transfer`. New default for `PrimitiveWetModel`.
- `TwoBandShortwaveTransmissivity`: visible and near-IR transmissivities (speedy.f90
  absorptivities, converted to kg/kg) in scratch arrays `a` and `b`. Returns
  `(; visible, near_infrared, zenith_factor)`.
- `SeasonalOzone` / `NoOzone` (new `ozone.jl`). The lower stratosphere absorption field
  `ozone_absorption_lower` is updated once per time step in the global parameterization step
  with a kernel. Absorption is distributed over layers by overlap with σ ∈ [0, 0.05] and
  [0.05, 0.14], so it works for any vertical resolution.
- `TwoBandShortwaveRadiativeTransfer`: reflection at the top of the cloud-top layer and the
  top of the surface layer. Visible-only surface reflection and upward beam.
- `DiagnosticClouds`: fixes 1–6 above, and `cloud_albedo` 0.6 → 0.43. This also affects
  `OneBandShortwave`. The layers are defined by σ (`σ_tropopause`, `σ_boundary_layer`) instead of
  speedy.f90's fixed layer indices, so the scheme works for any number of layers.

## Testing and verification

- `test/parameterizations/shortwave_radiation.jl` has new testsets for the two-band
  transmissivity (both bands in (0, 1], near-IR more strongly absorbed) and for ozone
  (non-negative, zero below σ_lower, column total = upper + lower, pole > equator), and for the
  ozone seasonal cycle following the orbit (NH > SH in January, symmetric in July, symmetric
  for zero tilt and flipped for negative tilt, shifted with the equinox).
- `free_troposphere_layers` is tested for 1–64 layers: the cloud layers stay within the σ-boundaries
  and never include the surface layer (for more than one layer), and 8 layers map to layers 2–7.
- Vertical resolution (T31, Δt = 20 min, same diagnostics as below): planetary albedo 0.334 / 0.314 /
  0.317, cloud cover 0.68 / 0.64 / 0.66, stratocumulus 0.029 / 0.027 / 0.018, atmospheric
  absorption 54 / 60 / 63 W/m² for 8 / 16 / 32 layers. Cloud tops stayed within layers 3–14 (L16)
  and 5–29 (L32). T21 with ≥16 layers blows up at time step 51, but it does so on `main` too,
  so that is unrelated to this PR.
- Global-mean shortwave budget: T31 with 8 layers from 1 January, area-weighted mean over
  days 35–40. speedy.f90 was built locally with an added diagnostic print.

| | Planetary albedo | Surface down [W/m²] | Atm. absorbed [W/m²] | Cloud cover |
|---|---|---|---|---|
| speedy.f90 | 0.306 | 208 | 62 | 0.59 |
| main `OneBandShortwave` | 0.238 | 233 | 69 | 0.31 |
| `TwoBandShortwave` | 0.342 | 209 | 53 | 0.70 |
| `OneBandShortwave` (fixed clouds) | 0.314 | 209 | 63 | 0.61 |

### Cloud comparison

Zonal means (area-weighted over latitude bands) of all cloud terms, T31 with 8 layers, days
35–40 from 1 January, shown as speedy.f90 / SpeedyWeather after fixes 1–6. `cltop` is the cloud
top layer index (N+1 = no cloud top), `GSE` is the stratocumulus stability Δs/ΔΦ, and
`P mm/d` is rain + snow.

| band | cloudc | clstr | Pterm | RHterm | cltop | RH(N-1) | P mm/d | GSE |
|---|---|---|---|---|---|---|---|---|
| -90..-60 | 0.60 / 0.51 | 0.06 / 0.03 | 0.08 / 0.10 | 0.54 / 0.42 | 6.22 / 5.45 | 0.77 / 0.67 | 0.51 / 0.54 | 0.62 / 0.14 |
| -60..-30 | 0.67 / 0.71 | 0.03 / 0.01 | 0.13 / 0.25 | 0.59 / 0.51 | 5.95 / 4.43 | 0.81 / 0.72 | 1.86 / 2.65 | 0.45 / 0.25 |
| -30..-10 | 0.53 / 0.66 | 0.06 / 0.04 | 0.11 / 0.13 | 0.44 / 0.56 | 4.85 / 3.48 | 0.74 / 0.68 | 3.84 / 1.17 | 0.35 / 0.27 |
| -10..10 | 0.57 / 0.77 | 0.06 / 0.04 | 0.10 / 0.18 | 0.51 / 0.66 | 5.38 / 3.19 | 0.76 / 0.72 | 5.18 / 1.79 | 0.37 / 0.32 |
| 10..30 | 0.48 / 0.67 | 0.05 / 0.03 | 0.07 / 0.18 | 0.43 / 0.52 | 6.06 / 3.64 | 0.67 / 0.61 | 1.89 / 1.76 | 0.32 / 0.25 |
| 30..60 | 0.66 / 0.54 | 0.08 / 0.04 | 0.15 / 0.16 | 0.59 / 0.44 | 5.86 / 5.18 | 0.77 / 0.62 | 2.17 / 1.62 | 0.47 / 0.12 |
| 60..90 | 0.63 / 0.66 | 0.10 / 0.05 | 0.09 / 0.06 | 0.57 / 0.62 | 6.06 / 5.41 | 0.80 / 0.81 | 0.59 / 0.23 | 0.85 / 0.42 |
| -90..90 | 0.59 / 0.66 | 0.06 / 0.03 | 0.11 / 0.16 | 0.52 / 0.54 | 5.70 / 4.21 | 0.76 / 0.68 | 2.68 / 1.62 | 0.44 / 0.25 |

With the algorithms now matching, the remaining differences come from the model climate:

- **Tropics and subtropics are too cloudy** (0.77 vs 0.57 at ±10°). The relative humidity
  maximum is in layer 3 (σ ≈ 0.31) in 70% of tropical columns, and cloud tops are much higher
  (3.2 vs 5.4). Our upper tropical troposphere is more humid than in speedy.f90, likely because
  Betts–Miller convection relaxes towards a humid reference profile, whereas speedy.f90 uses
  mass-flux convection.
- **Mid-latitudes are not cloudy enough** (0.54 vs 0.66 at 30–60°N). Relative humidity in the
  layer above the surface layer is lower (0.62 vs 0.77).
- **Stratocumulus is about half** (0.03 vs 0.06) because the lower troposphere is less stable
  (GSE 0.25 vs 0.44, and 0.42 vs 0.85 in the Arctic), which points to the boundary layer and
  vertical diffusion.
- **Tropical precipitation is low** (1.8 vs 5.2 mm/day at ±10°), which is worth a separate look.

## Documentation changes

- `docs/src/radiation.md`: new `TwoBandShortwave` section with the ozone formula. The cloud
  diagnosis description now covers max RH, the GSE definition and the land stratocumulus minimum.
- `docs/src/references.bib`: added Molteni (2003).

## Known limitations

- Cloud cover (0.70) is higher than in speedy.f90 (0.59) because of a different humidity
  climate (see "Cloud comparison"). The planetary albedo is slightly too high and atmospheric
  absorption too low.
- Instantaneous zenith angle with a diurnal cycle, whereas speedy.f90 uses daily means. Ozone
  is distributed by pressure instead of speedy.f90's pₛ/p₀ scaling (same heating rates, but the
  ozone column is kept over high orography).
- Our longwave cannot balance speedy.f90-strength ozone heating in the polar summer stratosphere
  (top layer up to 310 K), see revision log.
- The near-IR band is not reflected by the surface (as in speedy.f90), although real surfaces
  do reflect near-IR.

## Future work

- The σ-boundaries are fixed. A diagnosed boundary layer height (e.g. from `BulkRichardsonDiffusion`)
  or tropopause could replace them, but that would couple clouds to other components.
- T21 with ≥16 layers is unstable, which is unrelated to radiation (see testing).

- Tune the RH thresholds or cloud albedo against observed global-mean budgets.
- Longwave cloud effects (speedy.f90 `ablcl1`, `ablcl2`) are not implemented in our longwave
  schemes.
- Consider a prescribed ozone climatology instead of the analytic SPEEDY formula.
