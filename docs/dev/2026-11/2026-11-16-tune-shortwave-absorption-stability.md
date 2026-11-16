# Tuning shortwave water vapour and dry-air absorption for stable long T32L8 runs

> Status: **in progress**. The stratospheric water vapour feedback is identified and two fixes keep the
> stratosphere dry for 10 years. Combined candidates (15 years) and 1000-year GPU runs are running. No source code
> changed in the repository yet.

Date of initial draft: 2026-11-16

Base revision: `caed82d647c38b6fd793ff49ff674e62d437b8cc` (`main`, branch `mg/tune-water-vapour`)

## Originating prompt

> the model is currently unstable for long rollouts. this is most likely due to water vapour and dry absorption
> not being tuned correctly after the last two radiation related PRs were merged. We need to tune those for the
> model to run stable. In the end the T32L8 should run stable for millenia. Quick rollouts you can do here on the
> login nodes, long ones via slurm on cpu or gpu

## Revision log

- **2026-11-16, initial draft** (prompt: "make a markdown plan that we can continue to update while we wait").
  Written while the first 30-year diagnostic runs are in progress. Records the setup, first-year energy budgets
  and the planned experiments.
- **2026-11-16, scope decided** (prompt: "the observed one, yes cloud albedo, optimal depth can be tuned as well").
  The tuning targets the observed global energy budget (Wild et al. 2013). Cloud albedo and the Frierson longwave
  optical depths are in scope as well as the shortwave absorptivities.
- **2026-11-16, three sessions.** Sessions 7f, d3 and ed (forks of the same conversation) worked in parallel in the
  same workspace and split the work: 7f the baselines, `s1`/`s1b` and the combined candidates `c7`–`c13`; d3 the
  wv × dry grid and sweep2; ed the structural `water_vapor_pressure_top` variant and the long-run infrastructure.
  d3 and later 7f ended; ed continues with everything. No commits to the branch unless the user asks.
- **2026-11-16, long runs on GPU** (prompt: "submit a gpu job if you want it faster, not multithreaded"). Restartable
  `run_long.jl` with `ARCH=gpu` and `long_gpu.sh`, see Phase 3. First 1000-year runs submitted.

## Problem description

The default `PrimitiveWetModel` at T32L8 (`SpectralGrid(truncation = 32, nlayers = 8)`) is reported to be
unstable in long integrations after two radiation PRs:

- #1262 `TwoBandShortwave` with `SeasonalOzone` (new default), `DiagnosticClouds` fixes, cloud albedo 0.6 → 0.43,
  stratospheric emission in `OneBandLongwave`.
- #1371 Jin et al. ocean albedo (`JinOceanAlbedo`, new default ocean albedo).

The goal is to tune the shortwave water vapour and dry-air absorptivities of `TwoBandShortwaveTransmissivity` so
that the default T32L8 model runs stably for millennia, with a plausible global energy budget.

## Background

### Model setup

The default `PrimitiveWetModel` is fully coupled, so any persistent energy imbalance drifts the climate:

- `SlabOcean` with a 50 m mixed layer (heat capacity ≈ 2.1·10⁸ J/m²/K), no ocean heat transport
- `ThermodynamicSeaIce` (freezing at 271.35 K), ice albedo 0.6
- prognostic `LandModel`
- `Radiation(TwoBandShortwave, OneBandLongwave)`. The longwave uses `FriersonLongwaveTransmissivity`, a
  prescribed optical depth (τ₀ = 6 at the equator, 1.5 at the poles) that does not depend on humidity. Water
  vapour therefore only enters the radiation through the shortwave, so the model has no longwave water vapour
  feedback.

### Shortwave absorptivities (current defaults)

`TwoBandShortwaveTransmissivity` uses speedy.f90's values, converted from g/kg to kg/kg. Optical depth per layer
is `absorptivity · Δp/10⁵ Pa · zenith_factor`, with `zenith_factor = 1 + (1 − cos z)²`.

| Parameter | Value | speedy.f90 |
|---|---|---|
| `absorptivity_dry_air` (visible) | 0.033 per 10⁵ Pa | `absdry = 0.033` |
| `absorptivity_aerosol` (visible, × σ²) | 0.033 per 10⁵ Pa | `absaer = 0.033` |
| `absorptivity_water_vapor` (visible) | 22 per kg/kg per 10⁵ Pa | `abswv1 = 0.022` per g/kg |
| `absorptivity_water_vapor_near_infrared` | 15000 per kg/kg per 10⁵ Pa | `abswv2 = 15` per g/kg |
| `near_infrared_fraction` | 0.05 | `fband2 = 0.05` |
| `absorptivity_cloud_base` / `_limit` | 15 / 0.15 | `abscl1` / `abscl2` |

The previous one-band default (`BackgroundShortwaveTransmissivity`) had a visible water vapour absorptivity of 75
per kg/kg (tuned in #974 for ≈75 W/m² atmospheric absorption), about 3.4 times speedy.f90's visible value.

### Reference global-mean energy budget

Observational estimates (Wild et al. 2013, Trenberth et al. 2009, CERES), in W/m²: insolation 341, reflected
shortwave ≈100 (planetary albedo ≈0.29), OLR ≈239, atmospheric shortwave absorption ≈79 (range ≈74–91), surface
shortwave down ≈185, surface shortwave net ≈161, sensible heat ≈20, latent heat ≈84 (≈2.9 mm/day). speedy.f90 at
T31L8 (from the #1262 plan): planetary albedo 0.306, atmospheric absorption 62, surface down 208.

### Diagnostics

Script `run_diag.jl` (workspace `/p/tmp/maxgelbr/code/Speedy2/claude_tune_water_vapour_2026-11-16/`, outside the
repository) runs T32L8 with a callback that samples every 6 hours and writes 30-day global means (area-weighted
with the grid's quadrature weights) of the TOA and surface fluxes, precipitation, cloud cover, SST, land
temperature, sea ice, surface pressure, and layer means and extremes of temperature and humidity. It stops at the
first NaN/Inf. `annual.py` reduces this to annual means. Throughput is ≈12 s per simulated month on login02
(≈2.4 min per year, ≈600 SYPD including diagnostics).

## Findings so far

### First-year global means (instantaneous samples every 6 h)

| Run | Revision | Months | net TOA | OSR | OLR | albedo | atm SW | sfc SW down | SH | LH | P [mm/d] | cloud | SST [K] | sea ice |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `main_default` | `caed82d6` (HEAD) | 12 | −17.5 | 116.5 | 242.2 | 0.342 | 53.7 | 196.5 | 42.9 | 76.4 | 2.68 | 0.70 | 290.0 | 0.074 |
| `pre_albedo` | `3245c69f` (before #1371) | 5 | −17.0 | 114.5 | 243.7 | 0.336 | 54.4 | 196.8 | 41.2 | 78.0 | 2.78 | 0.70 | 291.1 | 0.046 |
| `pre_radiation` | `b43dc87c` (before #1262) | 5 | −42.2 | 153.8 | 229.6 | 0.451 | 49.8 | 159.2 | 37.8 | 70.1 | 2.54 | 0.69 | 290.1 | 0.047 |

Pending: the full 30 years of all three runs, the annual means for `pre_albedo` and `pre_radiation` (only
January–May so far, so not comparable with the 12-month `main_default` row), and the failure mode.

### Observations

1. **Large negative TOA imbalance on `main`** (−17.5 W/m² in year 1, −12.6 W/m² in January of year 2). SST falls
   from 290.0 to 288.7 K and sea ice rises from 0.07 to 0.10 between year 1 and January of year 2. A sea-ice albedo
   runaway towards a snowball, followed by a numerical blow-up, is the leading hypothesis for the instability, but
   is not yet confirmed.
2. **#1371 is not the main driver.** Over the same months, `pre_albedo` and `main` have nearly identical budgets.
3. **Before #1262 the budget was worse**, not better (−42 W/m², planetary albedo 0.45). This is the one-band scheme
   with cloud albedo 0.6 after the rain rate fix #1351 raised cloud cover. Whether that version was stable in long
   runs is pending.
4. **Atmospheric shortwave absorption is low**: 54 W/m² against ≈79 observed and 62 in speedy.f90. Surface
   shortwave down is high (196 vs ≈185).
5. **Sensible heat flux is about twice the observed value** (43 vs ≈20 W/m²). This is outside the scope of this
   tuning but is recorded here.
6. **The top layer warms and moistens year after year since #1262.** On `main`, the top-layer (σ = 0.0625,
   0–125 hPa) mean temperature rises 232 → 241 K and its maximum 264 → 283 K over 4 years. Its specific humidity
   rises 80-fold, 0.0018 → 0.15 g/kg (observed stratospheric water vapour is ≈0.003 g/kg). `pre_albedo` behaves the
   same. Before #1262 (`pre_radiation`) the top layer stays dry (0.001–0.004 g/kg), with maxima of 250–268 K. At
   0.15 g/kg, the near-IR absorptivity of 15000 per kg/kg gives the top layer a near-IR optical depth of ≈0.3–0.4.
   In polar summer that is ≈8 W/m², or ≈0.6 K/day, on top of ≈1.2 K/day of ozone heating, in a layer whose
   longwave cooling is weak (Frierson τ ≈ 0.02–0.08 plus the stratospheric emissivity of 0.05). This is a positive
   feedback: a warmer top layer has a much higher saturation humidity (≈3.8 g/kg at 240 K and 60 hPa), so moisture
   arriving from below is never condensed out, and it absorbs more near-IR. This is the leading candidate for the
   long-run instability (the "water vapour absorption" of the originating prompt).
7. **speedy.f90 has no water vapour absorption in its top layer.** Its layer 1 (0–50 hPa) only has dry-air
   absorption in the visible band (`tau2(:,:,1,1) = exp(-psaz*dhs(1)*absdry)`, and the near-IR loop is
   `do k = 2, kx`). Our 8 equally spaced σ-layers put 0–125 hPa into the top layer, which absorbs in both bands with
   its own humidity.
8. **The near-IR band is saturated in the troposphere.** With q ≈ 5 g/kg and Δp = 125 hPa, a lower-tropospheric
   layer has a near-IR optical depth of ≈9. Lowering `absorptivity_water_vapor_near_infrared` by a factor of ~10
   would therefore barely change the tropospheric absorption of that band, but would cut the stratospheric heating
   by the same factor.
9. **`ThermodynamicSeaIce` is an energy source, so net TOA can't reach zero.** On `main` in years 3–5, net TOA
   stays at −4.6 to −5.8 W/m² and the surface budget at ≈−6 W/m², while SST is flat (287.3 → 287.0 K). The
   atmosphere is in balance (+1 W/m²). In `filter!(vars, ::ThermodynamicSeaIce, model)`, SST below the freezing
   point is reset to 271.35 K after every time step, and the ice concentration is clamped at 1. An ice-covered
   ocean keeps losing heat to the atmosphere (in polar night ≈300 W/m² of longwave emission plus sensible heat into
   230–240 K air) but never cools. It isn't insulated by the ice. With ≈8% of the globe ice covered, this supplies
   the missing ≈5–6 W/m², and probably most of the doubled sensible heat flux (44 vs ≈20 W/m²). It also prevents a
   snowball. This is a structural issue of the sea ice model, not a tuning parameter, and it's out of scope here
   (see future work). For the tuning it means that the steady-state net TOA is about −(sea-ice heat source), not
   zero. The tuning therefore targets the individual budget terms (OSR, OLR, atmospheric and surface shortwave) and
   a steady SST, not net TOA = 0.
10. **Raising the absorptivities strengthens the stratospheric feedback** (session d3's wv × dry grid, 20 y runs).
   Visible water vapour absorptivity raises atmospheric shortwave absorption about linearly: 53.7 / 65.3 / 78.2 /
   92.6 W/m² for `wv` = 22 / 44 / 66 / 88 in year 1. `dry` = 0.066 adds about 12 W/m². However, OLR rises almost
   as much as OSR falls (`wv` = 66: OSR −11, OLR +9), so net TOA improves by only 1–3 W/m². Every increase of `wv`
   or `dry` heats the top layers faster. With `wv` = 88 and `dry` = 0.066, q₁ reaches 0.64 g/kg, q₂ 0.70 g/kg, and
   the top-layer maximum 307 K by May of year 3. The observed atmospheric absorption is therefore only reachable once
   the stratospheric water vapour feedback is cut (finding 6).
11. **A lower near-IR water vapour absorptivity removes the stratospheric feedback** (`s1_*`, year-2 annual means).
   The near-IR band heats layers 1 and 2 enough to let moisture in, and the moisture then absorbs more.

   | run | atm SW | T₁ / T₂ [K] | max T₁ | q₁ / q₂ [g/kg] | SST | net TOA | OSR | LW sfc down |
   |---|---|---|---|---|---|---|---|---|
   | default (`wv_nir` = 15000) | 50.7 | 237 / 239 | 269 | 0.090 / 0.300 | 288.0 | −10.6 | 120.9 | 307.9 |
   | `wv_nir` = 5000 | 51.6 | 230 / 220 | 265 | 0.0017 / 0.019 | 288.8 | −6.6 | 116.0 | 315.3 |
   | `wv_nir` = 1500 | 50.1 | 231 / 219 | 264 | 0.0008 / 0.011 | 289.2 | −6.1 | 114.8 | 318.9 |
   | `wv_nir` = 0 | 35.0 | 230 / 218 | 266 | 0.0009 / 0.014 | 289.5 | −5.6 | 117.4 | 317.4 |
   | `strat_emis` = 0.1 | 50.9 | 215 / 209 | 266 | 0.0059 / 0.009 | 288.1 | −9.6 | 118.3 | 308.8 |
   | `wv` = 60 | 66.2 | 243 / 246 | 283 | 0.166 / 0.446 | 288.1 | −12.0 | 115.9 | 311.1 |
   | `wv` = 100 | 84.8 | 255 / 254 | 322 | 0.337 / 0.652 | 288.0 | −13.6 | 109.5 | 314.7 |
   | `dry` = 0.066 | 62.3 | 243 / 247 | 279 | 0.111 / 0.384 | 288.0 | −11.9 | 116.1 | 309.3 |
   | `aer` = 0.066 | 54.0 | 236 / 240 | 264 | 0.084 / 0.306 | 288.1 | −10.9 | 119.6 | 309.1 |
   | `cloud_albedo` = 0.35 | 55.4 | 240 / 243 | 271 | 0.121 / 0.376 | 289.9 | −6.3 | 103.1 | 320.9 |
   | `cloud_albedo` = 0.30 | 58.9 | 241 / 245 | 276 | 0.118 / 0.413 | 291.1 | −3.1 | 91.8 | 328.8 |
   | `τ₀` = 7 / 1.75 | 52.1 | 235 / 234 | 271 | 0.097 / 0.236 | 289.2 | −5.0 | 117.4 | 323.2 |
   | `τ₀` = 8 / 2 | 53.9 | 232 / 230 | 264 | 0.075 / 0.195 | 290.3 | −1.4 | 114.4 | 336.8 |

   With `wv_nir` ≤ 5000 the top two layers stay dry, and layer 2 at 219 K is close to the observed ≈218 K at
   200 hPa. A value of 1500 leaves the tropospheric absorption unchanged, because the band stays saturated there;
   0 removes ≈16 W/m². Increasing the stratospheric emissivity also keeps the top dry, but cools layers 1–2 well
   below observed values. Increases of `wv` or `dry` with the default near-IR absorptivity strengthen the feedback.
   For the budget: `cloud_albedo` ≈ 0.34 gives OSR ≈ 100, `τ₀` ≈ 8 / 2 gives LW surface down ≈ 337 (obs ≈ 342),
   and `wv` ≈ 75 should give atmospheric absorption ≈ 79 W/m².
12. **Cutting water vapour absorption above a pressure level works too, but needs 250 hPa at L8** (`ed_pt*`,
   10-year runs with `water_vapor_pressure_top`, final-year means).

   | run | T₁ / T₂ [K] | max T₁ | q₁ / q₂ [g/kg] | atm SW | OSR | OLR | sfc SW down | SST |
   |---|---|---|---|---|---|---|---|---|
   | `main` (year 5) | 239 / 230, rising | 283 | 0.135 / 0.152 | 49.8 | 121.7 | 225.3 | 198.4 | 287.0 |
   | `p_top` = 50 hPa | 229 / 234 | ≈265 | 0.25 / 0.25 | 49.5 | 121 | 224 | 199 | 286.8 |
   | `p_top` = 100 hPa | 229 / 223 | ≈265 | 0.043 / 0.17 | 50.5 | 118 | 227.5 | 200 | 287.2 |
   | `p_top` = 140 hPa | 229 / 220 | ≈263 | 0.025 / 0.12 | 50.3 | 119 | 226.5 | 200.5 | 287.4 |
   | `p_top` = 250 hPa | 229 / 218 | ≈268 | 0.0022 / 0.013 | 51.5 | 117 | 228 | 200.6 | 287.9 |
   | 140 hPa, `wv` = 66, `dry` = 0.066 | 235 / 228 | 270–286 | 0.082 / 0.29 | 82.4 | 103 | 242.5 | 182 | 287.6 |
   | 140 hPa, `wv` = 88 | 230 / 223 | 268–278 | 0.0094 / 0.128 | 80.7 | 105 | 240.5 | 182 | 287.7 |
   | 250 hPa, `wv` = 66, `dry` = 0.066 | 235 / 225 | ≈275 | 0.0017 / 0.013 | 80 | 108 | 238 | 179 | 287.8 |

   50 hPa (speedy.f90's top-layer boundary) is too shallow for 8 equally spaced σ-layers: layer 2 (125–250 hPa)
   still heats and moistens. 100–140 hPa stops the warming but leaves layer 2 moist. Only 250 hPa, which excludes
   layer 2 as well, keeps the stratosphere as dry as before #1262, also with `wv` = 66 and `dry` = 0.066. This
   confirms that the moistening is driven by layer 2 absorbing with its own water vapour, which `wv_nir` = 1500
   (finding 11) prevents equally well (q₂ 0.011 g/kg). The tuning-only route needs no code change and keeps water
   vapour absorption in the tropical upper troposphere, so it is preferred if it stays dry with `wv` = 75–100
   (`s1b_*`). With the feedback cut, the observed shortwave budget is reachable: `wv` ≈ 66–88 gives atmospheric
   absorption ≈ 80 W/m², OSR ≈ 103–108 and surface shortwave down ≈ 180, with layers 1–2 flat over 10 years.
   The remaining gap is the surface longwave (LW down 315 vs ≈342 observed), which τ₀ addresses.
13. **`wv_nir` = 1500 keeps the stratosphere dry for any visible absorptivity tested** (`s1b_*`, 9 years, means
   over the last 3). With `wv` = 60 / 80 / 100, q₁ stays at 0.0010–0.0015 g/kg and q₂ at 0.007–0.020 g/kg, and
   T₁ / T₂ at 231–232 / 219–221 K, while atmospheric SW absorption rises 67 / 77 / 89 W/m². `wv_nir` 500 / 3000
   behave the same. This tuning-only route is preferred over `water_vapor_pressure_top`, which isn't needed.
   Sweep2 (`c1`–`c6`, default `wv_nir`) confirms the feedback: all six reach top-layer maxima of 297–321 K with
   q₁ = 0.2–0.36 g/kg. `τ₀` = 7.2 / 1.8 alone stayed dry for 15 years while neighbouring configurations didn't,
   so the default sits near a threshold and could flip at any time in a long run. That is a further reason to cut
   the feedback through `wv_nir` rather than rely on the climate staying below the threshold.
14. **First combined candidates are close to the observed budget** (`c7`–`c13`, all `wv_nir` = 1500, 11 years,
   means over the last 3).

   | run | OSR | OLR | atm SW | sfc SW down | LW down | LW up | SH | LH | SST | ice | net TOA |
   |---|---|---|---|---|---|---|---|---|---|---|---|
   | observed | 100 | 239 | 79 | 185 | 342 | 397 | 20 | 84 | | | |
   | `c7`: `wv` 75, `τ₀` 7.5 / 1.875 | 99.5 | 242.4 | 85.8 | 178.2 | 349.1 | 398.5 | 34.5 | 74.3 | 291.0 | 0.056 | −0.6 |
   | `c10`: `wv` 75, `τ₀` 8 / 2 | 98.7 | 242.5 | 88.2 | 175.5 | 355.4 | 402.3 | 33.7 | 75.3 | 291.6 | 0.046 | 0.0 |
   | `c13`: `wv` 85, `τ₀` 7.5 / 1.875 | 97.6 | 244.3 | 92.1 | 173.2 | 350.5 | 398.5 | 33.0 | 72.8 | 291.0 | 0.054 | −0.6 |
   | `c8`: `c7` + `cloud_albedo` 0.40 | 94.3 | 247.2 | 90.4 | 178.3 | 353.8 | 402.8 | 33.8 | 75.7 | 291.7 | 0.049 | −0.3 |

   Raising `τ₀` warms the climate and moistens the troposphere, so the shortwave absorption overshoots (85.8 instead
   of ≈79 W/m² for `wv` = 75). From `s1b`, absorption changes by ≈0.5 W/m² per unit of `wv`, so `wv` ≈ 60–65
   should give ≈79. The default `cloud_albedo` = 0.43 already gives OSR ≈ 100 once the atmosphere absorbs more.
   Net TOA approaches zero because the sea-ice heat source (finding 9) shrinks with the ice cover (0.12 → 0.05).
   Sensible heat (≈34 vs 20) and latent heat (≈74 vs 84) remain off, as expected from the out-of-scope surface
   fluxes. Round 2 (`r2_*`) refines `wv` 55–65 and `τ₀` 7.0–7.5.
15. **Caveat on the direction of the tuning.** With a grey longwave at fixed optical depth, moving shortwave
   absorption from the surface into the atmosphere cools the surface for a given TOA balance (anti-greenhouse
   effect). The only TOA effect of more absorption is a lower OSR, from absorption above the cloud tops and in the
   reflected beam. More absorption alone may therefore not stop a cold drift. The experiments below test this
   directly.

## Plan

### Phase 1: reproduce and characterise the instability (in progress)

- 30-year runs of `main_default`, `pre_albedo` and `pre_radiation` on login02.
- Record when and where the first NaN appears, and which precursor goes first (sea ice, SST, top-layer
  temperature, wind speed).
- If none of them blows up within 30 years, run longer (100+ years) on Slurm.

### Phase 2: sensitivity of the energy budget to the absorptivities

Short runs (2–5 years), many in parallel as Slurm jobs, varying one parameter at a time around the defaults:

- `absorptivity_water_vapor` (visible): 22 (default), 40, 60, 75, 100
- `absorptivity_dry_air` (visible): 0.033 (default), 0.05, 0.075, 0.1
- `absorptivity_aerosol`: 0.033 (default), 0.05, 0.075
- `absorptivity_water_vapor_near_infrared` and `near_infrared_fraction` only if the visible absorptivities can't
  reach the target

Diagnostics per run: annual net TOA, planetary albedo, atmospheric and surface shortwave, SST and sea-ice trend.
The target is atmospheric shortwave absorption of ≈70–80 W/m², a TOA imbalance that approaches zero without a
runaway, and global-mean SST and sea-ice extent in a plausible range.

### Phase 2b: cut the stratospheric water vapour feedback

Two alternatives, run in parallel:

- **Tuning only:** lower `absorptivity_water_vapor_near_infrared` (`s1_nir*`). The near-IR band stays saturated in the
  troposphere, so the cost is small.
- **Structural (session ed):** a new `TwoBandShortwaveTransmissivity` option `water_vapor_pressure_top` [Pa], with
  default 0 for the current behaviour. Water vapour (visible and near-IR) only absorbs in the part of each layer
  below it, using the same pressure overlap as `SeasonalOzone`. It is tested at 50 / 100 / 140 hPa with the default
  absorptivities and with `wv` = 66 / `dry` = 0.066 and `wv` = 88 / `dry` = 0.033. This follows speedy.f90, whose
  0–50 hPa layer has no water vapour absorption. It is developed in its own worktree (`tree_ed_strat`), not in the
  repository.

### Phase 3: validate candidates

- Long runs go on the GPU (prompt: "submit a gpu job if you want it faster, not multithreaded"). On the CPU,
  4 Julia threads were slightly slower than 1 (≈13 vs ≈12 s per simulated month). The CUDA run of T32L8 does
  ≈3600 SYPD on login02's A40 including the diagnostics, so a millennium takes a few hours on an H100.
- `run_long.jl` (workspace) runs T32L8 with the same monthly diagnostics as `run_diag.jl`, frozen as a copy, on
  CPU or GPU (`ARCH=gpu`). It runs in chunks and writes `runs/LABEL_restart/restart.jld2` after each one, using
  `WriteVariablesRestartFile`, written to a temporary file and then moved into place. If that file exists it
  resumes from it (`StartFromFile`, which also restores the clock, ocean, land and sea ice). Tested on CPU and on
  the A40: the resumed run continues the same climate (SST 291.1 vs 291.2 K, sea ice 0.049 vs 0.048 one month
  after the restart point).
- `long_gpu.sh` runs all lines of a config file concurrently on one GPU (`gpushort` QOS, 1 day), one process each.
  Resubmitting the same config resumes every run. The GPU environment (`tree/SpeedyWeather/test/GPU/CUDA`) must be
  precompiled on login02, never on a GPU-less node, because the depot is shared.
- 100-year runs of the 2–3 best candidates, then one 1000-year run of the chosen configuration.
- Acceptance: no NaN in 1000 years, global-mean SST and sea ice without a long-term trend after spin-up, top-layer
  temperature and maximum wind speed bounded.

### Phase 4: change the defaults

Only once the parameters are settled.

## Summary of changes

Planned, pending the experiments:

- New default values in `TwoBandShortwaveTransmissivity` (`shortwave_transmissivity.jl`), with a comment giving the
  tuning target and this document.
- `CHANGELOG.md` entry under `## Unreleased`.
- `SpeedyWeather` is already at `0.25.0-DEV`, so no version bump is needed for a change of default values.

## Testing and verification

- Existing unit tests in `test/parameterizations/shortwave_radiation.jl` and `radiation.jl` must still pass with
  the new defaults (for example, transmissivities in (0, 1] and near-IR absorbed more strongly than visible).
- `test/long_integrations/default_primitive_wet.jl` (1 year) must pass.
- The 1000-year T32L8 run (Phase 3) is the main verification. Its log and annual means will be summarised here.

## Documentation changes

- `docs/src/radiation.md`: the `TwoBandShortwave` section states that the absorptivities are speedy.f90's values.
  If they change, document the new values and the tuning target.

## Known limitations

- Global means are from instantaneous samples every 6 hours, not time-step averages. The diurnal cycle is sampled
  at 4 fixed UTC times, which is adequate for global means.
- The longwave has no water vapour dependence, so this tuning cannot introduce a water vapour longwave feedback.

## Future work

- The sensible heat flux is about twice the observed value. It is out of scope here and needs its own look at
  the surface fluxes and boundary layer.
- A longwave scheme with a water vapour dependence would replace the prescribed Frierson optical depth.
- Energy-conserving sea ice: an ice surface temperature that can drop below freezing (insulating ice with
  conductive flux), and ice growth that removes the latent heat of fusion from the ocean instead of resetting SST.
  Required for a closed TOA budget (finding 9).

## Decisions

1. **Tuning target:** the observed global-mean energy budget (Wild et al. 2013), not speedy.f90's.
2. **Scope:** the shortwave absorptivities, `cloud_albedo` of `DiagnosticClouds` (0.43) and the optical depths
   `τ₀_equator`, `τ₀_pole` of `FriersonLongwaveTransmissivity` may all be tuned in this PR. The surface fluxes
   stay out of scope.

## Experiment log

| Label | Revision / config | Length | Status | Notes |
|---|---|---|---|---|
| `smoke` | `main` | 60 days | completed | script check, net TOA −14 to −19 W/m² |
| `main_default` | `main` defaults | 30 years | running | year 4: net TOA −4.5, SST 287.2 K, ice 0.115, top layer 240 K / max 283 K, q₁ 0.15 g/kg |
| `pre_albedo` | `3245c69f` | 30 years | running | same as `main` within noise |
| `pre_radiation` | `b43dc87c` | 30 years | running | colder (SST 283 K), top layer dry (q₁ ≈ 0.002 g/kg) |
| `s1_*` (Slurm 3279280, 7f) | one-at-a-time: `wv_nir` 0/1500/5000, `wv` 60/100, `dry` 0.066, `aer` 0.066, `cloud_albedo` 0.35/0.30, `τ₀` ×7/6 and ×4/3, `strat_emis` 0.1, plus default | 6 years each | running | year-2 results in finding 11: `wv_nir` ≤ 5000 keeps the stratosphere dry |
| `s1b_*` (Slurm 3279354, 7f) | `wv_nir` 500/1500/3000 × `wv` 60/80/100, `dry` 0.05/0.066 combinations, and `wv_nir` = 1500 alone for 30 y | 10 years (one 30) | submitted | does the stratosphere stay dry with more visible absorption? |
| `wv{44,66,88}_dry{0.033,0.066}` (login02, d3) | 2D grid with default `wv_nir` | 20 years | running | finding 10 |
| sweep2 `cl*`, `tau*`, `fl*`, `c1`–`c6` (Slurm 3279283, d3) | cloud albedo, `τ₀`, `fₗ`, and combinations with default `wv_nir` | 15 years | running | session d3 ended. 7f reads the results |
| `c7`–`c13` (Slurm 3279366, 7f) | all `wv_nir` = 1500: `wv` 65/75/85, `cloud_albedo` 0.37/0.40/0.43, `τ₀` 7.5/1.875 and 8/2, one with `dry` = 0.045 | 15 years | submitted | first combined candidates |
| `ed_*` (Slurm 3279344, ed) | `water_vapor_pressure_top` 50/100/140/250 hPa, alone and with `wv` 66/88, `dry` 0.066 | 10 years | completed, all stable | finding 12, worktree `tree_ed_strat` |
| `r2_*` (Slurm 3279654, ed) | all `wv_nir` = 1500: `wv` 55/60/65 × `τ₀` 7.0/7.2/7.5, polar `τ₀` 2.0, `cloud_albedo` 0.45 | 20 years | running | round 2 around `c7` |
| `L_main`, `L_nir1500`, `L_pt250` (GPU Slurm 3279622, ed) | `main` defaults; `wv_nir` = 1500; `water_vapor_pressure_top` = 250 hPa | 1000 years, 25-year restart chunks | running (H100) | does `main` blow up, and when? Do both stratosphere fixes hold for a millennium? ≈2.1 s per simulated month per run, done ≈03:00 on 2026-11-17 |
| `L2_wv60_t75`, `L2_wv60_t72`, `L2_wv60_t70`, `L2_wv65_t72` (GPU Slurm 3279736, 5d) | all `wv_nir` = 1500: `wv` 60 with `τ₀` 7.5/1.875, 7.2/1.8, 7.0/1.75, and `wv` 65 with `τ₀` 7.2/1.8 (the `r2` grid's centre, see its year 1) | 1000 years, 25-year restart chunks | running (H100), expected ≈05:00–07:00 on 2026-11-17 | validation of the final candidates. `r2_*` (same parameters, 20 y on CPU) picks the best budget match |
