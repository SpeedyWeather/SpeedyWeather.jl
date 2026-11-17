# Tuning shortwave water vapour and dry-air absorption for stable long T32L8 runs

> Status: **in progress**. On `origin/main` the default Byrne–O'Gorman longwave runs away. The new defaults
> (decision 4: `b` = 750, `a` = 1.45 at the new default CO₂ of 360 ppm, shortwave `dry` 0.066, `wv_nir` 1500) are
> set in the worktree `tree_defaults` (branch `mg/tune-water-vapour-defaults`, uncommitted), and their 1000-year
> GPU validation `D_*` is running.

Date of initial draft: 2026-11-16

Base revision: `b06af911abf015a6cb2bcba991959908dbb87a03` (`origin/main` with #1279) since 2026-11-17. Before that,
`caed82d647c38b6fd793ff49ff674e62d437b8cc` (stale local `main`, see revision log).

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
- **2026-11-17, wrong base revision** (prompt: "yes, this is about using ByrneOGormanLongwaveTransmissivity. What were
  you doing before, why didn't you use the default?"). The branch was created from a local `main` at `caed82d6`
  and `origin/main` was never fetched. #1279 (`ByrneOGormanLongwaveTransmissivity` as the default, `dry_absorption`
  1.3, `water_vapor_absorption` 1997.9, Frierson defaults 8.4 / 2.1 / 0.25, default CO₂) had been merged at
  2026-11-16 15:48 CET, three hours before the originating prompt. All runs up to here therefore used the old
  default `FriersonLongwaveTransmissivity`, and "water vapour and dry absorption" was read as the shortwave
  `absorptivity_water_vapor` / `absorptivity_dry_air`. The prompt means the longwave `water_vapor_absorption` and
  `dry_absorption`. Consequences:
  - On `caed82d6`, `main` itself ran 1000 years on the GPU without a NaN, with a steady climate after ≈5 years
    (`L_main`, finding 16). The instability is therefore not in that code. The prime suspect is the longwave water
    vapour feedback of #1279, whose constants were calibrated with 2-year spin-ups.
  - Still valid: the shortwave stratospheric near-IR feedback (findings 6–8, 11–13), because the shortwave is
    unchanged in #1279. Its importance with Byrne–O'Gorman has to be re-tested: the dry longwave optical depth of
    the top layer is ≈0.16 instead of Frierson's ≈0.05, so the top layer emits more. Also still valid: the
    sea-ice energy source (finding 9) and the run infrastructure. The scripts now keep each tree's default
    longwave, and `lw_dry` / `lw_wv` / `lw_co2` set the Byrne–O'Gorman parameters.
  - Not transferable: the Frierson-based candidates (`c*`, `r2_*`, `L2_*`, `τ₀` 7.2 / 1.8 with `wv` 65).
  - New worktree `tree_1279` at `b06af911a`. Reproduction and sensitivity runs `M_*` were submitted.
- **2026-11-17, default CO₂ 360 ppm** (prompts in decision 3). `DEFAULT_CO2` = 360 and reused by
  `ByrneOGormanLongwaveTransmissivity.reference_co2`, so `dry_absorption` is the dry optical depth at the default
  CO₂. The change lives in the new worktree `tree_defaults`, which also gets the tuned defaults later; the tuning
  continues there (`W_*`), and earlier `a` values translate by −0.0508.
- **2026-11-17, new defaults chosen** (prompt: "and then continue the tuining of the water-vapour longwave
  absorption and dry_absorption"). After `W_*` (finding 22), decision 4 sets the defaults in `tree_defaults`. The
  1000-year validation `D_*` runs the code defaults, not a parameter list.
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
16. **On `caed82d6` (Frierson longwave), every configuration is stable for 1000 years** (GPU, H100, `run_long.jl`,
   century means; `L2_*` at year ≈835 when read).

   | run | net TOA | OSR | OLR | atm SW | sfc SW down | LW down | SST | ice | T₁ / T₂ | max T₁ | q₁ [g/kg] |
   |---|---|---|---|---|---|---|---|---|---|---|---|
   | `L_main` | −4.3 | 121.3 | 224.2 | 49.6 | 199.5 | 302.0 | 286.87 | 0.121 | 239 / 230 | 293–308 | 0.13 |
   | `L_nir1500` | −2.8 | 115.9 | 228.1 | 49.0 | 204.0 | 314.3 | 288.38 | 0.096 | 230 / 218 | 271–291 | 0.036 |
   | `L_pt250` | −4.0 | 117.9 | 227.3 | 51.2 | 201.3 | 309.3 | 287.85 | 0.118 | 229 / 218 | 272–281 | 0.002 |
   | `L2_wv65_t72` | −0.9 | 102.6 | 239.6 | 77.5 | 184.5 | 342.8 | 290.59 | 0.062 | 229 / 220 | 269–281 | 0.003 |

   After ≈5 years of spin-up all century means agree to within 0.05 K. The stratospheric feedback on `caed82d6`
   saturates at a warm, moist equilibrium (top layer up to 308 K) instead of running away. `L2_wv65_t72` (5d) matched
   the observed budget term by term on that base, but it is Frierson-based and so superseded by the
   Byrne–O'Gorman tuning.
17. **On `origin/main` (`b06af911a`, Byrne–O'Gorman longwave) the default runs into a runaway greenhouse**
   (`M_*`, CPU, annual means; τ from the logged humidity profile with equal Δσ, PW from the mass-weighted q).

   | run | year | SST | net TOA | OLR | OSR | PW [kg/m²] | τ_dry | τ_wv |
   |---|---|---|---|---|---|---|---|---|
   | default (`a` = 1.3, `b` = 1997.9) | 1 / 2 / 3 / 4 | 292.4 / 293.7 / 295.1 / 296.1 | +6.7 / +5.7 / +6.9 / +6.5 | 223.5 → 229.1 | 111.1 → 107.0 | 29 → 42 | 1.27 | 5.7 → 8.3 |
   | `a` = 1.5 | 1 / 2 / 3 | 292.9 / 295.2 / 297.3 | +13.1 / +12.3 / +11.8 | 219.0 → 228.2 | 109.1 → 101.5 | 32 → 50 | 1.46 | 6.2 → 9.8 |
   | `b` = 2500 | 1 / 2 / 3 | 292.8 / 295.1 / 297.3 | +12.6 / +12.5 / +11.5 | 219.5 → 227.5 | 109.2 → 102.2 | | | |
   | `a` = 1.1 | 1 / 2 / 3 / 4 | 291.9 / 292.3 / 292.7 / 293.0 | −0.1 / −1.5 / −1.5 / −1.6 | ≈229 | ≈113 | 27 → 30 | 1.07 | 5.3 → 5.9 |
   | `b` = 1500 | 1 / 2 / 3 / 4 | 291.8 / 292.2 / 292.6 / 293.0 | −0.2 / −1.5 / −1.8 / −2.0 | ≈230 | ≈112 | 27 → 30 | 1.27 | 4.0 → 4.4 |

   With the default the imbalance doesn't decrease as the planet warms: water vapour supplies ≈80% of the column
   optical depth, and precipitable water grows ≈12% per K, so OLR rises only ≈1.5 W/m² per K. Reflected
   shortwave meanwhile falls ≈1.8 W/m² per K as sea ice, snow and cloud cover retreat. The net feedback is ≈0 or
   positive, which is the long-run instability. Weaker water vapour absorption or a weaker greenhouse (`a` = 1.1)
   is close to balance, but still warms by ≈0.4 K per year, the steady net TOA of −1.5 to −2 W/m² being the
   sea-ice energy source. The tuning therefore needs a clearly negative feedback, not only a balanced mean
   state: shift optical depth from the water vapour term `b·q` to the dry term `a`. The total τ stays ≈5.3 (the
   balanced value), and dτ/dT shrinks. Each configuration is paired with a +2% solar constant run
   (1392.3 W/m², ≈+5 W/m² forcing) to measure the warming per W/m² directly.
18. **The default blows up after 24 years, and the climate has a hot second state** (`M_*`, `G_*`, `X_*`, 30 years,
   `_S` = +2% solar twin; sensitivity = ΔSST over the last 5 years per 4.8 W/m²).
   - **Blow-ups (NaN):** the default in year 24, `wv_nir` = 1500 in year 26, `a` = 1.5 in year 11, `b` = 2500 in
     year 9, `b` / `a` = 1000 / 2.5, 1000 / 3.0, 750 / 3.0, 750 / 3.5, 500 / 4.0, 500 / 4.5 and 500 / 4.0 with
     `cloud_albedo` 0.40 in years 7–16, and the +2% solar twins of the default and of `a` = 1.1 (years 11 and 21).
     All follow the same path. SST rises to 318–324 K, OSR collapses to 22–44 W/m² (planetary albedo 0.07–0.13) as
     cloud cover disappears and water vapour absorbs most of the shortwave, and net TOA reaches +25 to +54 W/m².
     The lower troposphere becomes almost isothermal (325–327 K) and precipitation drops below 1 mm/day before
     the NaN.
   - **Hot-house survivors:** with little water vapour longwave absorption, the climate settles at a hot second
     equilibrium instead of blowing up.

     | run | SST years 26–30 | trend [K/decade] | OSR | OLR | atm SW | ice | cloud | twin ΔSST | K per W/m² |
     |---|---|---|---|---|---|---|---|---|---|
     | `b` = 0, `a` = 4.8 | 310.1 | 0.0 | 68 | 272 | 136 | 0 | 0.43 | +0.8 | 0.17 |
     | `b` = 0, `a` = 6.0 | 314.0 | 0.0 | 64 | 277 | 158 | 0 | 0.42 | +0.7 | 0.15 |
     | `b` = 250, `a` = 5.0 | 315.6 | 0.0 | 68 | 273 | 164 | 0 | 0.45 | +1.8 | 0.37 |
     | `b` = 1000, `a` = 2.0 | 309.0 | +3.1 | 81 | 258 | 126 | 0 | 0.50 | +11.7 | 2.44 |
     | `b` = 1500, `a` = 1.3 | 297.7 | +1.8 | 102 | 240 | 78 | 0.037 | 0.60 | +3.6 | 0.74 |
     | `b` = 1998, `a` = 1.1 | 296.0 | −0.9 | 106 | 237 | 72 | 0.049 | 0.62 | +10.6, NaN in year 21 | 2.2 |

   Two corrections to finding 17. First, the assumption of equal column τ was wrong: a dry absorber uniform in
   pressure puts much more opacity into the upper troposphere than water vapour does, so it is a far stronger
   greenhouse for the same τ, and `a` = 4.5–6 makes a hot house. Second, the dominant positive feedback is in the
   shortwave. The water vapour shortwave absorptivity is linear in q without saturation (visible 22 per kg/kg),
   and the diagnostic cloud cover falls with warming. Atmospheric shortwave absorption therefore rises from ≈75 to
   130–165 W/m², and OSR from the normal ≈105 falls to 65–70 W/m² (hot house) or 22–44 W/m² (runaway). Near the
   observed climate (SST ≈ 291–297 K), every configuration so far is still drifting or extremely sensitive. Next
   (`Y_*`, 2 GPU jobs, 40 years with twins): much smaller `a` with weak `b` (0 / 250 / 500), and the paper value
   `a` = 0.87 or 0.7 with the default `b`, all with `wv_nir` = 1500. A water vapour shortwave absorption that
   saturates with q would be a structural fix of the shortwave feedback, if tuning can't make the normal state
   robust.
19. **With small `a`, all 24 `Y_*` runs are stable and equilibrated after 40 years** (GPU, `wv_nir` = 1500; means
   over years 35–39; drift = change since years 25–29, ≤ 0.06 K per decade in all runs; K per W/m² from the +2%
   solar twin, which is also stable in all cases).

   | `b` / `a` | SST | net TOA | OSR | OLR | atm SW | sfc SW down | LW down | LW up | SH | LH | P | cloud | ice | max T₁ | K per W/m² |
   |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
   | observed | ≈291 | 0 | 100 | 239 | 79 | 185 | 342 | 397 | 20 | 84 | 2.9 | | | | ≈0.8 |
   | 0 / 2.0 | 285.9 | −7.1 | 118 | 230 | 54 | 198 | 227 | 363 | 13 | 27 | 0.93 | 0.65 | 0.146 | 241 | 0.21 |
   | 0 / 2.5 | 291.8 | −2.2 | 107 | 236 | 66 | 191 | 279 | 399 | 13 | 39 | 1.32 | 0.62 | 0.058 | 242 | 0.22 |
   | 0 / 3.0 | 297.8 | +0.3 | 89 | 252 | 86 | 186 | 336 | 439 | 12 | 52 | 1.77 | 0.52 | 0.015 | 246 | 0.24 |
   | 0 / 3.5 | 302.3 | +0.9 | 82 | 259 | 101 | 176 | 377 | 469 | 11 | 56 | 1.90 | 0.49 | 0.001 | 245 | 0.23 |
   | 250 / 2.0 | 290.3 | −4.7 | 116 | 231 | 59 | 192 | 278 | 388 | 18 | 43 | 1.50 | 0.66 | 0.100 | 242 | 0.29 |
   | 250 / 2.5 | 298.4 | +0.2 | 93 | 248 | 84 | 185 | 368 | 445 | 18 | 70 | 2.38 | 0.56 | 0.021 | 251 | 0.34 |
   | 250 / 3.0 | 304.0 | +1.0 | 84 | 256 | 104 | 171 | 424 | 484 | 16 | 77 | 2.60 | 0.52 | 0.002 | 247 | 0.39 |
   | 500 / 1.5 | 287.5 | −9.8 | 120 | 231 | 53 | 199 | 254 | 370 | 21 | 41 | 1.42 | 0.66 | 0.172 | 245 | 0.24 |
   | 500 / 2.0 | 296.3 | −1.1 | 100 | 243 | 75 | 188 | 358 | 431 | 22 | 73 | 2.51 | 0.59 | 0.045 | 252 | 0.38 |
   | 500 / 2.5 | 303.5 | +1.1 | 87 | 253 | 101 | 171 | 432 | 483 | 19 | 84 | 2.84 | 0.53 | 0.003 | 249 | 0.81 |
   | 1997.9 / 0.7 | 285.1 | −17.5 | 126 | 233 | 46 | 208 | 254 | 353 | 34 | 54 | 1.87 | 0.65 | 0.255 | 258 | 0.60 |
   | 1997.9 / 0.87 (paper) | 289.6 | −7.8 | 116 | 233 | 55 | 200 | 310 | 385 | 34 | 70 | 2.42 | 0.65 | 0.142 | 252 | 0.70 |

   With little water vapour longwave absorption (`b` ≤ 250) the climate is very insensitive (0.2–0.4 K per W/m²,
   even at 300 K), but surface LW down and precipitation are far too low at Earth-like temperatures. The water vapour
   longwave emission near the surface is what drives them. The default `b` with the paper's `a` = 0.87 has an
   Earth-like sensitivity (0.70 K per W/m²), but it is only ≈3 K below the regime where `a` = 1.1–1.3 runs away
   (findings 17–18), and its shortwave budget is off (atm SW 55, OSR 116). `b` = 500 with `a` ≈ 1.75 interpolates
   to ≈292 K with a moderate sensitivity, and is refined next (`Z_*`), also with visible `wv` = 45 for more
   atmospheric absorption and `cloud_albedo` = 0.40 for OSR. With much sea ice, the steady net TOA becomes strongly
   negative (−17.5 W/m² at ice 0.255): the sea-ice energy source of finding 9 scales with the ice cover.
20. **`b` = 500–1000 is stable, the default `b` with more shortwave absorption tips over** (`Z_*`, GPU, 40 years,
   `wv_nir` = 1500; means over years 35–39; K per W/m² from the twin over years 20–39). 27 of 28 runs are stable.
   The exception is the +2% solar twin of `b` / `a` = 1997.9 / 0.87 with `wv` = 45, which runs away from 295.2 K
   to a NaN in year 15. Its neighbour, 1997.9 / 0.75 with `wv` = 45, warms by 10.5 K in its twin.

   | `b` / `a`, other | SST | net TOA | OSR | OLR | atm SW | sfc SW down | LW down | LW up | LH | P | cloud | ice | K per W/m² |
   |---|---|---|---|---|---|---|---|---|---|---|---|---|---|
   | observed (all-sky) | ≈291 | 0 | 100 | 239 | 79 | 185 | 342 | 397 | 84 | 2.9 | | | ≈0.8 |
   | 500 / 1.6 | 288.5 | −7.9 | 120 | 230 | 54 | 196 | 267 | 377 | 44 | 1.53 | 0.67 | 0.149 | 0.32 |
   | 500 / 1.75 | 291.1 | −5.2 | 114 | 232 | 60 | 193 | 298 | 394 | 54 | 1.87 | 0.66 | 0.105 | 0.70 |
   | 500 / 1.9 | 294.8 | −2.0 | 102 | 241 | 71 | 191 | 343 | 421 | 70 | 2.40 | 0.60 | 0.061 | 0.40 |
   | 500 / 1.6, `wv` 45 | 289.8 | −6.8 | 112 | 236 | 73 | 183 | 287 | 386 | 46 | 1.60 | 0.65 | 0.130 | 1.02 |
   | 500 / 1.75, `wv` 45 | 294.9 | −2.3 | 91 | 252 | 103 | 169 | 353 | 423 | 63 | 2.15 | 0.55 | 0.065 | 0.43 |
   | 500 / 1.75, `wv` 45, `cloud_albedo` 0.40 | 296.6 | −1.1 | 84 | 258 | 112 | 166 | 369 | 434 | 66 | 2.24 | 0.54 | 0.046 | 0.65 |
   | 1000 / 1.2 | 287.8 | −10.6 | 121 | 231 | 52 | 200 | 272 | 372 | 53 | 1.83 | 0.66 | 0.180 | 0.40 |
   | 1000 / 1.4 | 292.1 | −4.4 | 111 | 235 | 62 | 194 | 325 | 402 | 69 | 2.38 | 0.64 | 0.100 | 0.67 |
   | 1000 / 1.2, `wv` 45 | 289.3 | −8.3 | 111 | 239 | 71 | 189 | 299 | 384 | 59 | 2.03 | 0.62 | 0.153 | 1.08 |
   | 1000 / 1.2, `wv` 45, `cloud_albedo` 0.40 | 294.2 | −3.0 | 93 | 252 | 95 | 176 | 355 | 418 | 71 | 2.44 | 0.58 | 0.075 | 0.95 |
   | 1997.9 / 0.75, `wv` 45 | 288.7 | −10.0 | 114 | 238 | 66 | 192 | 303 | 379 | 65 | 2.25 | 0.63 | 0.166 | 2.18 |
   | 1997.9 / 0.87, `wv` 45 | 295.2 | −1.6 | 98 | 245 | 94 | 170 | 379 | 425 | 79 | 2.71 | 0.59 | 0.058 | twin NaN |

   - **The sensitivity depends on the state, not on noise.** The interannual standard deviation of the annual SST
     is 0.02–0.12 K. Along `b` = 500 the sensitivity is 0.3 K per W/m² at 288 K, 0.7 at 291 K and 0.4 at 295 K.
     It peaks where sea ice (0.10–0.15) and cloud cover (0.66 → 0.60) retreat with warming: from `a` = 1.75 to
     1.9, OSR falls by 12 W/m² for 3.7 K. Earth-like temperatures lie in this sensitive range.
   - **Visible water vapour absorption is the riskier lever.** `wv` = 45 adds 13–19 W/m² of atmospheric
     absorption, but triples the sensitivity in the cooler states (500 / 1.6: 0.32 → 1.02; 1000 / 1.2: 0.40 → 1.08)
     and tips the default `b` over. This is the shortwave absorption that is linear in q (finding 18). More
     atmospheric absorption should come from the dry-air absorptivity, which has no q dependence.
   - **The longwave has no cloud effect.** `ByrneOGormanLongwaveTransmissivity` is clear-sky, but the observed LW
     down of 342 W/m² and OLR of 239 W/m² are all-sky values. Clear-sky values are ≈314 and ≈266 W/m² (Wild et
     al. 2013). At an Earth-like surface temperature (LW up ≈ 395, `b` / `a` = 500 / 1.75) the greenhouse effect
     LW up − OLR is 162 W/m² against the observed 158. LW down (298) is 16 W/m² below the clear-sky value. The clear-sky longwave therefore
     emulates the cloud longwave effect at the TOA, but can't reach the all-sky LW down at the surface. The
     remaining biases at Earth-like temperatures are in the shortwave: OSR 114 instead of 100, atmospheric
     absorption 60 instead of 79.
   - Next: `V_*` (CPU, 40 years with +2% twins) adds `dry` = 0.066 / 0.1 and `cloud_albedo` = 0.38 to
     `b` / `a` = 500 / 1.6–1.8 and 750 / 1.4–1.5. `L3_*` (CPU, 1000 years) brackets the `b` = 500 family at
     `a` = 1.75 and 1.9, each with a +4% solar twin (≈+9.6 W/m²) as stress test.
21. **More dry-air shortwave absorption or a lower cloud albedo brings the budget close, but the climate jumps
   between ≈291 and ≈294 K** (`V_*`, CPU, 40 years, all 20 runs and twins stable; `a` given at 280 ppm and, in
   brackets, at the new default of 360 ppm, see decision 3).

   | `b` / `a` (360 ppm), `dry`, `cloud_albedo` | SST | OSR | OLR | atm SW | sfc SW down | LW down | LW up | LH | P | cloud | ice | K per W/m² |
   |---|---|---|---|---|---|---|---|---|---|---|---|---|
   | observed (all-sky) | ≈291 | 100 | 239 | 79 | 185 | 342 | 397 | 84 | 2.9 | | | ≈0.8 |
   | 500 / 1.7 (1.65), 0.066, 0.43 | 290.7 | 112 | 235 | 70 | 185 | 294 | 391 | 49 | 1.69 | 0.67 | 0.110 | 0.91 |
   | 500 / 1.8 (1.75), 0.066, 0.43 | 294.5¹ | 98 | 246 | 82 | 184 | 342 | 419 | 66 | 2.29 | 0.59 | 0.069 | 0.46 |
   | 500 / 1.6 (1.55), 0.033, 0.38 | 293.1 | 99 | 247 | 68 | 201 | 319 | 409 | 66 | 2.28 | 0.61 | 0.086 | 0.49 |
   | 500 / 1.7 (1.65), 0.033, 0.38 | 294.9 | 95 | 249 | 73 | 197 | 340 | 421 | 72 | 2.46 | 0.59 | 0.065 | 0.41 |
   | 500 / 1.6 (1.55), 0.066, 0.38 | 294.6 | 91 | 254 | 84 | 190 | 338 | 419 | 68 | 2.35 | 0.59 | 0.072 | 0.48 |
   | 500 / 1.7 (1.65), 0.1, 0.38 | 297.5 | 81 | 261 | 104 | 177 | 372 | 439 | 72 | 2.45 | 0.56 | 0.043 | 0.39 |
   | 750 / 1.4 (1.35), 0.066, 0.38 | 295.6 | 91 | 253 | 86 | 188 | 357 | 426 | 75 | 2.59 | 0.59 | 0.065 | 0.37 |

   ¹ still warming by 0.19 K per decade (decadal means 292.5, 293.2, 294.1, 294.5 K).

   - `dry` 0.033 → 0.066 adds 12–16 W/m² of atmospheric absorption and takes 10–15 W/m² from the surface, but
     also lowers OSR by ≈8 W/m² (absorption above the cloud tops), so the same `a` gives a warmer climate.
   - From `a` = 1.7 to 1.8 (`dry` 0.066, `cloud_albedo` 0.43) SST jumps by 3.8 K and OSR falls by 14 W/m², as
     cloud cover drops from 0.67 to 0.59 and sea ice from 0.11 to 0.07. The observed state (LW up 397, OSR 100)
     falls right into this jump. The diagnostic cloud cover, which falls by ≈2% per K, is the main shortwave
     feedback.
   - In the warm branch the budget is close to the observed one, except for the surface temperature. With 500 /
     1.8 / `dry` 0.066, OSR, atmospheric absorption, surface shortwave down and LW down are all within
     3 W/m² of the observations, but LW up is 22 W/m² (≈3 K) too high and the latent heat flux 18 W/m² too low.
   - All sensitivities are 0.37–0.59 K per W/m², except 0.91 in the cool branch.
   - Next: `W_*` (GPU, 360 ppm): `cloud_albedo` 0.40 with `dry` 0.066 and `a` = 1.40–1.70, the jump at
     `cloud_albedo` 0.43 (`a` = 1.68 / 1.72), `dry` 0.08, and `b` = 750.
22. **`b` = 750 with `a` = 1.45 (360 ppm), `dry` 0.066 and the default `cloud_albedo` 0.43 matches the observed
   budget best** (`W_*`, GPU, 40 years, all 25 runs and twins stable; means over years 35–39, K per W/m² from the
   +2% solar twin over years 20–39; all `wv_nir` = 1500). The check run with `a` = 1.6492 at 360 ppm reproduces
   `V_b500_a1.7_d066` (290.74 vs 290.7 K), so the 280 → 360 ppm translation of decision 3 is exact.

   | `b` / `a` (360 ppm), `dry`, `cloud_albedo` | SST | OSR | OLR | atm SW | sfc SW down | LW down | LW up | LH | P | cloud | ice | K per W/m² |
   |---|---|---|---|---|---|---|---|---|---|---|---|---|
   | observed | ≈291 | 100 | 239 | 79 | 185 | 314 clear / 342 all-sky | 397 | 84 | 2.9 | | | ≈0.8 |
   | 500 / 1.40, 0.066, 0.40 | 289.1 | 108.6 | 241.0 | 68.2 | 193.0 | 272.6 | 380.7 | 45.0 | 1.56 | 0.66 | 0.151 | 0.36 |
   | 500 / 1.47, 0.066, 0.40 | 290.1 | 107.3 | 240.8 | 69.9 | 190.8 | 284.2 | 387.1 | 48.0 | 1.66 | 0.66 | 0.127 | 0.54 |
   | 500 / 1.55, 0.066, 0.40 | 291.7 | 103.5 | 242.8 | 74.1 | 189.0 | 304.2 | 398.1 | 54.5 | 1.88 | 0.65 | 0.103 | 0.78 |
   | 500 / 1.62, 0.066, 0.40 | 294.3 | 93.7 | 250.4 | 82.9 | 188.0 | 337.5 | 417.7 | 66.9 | 2.31 | 0.59 | 0.072 | 0.49 |
   | 500 / 1.70, 0.066, 0.40 | 295.7 | 91.4 | 251.8 | 86.6 | 185.7 | 352.5 | 427.2 | 70.2 | 2.42 | 0.58 | 0.058 | 0.41 |
   | 500 / 1.68, 0.066, 0.43 | 291.3 | 110.2 | 235.9 | 71.3 | 184.4 | 301.5 | 395.6 | 50.7 | 1.75 | 0.66 | 0.100 | 0.89 |
   | 500 / 1.72, 0.066, 0.43 | 292.6 | 105.7 | 239.3 | 75.4 | 183.6 | 317.9 | 404.6 | 57.3 | 1.98 | 0.64 | 0.083 | 0.77 |
   | 500 / 1.60, 0.08, 0.43 | 290.4 | 110.0 | 237.1 | 73.7 | 182.5 | 291.0 | 389.1 | 46.3 | 1.60 | 0.67 | 0.114 | 0.69¹ |
   | 500 / 1.70, 0.08, 0.43 | 292.6 | 104.1 | 241.0 | 79.9 | 180.4 | 318.3 | 404.6 | 55.5 | 1.92 | 0.64 | 0.086 | 0.84 |
   | **750 / 1.45, 0.066, 0.43** | **291.5** | 108.6 | 237.9 | 71.8 | 185.8 | 313.7 | 397.3 | 59.8 | 2.07 | 0.65 | 0.105 | 0.92 |
   | 750 / 1.55, 0.066, 0.43 | 295.3 | 98.2 | 245.3 | 83.2 | 182.2 | 357.6 | 424.7 | 72.7 | 2.50 | 0.59 | 0.064 | 0.41 |
   | 750 / 1.35, 0.066, 0.40 | 293.1² | 98.0 | 247.4 | 78.1 | 190.0 | 331.1 | 409.1 | 68.1 | 2.35 | 0.61 | 0.092 | 0.71 |

   ¹ the twin is still warming (decadal means 292.4, 293.2, 293.4, 294.1 K). ² still warming by ≈0.2 K per decade.

   - **Two branches, nothing in between.** The cool branch (cloud cover 0.64–0.67, sea ice 0.08–0.15) reaches up
     to ≈292.6 K, the warm branch (cloud cover 0.58–0.59, sea ice ≤ 0.07) starts at ≈294 K. The jump lies at
     `a` ≈ 1.55–1.62 for `b` = 500 with `cloud_albedo` 0.40, at ≈1.72–1.75 with 0.43, and at ≈1.45–1.55 for
     `b` = 750. The observed surface temperature (LW up 397) is on the cool branch, about 1 K below its upper end.
   - **Every +2% solar twin of a cool-branch run ends on the warm branch** (295.5–296.2 K). That gives the cool
     branch its 0.7–0.9 K per W/m², while warm-branch runs warm by ≈2 K (≈0.4 K per W/m²). The observed climate
     sensitivity, ≈3 K for 3.7 W/m², is ≈0.8 K per W/m², so the cool branch is realistic in this respect too.
   - **750 / 1.45 / `dry` 0.066 / `cloud_albedo` 0.43** has LW up, LW down (clear sky), OLR and surface SW
     down all within 1–2 W/m² of the observations, steady over the 40 years (decadal SST 291.5, 291.5, 291.6,
     291.5 K). It is off in OSR (+9 W/m²), atmospheric SW (−7 W/m²) and the hydrological cycle (LH 60 vs 84 W/m²,
     P 2.07 vs 2.9 mm/day), like every run. Part of the OSR excess is surface reflection (SW up ≈30 vs 24 W/m²),
     and the net TOA of −5.3 W/m² is the sea-ice energy source (finding 9).
   - **500 / 1.55 / `cloud_albedo` 0.40** is the alternative. Its OSR is closer (103.5 W/m²), but LW down is
     10 W/m² too low and OLR 4 W/m² too high, P is lower (1.88 mm/day), and it changes one more default.
   - At the same SST a larger `b` gives more LW down and a stronger hydrological cycle (750 / 1.45 vs 500 / 1.68
     at 291.5 / 291.3 K: LW down 313.7 vs 301.5 W/m², P 2.07 vs 1.75 mm/day).
   - `L3_*` after 86 years (CPU, 280 ppm, `b` = 500): `a` = 1.75 (1.70 at 360 ppm) stays at 291.1 K with decadal
     means within 0.1 K, its +4% twin at 296.5 K; `a` = 1.9 stays at 294.8 K, its twin at 298.7 K. All stable.

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

Done (uncommitted, worktree `tree_defaults` on branch `mg/tune-water-vapour-defaults` from `origin/main`
`f6c4d99f9`, because `mg/tune-water-vapour` is still based on `caed82d6`):

- `DEFAULT_CO2` = 360 ppm in `greenhouse_gases.jl` (was 280), `TwoTimesCO2` / `FourTimesCO2` docstrings refer to
  the default instead of preindustrial CO₂. `ExponentialCO2` keeps its preindustrial base of 280 ppm in 1850,
  part of its fit to the Keeling curve (≈368 ppm in 2000).
- `ByrneOGormanLongwaveTransmissivity.reference_co2 = DEFAULT_CO2` (was the literal 360).
- NumericalRadiation extension: the `AnalyticBandLongwave` fallback without a CO₂ component is `DEFAULT_CO2`
  (was its own 280). `ClearSkyEcCKDRadiation` keeps its own `mole_fractions.CO₂` default (NumericalRadiation's).
- Tests: `numericalradiation.jl` uses `DEFAULT_CO2` instead of 280. `longwave_radiation.jl` checks that the dry
  optical depth equals `dry_absorption` at the default CO₂. `greenhouse_gases.jl`, `longwave_radiation.jl`
  and `numericalradiation.jl` pass with `--check-bounds=yes`.
- `docs/src/radiation.md` and `CHANGELOG.md` (PR number still to be filled in).
- Decision 4: `dry_absorption` = 1.45 and `water_vapor_absorption` = 750 in `longwave_transmissivity.jl`, with a
  comment on why Byrne and O'Gorman's values are not used. `absorptivity_dry_air` = 0.066 and
  `absorptivity_water_vapor_near_infrared` = 1500 in `TwoBandShortwaveTransmissivity`
  (`shortwave_transmissivity.jl`), with comments. The `ByrneOGorman` paragraph of `docs/src/radiation.md` now
  gives the tuning target and the resulting budget. Second `CHANGELOG.md` entry.
- `SpeedyWeather` is already at `0.25.0-DEV`, so no version bump is needed for a change of default values.

## Testing and verification

- Unit tests with the new defaults and `--check-bounds=yes` (2026-11-17, `tree_defaults`, identical to commit
  `17f3f7987` on `mg/tune-water-vapour`): `longwave_radiation` 21, `shortwave_radiation` 3267, `radiation` 82,
  `greenhouse_gases` 13, `albedo` 259, `numericalradiation` 26, `all_parametrizations` 61, all pass.
- `test/long_integrations/default_primitive_wet.jl` (1 year, `--check-bounds=yes`): passes.
- The 1000-year T32L8 runs `D_*` (experiment log) are the main verification. Their annual means will be
  summarised here.

## Documentation changes

- `docs/src/radiation.md`: the `TwoBandShortwave` section states that the absorptivities are speedy.f90's values.
  If they change, document the new values and the tuning target.

## Known limitations

- Global means are from instantaneous samples every 6 hours, not time-step averages. The diurnal cycle is sampled
  at 4 fixed UTC times, which is adequate for global means.
- The Byrne–O'Gorman longwave is clear-sky. It can match the all-sky OLR only by being more opaque than the clear
  sky, and it can't match the all-sky surface LW down at an Earth-like surface temperature (finding 20).

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
3. **Default CO₂ of 360 ppm** (2026-11-17, prompts: "we might also change the co2 ppm default throughout the
   model. By default we start in year 2000, so picking 360ppm seems right. Can you apply that througout the model
   and then incorporate that into the tuning", "we should reuse the DEFAULT_CO2 throughout the model and set this
   one to 360ppm"). `DEFAULT_CO2` = 360 (was preindustrial 280). `ByrneOGormanLongwaveTransmissivity.reference_co2`
   and the NumericalRadiation extension's fallback reuse it. `dry_absorption` is therefore the dry optical depth at
   the default CO₂. CO₂ enters only `ByrneOGormanLongwaveTransmissivity` (and NumericalRadiation's schemes), so
   all earlier runs translate exactly: `a`(360 ppm) = `a`(280 ppm) − 0.2023 ln(360/280) = `a`(280 ppm) − 0.0508.
   From `W_*` on, `a` is given at 360 ppm.
4. **New defaults** (2026-11-17, finding 22): `ByrneOGormanLongwaveTransmissivity` `dry_absorption` 1.3 → 1.45
   and `water_vapor_absorption` 1997.9 → 750; `TwoBandShortwaveTransmissivity` `absorptivity_dry_air`
   0.033 → 0.066 and `absorptivity_water_vapor_near_infrared` 15000 → 1500. `cloud_albedo` stays 0.43. Reasons:
   the closest match to the observed surface budget on the cool branch, steady over 40 years, the strongest
   hydrological cycle of the candidates, and no change to the clouds. Pending the 1000-year validation `D_*`.

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
| `M_*` (Slurm 3283731, CPU, `tree_1279` = `b06af911a`) | Byrne–O'Gorman default; `a` 1.1 / 1.5; `b` 1500 / 2500; `wv_nir` 1500; `a` 1.1 + `wv_nir` 1500; `strat_emis` 0.1 | 30 years | running | default runs away (finding 17) |
| `G_*` (GPU Slurm 3283750) and `M_b*`, `M_*_S` (CPU Slurm 3283751), `tree_1279` | `b` / `a` = 1000 / 2.5, 1000 / 3.0, 500 / 4.0, 500 / 4.5, 0 / 5.4, each with a +2% solar twin (`_S`); `_S` twins for the default, `a` = 1.1 and `b` = 1500 | 30 years | running | feedback strength, warming per W/m² |
| `Y_*` (GPU Slurm 3284305 + 3284306, `tree_1279`) | all `wv_nir` = 1500: `b` = 0 with `a` 2.0 / 2.5 / 3.0 / 3.5, `b` = 250 with `a` 2.0 / 2.5 / 3.0, `b` = 500 with `a` 1.5 / 2.0 / 2.5, `b` = 1997.9 with `a` 0.7 / 0.87; each with a +2% solar twin | 40 years | completed, all stable | finding 19 |
| `Z_*` (GPU Slurm 3289286 + 3289287, `tree_1279`) | all `wv_nir` = 1500: `b` = 500 with `a` 1.6 / 1.75 / 1.9; `b` = 1000 with `a` 1.0 / 1.2 / 1.4; with visible `wv` = 45: `b` / `a` 500 / 1.6, 500 / 1.75, 1000 / 1.0, 1000 / 1.2, 1997.9 / 0.75, 1997.9 / 0.87; `cloud_albedo` 0.40 with 500 / 1.75 / `wv` 45 and 1000 / 1.2 / `wv` 45; each with a +2% solar twin | 40 years | completed; the twin of 1997.9 / 0.87 / `wv` 45 NaN in year 15 | finding 20 |
| `V_*` (CPU Slurm 3289952, `tree_1279`, `long_cpu.sh`) | all `wv_nir` = 1500: `b` = 500 with `a` 1.6 / 1.7 / 1.8 and `dry` 0.066 + `cloud_albedo` 0.38; `a` 1.7 / 1.8 with `dry` 0.066; `a` 1.6 / 1.7 with `cloud_albedo` 0.38; `a` 1.7 with `dry` 0.1 + `cloud_albedo` 0.38; `b` = 750 with `a` 1.4 / 1.5, `dry` 0.066 + `cloud_albedo` 0.38; each with a +2% solar twin | 40 years | completed, all stable | finding 21 |
| `W_*` (GPU Slurm 3297308 + 3297309, `tree_defaults` = `f6c4d99f9` + `DEFAULT_CO2` = 360 ppm, `a` at 360 ppm) | all `wv_nir` = 1500: `b` = 500, `dry` 0.066, `cloud_albedo` 0.40 with `a` 1.40 / 1.47 / 1.55 / 1.62 / 1.70; `cloud_albedo` 0.43 with `a` 1.68 / 1.72 (`dry` 0.066) and 1.60 / 1.70 (`dry` 0.08); `b` = 750, `dry` 0.066 with `a` 1.45 / 1.55 (0.43) and 1.35 (0.40); each with a +2% solar twin; plus `a` = 1.6492 as a check of the 280 → 360 ppm translation (≡ `V_b500_a1.7_d066`) | 40 years | completed, all stable | finding 22 |
| `L3_*` (CPU Slurm 3289953, `tree_1279`) | `b` / `a` = 500 / 1.75 and 500 / 1.9 (`wv_nir` = 1500), each with a +4% solar twin (1419.6 W/m²) | 1000 years, 25-year restart chunks | running, ≈700 years per day; year 86 at 15:15, all stable | does the `b` = 500 family stay stable for centuries, also under strong forcing? Resubmit `configs/l3.txt` after the 1-day limit (≈12:15 on 2026-11-18) |
| `D_*` (GPU Slurm 3298103; `D_default_2xCO2` alone on GPU Slurm 3298179, `configs/d1_co2.txt`; `tree_defaults` with the decision-4 defaults) | code defaults; +4% solar (1419.6 W/m²); −4% solar (1310.4 W/m²); 2×CO₂ (`greenhouse_gases = (; co2 = CO2(spectral_grid, 720))`) | 1000 years, 25-year restart chunks | running since 15:21 (H100, ≈9 h); 2×CO₂ since 15:31 (H100, ≈5 h) | validation of the new defaults |
| `X_*` (CPU Slurm 3283769, `tree_1279`) | all `wv_nir` = 1500: `b` / `a` = 1000 / 2.0, 750 / 3.0, 750 / 3.5, 250 / 4.5, 250 / 5.0, 0 / 4.8, 0 / 6.0, and 500 / 4.0 with `cloud_albedo` 0.40, each with a +2% solar twin | 30 years | running, ready ≈09:00 | fills in `b` ≤ 1000. At year 7 the `M_*` runs show `b` ≥ 1500 still warming: default SST 298.9 K at +4.5 W/m²; `a` = 1.1 and `b` = 1500 at 294.0 / 293.6 K, still rising ≈0.3 K per year. With Byrne–O'Gorman the stratosphere stays cool (T₁ ≈ 221 K, max 252 K); `wv_nir` = 1500 halves q₂ |
| `r2_*` (Slurm 3279654, ed) | all `wv_nir` = 1500: `wv` 55/60/65 × `τ₀` 7.0/7.2/7.5, polar `τ₀` 2.0, `cloud_albedo` 0.45 | 20 years | running | round 2 around `c7` |
| `L_main`, `L_nir1500`, `L_pt250` (GPU Slurm 3279622, ed) | `main` defaults; `wv_nir` = 1500; `water_vapor_pressure_top` = 250 hPa | 1000 years, 25-year restart chunks | running (H100) | does `main` blow up, and when? Do both stratosphere fixes hold for a millennium? ≈2.1 s per simulated month per run, done ≈03:00 on 2026-11-17 |
| `L2_wv60_t75`, `L2_wv60_t72`, `L2_wv60_t70`, `L2_wv65_t72` (GPU Slurm 3279736, 5d) | all `wv_nir` = 1500: `wv` 60 with `τ₀` 7.5/1.875, 7.2/1.8, 7.0/1.75, and `wv` 65 with `τ₀` 7.2/1.8 (the `r2` grid's centre, see its year 1) | 1000 years, 25-year restart chunks | running (H100), expected ≈05:00–07:00 on 2026-11-17 | validation of the final candidates. `r2_*` (same parameters, 20 y on CPU) picks the best budget match |
