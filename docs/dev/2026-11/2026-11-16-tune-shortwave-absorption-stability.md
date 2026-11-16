# Tuning shortwave water vapour and dry-air absorption for stable long T32L8 runs

> Status: **in progress**. Reproducing the long-run instability with 30-year diagnostic runs on `main` and
> before the two radiation PRs. No source code changed yet.

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
6. **Caveat on the direction of the tuning.** With a grey longwave at fixed optical depth, moving shortwave
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

### Phase 3: validate candidates

- 100-year runs of the 2–3 best candidates (Slurm, `standard` partition).
- One 1000-year run of the chosen configuration. At ≈600 SYPD, a millennium takes ≈40 h, which exceeds the 1-day
  `priority`/`short` QOS, so it runs on the `medium` QOS (7 days) or is chained with restart files.
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

- If tuning the absorptivities alone can't remove the cold drift, the next candidates are cloud albedo (0.43), the
  Frierson optical depths, and the sensible heat flux (twice the observed value). These need a decision on scope
  (see open questions).

## Open questions

1. **Tuning target:** observed global budgets (atmospheric absorption ≈79 W/m²) or speedy.f90's budget (≈62 W/m²)?
2. **Scope if the absorptivities aren't enough:** may other parameters (cloud albedo, Frierson τ₀, surface fluxes)
   be tuned in the same PR, or should that be a follow-up?

## Experiment log

| Label | Revision / config | Length | Status | Notes |
|---|---|---|---|---|
| `smoke` | `main` | 60 days | completed | script check, net TOA −14 to −19 W/m² |
| `main_default` | `main` defaults | 30 years | running | |
| `pre_albedo` | `3245c69f` | 30 years | running | |
| `pre_radiation` | `b43dc87c` | 30 years | running | |
