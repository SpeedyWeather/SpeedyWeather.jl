# Humidity-dependent longwave transmissivity (Byrne and O'Gorman, 2013)

> Status: **in progress**. Implemented and unit-tested, fₗ re-evaluated on current `main`; awaiting review.

Date of initial draft: 2026-10-30

Base revision: 28069552f493ae7122c90588304b6b7ea427c743

## Originating prompt

> Have a look at PR1279 and compare to main, many PRs have now been merged. Can you remove all commits
> that qualitatively have been merged and reduce this PR to only commits that aren't in other PRs, so take
> out changes to the shortwave as that's in PR1262. But I think the changes you made here towards a better
> longwave transmissivity aren't elsewhere, then update PR description and title to reflect what this PR
> is about. It used to be more like an experimental draft to play with but now I don't want to lose the
> advances that are in here and want to evaluate them more in isolation

## Revision log

- 2026-10-30: #1279 started as an integration branch to pretest #1183, #1221, #1262 and #1277 together.
  #1183, #1221 and #1277 are now merged, and #1262 (shortwave, ozone, clouds, stratospheric longwave
  emission) is still open. The branch was rebuilt on `main` and keeps only the longwave transmissivity
  changes. The provisional shortwave absorptivities, ozone, ocean roughness and initial-humidity tweaks
  were dropped. They are covered by #1262 and #1284, or no longer apply after the surface-roughness
  changes on `main`.
- 2026-10-30: Byrne–O'Gorman CO₂ hooked up to `model.greenhouse_gases.co2`. Following Isca, the dry
  optical depth is now logarithmic in CO₂, `a + 0.2023·log(CO₂/360 ppm)`; this replaces the `co2_scaling`
  multiplier. The reference pressure is taken from `model.atmosphere`. The primitive equation models now
  default to `greenhouse_gases = (; co2 = CO2(spectral_grid))` at 280 ppm. Re-evaluated `fₗ` against the
  current `main` physics (see below) and kept 0.25.

## Problem description

`FriersonLongwaveTransmissivity` prescribes the optical depth as a function of latitude and σ only. With
its default `fₗ = 0.1`, `PrimitiveWetModel` has a weak greenhouse effect and a troposphere that is too
cold (~284 K near-surface air). The optical depth also does not respond to the model's humidity, so there
is no water vapor feedback.

## Background

Byrne and O'Gorman (2013, J. Climate, doi:10.1175/JCLI-D-12-00262.1) use a grey optical depth per layer
of `dτ = (a μ + b q) Δp / p₀`, as in Isca. Here `a` is the dry absorption (well-mixed greenhouse gases,
scaled by `μ` for CO₂ changes), `b` is the water vapor absorption and `q` is the specific humidity.

## Summary of changes

- New `ByrneOGormanLongwaveTransmissivity`, `dτ = (a + c·log(CO₂/CO₂ᵣ) + b·q)·Δp/p₀`, with
  `dry_absorption = 0.8678`, `co2_sensitivity = 0.2023`, `reference_co2 = 360` ppm and
  `water_vapor_absorption = 1997.9`. CO₂ comes from `model.greenhouse_gases.co2` (if there is none,
  `reference_co2` is used) and `p₀` from `model.atmosphere.reference_pressure`. Dry models use only the dry
  optical depth.
- Default `greenhouse_gases = (; co2 = CO2(spectral_grid))` (280 ppm) in `PrimitiveWetModel` and
  `PrimitiveDryModel`.
- Default for `PrimitiveWetModel`: `FriersonLongwaveTransmissivity(SG, fₗ = 0.25)`. Re-evaluated below.

### Re-evaluation of fₗ on current `main` (2026-10-30)

T31, 8 layers, default `PrimitiveWetModel` from `main` at 28069552 plus this PR. 360-day spin-up, then the
mean of 36 samples, one every 10 days, over the second year. Script: `benchmarks/pr1279/fl_sweep.jl`
(outside the repo). The temperature is the lowest model layer; the normalized greenhouse effect is
g = 1 − OLR/σTₛ⁴.

| case | T lowest layer [K] | T 500 hPa [K] | T top [K] | OLR [W/m²] | OSR [W/m²] | PW [kg/m²] | precip [mm/day] | g |
|---|---|---|---|---|---|---|---|---|
| Frierson fₗ=0.1 | 272.6 | 244.5 | 231.6 | 200.4 | 159.0 | 8.7 | 1.81 | 0.36 |
| Frierson fₗ=0.15 | 274.4 | 247.6 | 224.5 | 202.0 | 155.7 | 10.5 | 1.90 | 0.37 |
| Frierson fₗ=0.2 | 276.1 | 250.6 | 220.7 | 203.5 | 151.6 | 12.5 | 1.97 | 0.38 |
| Frierson fₗ=0.25 | 277.9 | 253.4 | 218.7 | 205.2 | 146.1 | 14.8 | 2.01 | 0.39 |
| Byrne–O'Gorman | 272.4 | 249.9 | 222.3 | 218.3 | 154.2 | 10.9 | 1.62 | 0.30 |

Conclusions:
- Since the earlier estimate, other processes have changed on `main`. All cases are now ~10 K colder than
  Earth and too dry (Earth: ~288 K, ~25 kg/m², ~2.7 mm/day).
- The cold bias is mainly a shortwave problem. Reflected shortwave is ~146–159 W/m², a planetary albedo of
  ~0.45 against Earth's ~0.3. OLR (~200–205 W/m²) is low because the atmosphere is cold, not because the
  greenhouse effect is too strong.
- With fₗ = 0.25, g = 0.39 is close to Earth's ≈ 0.40, while 0.1 gives 0.36. fₗ = 0.25 stays the default.
  Raising fₗ further to warm the model would compensate a shortwave error with a longwave one.
- The fₗ = 0.25 run is still ~11 W/m² short of TOA balance (an estimated ~340 W/m² incoming minus
  146 W/m² OSR minus 205 W/m² OLR). The 1-year spin-up may not be long enough to equilibrate.
- Byrne–O'Gorman has a weaker greenhouse effect here (g = 0.30) because the column is dry: its water vapor
  absorption is calibrated for a moister climate. It should be re-evaluated once the shortwave bias is
  addressed.

## Testing and verification

- `ByrneOGormanLongwaveTransmissivity` is added to the generic transmissivity testset.
- New testset: with zero humidity the column transmissivity is `exp(-dry_absorption)`. With humidity it is
  lower in every layer and matches the analytic column optical depth.
- Tests cover the logarithmic CO₂ dependence (2×CO₂ adds `0.2023·log 2`), the default CO₂ of 280 ppm and
  models without greenhouse gases.
- Climatology comparison: see the table above.

## Documentation changes

Docstrings; the greenhouse gas section in `docs/src/radiation.md` now says which schemes use CO₂.

## Known limitations

- `fₗ = 0.25` is chosen to match the normalized greenhouse effect, not calibrated against observations,
  and it should be revisited together with #1262.
- CO₂ is read from a `Ref` inside the column kernel, like the NumericalRadiation extension does; this is
  not tested on GPU yet.
- `FriersonLongwaveTransmissivity` does not use CO₂.
- `ByrneOGormanLongwaveTransmissivity` is not the default yet.

## Future work

- Decide on the default transmissivity after the isolated evaluation.
- Add a section to the radiation docs.
