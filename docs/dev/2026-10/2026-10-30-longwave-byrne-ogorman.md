# Humidity-dependent longwave transmissivity (Byrne and O'Gorman, 2013)

> Status: **in progress**. Implemented and unit-tested; the climate impact still has to be checked in isolation.

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

- New `ByrneOGormanLongwaveTransmissivity` with parameters `dry_absorption = 0.8678`,
  `water_vapor_absorption = 1997.9`, `co2_scaling = 1` and `reference_pressure = 1e5`. Dry models use only
  the dry optical depth.
- Provisional default for `PrimitiveWetModel`: `FriersonLongwaveTransmissivity(SG, fₗ = 0.25)`. This gives
  ~288 K and outgoing longwave of ~240 W/m², but it is marked TODO for proper calibration.

## Testing and verification

- `ByrneOGormanLongwaveTransmissivity` is added to the generic transmissivity testset.
- New testset: with zero humidity the column transmissivity is `exp(-dry_absorption)`. With humidity it is
  lower in every layer and matches the analytic column optical depth.
- Still to do: compare climatologies (temperature, OLR, humidity) for Frierson `fₗ = 0.1` vs `0.25` vs
  Byrne–O'Gorman, separately from #1262.

## Documentation changes

Docstrings only.

## Known limitations

- The `fₗ = 0.25` default is a provisional tuning, not a calibration.
- `ByrneOGormanLongwaveTransmissivity` is not the default yet.

## Future work

- Decide on the default transmissivity after the isolated evaluation.
- Add a section to the radiation docs.
