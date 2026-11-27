# Entrainment as a rate per unit σ in Betts-Miller convection

> Status: **completed**. Entrainment mixes a fraction per σ distance, making the parcel dilution independent of the number of layers.

Date of initial draft: 2026-11-22

Base revision: 32fbfa49d5c9e244340109cdb403e8642b3f90f2

## Originating prompt

> Make a separate PR for convection and clarify this in the issue.

## Revision log

## Problem description

`BettsMillerConvection` mixes the rising parcel with environmental air by a fixed fraction ε(σ) per
model layer (#1426). The fraction of the original parcel remaining is Π(1 - ε_k), so dilution grows
with the number of layers: ~40% remains at σ = 0.5 at L8, ~e⁻¹⁰ at L100, which switches off deep
convection at high vertical resolution.

## Background

Suggestion 1 of #1426: entrainment as a rate per unit σ, calibrated to the current scheme at Δσ_ref ≈ 1/8.

## Summary of changes

- New method `(E::AbstractEntrainment)(σ, Δσ) = 1 - (1 - E(σ))^(Δσ/Δσ_ref)`: `E(σ)` keeps its meaning,
  the fraction mixed over the reference σ distance. The remaining fraction multiplies to the same over
  many thin or few thick layers. It equals the old per-layer ε for Δσ = Δσ_ref.
- New parameter `reference_layer_thickness = 1/8` in `LinearEntrainment` and `ConstantEntrainment`.
- `pseudo_adiabat!` and `dry_adiabat!` call `entrainment(σ[k], σ[k+1] - σ[k])`, the σ distance the parcel
  rose. No new arguments are needed, σ is already passed.

## Testing and verification

Functor tests: reference distance gives ε, zero distance gives 0, n thin steps dilute as one thick step.
Existing convection tests still pass.

## Documentation changes

Docstrings of the entrainment types.

## Known limitations

ε is evaluated at the upper level of each step (σ[k]), a first-order quadrature that converges with
more layers. Results change at L8 since layers aren't exactly Δσ_ref thick. Parcel start from the lowest
layer (suggestion 3) is not addressed.

## Future work

Suggestions 3-5 of #1426.
