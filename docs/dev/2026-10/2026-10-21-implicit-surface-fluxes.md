# Implicit time integration of surface fluxes

> Status: **in progress**. Draft PR implements `ExplicitFluxes`, `ImplicitFluxes`, `ExponentialFluxes`, with `ImplicitFluxes` the default for surface momentum, heat and humidity fluxes.

Date of initial draft: 2026-10-21

Base revision: d6e3286d

## Originating prompt

> I would like to make FriersonSigmaCoordinates the default. It works on the dry dycore see vertical coordinate tests but fails for more complex models. Can you step by increase the model's complexity and identify problems when using Frierson spacing.

> Could you draft some example code here on how the implicit formulation, or even could we do exponential time differencing could be implemented right into the `surface_flux_to_tendency` function? [...] Maybe even with an trait to distinguish like ExplicitFluxes(), ImplicitFluxes() etc, that would then be chosen from the callsite, i.e. where the fluxes are computed?

> Yes please add this sketch as a comment in the issue, then draft a PR that addresses this. Apply to momentum, heat and humidity, use `ImplicitFluxes` as default for those.

## Revision log

- 2026-10-21: initial draft, following the analysis in #1335.

## Problem description

With `FriersonSigmaCoordinates`, the lowermost layer is thin (Δσ ≈ 0.04 at 8 layers vs 0.125 equi-spaced).
Surface fluxes relax the lowermost-layer u, v, T, q towards surface values on the timescale
τ = Δp / (g ρ C V), which drops below the time step Δt (2400 s at T31) with Frierson spacing.
The explicit treatment then becomes unstable, and the dry and the full wet model blow up within
a day (#1335). Halving Δt makes them stable.

## Background

Every bulk surface flux has the relaxation form F = -a (x - x_s), with conductance a = ρ C V
(times wind slowdown, sea-ice/snow insulation, soil moisture availability where applicable), so
dx/dt = -λ (x - x_s) with λ = g a / Δp. Over the step h (2Δt for Leapfrog, `default_time_step`)
all of the following only differ by a factor φ(z), z = λ h, applied to the explicit flux:

| Scheme | φ(z) |
|---|---|
| `ExplicitFluxes` | 1 |
| `ImplicitFluxes` (backward Euler) | 1/(1 + z) |
| `ExponentialFluxes` (exponential time differencing) | (1 - e^(-z))/z |

## Summary of changes

- `parameterizations/tendencies.jl`: `AbstractFluxIntegration` with `ExplicitFluxes`,
  `ImplicitFluxes`, `ExponentialFluxes`; `flux_damping(scheme, z)`;
  `surface_flux_damping(scheme, conductance, pₛ, model)` returning φ; a 5-argument
  `surface_flux_to_tendency(scheme, flux, conductance, pₛ, model)`.
- `SurfaceMomentumFlux`, `SurfaceOceanHeatFlux`, `SurfaceLandHeatFlux`,
  `SurfaceOceanHumidityFlux`, `SurfaceLandHumidityFlux` get a `time_integration` field
  (default `ImplicitFluxes()`). The damping factor is applied to the flux *before* it is
  stored, so the fluxes passed to ocean/land and to output are the ones the atmosphere actually
  receives (conservation in coupling).
- The conductance used for φ includes the land/ocean fraction weighting, as that is what enters
  the lowermost-layer tendency.
- Prescribed fluxes and radiation are unchanged (not relaxation fluxes).

## Testing and verification

- Unit tests for `flux_damping` (limits, monotonicity, z → 0) and for an idealised single
  relaxation step: implicit and exponential stay between x and x_s for large z, explicit overshoots.
- `time_integration` option on all five flux types; `ExplicitFluxes` reproduces the previous tendency.
- Short `PrimitiveDryModel` and `PrimitiveWetModel` runs with `FriersonSigmaCoordinates` and
  default physics at default Δt, which NaN'd before this change.

## Documentation changes

Docstrings of the new types; short section in the surface fluxes docs.

## Known limitations

- Ocean and land fluxes into the same layer are damped independently instead of with the summed
  conductance; stable but slightly less accurate.
- C and V are held fixed over the step (linearisation); x_s (SST, soil temperature) is treated as fixed.
- The first (Euler, Δt/2) step uses the 2Δt factor and is slightly over-damped.

## Future work

- Make `FriersonSigmaCoordinates` the default (#1335), once the near-surface warm bias in the
  wet model without convection is understood.
- Possibly treat vertical diffusion implicitly in the same way.
