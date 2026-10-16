# Euler forward time stepping for ocean, sea ice and land with leapfrog

> Status: **completed**. Ocean, sea ice and land variables are stepped with Euler forward over Δt
> (via a new `update_prognostic_surface!` hook) instead of being leapfrogged over 2Δt. This fixes
> the long-integration instability in #1183. The fixes from the review of #1183 are included.

Date of initial draft: 2026-09-19

Base revision: `8b7100d7` (`mk/landstepping`, PR #1183)

## Originating prompt

> Can you have a look at PR1183. I had generalised the time stepping for land. Review the pull
> request and flag if there's any bug or anything inconsistent and write those into the PR. I
> believe all tests are passing except that there is now some instability in the long integrations.
> If you could have a look at that.

> Previously, all land, ocean and sea ice variables simply used Euler forward. We can go back by
> setting the Leapfrog step to 1, which is a staggered Euler. If you want other ways to stabilise
> good, but we can also Euler these variables. The time stepping generalisation is mostly for the
> future when we want to use NCycleLorenz for all variables

> Yes please create a PR that points to this PR with your proposed changes

## Revision log

- **Review of #1183.** Reproduced the instability (forcing_drag 3D fails at step ~20, and the
  default PrimitiveDry year fails at step 6301). The cause is a 2Δt oscillation in the top soil
  layer where surface winds are strong: the leapfrog applies the surface-flux tendency, which is
  evaluated at t−Δt, over 2Δt. Prototyped two fixes (A: Euler over Δt; B: a linearised implicit
  surface flux) and posted the review:
  https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/1183#issuecomment-5743520635
- **User preferred Euler (Fix A).** A "staggered Euler" (only setting `which_prognostic_step = 1`)
  is still a 2Δt step and does not stabilise. So the surface update itself is changed to Euler
  over Δt, writing into both steps. The generic hook keeps other time steppers (NCycleLorenz) on
  their own `update_prognostic!`.
- **Included the review fixes** that are small and self-contained (see below). The `isfinite`
  guards removed in #1183 (they matter for `mask = false`) and the 285 K initialisation at
  ice-sheet land points are left to the author.

## Problem description

#1183 moved the soil temperature and moisture updates out of the kernels into the general
`update_prognostic!`, so with `Leapfrog` they became `T(t+Δt) = T(t−Δt) + 2Δt·f(T(t−Δt))`. That is
Euler with 2Δt, because the surface fluxes are parameterizations evaluated at step 1. The top soil
layer (z₁C₁ ≈ 2.3·10⁵ J/m²/K) relaxes towards the air temperature at a rate of about
k/(z₁C₁) with k ≈ 40–100 W/m²/K in strong winds. An explicit Euler step is only stable for
2Δt·κ < 2, and at T31 this is exceeded over Antarctica in winter and in the forcing_drag test.
`main` stepped with Δt from the current state, which gave twice the margin.

## Background

- The leapfrog in SpeedyWeather evaluates parameterizations at the previous step (i−1). Land,
  ocean and sea ice contain no advective/oscillatory terms that would benefit from centred
  (leapfrog) evaluation, only relaxation and damping terms.
- Snow melt cooled the soil with `latent_heat_sublimation` whereas `SnowModel` computes the melt
  rate with `latent_heat_fusion`, an 8.5× overestimate. This is pre-existing on `main`.

## Summary of changes

- `time_stepping/steppers/general.jl`: new `update_prognostic_surface!` (defaults to
  `update_prognostic!`) and `surface_time_step` (defaults to `default_time_step`).
- `time_stepping/time_integration.jl`: the `:ocean`/`:land` namespaces are stepped with
  `update_prognostic_surface!`.
- `time_stepping/steppers/leapfrog.jl`: `update_prognostic_surface!(…, ::Leapfrog, …)` does Euler
  forward over `surface_time_step` (Δt/2 on the first two steps, matching the clock, then Δt) and
  writes into both steps. All ocean, sea ice and land components read step 1.
- `ThermodynamicSeaIce` uses `surface_time_step` for its "restore to freezing within one time
  step" terms.
- `LandBucketMoisture`: the excess-water handling (infiltration into layer 2, river runoff, cap
  at field capacity) moved into `filter!`, where the excess of the new state is known. In #1183
  the clamp in `filter!` ran on every step, so at the next kernel call the excess was always
  zero: infiltration and runoff never happened and the excess water was deleted instead. River
  runoff is now accumulated in m (it was m·s).
- Snow melt cools the soil with `latent_heat_fusion` (new accessor with a zero fallback for dry
  atmospheres).
- Type hierarchy: `AbstractDynamic*`/`AbstractPrescribed*` are subtypes of their component type
  again (e.g. `AbstractDynamicSnow <: AbstractSnow`, `AbstractDynamicSeaIce <: AbstractSeaIce`).
  `AbstractDynamicLandComponent`/`AbstractPrescribedLandComponent` were removed; nothing
  dispatches on them any more.
- `variables(::AbstractOcean)` and `variables(::AbstractSeaIce)` are fallbacks again (no time
  dimension), so custom oceans as in `docs/src/custom_ocean.md` get their SST allocated.
- `SnowModel` no longer declares `soil_temperature` (it only reads it). This avoids a
  dimension clash with prescribed land temperatures.
- `@boundscheck` fixes: missing `|| throw(...)`, and `throw(BoundsError)` → `throw(BoundsError())`.

## Testing and verification

- New `test/parameterizations/surface_time_stepping.jl`: `surface_time_step`, the Euler update
  writing into both steps, water conservation in the `LandBucketMoisture` filter, the type
  hierarchy, the custom ocean SST fallback, a prescribed land temperature with snow, and
  `latent_heat_fusion`.
- Scripts from the review: forcing_drag 3D Dry/Wet and the year-long default Dry/Wet runs.
  With the same code as Fix A, forcing_drag and the Wet year pass. The Dry year ended in an
  atmospheric blow-up on 3 December over the North Pacific jet, which is unrelated to the surface
  (soil temperatures were normal); see the review comment.
- Existing tests: `parameterizations/{land,ocean_sea_ice,surface_fluxes,longwave_radiation,all_parametrizations}.jl`,
  `dynamics/{forcing_drag,set}.jl`, `variables/steps.jl`, `long_integrations/default_primitive_{dry,wet}.jl`.

## Documentation changes

- `docs/src/time_integration.md`: new subsection "Ocean, sea ice and land" explaining the Euler
  surface step with leapfrog and why.

## Known limitations

- The ocean, sea ice and land still need both leapfrog steps allocated even though they are
  identical after every step. Kept this way so that other time steppers can use the general path.
- The `isfinite` guards removed in #1183 (`SurfaceLandHeatFlux`, `land_snow_kernel!`) are not
  restored, so `mask = false` for land temperatures can still produce NaNs.

## Future work

- With NCycleLorenz (or longer time steps) the surface may again approach its explicit stability
  limit. The linearised implicit surface flux (Fix B in the review, as in SPEEDY's `dhfdt`) would
  remove that limit.
