# Sea ice freezing and snow melt cap as adjustments, not tendencies

> Status: **completed** (pending CI). Freezing of sea ice (with restoring SST to freezing) and the
> cap of the snow depth tendency to the available snow moved from tendencies into `filter!`.

Date of initial draft: 2026-09-24

Base revision: `2730c024` (`mk/landstepping`, PR #1183, after merging #1264)

## Originating prompt

> I merged this PR into the other PR but could you have look at what happens when you use
> NCycleLorenz for ocean but leapfrog for atmosphere? I see the sea ice appearing on step 1 but
> disappearing again on step 2?

> Then just do another PR, also include the snow melt cap then!

## Revision log

- Diagnosed the sea ice flip, proposed moving the freezing adjustment into `filter!`, which was
  asked to be done as a separate PR together with the snow melt cap.

## Problem description

With `Leapfrog(spectral_grid, ocean = NCycleLorenz(spectral_grid))` sea ice formed on one step,
disappeared on the next (with 85 ocean points below freezing at T21) and overshot on the one
after. `ThermodynamicSeaIce` restored SST below freezing within one time step
(`dsst -= dT/Δt`) and froze ice with `f/Δt·dT`. Such a term makes up the whole deficit in one
step, so the next tendency is ≈ 0. The Lorenz N-cycle accumulates tendencies
`G = w·F + (1−w)·G` with weights 1, 3/2, 3 (N = 3, variant A) assuming a smooth tendency, so
on the second substep `G = −½·G_prev` reverses half of the adjustment, and on the third it overshoots.
`filter!` clamping the ice to [0, 1] makes it worse, as G still remembers the unclamped tendency.
The snow model's cap of the melt to the available snow (`snow_depth/Δt`) is the same kind of term.

## Background

Terms that act "within one time step" are adjustments of the state, not physical rates. Like
the snow depth cap and the soil moisture excess (#1264) they belong into `filter!`, which
is applied after the time step independent of the time stepper.

## Summary of changes

- `sea_ice_kernel!` only computes the melt tendency `dℵ = −m·max(SST − T_freeze, 0)`; it no
  longer changes the SST tendency.
- `filter!(vars, ::ThermodynamicSeaIce, model)` freezes sea ice `f·max(T_freeze − SST, 0)` at
  ocean points after the time step and restores SST to `T_freeze` (only if the ocean has an SST
  tendency, i.e. SST is not prescribed), then clamps the ice to [0, 1].
- `SnowModel` writes the uncapped snow depth tendency (snowfall − melt), the melt rate passed
  to soil moisture stays capped to the available snow, and `filter!` clamps snow depth at 0.
  With Euler forward the water budget is unchanged.

## Testing and verification

- New testset in `test/parameterizations/ocean_land_time_stepping.jl`: 8 steps with
  `EulerForward` and `NCycleLorenz` for ocean and land under Leapfrog, no ocean point below
  freezing after any step, sea ice growing monotonically, snow depth within bounds.
- Before/after at T21 (total sea ice after steps 1–8): EulerForward 0.09…1.05 (unchanged),
  NCycleLorenz before 0.13, 0.00, 0.59, 0.54, 1.64 …; after identical to EulerForward.
- `test/parameterizations/ocean_sea_ice.jl` locally, rest on CI.

## Documentation changes

None, the docs describe the time steppers, not the sea ice formulation.

## Known limitations

- With a time stepper other than Euler forward, the snow melt passed to soil moisture is capped
  with an Euler estimate of the available snow, so water is only approximately conserved.
- The melt energy `E_avail = cₛ·δT·z₁/Δt` (soil warmth above the melting threshold used within
  one time step) is also an adjustment-type term and remains a tendency.

## Future work

- Move the snow melt energy term into an adjustment too.
