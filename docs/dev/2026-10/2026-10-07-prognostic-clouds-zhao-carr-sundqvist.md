# Prognostic clouds for SpeedyWeather: Sundqvist, Zhao–Carr, and how ICON does it

> Status: **in progress**. Stages 1, 2 and 3 with per-layer shortwave clouds are implemented and
> tested on CPU and GPU (see Summary of changes); all-sky ecCKD, the radiation call frequency and
> the tuning with SpeedyCalibration.jl are open.

Date of initial draft: 2026-10-07

Base revision: `3c033168` (`mg/clouds`). The initial draft was checked against `dbd4c661`
(`mg/version1-version023`), the 2026-10-07 re-audit against `3c033168`.

## Originating prompt

> What would SpeedyWeather need to have prognostic clouds? Specifically reference also how ICON is
> doing that and also research whether they are easier and simpler schemes avaliable

> Summarize this research in a markdown file, put particular focus on the Zhao-Carr and Sundqvist
> scheme and the recommendations

## Revision log

- 2026-10-08, tuning with SpeedyCalibration.jl, then the water loss, per-layer shortwave clouds,
  Stages 2 and 3.
  > Ok, before we continue note in the plan that our aim is to tune the cloud parameterization with
  > SpeedyCalibration.jl (~/SpeedyCalibration.jl) using Enzyme. For this purpose we might also make
  > changes to SpeedyCalibration. That's fine, we just want to reuse the core calibration loop.
  >
  > Then
  > * investigat ethe issue of the loss of water
  > * Implement per-layer shortwave clouds
  >
  > Then proceed with Stage 2 and 3

  - New section 7.8 on tuning with SpeedyCalibration.jl (its loop, the changes it needs, what
    SpeedyWeather needs); the generic Enzyme item of Future work replaced by a pointer.
  - Water loss of `ImplicitCondensation` investigated: two mechanisms, 5 % of large-scale
    precipitation globally at T32; fix proposed and tested in a script, not applied (appendix).
  - Per-layer shortwave clouds implemented (7.2.6): `CloudyShortwaveRadiativeTransfer` with
    absorbing two-stream cloud layers (single-scattering albedo 0.999) and the adding method;
    `OneBandCloudyShortwave` and `OneBandCloudyLongwave` as convenience constructors.
  - Stage 2 implemented as `SundqvistClosure`, a component of `PrognosticCloudCondensation`; the
    relaxation of Stage 1 moved into `RelaxationClosure` (default). Correction to 7.4: physics reads
    the lagged step, so consecutive calls are one step apart and alternate leapfrog parity; the
    supply needs the reference from two calls before, i.e. two copies of the reference state, as
    GFS's `tp` and `tp1`, not one. The reference includes only the scheme's condensation, not its
    precipitation processes, as GFS stores it after `gscond` and before `precpd`.
  - Stage 3 decided for option (a): `BettsMillerConvection(; detrainment)` detrains a fraction of
    the deep convective precipitation as condensate at the level of zero buoyancy, default 0.
  - Found, not fixed: Betts-Miller convective snow carries no latent heat of fusion.
  > just merge main back into this branch, this should also fix that

  - The full test suite failed to start because Terrarium 0.1.8, released overnight, broke
    `SpeedyWeatherTerrariumExt`. Merged `origin/main` (release 0.23.0 with the Terrarium 0.1.8 fix
    #1300) into `mg/clouds`; `SpeedyWeather` is now 0.24.0-DEV. Main also uses the matrix transform
    on GPU up to T100 now, which makes the larger transform scratch of 7.2.1 matter up to T100.
  - *State at the end of the session (2026-10-08):*
    - Before the merge, all new and changed cloud testsets passed on CPU (Two-stream, per-layer
      shortwave, Sundqvist closure, refactored Stage 1, detrainment), and both GPU testsets passed
      on an A40.
    - After the merge, the full CPU suite and the GPU testsets were started but had not finished;
      the suite got past the Terrarium precompilation that had failed before. Rerun both.
    - `SpeedyWeather/src/parameterizations/cloud_condensation.jl` was untracked when the merge was
      made, so commit `2feedc0a` does not build on its own; commit it with the remaining changes.
  - *Next steps:*
    1. Rerun the full test suite and the GPU tests on the merged branch.
    2. Commit the session's work (the scheme file, Stages 2 and 3, per-layer shortwave, docs).
    3. Enzyme differentiability tests of the cloud scheme and the cloudy radiation (7.8).
    4. Adapt SpeedyCalibration.jl (7.8) and tune: cloud cover and shortwave cloud effect are far too
       small (Summary of changes).
    5. Decide on the `ImplicitCondensation` water fix (appendix), as its own PR.
    6. Clear-sky fluxes for cloud radiative effects; all-sky ecCKD and the radiation call frequency.
- 2026-10-07, tuning by differentiation.
  > note in the plan that enzyme differentiabillity might be used to tune the cloud model

  - Added gradient-based tuning with Enzyme to Future work, with what it needs; pointers from the
    differentiability checklist (7.6) and the "not tuned" limitation.
- 2026-10-07, implementation of Stage 1 and the Stage 0 parts it needs.
  > okay start with the implementation

  - Implemented 7.2.1–7.2.4, the cloud state of 7.2.5, the one-band radiation of 7.2.6 at the
    "day one" level, 7.2.8, and Stage 1 (7.3). Not yet: per-layer shortwave reflection and
    absorption, all-sky ecCKD (needs the upstream solver), 7.2.7, Stages 2 and 3.
  - Hook 4 of 7.2.1 dropped: the condensate's grid copy is not clipped, the scheme fills
    negative condensate from vapour instead (clipping the grid copy would hide the negatives
    from the fill and add mass to the advected field). Vertical advection turned out to be a
    fifth hook: its call site names the variables explicitly.
  - The hooks are generated over a registry `ADVECTED_SCALARS` of `(name, u_flux, v_flux)`
    filtered by the names present in the `Variables` type; the humidity path emits the same
    calls as before.
  - New, not in the plan: the fused tendency batch with the condensate is 12L+1 = 97 layers,
    wider than the transform scratch of 9L+1. The matrix transform (the GPU default up to T64)
    throws for that, the FFT transform on GPU falls back to the serial path. `max_transform_batch`
    is now 12L+1; the eagerly planned GPU batches are unchanged (9L+1), wider ones are planned on
    first use. GPU scratch memory of the FFT transform grows by a third for all models.
  - New, not in the plan: SPPT read `vars.scratch.a_grid`, which does not exist, so SPPT failed
    for models without humidity. Fixed together with adding the condensate.
  - Stage 1 details decided during implementation: rain reevaporation capped at saturation;
    an extra latent heat term keeps the condensate at its ice fraction when liquid and ice
    autoconvert at different rates, which makes the enthalpy budget exact; effective radii as
    GFS (10 µm ocean, 5–10 µm land, 50 µm ice); the overlap parameter `α` is not part of the cloud
    state yet (the one-band schemes use maximum-random and random overlap).
  - The "bitwise identical" test of the plan does not hold: the wider batched transforms change
    rounding (Float32, relative differences of order 1e-7 after 20 steps). The test compares to
    `rtol = 1e-5` instead.
  - Found, not fixed: `ImplicitCondensation` does not conserve water with reevaporation on (the
    default). The evaporated rain is removed from the flux in full but added to humidity divided
    by the relaxation factor; in a test column only 43 % of the condensed water reached the
    surface or returned to vapour.
- 2026-10-07, review and revision.
  > Critically review and reaudit the plan and the respective SpeedyWeather functionality. Also be
  > aware that we could also set up a full ecCKD/ecRad model with NumericalRadiation that uses
  > prognostic clouds. The goal of every implementation should also be to have an implementation
  > that is efficient to compute on GPU

  > review how the implementation fits our fused variables approach

  > okay revise the plan based on this. Then re-audit the full plan again once more.

  - Base revision moved to `3c033168`; every path and line reference re-checked there (one
    citation corrected: the tracer hyperdiffusion of the primitive-equation models).
  - The Zhao–Carr formulas, bounds, phase rules and constants of sections 3 and 4 were checked
    against the v6.0.0 Fortran of `gscond`, `precpd` and `radiation_clouds.f` and hold. Two
    details added in 4.5: the `clwmin` offset inside the Xu–Randall exponent and the default
    ice effective radius.
  - Section 1 gained rows on variable fusion and on what column kernels can see (time step,
    `model` fields); section 2 two needs (SPPT, local mixed-phase saturation).
  - Stage 0 rewritten (7.2): the condensate is a dedicated prognostic variable fused with the
    atmospheric prognostics, not a `Tracer`; the four dynamical-core hooks it needs are listed
    (7.2.1). New: the time step physics sees versus the step the state is advanced with (7.2.2),
    SPPT (7.2.3), mixed-phase saturation kept local (7.2.4), a cloud state as the interface to
    radiation (7.2.5), the one-band schemes with layer clouds and an all-sky ecCKD scheme through
    NumericalRadiation with its GPU requirements (7.2.6), a radiation call frequency (7.2.7).
  - Stage 1 (7.3): cloud evaporation only in the clear fraction with its own time scale; sinks
    bounded with the prognostic step. Stage 2 (7.4): the forcing bookkeeping restated for
    tendency-based physics; the reference states are prognostic variables without tendencies.
  - Numerics checklist, tests, documentation, known limitations and future work updated.
  - Re-audit of the revised plan: every `file:line` citation resolved against `3c033168` by
    script (four path shorthands corrected), tables and cross-references checked, and the fused
    declaration of 7.2.1 verified empirically with a stub component in `Variables(model)`.
- 2026-10-07: initial draft.
  - Based on a code survey and literature research from 2026-09-28.
  - The SpeedyWeather facts were re-checked at the base revision.
  - The Sundqvist and Zhao–Carr formulas and constants were checked line by line against NCEP's
    GFS Fortran (`zhaocarr_gscond.f`, `zhaocarr_precpd.f` and `radiation_clouds.f` in
    ccpp-physics `v6.0.0`, the defaults in fv3atm `GFS_typedefs.F90`).
  - The ICON options were checked against the ICON namelist overview of 2025-04-22.

## Summary

- **Today:** SpeedyWeather has no condensate variable.
  - Large-scale condensation rains out in the same step and column.
  - Clouds are diagnosed SPEEDY-style: one cover and one cloud-top level per column.
  - They only affect shortwave radiation.
- **What already exists:** the variable system with fused parents (humidity is the template), the
  column-physics interface and the precipitation sweep of `ImplicitCondensation` cover most of
  what a one-condensate scheme needs.
- **What is missing:**
  - a condensate variable the physics writes a tendency for, fused with the other prognostics;
  - positivity of a spectrally transported condensate;
  - saturation over ice;
  - a cloud state (fraction, water, ice and effective radii in each layer) and radiation that
    uses it in shortwave *and* longwave, in the one-band schemes and in ecCKD;
  - a time step accessor for physics that returns the step the state is advanced with.
- **Sundqvist (1978; Sundqvist, Berge & Kristjánsson 1989)** has three ingredients:
  - one prognostic cloud condensate `m`;
  - a cloud fraction diagnosed from relative humidity, `b = 1 − √((1 − f)/(1 − u))`;
  - a condensation closure: the cloudy fraction `b` of the large-scale moisture supply condenses and
    the rest moistens the clear part.
  
  On top of that it uses a smooth precipitation-release law, `P = C₀ m [1 − exp(−(m/(b m_r))²)]`.
- **Zhao & Carr (1997)** turn Sundqvist into a complete NWP scheme:
  - one tracer for cloud water *and* ice, with the phase decided by temperature;
  - ice-phase processes;
  - precipitation diagnosed in one top-down sweep, with evaporation and melting;
  - in GFS, a separate Xu–Randall cloud fraction for radiation.
  
  It was operational in NCEP's *spectral* GFS from 2001 to 2019. Until 2015 that GFS used an Eulerian
  core with a three-time-level (leapfrog) step. This makes it the closest operational analogue to
  SpeedyWeather.
- **ICON** predicts the condensate and diagnoses the cloud fraction, but it is far heavier:
  - NWP and XPP: a one-moment scheme with 5–6 species and prognostic rain and snow, plus a
    PDF-based diagnostic cover;
  - Sapphire/AES: a graupel scheme with clouds either fully on or off in each cell;
  - the earlier ECHAM-physics ICON-A used Sundqvist cover, so it belongs to the same family as
    Zhao–Carr.
- **Recommendation** (section 7):
  - Stage 0: the fused condensate variable, the physics time step, the cloud state, layer clouds
    in the one-band radiation, and the prerequisites of an all-sky ecCKD scheme.
  - Stage 1: Zhao–Carr-type physics, keeping `ImplicitCondensation`'s relaxation for condensation.
  - Stage 2: the Sundqvist condensation closure.
  - Stage 3: convective detrainment.
  - Defer prognostic rain and snow, two-moment schemes and prognostic cloud fraction.

## Problem description

A prognostic cloud scheme carries cloud condensate from step to step, advects it, converts it to
precipitation and gives it to radiation as cloud water and ice. This changes cloud radiative effects
from a function of instantaneous relative humidity and precipitation to a function of the model's
condensate history. It also lets anvils and stratiform shields be advected away from where they form.

This document collects:

- what SpeedyWeather lacks;
- how the Sundqvist and Zhao–Carr schemes work;
- how ICON does it;
- which simpler schemes exist;
- what to implement first.

## Background

### 1. What SpeedyWeather has today

Paths are relative to `SpeedyWeather/src/` at the base revision.

| Component | State | Where |
|---|---|---|
| Water variables | Only specific humidity `q` is prognostic (spectral). There is no cloud water, ice, rain or snow variable. | `models/primitive_wet.jl` |
| Tracers | `add!(model, Tracer(:name))` creates a spectral prognostic variable, a grid copy, and spectral and grid tendencies (namespaces `:tracers` and `:grid_tracers`). Tracers are advected (flux form) and hyperdiffused. No parameterization writes a tracer tendency yet. Tracers are **not fused**: each costs one separate spectral→grid and three separate grid→spectral transforms per step through the shared `:a/:b` scratch slots, driven by a runtime `Dict` loop. There is no tracer test under `test/GPU/` or `test/reactant/`. | `dynamics/tracers.jl:9-15, 57-71`, `dynamics/tendencies.jl:805-819`, `dynamics/horizontal_diffusion.jl:357-361`, `time_stepping/transform.jl:181-185` |
| Positivity | `ClipNegatives` does `max(q, 0)`, but only on the **grid copy of humidity**. The spectral state is untouched and mass is not conserved. Tracers get nothing. Clipping another grid variable is one more line in the same place. | `dynamics/hole_filling.jl`, `time_stepping/transform.jl:153-156` |
| Large-scale condensation | `ImplicitCondensation` relaxes `q` to `0.95 q_s` over `3Δt`, with the implicit latent-heat factor `1 + Lᵥ/cₚ ∂q_s/∂T` (Frierson et al. 2006). The condensate rains out immediately. A top-down sweep handles melting (above 278 K), reevaporation (proportional to subsaturation) and snow (below 263 K), and sets `cloud_top`. | `parameterizations/large_scale_condensation.jl` |
| Convection | Betts–Miller (Frierson 2007) with linear entrainment. Rain is the net column drying. All of it falls as snow if the lowest layer is below 273.15 K. No condensate is detrained. `cloud_top` is set to the level of zero buoyancy. | `parameterizations/convection.jl:8-23, 141-179` |
| Clouds | `DiagnosticClouds` (SPEEDY): one cover per column, `min(1, 0.2√P + RH-term²)`. Cloud top is the level of maximum RH or the precipitation cloud top. Stratocumulus comes from static stability. | `parameterizations/radiation/clouds.jl:29-185` (cover at `:149-151`) |
| Radiation | `OneBandShortwave` reflects `cloud_albedo × cover` only at `k == cloud_top` and adds cloud absorptivity below it. The Frierson longwave has **no clouds**. The NumericalRadiation (ecCKD) extension is clear-sky; it declares its per-g-point work arrays as `Grid4D` parameterization variables in the `:ecckd` namespace and runs the clear-sky solvers allocation-free with caller-owned scratch. | `parameterizations/radiation/shortwave_radiation.jl:192-196, 233`, `parameterizations/radiation/shortwave_transmissivity.jl:109-131`, `parameterizations/radiation/longwave_transmissivity.jl:29-78`, `ext/SpeedyWeatherNumericalRadiationExt/ecckd_radiation.jl` |
| Saturation | Clausius–Clapeyron over liquid with one `Lᵥ`. `TetensEquation` has a Murray (1967) branch below freezing, but it is unused anywhere in `src/`. `saturation_humidity` is shared by convection, the diagnostic clouds, surface fluxes and condensation. `latent_heat_sublimation` is 2801 kJ/kg while `Lᵥ + L_f` is 2831 kJ/kg. | `dynamics/atmosphere.jl:46-52, 141-160, 162-216` |
| Physics infrastructure | Column kernels `parameterization!(ij, vars, scheme, model)` add into the grid tendencies and are fused into one GPU kernel. The call order is zenith, vertical diffusion, condensation, convection, albedo, radiation, boundary layer, surface fluxes, stochastic physics. Physics reads the *lagged* leapfrog step and runs before dynamics. On the device `model` is reduced to its `core_components`, so `model.tracers` and other `Dict`s are not available in the kernel. SPPT multiplies the u, v, T and q tendencies only. | `models/primitive_wet.jl:122-134, 230-235`, `parameterizations/tendencies.jl:32-88`, `time_stepping/steps.jl:187`, `time_stepping/time_integration.jl:83-90`, `parameterizations/stochastic_physics.jl:74-84` |
| Time stepping | Leapfrog with `Δt = 40 min` at T32; after the first two steps the state is advanced by `2Δt` (`default_time_step`), so physics tendencies act over `2Δt`. Column kernels only see `time_stepping.Δt`, and the GPU-adapted `LeapfrogCore` carries nothing else. `ImplicitCondensation` therefore removes two thirds, not one third, of the supersaturation per step. | `time_stepping/steppers/leapfrog.jl:171-179, 186, 246-256` |
| Post-step adjustments | `filter!` exists only for the grid-point components ocean, sea ice and land; there is none for spectral atmospheric variables. The precedent for adjustments instead of tendencies is `docs/dev/2026-09/adjustments-not-tendencies.md`. | `time_stepping/time_integration.jl:140-142` |
| Variable fusion | Variables sharing a `fuse` symbol within a namespace share one parent buffer, concatenated along the layer axis in declaration order (the model's own variables first, then the components in field order). The `:prognostic`↔`:grid` and `:spectral_tendencies`↔`:grid_tendencies` parents must declare the same members in the same order, asserted in `Variables(model)`. The batched transforms in both directions, tendency scaling, the leapfrog shift of the grid copies and the time stepping over all top-level tendency names iterate the parents or the names, so they pick up new members automatically. Hyperdiffusion, the horizontal flux-form advection and hole filling name their variables explicitly. In a parent with a step dimension a member without one collapses to a single slot. Several in-place broadcasts into views of one parent corrupt data under Reactant. | `variables/variables.jl:237-272, 398-465`, `variables/dimensions.jl:177-330`, `time_stepping/transform.jl:149-165`, `dynamics/tendencies.jl:160-175, 778-803`, `time_stepping/steppers/leapfrog.jl:125-144`, `dynamics/scaling.jl:15-60`, `time_stepping/time_integration.jl:152-165` |

### 2. What a prognostic cloud scheme needs

| # | Need | Status |
|---|---|---|
| 1 | A condensate variable `q_c` (cloud water + ice), advected | The variable system allocates it. Fused with humidity (7.2.1) it joins the batched transforms, tendency scaling, the leapfrog grid shift, vertical advection and the time stepping automatically; four hooks in the dynamical core remain. Physics writes `vars.tendencies.grid.cloud_condensate`, transformed together with the advection tendency. |
| 2 | Sources and sinks: condensation, cloud evaporation, autoconversion, accretion, phase | Missing. `ImplicitCondensation`'s relaxation can be reused as the condensation source. |
| 3 | Precipitation with evaporation and melting | The top-down sweep exists in `ImplicitCondensation`. |
| 4 | A cloud fraction consistent with `q_c` | Missing. The current cover is RH- and precipitation-based, one per column. |
| 5 | Radiation using cloud fraction, condensate path and effective radius in each layer, in SW **and** LW, with an overlap assumption | Missing, and the largest single work item. Needs a cloud state as the interface (7.2.5). NumericalRadiation has GPU-fit cloud optics, but its all-sky solvers allocate (7.2.6). |
| 6 | Convective condensate (anvils) | Missing. Betts–Miller has no updraft to detrain from. |
| 7 | Positivity and conservation of a spectrally transported, intermittent field | Missing for tracers. |
| 8 | Ice thermodynamics: `q_s` over ice or mixed phase, `L_s`, fusion heat | Partly present (Tetens–Murray). |
| 9 | Bounded fast sinks under leapfrog (`2Δt = 80 min` at T32) | A design rule: implicit or exponential forms, bounded with the step the state is advanced with, which physics cannot see today (7.2.2). |
| 10 | Output: liquid and ice water path, 3D cloud fraction | Missing. |
| 11 | SPPT consistent across `q`, `T` and `q_c` | Missing: SPPT perturbs u, v, T and q only (7.2.3). |
| 12 | Mixed-phase saturation without changing the rest of the model | Missing (7.2.4). |

### 3. The Sundqvist scheme

Sundqvist (1978) introduced it in QJRMS, and Sundqvist, Berge & Kristjánsson (1989, MWR) developed
it in a mesoscale NWP model. The idea:

- carry one prognostic cloud condensate `m`;
- diagnose cloud cover from relative humidity `f`;
- close the condensation rate with the large-scale moisture supply;
- release precipitation as a smooth function of the in-cloud water.

The equations below are in the form NCEP implemented them (GFS `zhaocarr_gscond.f`,
`zhaocarr_precpd.f`). The original papers could not be retrieved (AMS returned 403), so the constants
are GFS's tuning, which may differ from the 1989 values.

#### 3.1 Cloud fraction from relative humidity

```
b = 1 − √((1 − f)/(1 − u))   for u < f < 1,     b = 0 for f ≤ u,     b = 1 for f ≥ 1
inverse:  f = 1 − (1 − u)(1 − b)²
```

`u` is the critical relative humidity, the grid-mean RH at which subgrid condensation starts.

- **GFS:** `u` is linear in the Exner function between the surface, the PBL top and the model top. All
  three defaults are 0.90 (`crtrh`). It is blended towards 0.9999999 for grids finer than 192×94, with
  a weight from the log of the cell area, so at T31/T32 it is 0.90 everywhere.
- **ECHAM / ICON-A** (crs, crt and nex in the ECHAM-physics namelist):
  - `u(p) = u_top + (u_surf − u_top) exp(1 − (p_s/p)ⁿ)`;
  - plus a low-cloud adjustment below ocean inversions (Mauritsen et al. 2019).
  - Grundner et al. (2022) refitted the generalised form `b = 1 − √((min(f, f_sat) − f_sat)/(u − f_sat))`
    to coarse-grained storm-resolving ICON output. Their fitted values of
    `{f_sat, u_top, u_surf, n}` are `{1.12, 0.3, 0.92, 0.8}` over land and `{1.07, 0.42, 0.9, 1.1}`
    over sea.
  
  I could not find the ICON-A defaults.
- **AD caveat:** `∂b/∂f → ∞` as `f → 1`. A differentiable implementation needs `f ≤ 1 − ε`, or a
  smoothed form.

#### 3.2 Condensation closure

Let `A_T`, `A_q` and `A_p` be the tendencies of `T`, `q` and `p` from **all processes except** the
cloud scheme (dynamics plus other physics). The moisture supply relative to the moving saturation
value is

```
M = A_q − f (∂q_s/∂T) A_T − f (∂q_s/∂p) A_p,      ∂q_s/∂T = ε L q_s / (R_d T²),   ∂q_s/∂p = −q_s/p
```

Sundqvist's hypothesis:

- the part `bM` condenses in the already cloudy fraction;
- the part `(1 − b)M + E_c` raises the humidity of the clear part, and with it `b`.

Combined with `b(f)`, this gives the relative-humidity tendency and the condensation rate:

```
f_t = 2(1 − b)(1 − u) [(1 − b) M + E_c] / (2 q_s (1 − b)(1 − u) + m/b)
C   = (M − q_s f_t) / (1 + (L/cₚ) f ∂q_s/∂T)
```

For a cloudy cell (`E_c = 0`) this simplifies to the form below. I checked algebraically that it is
identical to the reorganised expression in `gscond`.

```
C = β M / (1 + (L/cₚ) f ∂q_s/∂T),     β = (b²(1 − b)(1 − u) q_s + m/2) / (b(1 − b)(1 − u) q_s + m/2)  ∈ [b, 1]
```

- With no cloud water, `β = b`: only the cloudy part of the supply condenses.
- With a lot of cloud water, `β → 1`.
- The denominator is the same implicit latent-heat factor that `ImplicitCondensation` already uses
  (Frierson et al. 2006, eq. 21), here weighted by `f`. `L` is `Lᵥ` or `L_s` depending on the phase.

GFS bounds the result:

- condensation only if `b > 10⁻³`;
- `0 ≤ C ≤ (q − u q_s)/Δt`, so the cell is never dried below the critical RH.

**Implementation consequence:** the closure needs the forcing `M`, that is, tendency bookkeeping.
Section 7.4 shows how GFS obtained it under leapfrog.

#### 3.3 Precipitation release (autoconversion)

```
P = C₀ m [1 − exp(−(m/(b m_r))²)]
```

- `m/b` is the in-cloud water and `m_r` a characteristic in-cloud value.
- For `m/b ≪ m_r`, `P ≈ C₀ m³/(b m_r)²`. This is a *smooth* threshold, with no `if` and differentiable.
- For `m/b ≫ m_r`, `P → C₀ m`.

Sundqvist added two enhancements: `C₀ → F C₀` and `m_r → m_r/F`, with `F = F_co F_BF`.

- Coalescence with precipitation falling from above: `F_co = 1 + c₁ √P_above`.
- Bergeron–Findeisen in supercooled cloud: `F_BF = 1 + c₂ √(min(max(268 K − T, 0), 20 K))`.

GFS values:

| Symbol | Value |
|---|---|
| `C₀` (`prautco`) | 1×10⁻⁴ s⁻¹ (e-folding about 2.8 h) |
| `m_r` | 3×10⁻⁴ kg kg⁻¹ |
| `c₁` | 300, with `P` in kg m⁻² s⁻¹; a commented-out line in the code has 100 |
| `c₂` | 0.5 |

Further GFS details:

- `m` is replaced by `m − w_min`, with `w_min = 10⁻⁵ kg kg⁻¹ × p/1000 hPa`.
- `b` is floored at 0.01.
- The argument of the exponential is capped at 50.
- The rate is limited to the available `m`.

#### 3.4 Evaporation of precipitation

```
E_r = K_E (u − f) √P_r        K_E = 2×10⁻⁵ (GFS evpco), E_r in kg kg⁻¹ s⁻¹, P_r in kg m⁻² s⁻¹
```

Evaporation only happens below the critical RH. It is limited so that the layer does not exceed `u`,
and by the flux itself. SpeedyWeather's current reevaporation is linear in `q_s − q` and removes a
fraction of the flux per layer. Either form fits the same sweep.

#### 3.5 Assessment

- **Pros:**
  - one tracer;
  - smooth, column-local formulas;
  - few parameters (`u`, `C₀`, `m_r`, `K_E`);
  - decades of use: the Norwegian LAM, the NRL model, NCEP Eta/GFS, and the ECHAM family's cover
    (via Xu & Krueger 1991 and Lohmann & Roeckner 1996).
- **Cons:**
  - RH-only cover can be cloudy without condensate, or the reverse. In ICON-A output about 7 % of
    the cloudy cells had no condensate because cover is diagnosed before the microphysics
    (Grundner et al. 2022).
  - `u` dominates the tuning and depends on resolution and layer thickness.
  - The closure needs the non-cloud forcing `M`.

### 4. The Zhao–Carr scheme

#### 4.1 History

- **Origin:** Zhao & Carr (1997, MWR) for NCEP's operational models. The Eta implementation is
  described in Zhao, Black & Baldwin (1997, *Wea. Forecasting*). The code was created by Q. Zhao in
  January 1995 and rewritten by S. Moorthi and H.-L. Pan in 1998–2000 for the MRF/AVN, which later
  became GFS.
- **GFS:** operational in May 2001. An August 2001 change to ice autoconversion fixed excessive light
  precipitation (CCPP Scientific Documentation). GFS v15 (June 2019) replaced it with GFDL
  microphysics on the FV3 core.
- **Why it is the closest analogue:** until the semi-Lagrangian T1534 upgrade in January 2015, GFS was
  an **Eulerian spectral** model, and `gscond` explicitly supports three-time-level (leapfrog)
  stepping. So Zhao–Carr ran for years in a spectral, leapfrog model with a spectrally transported
  condensate, like SpeedyWeather would.

#### 4.2 Prognostic variable and phase

- There is one condensate mixing ratio, `cwm`, for liquid *or* ice. The phase flag `IW` (0 = water,
  1 = ice) follows Table 2 of Zhao & Carr:
  - `T ≥ 0 °C`: water.
  - `T < −15 °C`: ice, if there is cloud or `q > u q_s`.
  - `−15 °C ≤ T < 0 °C`: water, unless the layer **above** is ice and there is cloud here (seeding).
    The paper also describes a memory of the previous time step's phase; the v6.0.0 code recomputes
    the flag each call.
- `L = Lᵥ` for water and `L_s` for ice. The code uses one saturation humidity from GFS's `fpvs` for
  both phases (an explicit ice-saturation formula is commented out).

#### 4.3 Condensation and cloud evaporation (`gscond`, column loop from the top down)

1. Compute the Sundqvist `b` from RH and `u(p)`.
2. **If `b ≤ 10⁻³`** and cloud exists, evaporate it towards `u q_s`. This is a saturation adjustment of
   three Newton iterations (the first increment halved), limited to `cwm/Δt`.
3. **If `b > 10⁻³`**, condense with the Sundqvist closure (section 3.2), with the bounds described
   there.
   - The forcing is `A_X = (X − X_ref)/dt`.
   - `X_ref` is the state stored at the end of the previous `gscond` call, so it includes everything
     except this scheme's own previous condensation.
   - **Under leapfrog** `dt = 2Δt` and the reference is from two calls back (`tp ← tp1, tp1 ← t`), so
     the difference is taken at the same leapfrog parity and avoids the computational mode.
4. Update `cwm += (C − E_c)Δt`, `q −= (C − E_c)Δt` and `T += L/cₚ (C − E_c)Δt`.

#### 4.4 Precipitation (`precpd`)

Precipitation is **diagnostic**. It is integrated from the top down in one sweep per time step: it
has no fall speed, no storage and no CFL limit. This is the same structure as SpeedyWeather's sweep.
The rain and snow fluxes accumulate:

```
P_r(k) = Σ_above (P_raut + P_racw + P_sacw + P_sm1 + P_sm2 − E_rr) Δp/g
P_s(k) = Σ_above (P_saut + P_saci − P_sm1 − P_sm2 − E_rs) Δp/g
```

| Process | Form (GFS v6.0.0) | Constants |
|---|---|---|
| Cloud water → rain, `P_raut` | Sundqvist (3.3) with `F_co`, `F_BF` | `C₀ = 1e-4 s⁻¹`, `m_r = 3e-4`, `c₁ = 300`, `c₂ = 0.5` |
| Cloud water collected by rain, `P_racw` | `C_r m P_r` | Switched off in operational GFS (commented out) |
| Cloud ice → snow, `P_saut` (Lin et al. 1983) | `a₁ (m − m_i0)`, `a₁ = c_saut exp(0.025 (T − T₀))` | `c_saut = 3e-4 s⁻¹` on coarse grids, `6e-4` on fine grids (Zhao & Carr: 1e-3). `m_i0 = 1e-5 × p/1000 hPa` |
| Ice accreted by snow, `P_saci` | `C_s m P_s`, `C_s ∝ exp(0.025 (T − T₀))`, zero above freezing | `1.25e-3` † |
| Rain evaporation, `E_rr` | `K_E (u − f) √P_r` | `K_E = 2e-5` |
| Snow sublimation, `E_rs` (T < 0 °C) | `[C_rs1 + C_rs2 (T − T₀)] (u − f)/u P_s` | `C_rs1 = 5e-6`, `C_rs2 = 6.67e-10` † |
| Snow melt, `P_sm1` (T > 0 °C) | `C_sm (T − T₀)² P_s` | `C_sm = 5e-8` †. Melts almost all snow before 5 °C |
| Melt by collecting cloud water, `P_sm2` | `C_ws C_r m P_s` | `C_ws = 0.025`, `C_r = 5e-4` † |

† In the code these coefficients multiply precipitation *per 800 s* (`zaodt = 800/dt`). For a flux in
kg m⁻² s⁻¹ the effective coefficient is 800 × the listed value. The Sundqvist terms (`P_raut`,
`E_rr`) carry no such factor.

Further details of `precpd`:

- Every process is limited to the condensate or flux available in that layer.
- `T` and `q` are updated with `Lᵥ`, `L_s` and `L_f`.
- Evaporation of melting snow is ignored.
- **Negative-condensate fix:** a negative `cwm` after the processes is refilled from vapour in the same
  cell, with the matching latent heating. This conserves total water and enthalpy. If there is not
  enough vapour, `q` is set to 0.
- Outputs are the surface rain + snow and the snow fraction.

#### 4.5 Clouds for radiation in GFS (`progcld_zhao_carr`)

The radiation does **not** use the Sundqvist `b`. Instead:

- **Cloud fraction** (Xu & Randall 1996 form, GFS constants), only where `q_c > 10⁻⁶ × p/1000 hPa`:

  ```
  C = f^¼ [1 − exp(−2000 q_c / ((1 − f) q_s)^¼)]
  ```

  - `((1 − f) q_s)^¼` is clamped to `[10⁻⁴, 1]`, `1 − f` to at least `10⁻¹⁰`, and the exponent to 50.
  - Inside the exponent the code uses `q_c − clwmin/(p/1000 hPa)` with `clwmin = 10⁻⁹`, negligible.
  - `C < 0.001` is set to 0.
  - Because `C → 0` as `q_c → 0`, there is no cloud without condensate.
- **Condensate path** `q_c Δp/g`:
  - The ice fraction is a linear ramp, `clamp((T₀ − T)/20 K, 0, 1)`, independent of the
    microphysics `IW`.
  - Optionally the path is normalised to in-cloud values (divided by `C`) and smoothed vertically
    1-2-1.
- **Effective radius:**
  - liquid: 10 µm over ocean, `5 + 5 f_ice` µm over land;
  - ice: a temperature-dependent power law of ice water content (Heymsfield & McFarquhar 1996),
    limited to 10–150 µm; 50 µm where the power law is not used.

#### 4.6 Known weaknesses

- **Phase:** liquid and ice share one variable, and the phase rule is discontinuous. When the phase
  flips the fusion heat of cloud ice is not accounted for. Cloud ice has no sedimentation (it only
  leaves via snow).
- **Two cloud fractions:** microphysics uses Sundqvist `b(f)` and radiation uses Xu–Randall
  `C(q_c, f)`.
- **Tuning:** constants are tied to time-step conventions (the 800 s factor) and to GFS, so they need
  retuning.
- **Early bias:** too much light precipitation after the first GFS implementation, fixed by changing
  ice autoconversion.

### 5. How ICON does it

ICON predicts the condensate and **diagnoses the cloud fraction**, in every configuration.

- **NWP physics** (DWD ICON; also ICON-XPP, the CMIP7 climate configuration, Müller et al. 2025).
  - **Microphysics** (`inwp_gscp`):
    - 1 = default "cloud-ice" scheme (COSMO-EU heritage; prognostic `qv, qc, qi, qr, qs`);
    - 2 = graupel scheme (adds `qg`);
    - 4–7 = Seifert–Beheng two-moment;
    - 8 = spectral bin;
    - 9 = Kessler.
    
    Rain and snow (and graupel) are prognostic, including sedimentation. Liquid condensation uses
    saturation adjustment (`satad`).
  - **Cloud cover** (`inwp_cldcover`): the default 1 is Martin Köhler's diagnostic scheme.
    - It combines grid-scale condensate with box/PDF subgrid estimates (`tune_box_liq` and
      `tune_box_ice` both 0.05, asymmetry 2.5) and convective contributions.
    - Other options: 3 = COSMO subgrid scheme, 4 = clouds from the turbulence scheme, 5 = grid-scale
      0/1.
    - Option 2, "prognostic total water variance", is listed as *not yet started*.
- **AES physics** (successor of ICON-A; ICON-Sapphire at km scale, Hohenegger et al. 2023).
  - Graupel one-moment microphysics (`mig`), two-moment optional.
  - Cloud cover is **all-or-nothing**: 1 if `q_c + q_i > cqx = 1e-8 kg kg⁻¹`, else 0.
  - Cloud optics use droplet number concentrations (20–180 × 10⁶ m⁻³ over land, 20–80 over sea),
    inhomogeneity factors (liquid 0.4, ice 0.8) and homogeneous freezing at `T₀ − 35 K`.
  - Its graupel kernel has been ported to GT4Py ("muphys" in icon4py), a reference for a GPU column
    kernel.
- **ECHAM-physics ICON-A** (Giorgetta et al. 2018).
  - Lohmann–Roeckner microphysics with prognostic `qv, qc, qi`; rain and snow are diagnosed in the
    column.
  - Sundqvist cover with the `crs, crt, nex` critical-RH profile (section 3.1) and an inversion
    adjustment.
  - The cover parameters moved to `echam_cov_nml` in 2019. The 2025 namelist overview only lists
    them in its change log; the current AES cover namelist has just `cqx`.
- **Radiation coupling:** RRTMG or ecRad (NWP), RTE+RRTMGP (AES). It receives layer cloud fraction,
  in-cloud `qc` and `qi`, effective radii and an overlap assumption.

**Lessons for SpeedyWeather:**

- Nobody advects a cloud fraction in ICON.
- The ICON configuration closest to SpeedyWeather's resolution and time step is the old ICON-A
  (Sundqvist cover with diagnostic precipitation), which is the same family as Zhao–Carr.
- Prognostic precipitation and instantaneous saturation adjustment presume time steps of minutes,
  many layers and a grid-point state.

### 6. Other schemes, simplest first

| Scheme | Prognostic | Cloud fraction | Fit for SpeedyWeather |
|---|---|---|---|
| SPEEDY, Isca SimCloud (Liu et al. 2021) | none | diagnostic (RH, inversion; Isca: liquid water and `r_e` as functions of `T`) | What exists now (SPEEDY) |
| CloudMicrophysics.jl 0-moment | none (condensate above a threshold removed over `τ`) | — | Equivalent to the current `ImplicitCondensation` |
| Kessler (1969) | `qc, qr` | — | Warm rain only; no ice |
| **Sundqvist et al. (1989)** | total condensate | RH (`b(f)`) | ⭐ the core of the recommendation |
| **Zhao & Carr (1997)** | one condensate, phase from `T` | Sundqvist in microphysics, Xu–Randall in radiation | ⭐ operational in an Eulerian spectral leapfrog model |
| Rasch & Kristjánsson (1998); Zhang et al. (2003) (CAM) | condensate (CAM3: `ql` + `qi`) | RH, with the condensation rate consistent with the fraction change | Good alternative, but also needs the forcing |
| CloudMicrophysics.jl 1-moment | `q_lcl, q_icl, q_rai, q_sno` | — | Julia functions of parameter structs, GPU/AD-friendly. Kessler-type autoconversion (`τ = 10³ s`, threshold `5e-4`) with an optional smooth threshold. Prognostic rain and snow are too much at 8 layers and 40 min |
| Tiedtke (1993) | `ql/qi` + cloud fraction | prognostic | Stiff, bounded variable, many sources; hard at 8 layers |
| Tompkins (2002); Smith (1990), PC2 (Wilson et al. 2008) | PDF moments, or fractions | statistical or prognostic | Overkill |
| IFS (Forbes et al. 2011), ICON graupel, two-moment | 5–6 species (+ numbers) | — | Overkill; the IFS implicit multi-species solver is a useful design reference |

Breeze.jl (Oceananigans-based) wraps CloudMicrophysics.jl. Its issues #1012 and #1013 document
pitfalls when saturation adjustment is combined with one-moment rain.

## 7. Recommendations

### 7.1 Choice

Adopt the **Zhao–Carr structure with Sundqvist physics**:

- **One condensate variable, `q_c`,** fused with the atmospheric prognostics (7.2.1). Its transforms
  ride in the batched calls that already exist, so it adds layers but no launches. Physics
  tendencies are transformed together with the advection tendency.
- **Diagnostic precipitation.** It needs no sedimentation stability (CFL) treatment, which matters
  because rain falls about 12 km in one 40 min step.
- **Column-local, smooth formulas**, which suit the fused GPU kernel and Enzyme/Reactant.
- **Proven numerics.** GFS ran exactly this kind of scheme in an Eulerian spectral leapfrog model.
- **Reuse.** It reuses the existing condensation and precipitation code.

### 7.2 Stage 0: prerequisites

#### 7.2.1 The condensate variable: fused with the atmospheric prognostics

*Decision:* `q_c` is a dedicated prognostic variable `vars.prognostic.cloud_condensate` next to
humidity, declared by the cloud scheme and fused into the parents that already exist. It is not
a `Tracer`.

Why not a tracer: tracers are not fused. Each costs one separate spectral→grid and three
separate grid→spectral transforms per step through the shared `:a/:b` scratch slots, is driven by
a runtime `Dict` loop, and has no GPU or Reactant test. A fused variable adds layers to the
batched calls but no launches, and runs on the code path every GPU test already exercises for
humidity. Transform work per step at 8 layers:

| Route | spectral→grid | grid→spectral |
|---|---|---|
| Today, fused atmosphere | 33 + 16 layers, 2 calls | 73 layers, 1 call |
| Condensate as `Tracer` | +8 layers, +1 call | +24 layers, +3 calls |
| Condensate fused with humidity | +8 layers, +0 calls | +24 layers, +0 calls |

What the scheme declares in `variables(::CloudScheme, model)`, mirroring the humidity lines in
`models/primitive_wet.jl:153-161`, with the step counts from `get_nsteps(model.time_stepping, model)`
as `variables(::Tracer, model)` does:

- `PrognosticVariable(:cloud_condensate, SpectralXYZT(ps), fuse = :prognostic)`;
- `GridVariable(:cloud_condensate, GridXYZT(pg), fuse = :grid)`;
- `TendencyVariable(:cloud_condensate, SpectralXYZT(ts), fuse = :spectral_tendencies)`;
- `TendencyVariable(:cloud_condensate, GridXYZT(tg), namespace = :grid, fuse = :grid_tendencies)`;
- `DynamicsVariable`s `uqc`, `vqc` on the grid (`GridXYZT(tg)`, namespace `:grid`, fuse
  `:grid_tendencies`) and in spectral space (`SpectralXYZT(ts)`, fuse `:spectral_tendencies`).

Because `all_variables` lists the model's variables first and then the components in field order,
these append one slot block to each parent in the same relative position, and the alignment
assertions in `Variables(model)` pass. The spectral and the grid variable must be declared in the
same order within the scheme; a mistake fails loudly at construction, not at run time.

Checked with a stub component at `3c033168` (T31, 8 layers): the `:prognostic` and `:grid`
parents grow from 33 to 41 slots and the tendency parents from 73 to 97; `cloud_condensate`
lands at slots 34:41 and 74:81 with `uqc`, `vqc` at 82:89 and 90:97; the grid copy has size
`(npoints, 8, 2)`; `tendency_names` includes `cloud_condensate`, so the time stepping and the
tendency scaling see it without any change.

What then works without any change: the batched transforms in both directions, tendency scaling,
the leapfrog shift of the grid copy, time stepping (it runs over all top-level tendency names),
vertical advection (generic over `vars.grid[name]`), `reset_tendencies!`, `copy!`, restarts
(`output/restart.jl:47` materializes the whole `prognostic` group), zero initial conditions, and the step selection: physics reads the lagged step through
`get_prognostic_step(vars.grid.cloud_condensate, …)` exactly as it reads humidity.

Four hooks in the dynamical core remain, each a copy of the humidity line:

1. the `+q_c D` term and the products `(u q_c, v q_c)` into the named grid slots, as
   `humidity_grid_tendency!` (`dynamics/tendencies.jl:778-789`);
2. the spectral flux divergence `−∇·(u q_c, v q_c)` after the batched transform, as
   `humidity_spectral_tendency!` (`dynamics/tendencies.jl:795-803`);
3. hyperdiffusion, which names vorticity, divergence, temperature and humidity explicitly
   (`dynamics/horizontal_diffusion.jl:331-355`);
4. hole filling of the grid copy next to humidity's (`time_stepping/transform.jl:153-156`).

*As implemented:* hook 4 is replaced by the fill from vapour inside the scheme (no clipping of the
grid copy), and the vertical advection call site (`dynamics/vertical_advection.jl:31-34`) is a
fifth hook. See the revision log.

Write these over a compile-time tuple of *advected scalar names* derived from the `Variables`
type, as `_tendency_names` is (`variables/variables.jl:641-648`), so that humidity and condensate
share one code path and a further species later costs nothing. Inside the fused column kernel
the scheme addresses the variable by its literal name; `model.tracers` is not available there
anyway.

Two rules of the fuse machinery the implementation must respect:

- **Shape.** In a parent with a step dimension, a member without one collapses to a single
  slot, so a grid copy declared `GridXYZ()` instead of `GridXYZT(pg)` becomes a 2D field with
  steps (`variables/dimensions.jl:252-262`). For the `:prognostic`↔`:grid` pair the alignment
  assertion catches this at construction (checked: slot range 34:41 against 34:34 is rejected);
  a fuse group without an aligned partner is not checked. Always use the `…XYZT` types with the
  step counts from `get_nsteps`, and pin the shape in a test (see Testing).
- **Reactant.** Several in-place broadcasts into different views of one fused buffer corrupt the
  data under Reactant (`dynamics/scaling.jl:20-35`). Hole filling and resets of fused members act
  on each member once per step, or on the parent.

What does not need fusion: the 3D cloud state (7.2.5) and the radiation work arrays. Column
kernels read them coalesced over `ij` either way, and per-g-point `Grid4D` arrays can only share
a parent when their trailing size agrees, which fails for ecCKD models whose longwave and
shortwave g-point counts differ.

Positivity, unchanged from the first draft:

- Physics reads `max(q_c, 0)`.
- GFS's **fill from vapour** as a tendency: a negative `q_c` is moved to 0 from `q` with latent
  heating, which conserves total water and enthalpy.
- Optionally the upwind or WENO vertical advection (`dynamics/vertical_advection.jl:8-11`). It is
  a model-wide choice and reduces overshoots at cloud edges.

#### 7.2.2 The time step physics sees

Column kernels get `time_stepping.Δt`, but leapfrog advances the state by `2Δt` after the first
two steps (`default_time_step`, `time_stepping/steppers/leapfrog.jl:246-256`). `ImplicitCondensation`
therefore removes two thirds of the supersaturation per step, not one third, and any condensate
sink bounded by `q_c/Δt` makes `q_c` negative. The GPU-adapted stepper `LeapfrogCore` carries
only `Δt` (`time_stepping/steppers/leapfrog.jl:171-179`).

- Add the prognostic step to `LeapfrogCore` and an accessor, e.g. `physics_time_step(time_stepping)`,
  returning `default_time_step` (`2Δt` for leapfrog, `Δt` for every other stepper,
  `time_stepping/steppers/general.jl:5-9`), following `dynamics/horizontal_diffusion.jl:430-432`,
  which already uses the prognostic step for the implicit diffusion.
- Every bounded sink of the cloud scheme uses that step, written `Δt_p` below. The first two
  steps advance by `Δt/2` and `Δt`, so the bound is conservative there.
- Whether `ImplicitCondensation` and the other schemes switch too is a separate decision; it
  changes their results.

#### 7.2.3 SPPT

`StochasticallyPerturbedParameterizationTendencies` multiplies the u, v, T and q tendencies
(`parameterizations/stochastic_physics.jl:74-84`). With a condensate tendency it has to multiply
that one with the same pattern, otherwise water and enthalpy are no longer conserved between `q`,
`T` and `q_c`. Add the condensate, and in general every advected scalar, to the SPPT loop.

#### 7.2.4 Ice thermodynamics, kept local

- A mixed-phase `q_s` as a blend of liquid and ice saturation over 0 to −20 °C, with the matching
  `∂q_s/∂T` for the implicit factor, inside the cloud scheme. `saturation_humidity` is shared by
  convection, the diagnostic clouds, surface fluxes and condensation; changing it globally
  changes the whole model. A global switch, if wanted later, is an `EarthAtmosphere` option with
  the current behaviour as default.
- `L = Lᵥ + f_ice L_f`, never `latent_heat_sublimation`: that constant (2801 kJ/kg) is not
  `Lᵥ + L_f` (2831 kJ/kg). Fix the constant in a separate PR.
- Keep the column enthalpy budget closed, as `test/parameterizations/large_scale_condensation.jl`
  tests today.

#### 7.2.5 The cloud state as the interface to radiation

The cloud scheme writes, per layer, as 3D parameterization variables:

- `cloud_fraction`;
- in-cloud liquid and ice water (mixing ratio, or path per layer);
- effective radii of liquid and ice (GFS rules: 10 µm over ocean, `5 + 5 f_ice` µm over land;
  ice from temperature and ice water content, or a fixed value to start);
- the overlap parameter `α` between adjacent layers from a decorrelation length `L`,
  `α = exp(−Δz/L)` (ecRad's default `L` is about 2 km), and optionally the fractional standard
  deviation of in-cloud condensate (ecRad uses 1).

Every radiation scheme reads this state and nothing else from the cloud scheme. The call order
already runs condensation before radiation (`models/primitive_wet.jl:122-134`), so radiation
sees this step's clouds. Column cover by maximum overlap and the highest cloudy layer are
derived from the state for the existing `cloud_cover`, `cloud_top` and `cloud_top_height`
outputs.

#### 7.2.6 Radiation with layer clouds

*One-band schemes, first.* A `PrognosticClouds <: AbstractShortwaveClouds` returns the
NamedTuple that `clouds!` returns today (cover, top, albedos), so `OneBandShortwave` works
unchanged on day one. Then:

- **Shortwave:** reflection and absorption in every layer from the optical depth. This needs an
  adding method, since the current code reflects once at `cloud_top`. *Implemented 2026-10-08:*
  `CloudyShortwaveRadiativeTransfer`, absorbing two-stream cloud layers for diffuse light, random
  overlap, the adding method; `OneBandCloudyShortwave` drops the SPEEDY cloud absorption of the
  background transmissivity
  (`parameterizations/radiation/shortwave_radiation.jl:192-196, 233`). For liquid, `τ = 3 LWP/(2 ρ_w r_e)`.
- **Longwave:** cloud emissivity `ε = 1 − exp(−D κ LWP)`, `D ≈ 1.66` (CCM3 style), multiplied into
  the Frierson transmissivity.
- **Overlap:** maximum-random with the layer cloud fractions.

*All-sky ecCKD through NumericalRadiation, second.* State of NumericalRadiation 0.1.1:

- Cloud optics are ready and GPU-fit: `SpectralCloudOptics` maps ecRad's Mie droplet and Baum
  ice tables onto the ecCKD g points (`:ecrad` averaging, delta-Eddington) on effective-radius
  nodes; at run time `effective_radius_bracket`, `cloud_layer_optics` and `add_scattering_layer`
  are allocation-free and branchless, and the tables are `Adapt`-able.
- The all-sky solvers are not: `CloudOverlapShortwave` and `CloudOverlapLongwave` allocate per
  column and per g point and are documented as diagnostic solvers, not Tripleclouds or McICA.
  The clear-sky extension only became allocation-free through `streaming_longwave_fluxes!` and
  the caller-owned `ShortwaveColumnScratch`.

An `AllSkyEcCKDRadiation` therefore needs, upstream in NumericalRadiation, a streaming two-region
(clear and cloudy) or three-region (Tripleclouds) adding solver with caller-owned scratch for
both streams. In SpeedyWeather it then follows the clear-sky extension exactly: per-g-point
optical depths for the clear and the cloudy region as `Grid4D` parameterization variables (twice
the clear-sky arrays), the cloud optics folded in per layer from the cloud state of 7.2.5 with one
radius bracket per layer and phase, and the fluxes blended by the solver. Cost: clear-sky ecCKD
is about 80× the one-band pair per column; all-sky adds a factor two to three.

#### 7.2.7 Radiation call frequency

Run radiation every `N` steps and hold its heating in a 3D parameterization variable in between.
This was future work in `docs/dev/2026-09/numericalradiation-extension.md`; with all-sky ecCKD as
a target it belongs here. The fused column kernel keeps a cheap "add the stored heating" branch
on the other steps.

#### 7.2.8 Output

Liquid and ice water path, 3D cloud fraction, condensate, and the process rates for debugging,
as output variable definitions like those in `output/variables/precipitation.jl`.

### 7.3 Stage 1: Zhao–Carr with relaxation condensation (recommended first implementation)

This is one column kernel that replaces or extends `ImplicitCondensation`, working from the top down:

1. **Phase.** A smooth ice fraction `f_ice(T)`, either GFS radiation's ramp (0 to −20 °C) or
   Zhao–Carr's thresholds at 0 and −15 °C written as a ramp. Avoid the branchy `IW` seeding logic
   for GPU and AD.
2. **Condensation and cloud evaporation.** Reuse the implicit relaxation of `ImplicitCondensation`:
   `δq = (q − RH_c q_s) / ((1 + L/cₚ RH_c ∂q_s/∂T) τ)`.
   - Supersaturated (`δq > 0`): condense into `q_c` instead of raining out.
   - Subsaturated with `q_c > 0`: evaporate cloud only in the clear fraction `1 − C` of the cell,
     with its own time scale `τ_evap ≥ τ`, limited by `q_c`. Relaxing the whole cell towards
     `RH_c q_s` would clear advected cloud within one or two steps at 70 % RH and no anvil would
     survive; GFS keeps the cloud unless RH falls below the critical value.
   - This is a saturation adjustment towards `RH_c` with timescale `τ`. It needs no forcing `M`, is
     already energy-checked and is stable under leapfrog.
3. **Cloud fraction:** Xu–Randall in GFS form, `C(q_c, f, q_s)` (section 4.5). Use it for both
   autoconversion and radiation. Do *not* use Sundqvist `b(f)` here: with the relaxation, RH stays
   near `RH_c` in cloudy air, so `b(f)` would be degenerate.
4. **Conversion to precipitation:**
   - liquid: Sundqvist with `F_co` from the flux above and `F_BF` (section 3.3);
   - ice: `P_saut` from section 4.4.
   - Apply both in **exponential form**, `Δq_c = q_c (1 − exp(−k Δt_p))` with `k = P/q_c` and `Δt_p`
     the prognostic step of 7.2.2 (`2Δt`, 80 min at T32), so `q_c` stays non-negative even when
     `F` is 5 to 10. The tendency is `Δq_c/Δt_p`.
5. **Precipitation sweep:** keep the existing melting, reevaporation and snow code. Optionally add
   collection of cloud by falling precipitation (`P_racw`, `P_saci`).
6. **Outputs:** tendencies for `q`, `T` and `q_c`, rain and snow rates, and the cloud state of 7.2.5
   (fraction, in-cloud water and ice, effective radii). Keep `cloud_top` and the column cover,
   derived from the cloud state, until no radiation scheme needs them.

**Parameters**, starting from GFS or the current SpeedyWeather values:

| Parameter | Start value |
|---|---|
| `RH_c` | 0.95 |
| `τ` | `3Δt` |
| `τ_evap` | `2τ` |
| `C₀` | `1e-4 s⁻¹` |
| `m_r` | `3e-4` |
| `c₁`, `c₂` | 300, 0.5 |
| `c_saut` | `3e-4 s⁻¹` |
| `w_min` | `1e-5 × p/1000 hPa` |
| Xu–Randall `α` | 2000 |

Stage 1 is *not* a published scheme. It is Zhao–Carr with its condensation closure replaced by
SpeedyWeather's relaxation. Call it that in the docs.

### 7.4 Stage 2: the Sundqvist condensation closure

This gives partial cloudiness before saturation and condensation consistent with `b(f)`. It needs the
non-cloud forcing `A_T`, `A_q`, `A_p` (section 3.2). Physics runs *before* dynamics and the dynamical
tendencies only exist in spectral space, so there are three options:

1. **The GFS approach (recommended), restated for tendency-based physics.** GFS stores the state
   after its own update; SpeedyWeather's physics adds tendencies. The equivalent is:
   - at each call store `X_ref = X_lagged + Δt_p F_cond` for `T`, `q` and `pₛ`, where `F_cond` is
     this scheme's own condensation tendency of that call (not its precipitation processes, as GFS
     stores the state after `gscond`, before `precpd`) and `Δt_p` the prognostic step of 7.2.2;
   - use `A_X = (X_lagged − X_ref)/Δt_p` with the reference from *two calls before*. Physics reads the
     lagged step, so consecutive calls are one step apart and alternate leapfrog parity; two calls
     before has the same parity, as GFS's `tp ← tp1 ← t`. *(Corrected 2026-10-08; the first draft
     said consecutive calls.)*
   - Cost: two copies (parities) of two 3D and one 2D field. Declare them as `PrognosticVariable`s
     *without* a tendency: they are allocated, saved in restarts and copied, but never stepped,
     since stepping runs only over names that have a tendency. This is the home for a state with
     memory; the `parameterizations` group is declared memoryless.

   *As implemented:* `SundqvistClosure` with the references `temperature_reference`,
   `humidity_reference` (`GridXYZT(2)`) and `surface_pressure_reference` (`GridXYT(2)`) in the
   `clouds` namespace of the prognostic variables, step 1 = two calls before, step 2 = last call,
   shifted per column in the kernel. The supply is zero until two calls have stored a reference.
2. **Reuse the two grid time levels SpeedyWeather already keeps**, minus the scheme's own previous
   tendency. No new `T`/`q` state is needed, but the difference includes the leapfrog computational
   mode. Experimental.
3. **Transform the dynamics tendencies to grid space before physics.** Extra transforms and a
   reordering; not recommended.

With the forcing available:

- use `u(σ)` = 0.9 (GFS) or the ECHAM form, and retune for 8 thick layers;
- use `b(f)` in the closure and in autoconversion;
- radiation can stay on Xu–Randall, as in GFS.

### 7.5 Stage 3: convective condensate

Betts–Miller relaxes towards a reference profile and has no updraft to detrain from. The options are:

- (a) detrain a fraction `f_det` of the convective condensate as `q_c` in the top one or two
  convective layers, below the level of zero buoyancy (anvil source);
- (b) keep the SPEEDY precipitation-based convective cover for radiation alongside the stratiform
  cloud.

Mass-flux schemes (SAS in GFS, Tiedtke in IFS and ICON) detrain updraft condensate into the cloud
tracer. Either SpeedyWeather option is an ad-hoc substitute, so treat it as a *decision*.

*Decided and implemented (2026-10-08):* option (a), `BettsMillerConvection(; detrainment)`, a fraction
of the deep convective precipitation detrained as condensate into the layer of zero buoyancy,
default 0 so the default model is unchanged. The condensate gets no freezing heat in cold layers,
like condensate advected into them. The detrainment fraction is a tuning parameter (7.8).

### 7.6 Numerics checklist

- **Bounded sinks:** every term that can remove `q_c` within a step is implicit or exponential, or is
  capped at `q_c/Δt_p` with the prognostic step of 7.2.2, never `q_c/Δt`. There is no atmospheric
  `filter!` to fix overshoots afterwards; `adjustments-not-tendencies.md` is the precedent if one
  is added.
- **SPPT:** the condensate tendency is perturbed with the same pattern as `q` and `T` (7.2.3).
- **Conservation tests:** column total water (`q + q_c` + precipitation) and enthalpy, extending the
  existing condensation budget tests.
- **Differentiability:**
  - clamp `f ≤ 1 − ε` in `b(f)`;
  - clamp the Xu–Randall denominator (GFS uses `≥ 1e-4`);
  - prefer the smooth Sundqvist threshold to hard thresholds;
  - use a phase ramp instead of the `IW` branches.
  - These choices also keep the parameters tunable by gradients, see 7.8.
- **GPU:** one fused column kernel, fixed loop lengths, no data-dependent early exits beyond what
  `ImplicitCondensation` already has. The condensate loop runs over all layers (`q_c ≥ 0`
  everywhere) rather than branching on cloud presence. Work arrays are parameterization
  variables, nothing is allocated in the kernel. The effective radius is bracketed once per
  layer and phase, not per g point. `model.tracers` and other `Dict`s are not available in the
  kernel; variables are addressed by literal name.
- **Reactant:** no repeated in-place broadcasts into views of one fused parent (7.2.1).
- **Spectral ringing:** `q_c` is intermittent. Hyperdiffusion damps the smallest scales. GFS shows
  that clipping on read plus fill-from-vapour is enough for an operational spectral model.

### 7.7 What not to do first, and why

| Approach | Why not now |
|---|---|
| Prognostic rain and snow (ICON, CloudMicrophysics 1-moment) | Sedimentation at about 5 m s⁻¹ with 40 min steps and 8 layers needs implicit or semi-Lagrangian fall. Little gain at T31. |
| Two-moment schemes | Need aerosol, which SpeedyWeather does not have. |
| Prognostic cloud fraction (Tiedtke, PC2) | Stiff, bounded variable with many source terms. |
| Instantaneous saturation adjustment (ICON `satad`) | Incompatible with a spectral state updated by leapfrog tendencies. The relaxation with `τ ≥ Δt` is SpeedyWeather's equivalent. |

### 7.8 Tuning with SpeedyCalibration.jl

The aim is to tune the cloud parameterization with SpeedyCalibration.jl (`~/SpeedyCalibration.jl`,
[github.com/SpeedyWeather/SpeedyCalibration.jl](https://github.com/SpeedyWeather/SpeedyCalibration.jl))
using Enzyme. We reuse its core calibration loop and change SpeedyCalibration where it does not
fit the cloud problem yet.

**What SpeedyCalibration does** (its `main` branch, read on 2026-10-08):

- *Online statistical gradient estimation:* the model runs continuously; every `steps_per_sample`
  steps the gradient of one `SpeedyWeather.timestep!` is taken with Enzyme reverse mode
  (`Duplicated(variables)`, `Duplicated(model)`), averaged over `samples_per_batch` samples per
  batch, scaled, clipped and applied with an Optimisers.jl optimizer. Single-step gradients avoid
  differentiating through the chaotic trajectory.
- *Parameters:* a `ParamSpec` per parameter with a property path into the model, bounds (sigmoid
  reparameterization) and a gradient scale; values are written back with `set_by_path!` and
  `reconstruct`.
- *Loss:* `LossConfig`, a weighted mean squared error of global-mean fluxes from
  `variables.parameterizations` (outgoing shortwave and longwave, surface shortwave and longwave up
  and down) against Trenberth-type targets (`TRENBERTH_LOSS`).
- *Model:* it builds `PrimitiveWetModel(spectral_grid; planet)` itself and pins SpeedyWeather 0.21.

**Changes to SpeedyCalibration** (fine to make; the calibration loop stays):

1. SpeedyWeather 0.23: the component paths changed with the `Radiation` bundle, e.g.
   `[:radiation, :shortwave, :clouds, …]` instead of `[:shortwave_radiation, :clouds, …]`, and
   `trunc` became `truncation`.
2. Model components passed in, instead of the fixed constructor: the cloud scheme, the cloudy
   radiation and the convective detrainment.
3. Cloud loss terms: global and zonal-mean cloud cover, liquid and ice water path, precipitation.
   Cloud radiative effects need clear-sky fluxes, which SpeedyWeather does not compute yet
   (Future work).
4. *Multi-step windows, possibly.* A single-step gradient only sees a parameter's effect within that
   step. Radiation reads the cloud state that the condensation scheme writes in the same step, so
   the microphysics parameters do reach the fluxes within one step. Their effect through the
   prognostic condensate, whose lifetime is hours, is missed, and the gradient is biased. A window
   of a few differentiated steps as an option of the loop may be needed. Check the single-step
   gradient against finite differences of long-run statistics first.

**What SpeedyWeather needs:**

- Enzyme differentiability tests: reverse mode through one `timestep!` with `Duplicated(model)`
  gives finite gradients for every candidate parameter, checked against finite differences, with
  the cloud scheme, the cloudy radiation and (Stage 2) the reference-state bookkeeping.
- *Candidate parameters:* `relative_humidity_threshold` (Stage 1) or the critical relative
  humidity (Stage 2), `time_scale`, `evaporation_time_scale`, `autoconversion_rate`,
  `autoconversion_water`, `ice_autoconversion_rate`, `cloud_fraction_coefficient`, the effective
  radii, the convective detrainment, and in radiation the asymmetry factor, the single-scattering
  albedo and the mass absorption coefficients. All are `@param` with bounds.
- *Non-smooth points:* the formulations were chosen smooth (7.6), but the cloud-presence thresholds
  (`q_c > 10⁻⁶ p/1000 hPa`, cloud fraction ≥ 0.001), the `min`/`max` limiters of the sinks and the
  clamps give zero or one-sided gradients where they are active.
- *Cost:* Enzyme compile time on the full `PrimitiveWetModel` is substantial.

## Summary of changes (Stage 1 and its Stage 0 parts)

- `parameterizations/cloud_condensation.jl` (new): `PrognosticCloudCondensation <: AbstractCondensation`
  with the fused condensate variables (`cloud_condensate_variables`), the cloud state
  (`cloud_state_variables`), the mixed-phase saturation, Xu–Randall cloud fraction, Sundqvist and
  Lin autoconversion, and the column sweep `cloud_condensation!`.
- `variables/variables.jl`: `ADVECTED_SCALARS` and `_advected_scalars`.
- `dynamics/tendencies.jl`, `dynamics/horizontal_diffusion.jl`, `dynamics/vertical_advection.jl`:
  generated drivers over the advected scalars for the grid and spectral halves of the flux-form
  advection, hyperdiffusion and vertical advection; the humidity functions are thin wrappers.
- `time_stepping/steppers/leapfrog.jl`: `default_time_step(::LeapfrogCore) = 2Δt`.
- `parameterizations/stochastic_physics.jl`: SPPT perturbs the condensate; scratch path fixed.
- `parameterizations/large_scale_condensation.jl`: precipitation variables shared by both schemes.
- `parameterizations/radiation/clouds.jl`: `PrognosticClouds` (maximum-random column cover,
  two-stream cloud albedo from the in-cloud optical depth).
- `parameterizations/radiation/longwave_transmissivity.jl`: `CloudyLongwaveTransmissivity`
  wrapping a clear-sky transmissivity; `OneBandLongwave` now declares its transmissivity's variables.
- `dynamics/spectral_grid.jl`: `max_transform_batch` 12L+1, `primitive_wet_tendency_batch` 9L+1.
- `output/variables/clouds.jl` (new): condensate, cloud fraction, liquid and ice water path.

Ten-day global means at T32 L8 (CPU, default initial conditions, `OneBandCloudyShortwave` and
`OneBandCloudyLongwave` for all prognostic configurations):

| Configuration | Cloud cover | Outgoing SW [W/m²] | OLR [W/m²] | LWP [g/m²] | IWP [g/m²] |
|---|---|---|---|---|---|
| Default (`ImplicitCondensation`, `DiagnosticClouds`) | 0.60 | 126 | 252 | – | – |
| Stage 1, relaxation closure | 0.26 | 38 | 250 | 5.8 | 21.7 |
| Stage 2, Sundqvist closure, u = 0.9 | 0.27 | 37 | 250 | 5.3 | 21.4 |
| Stages 2 + 3, detrainment 0.2 | 0.35 | 45 | 244 | 5.7 | 30.5 |
| Stages 2 + 3, detrainment 0.2, u = 0.8 | 0.38 | 54 | 241 | 9.1 | 34.4 |

The per-layer shortwave alone changes little against the reflection at the cloud top (38 against
40 W/m² with Stage 1): the cloud amount, not the transfer, limits the shortwave cloud effect.

Added on 2026-10-08:

- `parameterizations/cloud_condensation.jl`: the closure as a component, `RelaxationClosure` and
  `SundqvistClosure` with its reference state, critical relative humidity and Sundqvist cloud
  fraction.
- `parameterizations/radiation/shortwave_radiation.jl`: `CloudyShortwaveRadiativeTransfer`,
  `two_stream_diffuse_layer`, `OneBandCloudyShortwave`; `OneBandShortwave` declares its radiative
  transfer's variables.
- `parameterizations/radiation/longwave_radiation.jl`: `OneBandCloudyLongwave`.
- `parameterizations/convection.jl`: `detrainment` of `BettsMillerConvection`.

First 10-day run at T31 L8 (CPU, from the default initial conditions): global mean cloud cover
0.26, liquid water path 5.7 g/m², ice water path 22 g/m², large-scale rain 0.21 mm/day and snow
0.16 mm/day, OLR 249 W/m², outgoing shortwave 40 W/m². The cloud cover and the liquid water path
are far below observed values; tuning is open. For comparison, the default model (ImplicitCondensation,
DiagnosticClouds) gives cloud cover 0.60 and outgoing shortwave 126 W/m² over the same 10 days.

Cost on an NVIDIA A40 (8 layers, prognostic clouds including cloud radiation against the default
model, 200 steps after spin-up):

| Truncation | Transform | Default [ms/step] | Prognostic clouds [ms/step] |
|---|---|---|---|
| T32 | matrix | 1.42 | 1.63 (+15 %) |
| T128 | FFT | 5.71 | 6.09 (+7 %) |

At T128 the new batch sizes 41 (prognostic) and 97 (tendencies) are planned on first use and run
batched, not serially. GPU and CPU agree in global means to about 10⁻³ after two days at T32.

## Testing and verification

*Status:* `test/parameterizations/cloud_condensation.jl` implements all items below except the
long runs and the differentiability test, plus tests of the Sundqvist closure (cloud fraction,
critical humidity, the reference bookkeeping, the closure formula against a prescribed supply,
budgets), the two-stream cloud layer, the per-layer shortwave (energy conservation, reduction to
the one-band transfer in clear sky) and the convective detrainment (water budget per column). The
GPU tests are in `test/GPU/primitive_wet.jl`, one with the relaxation closure and one with the
Sundqvist closure, per-layer shortwave and detrainment.

- **Unit tests on prescribed columns:**
  - a supersaturated layer produces `q_c`;
  - `q_c` decays by autoconversion with the expected rate;
  - cloud evaporates in subsaturated air;
  - the water and enthalpy budgets close;
  - `q_c ≥ 0` after one step for extreme rates.
- **Limit test:** with instantaneous autoconversion the new scheme reproduces `ImplicitCondensation`'s
  precipitation and `q`, `T` tendencies on a prescribed column.
- **Variable layout:** `Variables(model)` with the cloud scheme passes the fuse alignment
  assertions; the condensate's grid copy has `nlayers` layers and the step dimension; the
  spectral and grid slots of `cloud_condensate` agree, and a grid copy without the step
  dimension is rejected (pins the behaviour observed in the stub check of 7.2.1). A run of a few
  steps with the scheme's tendencies set to zero agrees in all other variables with a run without
  the scheme (to `rtol = 1e-5`, not bitwise: the wider batched transforms change rounding).
- **Time step:** a sink at the maximum rate leaves `q_c ≥ 0` after a full leapfrog step (`2Δt`),
  not only after `Δt`.
- **GPU and differentiability tests** like those of the existing parameterizations: a fused
  column-kernel run in `test/GPU/` and a parameter-AD test as in
  `test/differentiability/primitivewet.jl`.
- **Long runs at T31 and T63.** Compare global means against observations: total cloud fraction of
  about 0.6–0.7, and CERES-EBAF cloud radiative effects of about −45 to −47 W m⁻² (SW),
  +26 to +28 W m⁻² (LW) and about −20 W m⁻² (net). Check that precipitation stays close to the
  current scheme. Also check plausible liquid and ice water paths, and the zonal-mean `q_c`
  structure (storm tracks, tropical upper troposphere).

## Documentation changes (planned)

A new docs page for the cloud scheme, following `docs/src/large_scale_condensation.md`. Document
the physics time step accessor and the cloud state in `docs/src/parameterizations.md`, and that
components may declare fused members in `docs/src/variable_system.md`. Update
`docs/src/radiation.md` once clouds enter the longwave and again for the all-sky ecCKD scheme.

## Known limitations of this review

- **Primary papers not read:** the Sundqvist et al. (1989) and Zhao & Carr (1997) PDFs were not
  accessible (AMS returned 403). Formulas and constants come from the GFS code and its embedded
  documentation, which follow the papers but carry GFS tuning (e.g. `c₁`, ice autoconversion).
- **Unverified constants:** the original Xu & Randall (1996) constants were not checked against the
  paper; from memory they are `α₀ = 100`, `γ = 0.49`, `p = 0.25`, and GFS's 2000 and ¼ are a
  retune. GFS's constants are verified. The ICON-A default critical-RH values were not found,
  only the parameter names and the refits by Grundner et al. (2022). The GFS Exner-interpolated
  critical RH profile was not checked in `GFS_suite_interstitial_3.F90`.
- **Time step convention:** the existing parameterizations use `Δt` while the state advances by
  `2Δt`. Section 7.2.2 adds the accessor; whether the existing schemes switch is undecided.
- **NumericalRadiation:** the all-sky solvers allocate and are diagnostic. The streaming all-sky
  solver is upstream work and not scheduled.
- **Not tuned:** only a 10-day run exists (see Summary of changes); cloud cover and liquid water
  path are far too low. The parameter values are starting points. Tuning is planned with
  SpeedyCalibration.jl, see 7.8.
- **No stratocumulus cloud;** convective cloud only through the detrainment (Stage 3).
- **Betts-Miller convective snow** carries no latent heat of fusion; found, not fixed.
- **No differentiability test yet** for the new scheme.
- **`ImplicitCondensation` loses water** with reevaporation on, see the appendix; a fix is proposed,
  not applied.

## Future work

- Two condensate variables (`q_l`, `q_i`) with sedimentation of cloud ice, as two more advected
  scalars.
- A convective cloud scheme tied to Betts–Miller.
- `AllSkyEcCKDRadiation` once NumericalRadiation has a streaming all-sky solver (7.2.6).
- Fusing tracers into their own parents so that passive tracers get batched transforms too.
- Prognostic cloud fraction, only if Stage 2 shows systematic cover errors that RH and condensate
  cannot fix.
- Tuning with SpeedyCalibration.jl, see 7.8.
- Clear-sky fluxes of the one-band schemes (a second pass without clouds) for cloud radiative
  effects as tuning targets and output.

## Appendix: water loss in `ImplicitCondensation`

Investigated on 2026-10-08 at `3c033168` plus the changes of this plan (the scheme itself is
unchanged since the base revision).

**Mechanism.** In a layer below precipitation the scheme removes the evaporated rain from the
downward flux in full, `rain_flux_down -= rain_evaporated`, but adds the corresponding humidity
through the condensation relaxation: `δq = min(0, δq_cond) + δq_evap` is divided by
`(1 + Lᵥ/cₚ ∂q*/∂T) × time_scale × Δt`. The humidity gained is therefore the evaporated rain
divided by `time_scale × (1 + Lᵥ/cₚ ∂q*/∂T)`, about 1/3 to 1/9 of it; the rest vanishes. A second,
smaller path: rain from a layer is `-min(0, δq)` of the *net* tendency, so in a layer that both
condenses and reevaporates (relative humidity between the threshold and 100 %) the evaporated rain
is subtracted twice. The column enthalpy test passes because the heating is consistent with the
humidity tendency; the lost water's latent heat stays in the atmosphere as heating.

**History.** The structure dates from the introduction of reevaporation (`a7c23ab7`, 2025-08-25,
"Reevaporation, sublimation and snow fall") and was kept by the snow PR #817. The docs
(`docs/src/large_scale_condensation.md`, "Re-evaporation") state the intent: subtract the
reevaporation from the condensation tendency and use the implicit time stepping as before. That
works for the implicit factor alone but not with the relaxation time scale, while the flux is
reduced in full.

**Magnitude.**

- A prescribed column (rain from two supersaturated layers falling into a 60 % relative humidity
  layer below, 285 K, Float64): surface rain is 43 % of the vapour removed; with `reevaporation = 0`
  100 %.
- The global state after 20 days at T32 L8 (default model, `ImplicitCondensation` called alone on
  every column): vapour removed 0.0801 mm/day, large-scale precipitation 0.0759 mm/day, so
  0.0042 mm/day (5 % of large-scale precipitation, 0.25 % of total precipitation) is lost. That
  is a spurious heating of about 0.12 W/m². Convective rain is not reevaporated and not affected.
  The share grows with the share of large-scale precipitation, e.g. at higher resolution.

**Proposed fix** (tested in a script, not applied): evaporated rain enters humidity in full, as
`rain_evaporated × gρ/Δp`, not relaxed, and capped at saturation over the prognostic time step;
rain from a layer is its gross condensation, `-min(0, δq_cond)` after the implicit relaxation, not
the net tendency; the freezing heat of snow uses the gross condensation too. With this the global
vapour removed equals the large-scale precipitation to 10⁻⁸ and the test column to round-off.
It changes the default model: more water is returned to the lower layers. This needs its own PR
with a test of the column water budget.

`PrognosticCloudCondensation` does not have the problem: rain is its gross autoconversion,
reevaporation enters humidity in full, and its water budget is tested to round-off.

## References

**Sundqvist, Zhao–Carr and GFS**

- Sundqvist, H. (1978): A parameterization scheme for non-convective condensation including
  prediction of cloud water content. *Q. J. R. Meteorol. Soc.*, 104, 677–690.
- Sundqvist, H., E. Berge, J. E. Kristjánsson (1989): Condensation and cloud parameterization studies
  with a mesoscale numerical weather prediction model. *Mon. Wea. Rev.*, 117, 1641–1657.
  [AMS](https://journals.ametsoc.org/view/journals/mwre/117/8/1520-0493_1989_117_1641_cacpsw_2_0_co_2.xml)
- Zhao, Q., F. H. Carr (1997): A prognostic cloud scheme for operational NWP models.
  *Mon. Wea. Rev.*, 125, 1931–1953.
  [AMS](https://journals.ametsoc.org/view/journals/mwre/125/8/1520-0493_1997_125_1931_apcsfo_2.0.co_2.xml)
- Zhao, Q., T. L. Black, M. E. Baldwin (1997): Implementation of the cloud prediction scheme in the
  Eta model at NCEP. *Wea. Forecasting*, 12, 697–711.
- GFS Zhao–Carr code, ccpp-physics `v6.0.0`:
  - [zhaocarr_gscond.f](https://github.com/NCAR/ccpp-physics/blob/v6.0.0/physics/zhaocarr_gscond.f)
  - [zhaocarr_precpd.f](https://github.com/NCAR/ccpp-physics/blob/v6.0.0/physics/zhaocarr_precpd.f)
  - [radiation_clouds.f](https://github.com/NCAR/ccpp-physics/blob/v6.0.0/physics/radiation_clouds.f)
  - [GFS_suite_interstitial_3.F90](https://github.com/NCAR/ccpp-physics/blob/v6.0.0/physics/GFS_suite_interstitial_3.F90)
    (critical RH)
  - default parameters in fv3atm
    [GFS_typedefs.F90](https://github.com/NOAA-EMC/fv3atm/blob/develop/ccpp/data/GFS_typedefs.F90)
- [CCPP Scientific Documentation: GFS Zhao–Carr microphysics](https://dtcenter.ucar.edu/gmtb/users/ccpp/docs/sci_doc_v2/GFS_ZHAOC.html)
- [EMC GFS documentation and implementation history](https://emc.ncep.noaa.gov/emc/pages/numerical_forecast_systems/gfs/documentation.php)
- Xu, K.-M., D. A. Randall (1996): A semiempirical cloudiness parameterization for use in climate
  models. *J. Atmos. Sci.*, 53, 3084–3102.
- Lin, Y.-L., R. D. Farley, H. D. Orville (1983): Bulk parameterization of the snow field in a cloud
  model. *J. Climate Appl. Meteor.*, 22, 1065–1092.

**ICON and ECHAM**

- [ICON Namelist Overview (2025-04-22)](https://icon-training-2025-scripts-rendering-cc74a6.gitlab-pages.dkrz.de/_downloads/6bd9f7b7262f616f699617df3098cc24/Namelist_overview.pdf)
- Giorgetta, M. A., et al. (2018): ICON-A, the atmosphere component of the ICON Earth system model:
  I. Model description. *JAMES*, 10, 1613–1637.
  [doi](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2017MS001242)
- Hohenegger, C., et al. (2023): ICON-Sapphire. *GMD*, 16, 779–811.
  [GMD](https://gmd.copernicus.org/articles/16/779/2023/)
- Müller, W. A., et al. (2025): ICON-XPP. *GMD*, 18, 9385.
  [GMD](https://gmd.copernicus.org/articles/18/9385/2025/)
- Grundner, A., et al. (2022): Deep learning based cloud cover parameterization for ICON. *JAMES*.
  [arXiv:2112.11317](https://arxiv.org/abs/2112.11317)
- Lohmann, U., E. Roeckner (1996): Design and performance of a new cloud microphysics scheme
  developed for the ECHAM general circulation model. *Clim. Dyn.*, 12, 557–572.
- Xu, K.-M., S. K. Krueger (1991): Evaluation of cloudiness parameterizations using a cumulus
  ensemble model. *Mon. Wea. Rev.*, 119, 342–367.
- Mauritsen, T., et al. (2019): Developments in the MPI-M Earth System Model version 1.2 (MPI-ESM1.2)
  and its response to increasing CO₂. *JAMES*, 11, 998–1038.
- [icon4py muphys (GT4Py graupel port)](https://github.com/C2SM/icon4py/pull/1454)

**Other schemes**

- Tiedtke, M. (1993): Representation of clouds in large-scale models. *Mon. Wea. Rev.*, 121,
  3040–3061.
  [AMS](https://journals.ametsoc.org/view/journals/mwre/121/11/1520-0493_1993_121_3040_rocils_2_0_co_2.xml)
- Rasch, P. J., J. E. Kristjánsson (1998): A comparison of the CCM3 model climate using diagnosed and
  predicted condensate parameterizations. *J. Climate*, 11, 1587–1614.
- Zhang, M., W. Lin, C. S. Bretherton, J. J. Hack, P. J. Rasch (2003): A modified formulation of
  fractional stratiform condensation rate in the NCAR Community Atmospheric Model (CAM2).
  *J. Geophys. Res.*, 108(D1), 4035.
- Tompkins, A. M. (2002): A prognostic parameterization for the subgrid-scale variability of water
  vapor and clouds in large-scale models and its use to diagnose cloud cover. *J. Atmos. Sci.*, 59,
  1917–1942.
  [AMS](https://journals.ametsoc.org/view/journals/atsc/59/12/1520-0469_2002_059_1917_appfts_2.0.co_2.xml)
- Smith, R. N. B. (1990): A scheme for predicting layer clouds and their water content in a general
  circulation model. *Q. J. R. Meteorol. Soc.*, 116, 435–460.
- Wilson, D. R., et al. (2008): PC2: A prognostic cloud fraction and condensation scheme.
  *Q. J. R. Meteorol. Soc.*, 134, 2093–2107. [doi](https://rmets.onlinelibrary.wiley.com/doi/10.1002/qj.333)
- Forbes, R. M., A. M. Tompkins, A. Untch (2011): A new prognostic bulk microphysics scheme for the
  IFS. ECMWF Tech. Memo. 649.
  [PDF](https://www.ecmwf.int/sites/default/files/elibrary/2011/9441-new-prognostic-bulk-microphysics-scheme-ifs.pdf)
- Kessler, E. (1969): On the distribution and continuity of water substance in atmospheric
  circulations. *Meteor. Monogr.*, 10(32).
- Liu, Q., et al. (2021): SimCloud version 1.0: a simple diagnostic cloud scheme for idealized
  climate models (in Isca). *GMD*, 14, 2801.
  [GMD](https://gmd.copernicus.org/articles/14/2801/2021/)
- [CloudMicrophysics.jl](https://github.com/CliMA/CloudMicrophysics.jl)
  ([0-moment](https://clima.github.io/CloudMicrophysics.jl/dev/Microphysics0M/),
  [1-moment](https://clima.github.io/CloudMicrophysics.jl/dev/Microphysics1M/))
- [Breeze.jl](https://github.com/NumericalEarth/Breeze.jl)

**Radiation**

- Heymsfield, A. J., G. M. McFarquhar (1996): High albedos of cirrus in the tropical Pacific warm
  pool. *J. Atmos. Sci.*, 53, 2424–2451. (GFS ice effective radius.)
- Hogan, R. J., A. J. Illingworth (2000): Deriving cloud overlap statistics from radar.
  *Q. J. R. Meteorol. Soc.*, 126, 2903–2909. (Overlap parameter `α` and decorrelation length.)
- Hogan, R. J., A. Bozzo (2018): A flexible and efficient radiation scheme for the ECMWF model.
  *JAMES*, 10, 1990–2008. (ecRad.)
- [NumericalRadiation.jl documentation](https://NumericalEarth.github.io/NumericalRadiation.jl/dev/):
  cloud optics and radiative transfer pages (`SpectralCloudOptics`, `CloudOverlapShortwave`,
  `CloudOverlapLongwave`).
- [NWS Technical Implementation Notice 14-46](https://www.weather.gov/media/notification/tins/tin14-46gfs_cca.pdf):
  the January 2015 GFS upgrade from Eulerian T574 to semi-Lagrangian T1534.

**SpeedyWeather development plans referenced**

- `docs/dev/2026-09/adjustments-not-tendencies.md` (state adjustments in `filter!`)
- `docs/dev/2026-09/condensation-energy-budget.md` (the enthalpy budget test)
- `docs/dev/2026-09/numericalradiation-extension.md` (the clear-sky ecCKD extension, its work
  arrays and future work)
- `docs/dev/2026-08/gpu-primitive-wet-model-profiling.md` (per-phase GPU cost)

**SpeedyWeather's current physics**

- Frierson, D. M. W., I. M. Held, P. Zurita-Gotor (2006): A gray-radiation aquaplanet moist GCM.
  Part I. *J. Atmos. Sci.*, 63, 2548–2566. (Implicit condensation factor.)
- Frierson, D. M. W. (2007): The dynamics of idealized convection schemes and their effect on the
  zonally averaged tropical circulation. *J. Atmos. Sci.*, 64, 1959–1976. (Simplified Betts–Miller.)
- Kiehl, J. T., et al. (1998): The National Center for Atmospheric Research Community Climate Model:
  CCM3. *J. Climate*, 11, 1131–1149. (Longwave cloud emissivity.)
