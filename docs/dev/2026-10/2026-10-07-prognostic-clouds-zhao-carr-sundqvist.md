# Prognostic clouds for SpeedyWeather: Sundqvist, Zhao–Carr, and how ICON does it

> Status: **planned**. Literature review and recommendation, no code changed yet. The recommendation
> is a Zhao–Carr-type scheme with one prognostic condensate tracer, introduced in stages (section 7).

Date of initial draft: 2026-10-07

Base revision: `dbd4c661` (`mg/version1-version023`)

## Originating prompt

> What would SpeedyWeather need to have prognostic clouds? Specifically reference also how ICON is
> doing that and also research whether they are easier and simpler schemes avaliable

> Summarize this research in a markdown file, put particular focus on the Zhao-Carr and Sundqvist
> scheme and the recommendations

## Revision log

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
- **What already exists:** the tracer framework, the column-physics interface and the precipitation
  sweep of `ImplicitCondensation` cover most of what a one-tracer scheme needs.
- **What is missing:**
  - physics writing tracer tendencies;
  - positivity of a spectrally transported condensate;
  - saturation over ice;
  - clouds in each layer of shortwave *and* longwave radiation.
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
| Tracers | `add!(model, Tracer(:name))` creates a spectral prognostic variable, a grid copy, and spectral and grid tendencies (namespaces `:tracers` and `:grid_tracers`). Tracers are advected (flux form) and hyperdiffused. No parameterization writes a tracer tendency yet. | `dynamics/tracers.jl:9-15, 57-71`, `dynamics/tendencies.jl:805-819`, `dynamics/horizontal_diffusion.jl:271-274` |
| Positivity | `ClipNegatives` does `max(q, 0)`, but only on the **grid copy of humidity**. The spectral state is untouched and mass is not conserved. Tracers get nothing. | `dynamics/hole_filling.jl`, `time_stepping/transform.jl:153-156` |
| Large-scale condensation | `ImplicitCondensation` relaxes `q` to `0.95 q_s` over `3Δt`, with the implicit latent-heat factor `1 + Lᵥ/cₚ ∂q_s/∂T` (Frierson et al. 2006). The condensate rains out immediately. A top-down sweep handles melting (above 278 K), reevaporation (proportional to subsaturation) and snow (below 263 K), and sets `cloud_top`. | `parameterizations/large_scale_condensation.jl` |
| Convection | Betts–Miller (Frierson 2007) with linear entrainment. Rain is the net column drying. All of it falls as snow if the lowest layer is below 273.15 K. No condensate is detrained. `cloud_top` is set to the level of zero buoyancy. | `parameterizations/convection.jl:8-23, 141-179` |
| Clouds | `DiagnosticClouds` (SPEEDY): one cover per column, `min(1, 0.2√P + RH-term²)`. Cloud top is the level of maximum RH or the precipitation cloud top. Stratocumulus comes from static stability. | `parameterizations/radiation/clouds.jl:29-185` (cover at `:149-151`) |
| Radiation | `OneBandShortwave` reflects `cloud_albedo × cover` only at `k == cloud_top` and adds cloud absorptivity below it. The Frierson longwave has **no clouds**. The NumericalRadiation (ecCKD) extension is clear-sky. | `radiation/shortwave_radiation.jl:125-127, 194`, `radiation/shortwave_transmissivity.jl:109-131`, `radiation/longwave_transmissivity.jl:29-78` |
| Saturation | Clausius–Clapeyron over liquid with one `Lᵥ`. `TetensEquation` has a Murray (1967) branch below freezing, but `saturation_humidity` does not use it. | `dynamics/atmosphere.jl:141-160, 162-216` |
| Physics infrastructure | Column kernels `parameterization!(ij, vars, scheme, model)` add into the grid tendencies and are fused into one GPU kernel. The call order is zenith, vertical diffusion, condensation, convection, albedo, radiation, boundary layer, surface fluxes, stochastic physics. Physics reads the *lagged* leapfrog step and runs before dynamics. | `models/primitive_wet.jl:122`, `parameterizations/tendencies.jl:32-88`, `time_stepping/steps.jl:187`, `time_stepping/time_integration.jl:83-90` |
| Time stepping | Leapfrog with `Δt = 40 min` at T32, so physics tendencies act over `2Δt`. | `time_stepping/steppers/leapfrog.jl:186` |
| Post-step adjustments | `filter!` exists only for the grid-point components ocean, sea ice and land; there is none for spectral atmospheric variables. | `time_stepping/time_integration.jl:140-142` |

### 2. What a prognostic cloud scheme needs

| # | Need | Status |
|---|---|---|
| 1 | A condensate variable `q_c` (cloud water + ice), advected | Tracer framework exists. Physics must write `vars.tendencies.grid_tracers[:q_c]`, which is then transformed together with the advection tendency at no extra cost. |
| 2 | Sources and sinks: condensation, cloud evaporation, autoconversion, accretion, phase | Missing. `ImplicitCondensation`'s relaxation can be reused as the condensation source. |
| 3 | Precipitation with evaporation and melting | The top-down sweep exists in `ImplicitCondensation`. |
| 4 | A cloud fraction consistent with `q_c` | Missing. The current cover is RH- and precipitation-based, one per column. |
| 5 | Radiation using cloud fraction, condensate path and effective radius in each layer, in SW **and** LW, with an overlap assumption | Missing, and the largest single work item. |
| 6 | Convective condensate (anvils) | Missing. Betts–Miller has no updraft to detrain from. |
| 7 | Positivity and conservation of a spectrally transported, intermittent field | Missing for tracers. |
| 8 | Ice thermodynamics: `q_s` over ice or mixed phase, `L_s`, fusion heat | Partly present (Tetens–Murray). |
| 9 | Bounded fast sinks under leapfrog (`2Δt = 80 min` at T32) | A design rule: implicit or exponential forms. |
| 10 | Output: liquid and ice water path, 3D cloud fraction | Missing. |

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

  - `((1 − f) q_s)^¼` is clamped to `[10⁻⁴, 1]` and the exponent to 50.
  - `C < 0.001` is set to 0.
  - Because `C → 0` as `q_c → 0`, there is no cloud without condensate.
- **Condensate path** `q_c Δp/g`:
  - The ice fraction is a linear ramp, `clamp((T₀ − T)/20 K, 0, 1)`, independent of the
    microphysics `IW`.
  - Optionally the path is normalised to in-cloud values (divided by `C`) and smoothed vertically
    1-2-1.
- **Effective radius:**
  - liquid: 10 µm over ocean, `5 + 5 f_ice` µm over land;
  - ice: a temperature-dependent power law of ice water content, limited to 10–150 µm.

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

- **One tracer, `q_c`.** It costs one spectral→grid transform (`transform.jl:181-185`) and three
  grid→spectral transforms (`tendencies_sequential.jl:151-175`) per step, the same as any tracer.
  Physics tendencies are transformed together with the advection tendency.
- **Diagnostic precipitation.** It needs no sedimentation stability (CFL) treatment, which matters
  because rain falls about 12 km in one 40 min step.
- **Column-local, smooth formulas**, which suit the fused GPU kernel and Enzyme/Reactant.
- **Proven numerics.** GFS ran exactly this kind of scheme in an Eulerian spectral leapfrog model.
- **Reuse.** It reuses the existing condensation and precipitation code.

### 7.2 Stage 0: prerequisites

1. **Condensate variable.** Either a `Tracer(:cloud_condensate)` that the parameterization owns, or a
   dedicated prognostic variable like humidity, fused into the main batched transform. *Decision
   needed.* The tracer is less work; the dedicated variable may be faster.
2. **Positivity.**
   - Physics reads `max(q_c, 0)`.
   - Add GFS's **fill from vapour** as a tendency: a negative `q_c` is moved to 0 from `q` with
     latent heating, which conserves total water and enthalpy.
   - Optionally use the upwind or WENO vertical advection (`dynamics/vertical_advection.jl:8-11`).
     It is a model-wide choice and reduces overshoots at cloud edges.
3. **Ice thermodynamics.**
   - A mixed-phase `q_s`, e.g. the existing Tetens–Murray branch, or a blend of liquid and ice over
     0 to −20 °C.
   - `L = Lᵥ + f_ice L_f`.
   - Keep the column enthalpy budget closed, as `test/parameterizations/large_scale_condensation.jl`
     tests today.
4. **Radiation with layer clouds.** This is the largest item and can be developed against Stage 1
   output.
   - **Shortwave:** reflection and absorption in every layer from the optical depth, instead of only
     at `cloud_top`. For liquid, `τ = 3 LWP/(2 ρ_w r_e)`.
   - **Longwave:** cloud emissivity `ε = 1 − exp(−D κ LWP)`, with `D ≈ 1.66` (CCM3 style), multiplied
     into the Frierson transmissivity.
   - **Overlap:** random or maximum-random, weighting each layer by its cloud fraction.
   - The NumericalRadiation extension would need an all-sky path.
5. **Output:** liquid and ice water path, 3D cloud fraction, and the process rates for debugging.

### 7.3 Stage 1: Zhao–Carr with relaxation condensation (recommended first implementation)

This is one column kernel that replaces or extends `ImplicitCondensation`, working from the top down:

1. **Phase.** A smooth ice fraction `f_ice(T)`, either GFS radiation's ramp (0 to −20 °C) or
   Zhao–Carr's thresholds at 0 and −15 °C written as a ramp. Avoid the branchy `IW` seeding logic
   for GPU and AD.
2. **Condensation and cloud evaporation.** Reuse the implicit relaxation of `ImplicitCondensation`:
   `δq = (q − RH_c q_s) / ((1 + L/cₚ RH_c ∂q_s/∂T) τ)`.
   - Supersaturated (`δq > 0`): condense into `q_c` instead of raining out.
   - Subsaturated with `q_c > 0`: evaporate cloud, limited by `q_c`.
   - This is a saturation adjustment towards `RH_c` with timescale `τ`. It needs no forcing `M`, is
     already energy-checked and is stable under leapfrog.
3. **Cloud fraction:** Xu–Randall in GFS form, `C(q_c, f, q_s)` (section 4.5). Use it for both
   autoconversion and radiation. Do *not* use Sundqvist `b(f)` here: with the relaxation, RH stays
   near `RH_c` in cloudy air, so `b(f)` would be degenerate.
4. **Conversion to precipitation:**
   - liquid: Sundqvist with `F_co` from the flux above and `F_BF` (section 3.3);
   - ice: `P_saut` from section 4.4.
   - Apply both in **exponential form**, `Δq_c = q_c (1 − exp(−k · 2Δt))` with `k = P/q_c`, so `q_c`
     stays non-negative at `2Δt = 80 min` even when `F` is 5 to 10.
5. **Precipitation sweep:** keep the existing melting, reevaporation and snow code. Optionally add
   collection of cloud by falling precipitation (`P_racw`, `P_saci`).
6. **Output tendencies** for `q`, `T` and `q_c`, plus rain and snow rates. Keep `cloud_top` until the
   radiation no longer needs it.

**Parameters**, starting from GFS or the current SpeedyWeather values:

| Parameter | Start value |
|---|---|
| `RH_c` | 0.95 |
| `τ` | `3Δt` |
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

1. **The GFS approach (recommended).**
   - Store `T` and `q` (and `p_s`) as they are after each cloud-scheme call, keeping two copies for
     leapfrog parity.
   - Use `A_X = (X_lagged − X_ref)/(2Δt)`.
   - Cost: about four extra 3D grid fields. This is proven under leapfrog in GFS.
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

### 7.6 Numerics checklist

- **Bounded sinks:** every term that can remove `q_c` within a step is implicit or exponential, or is
  capped at `q_c/(2Δt)`. There is no atmospheric `filter!` to fix overshoots afterwards.
- **Conservation tests:** column total water (`q + q_c` + precipitation) and enthalpy, extending the
  existing condensation budget tests.
- **Differentiability:**
  - clamp `f ≤ 1 − ε` in `b(f)`;
  - clamp the Xu–Randall denominator (GFS uses `≥ 1e-4`);
  - prefer the smooth Sundqvist threshold to hard thresholds;
  - use a phase ramp instead of the `IW` branches.
- **GPU:** one fused column kernel, fixed loop lengths, no data-dependent early exits beyond what
  `ImplicitCondensation` already has.
- **Spectral ringing:** `q_c` is intermittent. Hyperdiffusion damps the smallest scales. GFS shows
  that clipping on read plus fill-from-vapour is enough for an operational spectral model.

### 7.7 What not to do first, and why

| Approach | Why not now |
|---|---|
| Prognostic rain and snow (ICON, CloudMicrophysics 1-moment) | Sedimentation at about 5 m s⁻¹ with 40 min steps and 8 layers needs implicit or semi-Lagrangian fall. Little gain at T31. |
| Two-moment schemes | Need aerosol, which SpeedyWeather does not have. |
| Prognostic cloud fraction (Tiedtke, PC2) | Stiff, bounded variable with many source terms. |
| Instantaneous saturation adjustment (ICON `satad`) | Incompatible with a spectral state updated by leapfrog tendencies. The relaxation with `τ ≥ Δt` is SpeedyWeather's equivalent. |

## Testing and verification (planned)

- **Unit tests on prescribed columns:**
  - a supersaturated layer produces `q_c`;
  - `q_c` decays by autoconversion with the expected rate;
  - cloud evaporates in subsaturated air;
  - the water and enthalpy budgets close;
  - `q_c ≥ 0` after one step for extreme rates.
- **GPU and differentiability tests** like those of the existing parameterizations.
- **Long runs at T31 and T63.** Compare global means against observations: total cloud fraction of
  about 0.6–0.7, and CERES-EBAF cloud radiative effects of about −45 to −47 W m⁻² (SW),
  +26 to +28 W m⁻² (LW) and about −20 W m⁻² (net). Check that precipitation stays close to the
  current scheme. Also check plausible liquid and ice water paths, and the zonal-mean `q_c`
  structure (storm tracks, tropical upper troposphere).

## Documentation changes (planned)

A new docs page for the cloud scheme, following `docs/src/large_scale_condensation.md`. Update the
radiation docs once clouds enter the longwave.

## Known limitations of this review

- **Primary papers not read:** the Sundqvist et al. (1989) and Zhao & Carr (1997) PDFs were not
  accessible (AMS returned 403). Formulas and constants come from the GFS code and its embedded
  documentation, which follow the papers but carry GFS tuning (e.g. `c₁`, ice autoconversion).
- **Unverified constants:** the original Xu & Randall (1996) constants were not checked; GFS's are.
  The ICON-A default critical-RH values were not found, only the parameter names and the refits by
  Grundner et al. (2022).
- **No experiments:** nothing has been run. The parameter values are starting points.

## Future work

- Two tracers (`q_l`, `q_i`) with sedimentation of cloud ice.
- A convective cloud scheme tied to Betts–Miller.
- An all-sky NumericalRadiation path.
- Prognostic cloud fraction, only if Stage 2 shows systematic cover errors that RH and condensate
  cannot fix.

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

**SpeedyWeather's current physics**

- Frierson, D. M. W., I. M. Held, P. Zurita-Gotor (2006): A gray-radiation aquaplanet moist GCM.
  Part I. *J. Atmos. Sci.*, 63, 2548–2566. (Implicit condensation factor.)
- Frierson, D. M. W. (2007): The dynamics of idealized convection schemes and their effect on the
  zonally averaged tropical circulation. *J. Atmos. Sci.*, 64, 1959–1976. (Simplified Betts–Miller.)
- Kiehl, J. T., et al. (1998): The National Center for Atmospheric Research Community Climate Model:
  CCM3. *J. Climate*, 11, 1131–1149. (Longwave cloud emissivity.)
