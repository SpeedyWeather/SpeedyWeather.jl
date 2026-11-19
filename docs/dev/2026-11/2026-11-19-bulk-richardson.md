# Overhaul of `BulkRichardsonDiffusion` (vertical diffusion)

> Status: **in progress**. Implementation started 2026-11-19 on `mg/bulk-richardson-reimplemented`, steps 1–6
> of section 5 being implemented; step 7 (GPU performance, on a GPU node), step 8 (default on) and the model
> validation of section 6 are open.

Date of initial draft: 2026-11-19

Base revision: `d8ceb78f` (`mg/bulk-richardson-reimplemented`, = main `2c729bee` + this plan). The plan itself was
written from the T512L24 stability investigation on `mg/t512l24-stability` at `81ad7df44`
(= `mg/reactant-output` `6ec9ce35c` + advection option, `STABILITY_NOTES.md` there).
Tracking issue in the repo: #1282 (vertical diffusion), related #1283, #1335/#1336, #1346, #1363.

## Originating prompt

> Start the implementation

## Revision log

- 2026-11-19: sign convention of 4.1 fixed. With F = ρK ∂X/∂z (positive for X increasing upwards, i.e. a
  downward turbulent flux) the tendency is −g ∂F/∂p = g (F_{k−½} − F_{k+½})/Δp_k, not g (F_{k+½} − F_{k−½})/Δp_k.
  The implicit equation of 4.2 already had the diffusive sign and is unchanged.
- 2026-11-19: new step 7 in section 5, GPU performance optimisation, to be done later in a session on a GPU node
  (user request). Initial GPU correctness checks on the local Metal GPU. The former step 7 (default on) is now step 8.
- 2026-11-19, decisions while implementing steps 1–6 (one rewrite of `vertical_diffusion.jl` rather than one commit
  per step, as the explicit operator of step 1 is the Δt = 0 limit of the implicit solve of step 2):
  - Formulas checked against the literature: Holtslag & Boville (1993) read from the paper (K = l²S F(Ri),
    1/l = 1/κz + 1/λ, λ = 300 m below 1 km and λ = 30 + 270 exp(1 − z/1000) above, unstable F = √(1 − 18 Ri),
    bulk Ri with |V|² ≥ 1 m²/s², Ri_cr = 0.5 for h; stable F = 1/(1 + 10Ri(1 + 8Ri)), a short tail). Sandu et al.
    (2013) read from ECMWF Tech. Memo 684: IFS uses λ = 150 m with long tails near the surface, UKMO 40 m, NCEP 30 m
    plus a background diffusivity (heat 1 m²/s decaying with height, momentum 3 m²/s). The Louis, Tiedtke & Geleyn
    (1982) long-tail functions could not be checked against the original (workshop proceedings); the forms of
    section 3 are used as quoted in secondary sources.
  - Separate diffusivities for momentum (f_m) and heat/humidity (f_h); profile and background are shared.
  - Asymptotic mixing length: 150 m (`mixing_length`, IFS/plan) below 1 km, decaying with the HB93 eq. 3.6 shape to
    30 m (`free_mixing_length`) above, instead of a single λ.
  - Background diffusivity K₀ exp(−z/H), K₀ = 0.1 m²/s, H = 1 km instead of a step at 1 km.
  - Surface Ri: land and ocean values of `bulk_richardson_surface` weighted with the land fraction. The function
    was already standalone in `drag.jl` and is called from there, not moved. ln(z/z₀) of the stability correction
    uses the grid-box `momentum_roughness` of the surface roughness scheme (fallback: option `roughness_length`).
  - Friction velocity u⋆ = √C_M Vₛ with the surface wind speed (incl. gustiness) of the boundary layer.
  - `boundary_layer_height` output is now in metres (was a layer index).
  - Source coupling (step 5) is the option `couple_surface_fluxes = true`: the solve is (1 − ΔtL) y = L X⁻ + S
    for the tendency y, replacing the accumulated tendency S.
  - Memory: two `ParameterizationVariable`s (K_m, K_h, `GridXYZ`) and three `ScratchVariable`s (ρ/Δz and the two
    Thomas vectors).
  - Unstable diffusivity: with HB93's λ = 30 m in the free atmosphere, a −30 K superadiabatic layer at 5–6 km
    (T32L8, steps of 2Δt = 4800 s) was mixed halfway to neutral only after ~20 steps (≈ 1 day; 75 m: ~3 steps,
    150 m: 1 step). The unstable K is therefore split into a shear part l²S with the HB93 mixing length and a
    buoyancy part l_c²(√(S² − 18N²) − S) with λ_c = max(λ(z), `convective_mixing_length` = 150 m). Continuous at
    N² = 0; `convective_mixing_length = 0` recovers HB93 eq. 3.7. Superadiabatic layers are now mixed halfway
    within 1–2 steps (test T8) without raising the stable free-atmosphere (jet) mixing. This answers the first
    open question of section 7 in favour of the hybrid K for now; the SPEEDY-style relaxation is not implemented.
- 2026-11-19: long GPU runs for stability added to the model validation of section 6 as the final validation
  (user request), required before step 8 (default on).

## 1. Why this matters

`vertical_diffusion = nothing` is the default of both primitive models since #1346. Two instabilities trace back
to the missing vertical mixing:

- **T512L24/L16 cold pools** (this investigation): the thin lowest layer (Δz ≈ 80 m) decouples in polar-night
  inversions, becomes sharp at the grid scale, and aliasing on the octahedral grid then drives a runaway.
  T512L8 (lowest layer ~500 m) has no cold pools. Vertical mixing removes the precondition.
- **T256L8 Antarctic cold spots** (`docs/dev/2026-10/2026-10-19-t256-instability.md`, finding 3): the
  mid-troposphere becomes strongly superadiabatic (θ_k4 − θ_k5 = −41 K). The model has no dry convective
  adjustment: Betts–Miller needs moisture and the present scheme only acts in the boundary layer.

So the overhaul should deliver (a) boundary-layer mixing that does not switch off completely in stable
conditions and (b) mixing of superadiabatic layers anywhere in the column.

## 2. What is wrong in the current code

File: `SpeedyWeather/src/parameterizations/vertical_diffusion.jl` (323 lines). Docs: `docs/src/vertical_diffusion.md`.

| # | Defect | Where | Evidence |
|---|---|---|---|
| D1 | Neighbour indices are `k₋ = max(k, 1)`, `k₊ = min(k, nlayers)`, so every vertical gradient is zero. The scheme is a no-op | `_vertical_diffusion!`, lines 263–264 | T512L24 run with the scheme is bit-identical to the control (this investigation); already noted in #1282 |
| D2 | Units: K is computed in m²/s (`κ·u·√C·z`) but applied with a Laplacian in σ (`∇²_above/below` ∝ 1/Δσ²). The missing factor is (∂σ/∂z)² = (gσ/(R_d T_v))² ≈ 2·10⁻⁸ m⁻² near the surface | `initialize!` (lines 46–74) and `_vertical_diffusion!` | With D1 fixed, T63L8 gives NaN after 17 h (tested 2026-11-19, change reverted) |
| D3 | Explicit time stepping. With the L24 levels (see table below) and K = 10–50 m²/s the limit Δz²/(2K) is 60–300 s; the leapfrog parameterization step is 2Δt = 300 s at T512 | whole scheme | 24-layer run blew up over Tibet in the #1283 attempt |
| D4 | Surface Richardson number uses the air temperature for the surface value: Θ₀ = cₚT_v(T_N), Θ₁ = Θ₀ + Φ, so Ri_N = Φ²/(Θ₀V²) > 0 always and independent of the actual surface temperature | `bulk_richardson!`, lines 309–313 | Same bug was fixed for the drag in #1363 (`bulk_richardson_surface` in `surface_fluxes/drag.jl`), deliberately not here |
| D5 | Ri uses the absolute geopotential Φ (includes orography) where the height above ground is meant. The K profile below does subtract the orography, the Ri does not | `bulk_richardson!`, lines 313 and 319 | `vars.dynamics.geopotential` includes orography: the drag computes `z = Φ/g − orography` |
| D6 | √C is recomputed from a reference height and a constant roughness length `z₀ = 3.21e-5`, inconsistent with the drag coefficients that the boundary layer scheme already provides per surface type | `get_diffusion_coefficients!`, lines 160–206 | `vars.parameterizations.boundary_layer_drag_momentum` exists |
| D7 | K = 0 above the diagnosed boundary layer and K → 0 as Ri_N → Ri_c: no mixing at all in strongly stable conditions and none for elevated unstable layers | lines 179–228 | T256 finding 3; T512 cold pools |
| D8 | Called *before* the drag and the surface fluxes (`parameterizations` tuple in `models/primitive_*.jl`: `:vertical_diffusion` is second), so it sees the previous step's drag and cannot spread the surface flux of the current step | model structs | – |
| D9 | The docs say Ri_c = 1, the code default is 10 | `docs/src/vertical_diffusion.md`, struct default | #1363 decided Ri_c = 10 for the drag |
| D10 | No tests | `test/parameterizations/` has `drag.jl` only | – |

L24 default levels (`FriersonSigmaCoordinates`), heights for T ≈ 260 K:

| k | σ | Δσ | z (m) | Δz (m) |
|---|---|---|---|---|
| 24 | 0.9946 | 0.0107 | 41 | 82 |
| 23 | 0.9830 | 0.0126 | 130 | 98 |
| 22 | 0.9685 | 0.0164 | 244 | 129 |
| 21 | 0.9493 | 0.0220 | 396 | 176 |
| 20 | 0.9238 | 0.0290 | 603 | 239 |
| 19 | 0.8908 | 0.0371 | 880 | 317 |

## 3. Literature

Equation numbers for Frierson (2006) are the ones already cited in the code and docs. The other references are
quoted from memory of the papers; check the formulas against the originals before coding them.

- **Frierson, Held & Zurita-Gotor (2006), JAS 63, "A gray-radiation aquaplanet moist GCM. Part I"**. The scheme
  the current code follows. Surface drag eqs. 12–14, bulk Ri eq. 15, boundary layer depth h = height where
  Ri(z) = Ri_c with **Ri_c = 1**, K profile eqs. 16–20 (K_b = κ u_a √C z in the surface layer z < f_b h,
  f_b = 0.1, cubic-like decay to zero at h; stable factor `[1 + (Ri/Ri_c) ln(z/z₀)/(1 − Ri/Ri_c)]⁻¹`).
  Dry static energy, humidity and momentum are diffused. Their model (GFDL FMS) solves the diffusion
  **implicitly, together with the surface fluxes**. Already in `docs/src/references.bib` as `Frierson2006`.
- **Troen & Mahrt (1986), BLM 37** and **Holtslag & Boville (1993), J. Climate 6**: the K-profile family that
  Frierson simplifies. h from a bulk Richardson criterion (Ri_c = 0.25–0.5), K ∝ κ w_s z (1 − z/h)². Holtslag &
  Boville combine it with a local Ri-dependent diffusion above the boundary layer, which is the hybrid
  recommended below.
- **Louis (1979), BLM 17; Louis, Tiedtke & Geleyn (1982), ECMWF workshop**: local first-order closure,
  K = l² |∂V/∂z| f(Ri) with Blackadar mixing length l = κz/(1 + κz/λ). The 1982 stability functions have
  "long tails" (f ∝ 1/Ri for large Ri instead of a cutoff): f_m = 1/(1 + 2b Ri/√(1 + d Ri)),
  f_h = 1/(1 + 3b Ri √(1 + d Ri)), b = d = 5; for Ri < 0 f grows, which mixes unstable layers.
- **Derbyshire (1999), BLM 90, "Boundary-layer decoupling over cold surfaces as a physical boundary
  instability"**: with a cutoff stability function, surface cooling reduces the turbulent heat flux, which
  cools the surface further: a positive feedback that decouples surface and atmosphere. This is the
  precondition seen at T512L24.
- **Viterbo, Beljaars, Mahfouf & Teixeira (1999), QJRMS 125**: enhanced stable-regime diffusion (long tails)
  in the ECMWF model removed a winter cold bias and runaway surface cooling over land. **Sandu, Beljaars,
  Bechtold, Mauritsen & Balsamo (2013), JAMES 5**: review of why NWP models keep more stable mixing than
  observations support (synoptic activity, 2 m temperature, numerical robustness).
- **Kalnay & Kanamitsu (1988), MWR 116, "Time schemes for strongly nonlinear damping equations"**: with K
  depending on the state, evaluate K explicitly and the operand implicitly; backward Euler is unconditionally
  stable but can oscillate when K ∝ gradient, an "over-implicit" weight (1.5) damps that. Standard practice
  (ECMWF IFS documentation Part IV, turbulent transport; **Richtmyer & Morton 1967** for the tridiagonal solve).
- **Molteni (2003), Clim. Dyn. 20** (SPEEDY): its vertical diffusion has three parts: shallow-convection-like
  redistribution of dry static energy and moisture between the two lowest layers, moisture diffusion in stable
  conditions, and a **relaxation of superadiabatic lapse rates** anywhere in the column. The last is the cheap
  alternative for requirement (b).

## 4. Proposed design

### 4.1 Formulation in physical units, on pressure thicknesses

Diffusion of X ∈ {u, v, s = cₚT + Φ, q} as a flux divergence:

    ∂X/∂t = (1/ρ) ∂/∂z (ρ K ∂X/∂z) = −g ∂F/∂p,     F = ρ K ∂X/∂z   [X · kg m⁻² s⁻¹]

Discretised with layers k = 1 (top) … N (surface) and interfaces k+½ between k and k+1:

    F_{k+½} = ρ_{k+½} K_{k+½} (X_k − X_{k+1}) / (z_k − z_{k+1}),   z_k − z_{k+1} = (Φ_k − Φ_{k+1})/g
    ρ_{k+½} = p_{k+½} / (R_d T_v,k+½)
    (∂X_k/∂t) = g (F_{k−½} − F_{k+½}) / Δp_k,      F_{½} = 0,   F_{N+½} = 0 (surface flux handled separately, see 4.4)

Reasons: correct units by construction (fixes D2); conserves the mass-weighted column integral Σ X_k Δp_k
exactly; uses `pressure_thickness(k, pₛ, vertical_coordinates)` and `pressure_half`, so it works for
`SigmaCoordinates` and `SigmaPressureCoordinates` (the current σ-Laplacian in `initialize!` does not, see the
TODO there). The precomputed `∇²_above/∇²_below` vectors are dropped. K is defined at interfaces directly
(no averaging of full-level K).

### 4.2 Implicit solve

Backward Euler per column, K evaluated at the old state (Kalnay & Kanamitsu):

    X⁺_k − (Δt g/Δp_k) [ a_{k+½}(X⁺_{k+1} − X⁺_k) − a_{k−½}(X⁺_k − X⁺_{k−1}) ] = X⁻_k + Δt S_k
    a_{k+½} = ρ_{k+½} K_{k+½} g / (Φ_k − Φ_{k+1})

Tridiagonal, diagonally dominant, solved with the Thomas algorithm (two sweeps over k, fixed trip count).
The scheme returns a tendency (X⁺ − X⁻)/Δt − S so it fits the parameterization interface.

- Δt = `default_time_step(model.time_stepping)` (2Δt for leapfrog), as in `surface_flux_damping`
  (`parameterizations/tendencies.jl` line 175). Parameterizations read the previous leapfrog level
  (`which_prognostic_step(…, ::AbstractParameterization) = 1`).
- Dry static energy: solve for s with Φ fixed, then T tendency = (s⁺ − s⁻)/(cₚ Δt).
- Work arrays: the Thomas algorithm needs two per-column vectors in addition to K. Declare them through
  `variables(::BulkRichardsonDiffusion)` as `ScratchVariable(…, GridXYZ(), namespace = :grid)` instead of
  reusing `vars.scratch.grid.a/b` (those are already used for Ri and K, and by the dynamical core).
- GPU/Reactant: the column code runs inside the fused column kernel (`column_parameterizations_kernel!`).
  Keep loops over `1:nlayers` with masks (`ifelse`) as the current `boundary_layers(kₕ, nlayers, ::ReactantDevice)`
  does; no data-dependent loop bounds, no early exits. Enzyme: no allocations, branch-free `ifelse`, avoid
  `0·∞` (see the comment in `neutral_transfer_coefficient`).

### 4.3 Diffusion coefficient

Hybrid as in Holtslag & Boville: K_{k+½} = max(K_profile, K_local, K_min).

1. **Boundary-layer K profile (Frierson eqs. 16–20), corrected**:
   - height above ground z = Φ/g − orography everywhere (D5);
   - surface Ri from `bulk_richardson_surface` (land/ocean, soil temperature and SST), combined with the land
     fraction (D4); factor this function out of `drag.jl` so drag and diffusion share it;
   - √C from `vars.parameterizations.boundary_layer_drag_momentum` (D6), which requires the drag to run first (4.4);
   - h by linear interpolation between the two layers that bracket Ri = Ri_c, not the layer index alone
     (at L8 the index jumps by hundreds of metres).
2. **Local Louis-type term** for everything the profile misses (D7): K_local = l² |ΔV/Δz| f(Ri_loc) with
   Ri_loc = g Δz Δθ_v / (θ_v |ΔV|²) between adjacent layers, l = κz/(1 + κz/λ), λ ≈ 150 m (30 m above the
   boundary layer in Holtslag & Boville), long-tail f for Ri > 0 and growing f for Ri < 0. This gives dry
   convective adjustment of superadiabatic layers (T256 finding) and residual mixing in inversions (T512).
   Cheaper alternative if this proves hard to tune: a SPEEDY-style relaxation of superadiabatic lapse rates
   plus item 3 only.
3. **Background K_min** (0.1–1 m²/s below ~1 km, option, default to be tuned): numerical floor against complete
   decoupling.

Options to expose: `critical_Richardson`, `surface_layer_fraction`, `mixing_length` (λ), `background_diffusivity`,
`diffuse_momentum/static_energy/humidity` (existing), and an upper bound `K_max` (e.g. 1000 m²/s) as a safety net.

Ri_c: Frierson uses 1. #1363 kept 10 for the drag (reasoning in `docs/dev/2026-10/bulk-richardson-surface-temperature.md`).
For the boundary-layer *height* 10 is far outside the literature range (0.25–1) and makes h very deep. Use a separate
parameter for the diffusion with default 1 and fix the docs (D9).

### 4.4 Order and coupling with the surface fluxes

Move `:vertical_diffusion` after the surface fluxes in the `parameterizations` tuple of `PrimitiveDryModel` and
`PrimitiveWetModel` (D8): `… :boundary_layer, :surface_momentum_flux, :surface_heat_flux, :surface_humidity_flux,
:vertical_diffusion, …`.

Then pass the tendencies already accumulated in the column for this step as the source S in 4.2 (at minimum the
surface-flux tendency of layer N). The surface flux is then spread over the boundary layer within the step instead
of being dumped into an 80 m layer. This is the process-split equivalent of Frierson's coupled implicit solve.
It should also make the surface flux damping of #1336 (`RationalFluxes`) less critical, but keep that unchanged
in the first PR. Known risk from #1283: explicit land–atmosphere coupling with thin lowest layers (Tibet).

## 5. Implementation steps

Each step is a separately testable commit; 1–3 are the minimum for a working scheme.

1. **Operator**: new `_vertical_diffusion!` on pressure thicknesses with interface K (4.1), still explicit, with a
   prescribed constant K for testing. Remove `∇²_above/below`. Tests T1–T3.
2. **Implicit solve** (4.2) with scratch variables. Tests T4–T5.
3. **K profile fixes**: z above ground, shared `bulk_richardson_surface`, drag from the boundary layer scheme,
   interpolated h, reorder parameterizations (4.3 item 1, 4.4). Tests T6–T7.
4. **Local term + background** (4.3 items 2–3). Tests T8.
5. **Source coupling** with surface-flux tendencies (4.4). Test T9.
6. **Docs** (`docs/src/vertical_diffusion.md`: operator, implicit solve, hybrid K, Ri_c), CHANGELOG, references
   (`Louis1979`, `HoltslagBoville1993`, `Viterbo1999`, `KalnayKanamitsu1988`).
7. **GPU performance** (to be done in a session on a GPU node, after steps 1–6 work on CPU and Metal): the
   scheme runs per column inside the fused `column_parameterizations_kernel!`. The first implementation is
   written for correctness and keeps three per-column scratch arrays in global memory (ρ/Δz and the two
   Thomas-algorithm vectors) plus the two diffusivities, re-reads u, v, T, q, Φ in several loops, and calls
   `exp`/`log`/`sqrt` per interface. To check and optimise there: cost of `parameterization_tendencies!` with
   and without the scheme at T256/T512 L24 on an H100 (benchmark suite, `--suites` for the PrimitiveWet sweep),
   register pressure and occupancy of the fused kernel, whether the Thomas algorithm should run for all
   variables in one sweep (shared pivots of the momentum and heat systems), and the memory of the five
   additional `GridXYZ` arrays at T512L24 (≈ 130 MB each in Float32). CUDA and Reactant correctness (section 6)
   are tested in the same session.
8. **Default on?** Only after the validation runs in section 6, including the long GPU stability runs; changes
   the climate of every multi-layer run.

## 6. Tests and validation

Unit tests (`SpeedyWeather/test/parameterizations/vertical_diffusion.jl`, add to the test runner):

- T1 conservation: Σ X_k Δp_k unchanged to round-off for random columns and random K, for σ and σ-pressure coordinates.
- T2 well-mixed column (constant X): zero tendency. For s: an adiabatic column (constant cₚT + Φ) gives zero T tendency.
- T3 units: constant K on equally spaced levels reproduces the analytic decay rate K(π/H)² of a cosine mode.
- T4 unconditional stability: K = 10⁴ m²/s, Δt = 3600 s on the L24 levels stays bounded and monotone (no new extrema).
- T5 implicit = explicit in the limit of small Δt.
- T6 surface 20 K warmer / 10 K colder / 20 K colder with weak wind: h and K respond as in `test/parameterizations/drag.jl`.
- T7 orography independence: same column shifted to 4 km surface height gives the same K.
- T8 superadiabatic interior layer is mixed towards neutral within a few steps; strong inversion keeps K ≥ K_min.
- T9 with a surface heat flux into a thin lowest layer, the warming is spread over the boundary layer.
- CPU, CUDA (`using CUDA` first, do not load system cuda modules) and Reactant; Enzyme via the existing
  differentiability tests (`test/differentiability/`).

Model validation:

- T31L8 and T31L24 default wet model, 1 year: stable, global means and precipitation vs `vertical_diffusion = nothing`.
- **Cold-pool restart tests from this investigation**: `run_t512l24_restart_test.jl` with `CONFIG=vdiff`
  (`submit_t512l24_restart_test.sh`), restarts in `output/run_t512l24_gpu_diag_0001/restarts/`. The control
  blows up at Lake Baikal on 02-06 17:55 and reaches −267 °C at Committee Bay on 01-29; compare with
  `python3 compare_restart_tests.py`. About 15 GPU-minutes per configuration on one H100. Also the L16 restarts in
  `output/run_t512l16_gpu_flux_0001/restarts/` (blow-up 01-14 02 UTC).
- T512L24 45 days from 2000-01-01 with `run_t512l24_gpu_diag.jl` (monitor output for bottom-layer extremes and inversions).
- T256L8 Antarctic winter (drivers in `/p/tmp/maxgelbr/code/claude_t256_instability_2026-10-19/`): number of
  superadiabatic columns and cold points < 170 K.
- Over Tibet at L24 (the #1283 failure): check soil temperature and lowest-layer temperature for 2Δt oscillations.
- **Long GPU runs for stability (final validation, user request)**: multi-year integrations of the default wet
  model with `BulkRichardsonDiffusion` on GPU, at several resolutions and layer numbers (L8 to L24; resolutions
  and lengths to be decided, as long as the budget allows). Pass criteria: no blow-up or NaNs, no drift of the
  global means (temperature, precipitation, surface pressure) beyond that of `vertical_diffusion = nothing`, and
  no growing extremes in the lowest layers (cold pools, polar inversions). Run after the GPU performance work
  (step 7) so that the runs use the final kernel.

## 7. Open questions for whoever picks this up

- Hybrid K (4.3 item 2) or the simpler SPEEDY-style lapse-rate relaxation for superadiabatic layers?
- Should the first PR already reorder the parameterizations (changes results of runs that enable the scheme only),
  and when should the scheme become the default?
- Momentum diffusion and the `SpeedLimitDrag`/jet strength interact; check maximum winds (L24 sits at 75–83 m/s,
  the advective limit at Δt = 150 s is ≈ 83 m/s).
- Vertical diffusion alone does not remove the aliasing on the octahedral grid (see `STABILITY_NOTES.md`,
  finding 4). The `dealiasing = 3` runs test that separately.
