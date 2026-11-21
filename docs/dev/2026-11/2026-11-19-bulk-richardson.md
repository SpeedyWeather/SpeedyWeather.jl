# Overhaul of `BulkRichardsonDiffusion` (vertical diffusion)

> Status: **in progress**. Steps 1–8 are implemented: the scheme works on CPU, CUDA and Reactant, costs 1.4–2.8 %
> of a time step on the H100 and about 21 % on the CPU (8.17), and is the default of both primitive models
> (committed in `b7bdfdd1` with the Reactant workaround).
> Tests: see 8.11. Validation: T32, T256L8 and the T512L24 cold-pool restart tests pass, but T32L24 with the scheme
> blows up in the lowest layers at the Antarctic coast in 2 of 4 years with main's drag (0 of 4 with the new damped
> drag, 8.13.2). T512L24 needs
> a stronger speed-limit drag against the top-layer jet failure (which it also has without the scheme), and in the
> multi-year runs a rare cold event at one Antarctic coastal point remains, which a shorter hyperdiffusion time
> scale removes in restart tests. The final verification, T512L24 over three years, is running in four setups
> (8.8, 8.10). Since `1cdc3b1a` the boundary layer height uses Ri_c = 10, which improves the tropics but has not
> been through these stability checks (8.16). A performance audit on 2026-11-21 (8.17) found bounds checks
> left active on the CPU and in the scheme (they are also what made the Reactant workaround of 8.3 necessary) and
> a form of the scheme with the layers in the outer loop that is 2 to 15 times faster on the CPU and 40–45 %
> faster on the GPU; none of it is implemented yet. Details in section 8.

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

- 2026-11-19, user requests during the GPU session: (1) the first goal is a fast, efficient and
  Reactant-compatible implementation of the scheme; (2) the goal remains a stable T512L24, to be reached by
  parameter choices of this and other parameterizations only (structural changes elsewhere are out of scope),
  trying setups with a shorter time step and a stronger speed-limit drag against the top-layer jet failure;
  (3) the last verification is a T512L24 run over several years. Section 6 is extended accordingly.
- 2026-11-20, user request: the scheme becomes the default of both primitive models, and all tests have to pass.
  Done in the working tree (8.11). Two things followed from it: the scheme no longer throws without a boundary
  layer (it then has no boundary-layer profile), and `reset_variable!` takes the variable name as a `Val`, because
  with the scheme's three additional variables the lookup by `Symbol` was no longer inferred (runtime dispatch).
- 2026-11-20, user requests: audit whether all variables of the scheme are necessary; use the model's general
  scratch arrays for the solver; store the conductances ρK/Δz instead of the diffusivities and ρ/Δz; keep two
  diffusivities (one K for momentum and heat is not wanted). Done and checked on Reactant (8.15).
- 2026-11-20, user request: a draft of the PR description, comparing the implementation with main and with the
  draft of #1283 (`milankl-claude:mc/richardson-surface-temperature`), with the GPU measurements already made (no
  new benchmarks). It is in `PR_bulk_richardson_diffusion.md` in the repository root, untracked.
- 2026-11-20, user decisions on the two parameters that T512L24 needs (8.8, 8.10): (1) the stronger drag becomes
  the default; (2) the shorter hyperdiffusion time scale is **not** changed and stays a to-do: check whether it
  is really needed (8.12).
- 2026-11-20, how the stronger drag became the default (8.13): a fixed `drag = 1e-5` (committed in `89be0251`) made
  T32L24 with the scheme blow up within a year and one test fail. I first added a limiter to the drag (at most
  excess/Δt); on the user's question whether the drag could instead just be larger at high resolution, that was
  replaced by a resolution-scaled coefficient, `drag * (truncation/32)^resolution_scaling` with `drag = 3e-6` and
  `resolution_scaling = 0.5`. The limiter is not in the code.
- 2026-11-20, user request: replace the resolution scaling of `SpeedLimitDrag` by a damped time integration as for
  the surface fluxes in #1336 (8.13.1). The drag is multiplied by `flux_damping(time_integration, 2Δt·c·excess)`,
  default `ImplicitFluxes()` (exact for the quadratic drag), with `drag = 1e-5` at all resolutions (`41854301`).
  I also moved it to the previous leapfrog step, the forward step over 2Δt that φ is derived for. With that,
  T32L24 with the scheme blew up in the lowest layers at the Antarctic coast in the first validation runs and in
  every restart variant (2 of 4 years in the ensemble); on the user's request it reads the current step again
  (`b7156e83`) and completes 4 of 4 years there (8.13.2). The same runs show that T32L24 with the scheme also
  blows up with main's drag (2 of 4 years).
- 2026-11-20, user decision: the Richardson number at which the stable surface correction reaches zero is a
  cutoff, not a critical Richardson number in the sense of a turbulence threshold (0.25–1). With a cutoff of 1
  the drag vanishes where buoyancy and shear are still comparable; 10 means that buoyancy dominates by an order
  of magnitude. Hence `BulkRichardsonDrag.critical_Richardson` is renamed to `Richardson_cutoff` (default 10),
  and `BulkRichardsonDiffusion` gets a separate `Richardson_cutoff = 10` for the stability correction of the
  profile, while `critical_Richardson = 1` keeps defining the boundary layer height. Before the split the
  diffusion used Ri_c = 1 for both. The validation runs of section 8 used that, so the profile is now somewhat
  stronger in stable surface conditions than in those runs.
- 2026-11-20, user decision: `critical_Richardson = 10` also for the boundary layer height (committed in
  `1cdc3b1a`, together with the new output `OutputVariables.BoundaryLayerHeight()`, `blh` in m). In the user's
  simulations this clearly improves the tropics: tropical cyclones at T86, patchier convective rain instead of
  large areas of weak rain, and much more small-scale dynamics. The tropical boundary layer is 1–2 km deep with
  Ri_c = 10 against below 1 km with 1, closer to observations (section 8.16). The cost is in the mid-latitudes,
  where h reaches the tropopause (> 10 km) under the jet, so the profile diffusivity extends through the
  troposphere there. Ri_c = 10 is therefore an empirical choice for this model, not a critical Richardson number
  in the sense of the literature (0.25–1), against the recommendation of 4.3.
- 2026-11-20, tried and reverted: surface excess terms in the bulk Richardson number of the boundary layer height
  (Troen and Mahrt 1986, Holtslag and Boville 1993, Vogelezang and Holtslag 1996), meant to deepen the tropical
  boundary layer without the mid-latitude overshoot so that Ri_c could return to 0.25–1. See section 8.16.
- 2026-11-21, user decision: the time integration schemes are renamed from fluxes to relaxation, as they are
  used for any relaxation dx/dt = −(x − x₀)/τ (now also the drag), not only for fluxes: `AbstractFluxIntegration`
  → `AbstractRelaxationIntegration`, `ExplicitFluxes`, `ImplicitFluxes`, `ExponentialFluxes`, `RationalFluxes` →
  `ExplicitRelaxation`, `ImplicitRelaxation`, `ExponentialRelaxation`, `RationalRelaxation` (old names removed,
  no deprecations, although exported since v0.24), `flux_damping` → `relaxation`, and
  `time_stepping/flux_integration.jl` → `time_stepping/relaxation_integration.jl`. `surface_flux_damping` and
  `DEFAULT_FLUX_INTEGRATION` (the default of the surface fluxes) keep their names. Older names in this plan
  refer to the code before the rename.
- 2026-11-21, user request: audit the performance of the scheme for further optimisations on CPU and GPU, and
  protocol the progress here. No source change; prototypes and measurements outside the repository. See 8.17.
- 2026-11-19 (GPU session, A40 on login02 and H100 via Slurm), step 7, see section 8 for the numbers:
  - No GPU optimisation is kept. The scheme doubles the cost of the fused column kernel, but that kernel is only
    about 4 % of a time step, so the scheme costs 2–5 % of a time step at T256/T512.
  - An in-place tridiagonal solve (right-hand side stored in the tendency, one scratch array fewer) was implemented
    and measured: 14–17 % less added kernel cost. It was reverted because the Reactant compilation of the column
    kernel failed with it. The patch is kept outside the repository (section 8.5).
  - Solving u and v (and T and q) in one sweep with shared pivots was not implemented: estimated gain below 1 % of
    a time step.
  - The model validation runs use their own environment and copies of the T512L24 driver scripts in
    `/p/tmp/maxgelbr/code/Speedy2/claude_bulk_richardson_2026-11-19/validation/`, so that nothing in
    `/p/tmp/maxgelbr/t512l24-reactant-10y/` (another session's work) is modified. Its restart files are read only.

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
  (step 7) so that the runs use the final kernel. **The last verification is a T512L24 run over several years**
  (user request 2026-11-19), with whatever time step and parameters of the other parameterizations make T512L24
  stable (8.8); structural changes outside the scheme are out of scope.

## 7. Open questions for whoever picks this up

- Hybrid K (4.3 item 2) or the simpler SPEEDY-style lapse-rate relaxation for superadiabatic layers?
- Should the first PR already reorder the parameterizations (changes results of runs that enable the scheme only),
  and when should the scheme become the default?
- Momentum diffusion and the `SpeedLimitDrag`/jet strength interact; check maximum winds (L24 sits at 75–83 m/s,
  the advective limit at Δt = 150 s is ≈ 83 m/s).
- Vertical diffusion alone does not remove the aliasing on the octahedral grid (see `STABILITY_NOTES.md`,
  finding 4). The `dealiasing = 3` runs test that separately.

## 8. Progress protocol

State on 2026-11-20, 11:25 CET (sections 8.10–8.14 carry their own times). Everything below was run on the committed code `1179b149` unless stated otherwise.
Scripts and logs are in `/p/tmp/maxgelbr/code/Speedy2/claude_bulk_richardson_2026-11-19/` (called `$D` below),
not in the repository.

### 8.1 Where the plan stands

| Step | State |
|---|---|
| 1–6 operator, implicit solve, K profile, local term, source coupling, docs | done, committed in `1179b149` |
| 7 GPU correctness, CUDA | done (8.2) |
| 7 GPU correctness, Reactant | dry and wet model pass with the workaround in the working tree; on the GPU backend the scheme works and changes nothing about the 3 comparisons that fail there without it (8.3) |
| 7 GPU performance | measured, three optimisations tried, none kept (8.4, 8.5) |
| 6 validation: T32L8 and T32L24, one year | done (8.6) |
| 6 validation: T512L24 cold-pool restart tests | done (8.7) |
| 6 validation: T512L24 45 days | no cold pools; NaN on 02-08 from the top-layer jet (8.8) |
| 6 validation: T256L8, 300 days through the Antarctic winter | done (8.9) |
| 6 validation: T512L24 with shorter time step / stronger speed-limit drag | all 4 runs completed 60 days (8.8) |
| 6 final verification: T512L24 over several years | four runs of 3 years running, the first three at 1.4–1.9 years without NaN at 09:20; rare cold event at one Antarctic point (8.10) |
| 6 Enzyme differentiability tests | pass with the scheme as default (8.11) |
| 8 default on | done in the working tree (8.11) |

The source and test changes of these sessions are committed in `b7bdfdd1` ("make it default", 2026-11-20 08:06):
the Reactant workaround of 8.3, no error without a boundary layer, the defaults of both primitive models,
`reset_variable!` in `parameterizations/tendencies.jl`, the tests of the scheme and one GPU test (8.11). Not
committed yet: `CHANGELOG.md`, `docs/src/vertical_diffusion.md` (both say that the scheme is the default) and
this document.

### 8.2 Unit tests and CUDA

- CPU unit tests (`test/parameterizations/vertical_diffusion.jl`, `--check-bounds=yes`): 67 of 67 pass. The file
  takes 7.5 minutes, almost all of it compilation of four model configurations.
- CUDA on the A40 (`$D/cuda_correctness.jl`, log next to it; repeated with the workaround of 8.3 for 2 and 20
  steps, same numbers): PrimitiveWet and PrimitiveDry at T32 with 8 and 24
  layers, CPU and GPU from the same initial conditions, compared after 2, 20 and 72 steps. No NaN in any of the
  12 runs. After 2 steps the diffusivities and the boundary layer height agree to 1–4·10⁻⁴ of their maximum and
  temperature to 2–5·10⁻⁶. Later the pointwise differences grow (up to 0.9 of the maximum diffusivity after 72
  steps) while the global means stay within 1 %, i.e. the two runs diverge chaotically, as they do without the
  scheme. This is a check of the whole model run, not of the kernel in isolation.

### 8.3 Reactant

`$D/reactant_correctness.jl` runs the tendency and time-stepping comparison of `test/reactant/test_correctness.jl`
(T32L8, CPU backend of Reactant, 100 steps) with `vertical_diffusion = BulkRichardsonDiffusion{Float32}()`.

- **PrimitiveDryModel: 20 of 20 pass.** Mean diffusivities agree (momentum 4.4263 on CPU, 4.4267 on Reactant).
- **PrimitiveWetModel: 10 of 23 fail** (tendencies and state outside the tolerances), and the Reactant model has
  `vertical_diffusivity_momentum`, `vertical_diffusivity_heat` and `boundary_layer_height` exactly zero
  everywhere, so the scheme writes nothing there. The cause is not known. What is established, from 5 time steps
  at T32L8 (`$D/reactant_probe.jl`, `$D/reactant_variant.jl`, CPU reference `$D/reactant_probe_cpu.jl`):
  - Without the scheme the wet model on Reactant agrees with the CPU (mean drag coefficient 5.3997·10⁻⁴ on both,
    and switching the large-scale condensation off changes the cloud top on both). So the scheme is the cause.
  - With the scheme, the Reactant result is the same to all printed digits whether the humidity diffusion, the
    large-scale condensation or the default instead of the test physics is used (mean drag 7.13·10⁻⁴, CPU
    5.39·10⁻⁴). Parameterizations other than the scheme therefore also seem to have no effect, i.e. the compiled
    column kernel is wrong as a whole, not only the scheme's part of it.
  - On the CPU the wet model has a mean diffusivity of momentum of 0.22 m²/s and a mean boundary layer height of
    222 m after these 5 steps, so zero is not a valid result.
  - Reproducer (`$D/reactant_clean.jl`): a kernel that only calls `vertical_diffusivity!` for every column of the
    wet model, compiled twice on the same state. With Reactant's default options (kernels not raised to StableHLO)
    it gives the expected diffusivities (mean of momentum 0.203 m²/s, boundary layer height 196 m). With
    `reactant_compile_options()` (kernels raised, as in every model run) the kernel writes nothing: both
    diffusivities, the boundary layer height and the scratch array stay zero. So the fault is in the raising of
    this kernel or a pass after it, not in the tracing, and it needs only the scheme's own function to show.
  - The whole `time_step!(vars, time_stepping, model)` compiled without raising also gives correct
    diffusivities (logs in `$D/default_options/`).
  - Variants of the function on a valid state (`$D/reactant_clean_[ab].log`): the raised kernel is correct if
    humidity is not read at all, if the virtual temperature is computed from humidity only in the loop for the
    boundary layer height, or only in the loop for the diffusivities. It writes nothing if both loops compute
    it. It also writes nothing if the first loop computes and stores it and the second reads the stored values
    (`$D/reactant_diffusion_check_folded.log`). Removing the surface Richardson number or the roughness length
    changes nothing. Simple kernels that read humidity at layers k and k+1 are raised correctly
    (`$D/reactant_minimal.jl`). What exactly breaks the raising is not known; it has not been reduced to a small
    example for an upstream issue.
  - An earlier round of variants (logs in `$D/bisect1/`, `reactant_bisect_[c-f].log`) is not usable: its state had
    zero surface wind (no boundary layer step had run), for which the surface Richardson number is NaN on any
    backend.
- **Workaround, in the working tree, not committed**: `vertical_diffusivity!` first writes the virtual temperature
  of all layers into the scratch array `tridiagonal_rhs` (free until the solves) in a loop of its own, and the
  two other loops read that instead of temperature and humidity. With it:
  - the reproducer and the whole `vertical_diffusion!` as one kernel (`$D/reactant_diffusion_check.jl`, wet and
    dry) give the same result raised and not raised (relative differences below 10⁻⁴);
  - the full comparison above passes for both models: PrimitiveDryModel 20 of 20, PrimitiveWetModel 23 of 23
    (`$D/reactant_correctness_workaround_only.log`; mean diffusivity of momentum 2.698 m²/s on the CPU and
    2.696 m²/s on Reactant after 100 steps), and again on the final text of the file
    (`$D/reactant_correctness_final_cpu.log`, same counts);
  - CPU unit tests: 67 of 67 pass.
- Reactant on the GPU backend (A40, `reactant_correctness.jl gpu …`, `$D/reactant_correctness_gpu.log` and
  `reactant_correctness_gpu_dry_vdiff.log`): the scheme runs and has non-zero diffusivities close to the CPU (wet
  model after 100 steps: mean of momentum 3.03 against 2.96 m²/s, boundary layer height 804 against 798 m). In
  every configuration 3 comparisons are outside the tolerances, with and without the scheme:

  | Model | Without the scheme | With the scheme |
  |---|---|---|
  | PrimitiveDryModel | 14 pass, 3 fail | 17 pass, 3 fail |
  | PrimitiveWetModel | 17 pass, 3 fail | 20 pass, 3 fail |

  (The three additional passes are the checks that the diffusivities are finite.) The differences between CPU and
  Reactant are of the same size with and without the scheme (wet model, largest difference of the divergence
  tendency 5.1·10⁻³ and 5.4·10⁻³, of temperature after 100 steps 1.5 K and 2.8 K). So the scheme does not make
  the agreement on the GPU backend worse. Why three comparisons fail there at all is a separate question that
  was not looked into; the tolerances are those of the CPU backend, and matrix products on the GPU run in
  reduced precision (TF32) by default.

### 8.4 GPU performance

H100, CUDA.jl, Float32, default `PrimitiveWetModel`, `$D/gpu_perf.jl` and `$D/gpu_step_profile.jl`
(Slurm jobs 3342330 and 3342587). Column kernel: minimum of 30 calls of `column_parameterizations!`.

| | Column kernel without [ms] | with scheme [ms] | added [ms] | Time step without [ms] | with [ms] |
|---|---|---|---|---|---|
| T256 L24 | 0.65 | 1.26 | 0.62 | 17.8 | 18.3 |
| T512 L24 | 2.21 | 4.45 | 2.25 | 72.6 – 74.7 | 76.4 – 76.6 |
| T512 L8 | 0.97 | 1.77 | 0.80 | 38.9 | 39.7 |

- The scheme about doubles the column kernel. A time step is dominated by the dynamics and the transforms
  (T256 L24: dynamics 7.6 ms, transform 8.9 ms, parameterizations 0.8 ms without and 1.8 ms with the scheme when
  called inside a time step), so the scheme costs about 3 % of a time step at T256 L24, 2 % at T512 L8 and
  2–5 % at T512 L24.
- The time step numbers are noisier than the kernel numbers. The T512 L24 range is from two jobs; the first
  model built in a Julia session ran up to 8 % faster than later ones (T256 L24: 16.5 ms as first, 17.8 ms when
  alternating), so only runs in alternating order should be compared.
- Parts of the added cost at T256 L24, from runs with only momentum or only dry static energy diffused: the
  diffusivities about 0.16 ms and each tridiagonal solve about 0.15 ms. `couple_surface_fluxes` makes no
  difference.
- The added cost scales with the number of grid points times layers. My reading is that it is limited by memory
  traffic (I count about 47 reads or writes of a 3D array per time step for the scheme), not by arithmetic. The
  register use and occupancy of the kernel were not measured.
- Memory of the five additional 3D arrays: not measured (the script's numbers accumulate over the models built in
  one session). By size it is 5 × 58 MB at T512 L24 in Float32, about 3 % of the 9.7 GB the model uses there.

### 8.5 Optimisations tried or considered

Column kernel on the H100 at T256 L24 (0.65 ms without the scheme), same script as 8.4. The same code measured
in two jobs differs by up to 0.03 ms (the workaround: 1.27 and 1.30 ms), so differences of that size are noise.

| Variant | Kernel with scheme [ms] | Reactant |
|---|---|---|
| committed (`1179b149`) | 1.26 | dry model correct, wet model writes nothing (8.3) |
| + workaround of 8.3 (virtual temperature in a loop of its own), **the code in the working tree** | 1.27, 1.30 | correct |
| in-place solve, first version (right-hand side in the tendency for both couplings, layer values carried through the loop) | 1.16 | does not compile ("cannot raise op to stablehlo") |
| workaround + right-hand side in the tendency for the coupled solve only, nothing carried | 1.26 | correct |
| workaround + layer values carried through the loops of solver and diffusivities, each layer read once | 1.30 | does not compile |
| virtual temperature computed and stored inside the boundary layer height loop | 1.30 | wet model writes nothing |

- None of the variants that Reactant compiles correctly is faster than the committed code, so none is kept. The
  only faster one (8 %, i.e. 0.5 % of a time step) is the first in-place solve, which Reactant cannot raise. Why
  it is faster is not understood: neither of its two ingredients alone (in-place, carried values) gives the gain.
- Carried values (`x = x_below` at the end of the forward sweep, and the same for temperature, geopotential and
  wind in the diffusivity loop) are what Reactant cannot raise. The other carried values of these loops (pivots,
  conductance, height and Richardson number of the layer below) are raised without problems.
- The workaround of 8.3 costs little: within the noise at T256 L24, and 0.12–0.16 ms at T512 L24 (kernel 4.57 and
  4.61 ms in two jobs against 4.45 ms), i.e. 0.2 % of a time step.
- **Not tried: one sweep for u and v, and for T and q, with shared pivots.** It would save the second read of
  diffusivity and density and the second upper diagonal, about a tenth of the scheme's memory accesses by my
  count, for one more scratch array. Expected gain below 0.5 % of a time step.
- Conclusion: the scheme stays at about 0.6 ms at T256 L24 and 2.4 ms at T512 L24 on the H100, 2–5 % of a time
  step. Making it substantially cheaper would need a different structure (all four variables and the
  diffusivities in one sweep), which is larger than what was tried here and would have to be checked against
  Reactant's raising at every step. Patches of the first variant: `$D/inplace_solve.patch`.

### 8.6 One-year runs at T32

`$D/validation/run_means.jl`, H100, default `PrimitiveWetModel` from 2000-01-01 for 360 days, statistics every
10 days (`$D/validation/means/*.csv`, job 3342771). All four runs finish without NaN. Means over days 190–360,
extremes over the same period:

| | L8 without | L8 with | L24 without | L24 with |
|---|---|---|---|---|
| Lowest-layer temperature, global mean [K] | 283.0 | 283.3 | 288.9 | 285.5 |
| Lowest-layer temperature, maximum [K] | 308.9 | 303.4 | 325.9 | 309.5 |
| Lowest-layer temperature, minimum [K] | 215.5 | 213.0 | 207.5 | 210.1 |
| Column-mean temperature [K] | 252.4 | 252.0 | 252.1 | 252.3 |
| Lowest-layer humidity [g/kg] | 7.9 | 6.9 | 11.9 | 7.8 |
| Soil temperature over land, mean [K] | 279.6 | 277.8 | 286.1 | 279.0 |
| Precipitation [mm/day] | 1.95 | 2.06 | 1.70 | 2.17 |
| Surface pressure [hPa] | 974.39 | 974.38 | 974.38 | 974.31 |
| Maximum wind speed [m/s] | 96 | 101 | 96 | 110 |
| Layer pairs more than 1 K superadiabatic (of 3168 columns) | 289 | 0.1 | 3202 | 3.3 |
| Strongest inversion in the lowest 6 layers [K] | 22.5 | 18.9 | 38.2 | 18.7 |
| Boundary layer height, global mean [m] | – | 884 | – | 1158 |
| Diffusivity of heat, mean / maximum [m²/s] | – | 7.3 / 350 | – | 8.5 / 542 |

- The scheme removes the superadiabatic layers and limits the inversions, which is what it was written for.
- With the scheme the 8- and 24-layer climates are much closer to each other near the surface (lowest-layer
  temperature, humidity, soil temperature, precipitation). Without it the thin lowest layer at L24 is 6 K warmer
  and 4 g/kg moister than at L8.
- The maximum wind speed rises by 5 m/s (L8) and 14 m/s (L24). This is the interaction with the jet raised in
  section 7 and has to be looked at at high resolution, where the advective limit is close.
- One year from a cold start says nothing about drift. That is left to the long runs.
- The L8 run with the scheme used the in-place variant of 8.5 (it started before the revert); the L24 runs used
  the committed code. The two differ only by rounding.

### 8.7 T512L24 cold-pool restart tests

`$D/validation/run_t512l24_restart_test.jl` (copy of the driver of the stability investigation) with
`CONFIG=vdiff`, H100, job 3342860: restart of the diagnostic run before the Lake Baikal blow-up (8 days) and
before the Committee Bay event (4 days). Lowest layer is layer 24. The other rows are the earlier runs in
`/p/tmp/maxgelbr/t512l24-reactant-10y/output/restart_tests`, made with the branch `mg/t512l24-stability`
(SpeedyWeather 0.25.0-DEV); the new run uses this branch (0.26.0-DEV), so differences between the branches other
than the scheme could play a role. They do not: the control repeated on this branch (job 3343104) reproduces
the earlier control (Baikal: NaN at the same time, coldest point −334.0 °C against −333.9 °C).

| Segment | Configuration | Ends | Coldest lowest-layer T [°C] | Strongest inversion [K] | Most points below −60 °C |
|---|---|---|---|---|---|
| Baikal | control (other branch) | NaN on 02-06 17:55 | −334 | 405 | 19 |
| Baikal | advective form (other branch) | 02-13, completed | −79 | 61 | 66 |
| Baikal | 4× hyperdiffusion (other branch) | 02-13, completed | −76 | 58 | 50 |
| Baikal | control (this branch) | NaN on 02-06 17:55 | −334 | 405 | 19 |
| Baikal | **new scheme** | 02-13, completed | **−53** | **27** | **0** |
| Committee Bay | control (other branch) | 02-01, completed | −267 | 259 | 44 |
| Committee Bay | advective form (other branch) | 02-01, completed | −65 | 51 | 2 |
| Committee Bay | 4× hyperdiffusion (other branch) | 02-01, completed | −60 | 54 | 0 |
| Committee Bay | control (this branch) | 02-01, completed | −267 | 259 | 44 |
| Committee Bay | **new scheme** | 02-01, completed | **−51** | **28** | **0** |

- With the scheme neither event happens: no NaN, no monitor event (lowest layer below −90 °C), and the coldest
  point of each run is at its start, i.e. inherited from the restart file.
- The largest superadiabatic difference between the lowest layers falls from 40–66 K in all earlier runs to 2 K.
- The restart files of the other branch load on this branch (with a version warning).
- 14 minutes for 8 days at T512L24 with 30-minute monitoring and 6-hourly output.

### 8.8 T512L24, 45 days from 2000-01-01

`$D/validation/run_t512l24_gpu_diag.jl` with `VDIFF=1` (copy of the diagnostic driver with the cold pool monitor),
H100, job 3343105, output in `$D/validation/output/run_t512l24_vdiff_0001`.

- No cold pools in 38.8 days: the coldest lowest-layer temperature of the run stays near −47 °C, the strongest
  inversion below 15 K, no point below −60 °C. The run passes 02-06, when the run without the scheme blew up over
  Lake Baikal.
- **The run ends with NaN on 02-08 20:27**, from a different failure: the wind in the top layer (layer 1) at
  49°N 270°E goes from 104 m/s to 2400 m/s within 30 minutes. The lowest layer is unremarkable until the last
  monitor step.
- The maximum wind speed grows through the run, from 37 m/s to about 100 m/s in layers 1–2 from late January on,
  above the advective limit of about 83 m/s at Δt = 150 s. The original run without the scheme shows the same
  growth (90 m/s on 01-21, 101 m/s on 01-26, 110 m/s on 01-31, all in layer 2), so the scheme does not cause
  it. This is the "top-layer jet" failure of the stability investigation (`STABILITY_NOTES.md`, finding 10:
  runs without the scheme and without cold pools die of it on 02-11), which is being worked on there with a
  shorter time step and a stronger speed-limit drag.
- Whether the scheme makes this failure more likely cannot be said from one run. It ended three days earlier than
  the comparable runs without the scheme, and the T32 runs have higher maximum winds with the scheme (8.6).
- Consequence: a 45-day or longer T512L24 validation of the scheme needs a remedy for the jet failure.
- The stability investigation tried two remedies without the scheme (with `dealiasing = 3`): the time step
  reduced from 150 s to 112.5 s (`Leapfrog(Δt_at_T32 = Minute(30))`) completed 60 days; `SpeedLimitDrag(speed_limit
  = 70, drag = 1e-4)` failed earlier, on 02-05. A drag that strong is probably unstable by itself: it is applied
  explicitly and its rate 2c(|V| − limit) exceeds 1/(2Δt) for an excess of 17 m/s (my estimate, not tested).
- Runs with the scheme, default grid, 60 days each (`$D/validation/output/run_t512l24_vdiff_*`). None has a
  NaN or a point below −60 °C in the lowest layer (coldest −53 to −55 °C, strongest inversion 20–21 K):

  | Run | Time step | Speed-limit drag | Result | Highest wind [m/s] | Speed with 6-hourly output |
  |---|---|---|---|---|---|
  | (above) | 150 s | default (80 m/s, 3·10⁻⁶ /m) | NaN on day 39 | 104 before the blow-up | 2.55 SYPD |
  | `drag1e-5` | 150 s | 80 m/s, 10⁻⁵ /m | completed 60 days | 109 | 2.28 SYPD |
  | `limit60_drag1e-5` | 150 s | 60 m/s, 10⁻⁵ /m | completed 60 days | 95 | 2.27 SYPD |
  | `dt112` | 112.5 s | default | completed 60 days | 122 | 1.91 SYPD |
  | `dt131` | 131 s | default | completed 60 days | 107 | 1.99 SYPD |

- The stronger drag at the standard time step costs no speed and keeps the run alive; the lower limit of 60 m/s
  also lowers the highest wind by 14 m/s. The shorter time steps cost an eighth to a quarter of the speed and
  let the wind go as high as or higher than in the run that failed. The global-mean temperatures of the runs with
  150 s and 112.5 s agree within 0.5 K at day 60 (layers 1, 12 and 24).
- 60 days is weak evidence for a failure that took 39 days to appear. That is what the multi-year runs are for.
- **Final verification, running**: `$D/validation/run_t512l24_years.jl` (the same driver in segments of 183 days,
  each with its own run folder and restart file, monthly output, monitor every 2 hours; `START_RUN_NUMBER=n`
  continues after a job's time limit; segments and continuation tested at T32). Three runs of 6 segments
  (3 years) with the scheme, output in `$D/validation/output/run_t512l24_years_<setup>_000n`:

  | Setup | Time step | Speed-limit drag | Job | QOS |
  |---|---|---|---|---|
  | `drag1e-5` | 150 s | 80 m/s, 10⁻⁵ /m | 3345833 | gpushort (24 h limit) |
  | `limit60_drag1e-5` | 150 s | 60 m/s, 10⁻⁵ /m | 3346081 | gpupreempt |
  | `dt112` | 112.5 s | default | 3346082 | gpupreempt |

  In their first half hour the runs with 150 s made about 1.2 simulated days per minute including compilation
  (roughly 4.8 SYPD, 15 hours for the 3 years), the one with 112.5 s about 0.8 (22 hours).

  A preempted or timed-out job is continued with `sbatch --export=ALL,RUN_ID=…,START_RUN_NUMBER=<last complete
  folder>,SEGMENTS=<remaining>,… submit_years.sh` (same options as the original, see the first line of its log
  in `$D/validation/logs/vdiff-years-<job>.log`). To evaluate: `monitor_global.csv` (extremes of all layers and
  layer-mean temperatures every 2 hours) and `monitor_stats.csv` (lowest layer) in every run folder.
- Disk: each 60-day run with 6-hourly output is 43 GB; about 250 GB are in `$D/validation/output`. The
  multi-year runs write monthly, about 0.7 GB per output time.

### 8.9 T256L8, 300 days from 2000-01-01

`run_means.jl` as in 8.6, job 3343106, `means/means_T256L8_*.csv`.

- Without the scheme the run ends with NaN between day 60 and 70 (2000-03-10). Until then it has 5800–8800 layer
  pairs that are more than 1 K superadiabatic and inversions up to 37 K in the lowest layers.
- With the scheme the run completes the 300 days, through the Antarctic winter: 6–36 superadiabatic layer pairs,
  inversions of 9–21 K. Coldest point 169.0 K, 12 points below 170 K, in the top layer on day 240. (Corrected: an
  earlier version of this line said no point below 170 K, from reading only every third output row.) Maximum wind
  speeds are 80–95 m/s, as in the run without the scheme before it failed (84–96 m/s).
- 17 minutes on the H100 for 300 days with the scheme.

### 8.10 Multi-year T512L24 runs: state and the Antarctic cold event

State at 08:15 on 2026-11-20 (`python3 $D/validation/evaluate_years.py` prints the table below for all run
folders). No run has a NaN.

| Setup | Simulated so far | Speed | Highest wind [m/s] | Coldest lowest-layer T [°C] |
|---|---|---|---|---|
| `drag1e-5` (150 s) | 1.6 years | 5.0 SYPD | 125 | −67 (Antarctic winter) |
| `limit60_drag1e-5` (150 s) | 1.6 years | 5.0 SYPD | 110 | **−93** on 2001-05-15, otherwise −74 |
| `dt112` (112.5 s) | 1.2 years | 3.8 SYPD | 132 | **−155** on 2000-06-26, otherwise −71 |
| `drag1e-5_hyper2` (150 s, hyperdiffusion 2 h), job 3348592 | started 08:05 | – | – | – |

- Layer-mean temperatures of the second year are within 1 K of the first in the same season (layers 1, 12, 24);
  the first four months include the spin-up. Drift can only be judged at the end.
- Lowest-layer temperatures of −60 to −74 °C with several hundred points below −60 °C occur in the Antarctic
  winter in all runs, on the plateau. I take these to be the model's Antarctic winter, not cold pools; this has
  not been compared with observations.
- A third event of the same kind occurred later in the `dt112` run, on 2001-08-12 at the same grid point: lowest
  layer −264 °C, inversion 159 K, again recovering.
- **Two cold events remain** (state at 08:15), both at the same grid point on the Antarctic coast (79°S, 81°W, land fraction
  0.12–0.25, next to steep orography): −155 °C with an inversion of 133 K in the `dt112` run, −93 °C in the
  `limit60_drag1e-5` run. Both recover within two days and neither ends the run. The `drag1e-5` run has had none
  in 1.6 years. The monitor's budget term shows the vertical advection warming the point by 20–40 K per hour, as
  in the events of the stability investigation, which it traced to aliasing of the flux-form horizontal
  advection on the octahedral grid (`STABILITY_NOTES.md`, finding 4). The scheme removes the decoupled polar-night
  cold pools that used to trigger this, but not this point.
- Restart tests from 2 days before each event (`run_t512l24_restart_test.jl` with the options of
  `config_components`, output in `$D/validation/output/event_tests`), 4–5 days each:

  | From the event of | Setup | Coldest lowest-layer T [°C] | Strongest inversion [K] |
  |---|---|---|---|
  | 2000-06-26 | 112.5 s (as the run) | −161 | 139 |
  | 2000-06-26 | 150 s, drag 10⁻⁵ | −140 | 120 |
  | 2000-06-26 | 112.5 s, background diffusivity 1 m²/s instead of 0.1 | −130 | 79 |
  | 2000-06-26 | 150 s, drag 10⁻⁵, hyperdiffusion 3 h | −78 | 62 |
  | 2000-06-26 | 112.5 s, hyperdiffusion 2 h | −60 (elsewhere) | 29 |
  | 2000-06-26 | 112.5 s, hyperdiffusion 1 h | −60 (elsewhere) | 25 |
  | 2000-06-26 | **150 s, drag 10⁻⁵, hyperdiffusion 2 h** | −60 (elsewhere) | 28 |
  | 2001-05-15 | 150 s, limit 60 m/s, drag 10⁻⁵ (as the run) | −94 | 51 |
  | 2001-05-15 | **150 s, drag 10⁻⁵, hyperdiffusion 2 h** | −63 (plateau) | 14 |

  The event does not depend on the time step or the drag, more vertical background mixing reduces it only a
  little, and `HyperDiffusion(time_scale = Hour(2))` instead of the default 4 hours removes it in both cases.
- Hence the fourth multi-year run: scheme, Δt = 150 s, `SpeedLimitDrag(drag = 1e-5)`,
  `HyperDiffusion(time_scale = Hour(2))`. It costs no speed. This is the setup I would propose for T512L24 if it
  completes the three years without an event.

**Update 14:45 on 2026-11-20.**

- `drag1e-5` (job 3345833) **completed the three years**: no NaN, no cold event (coldest lowest-layer temperature
  −67.5 °C, never below −75 °C in 13176 monitor times), highest wind 128 m/s, 5.0 SYPD. Half-year means of the
  layer-mean temperatures [°C], layers 1 / 12 / 24:

  | | Year 1 | Year 2 | Year 3 |
  |---|---|---|---|
  | January–June | −50.7 / −19.7 / 13.9 | −51.5 / −18.5 / 13.8 | −51.4 / −18.0 / 13.8 |
  | July–December | −50.8 / −18.1 / 14.5 | −50.8 / −17.6 / 14.5 | −50.7 / −17.3 / 14.6 |

  Layer 12 warms by about 0.4 K per year after the first year; layers 1 and 24 are steady within 0.15 K per year.
  So the run is not equilibrated in the mid-troposphere. There is no T512L24 run without the scheme to compare
  with.
- The other three jobs were **preempted** at 14:31–14:34 (`limit60_drag1e-5` 16 days before its end, `dt112` in
  segment 5, `drag1e-5_hyper2` in segment 3). Their partial run folders were moved to
  `$D/validation/output/preempted_partial/` and they are continued from the last complete segment: jobs 3372065
  (`drag1e-5_hyper2`, 4 segments, gpushort), 3372066 (`limit60_drag1e-5`, 1 segment, gpushort), 3372067 (`dt112`,
  2 segments, gpupreempt). All three were pending at 14:45: both gpushort slots are taken by other jobs of this account until about 18:00, the gpupreempt job is estimated to start at 15:40.
- The continuations run on a frozen copy of the code the runs were started with: worktree
  `$D/validation/tree_frozen` at `372eb321` with the environment `$D/validation/env_frozen`
  (`submit_years_frozen.sh`). The working tree has since changed `SpeedLimitDrag` (scaling, then damping), and
  `SPEED_DRAG=1e-5` would mean something else there.

### 8.11 Scheme as the default (step 8) and the state of the tests

Changes (committed in `b7bdfdd1`, except docs page and changelog):

- `vertical_diffusion = BulkRichardsonDiffusion(spectral_grid)` is the default of `PrimitiveDryModel` and
  `PrimitiveWetModel`; docs page and changelog say so.
- The scheme does not throw any more in a model without a boundary layer: `surface_layer` returns a zero friction
  velocity if the drag or the surface wind are missing, so only the local and the background diffusivity act.
  (A model with `boundary_layer = nothing` does not run on main for other reasons: surface fluxes, ocean albedo
  and slab ocean read variables of the boundary layer and the surface fluxes. Not touched.)
- `reset_variable!(vars, ::Val{name}, value)` instead of a `Symbol` argument. With 43 instead of 40
  parameterization variables the compiler no longer inferred the lookup by `Symbol`, which gave a runtime
  dispatch in every time step and failed `test/dynamics/dispatch.jl`.
- `test/GPU/set.jl` used `SpeedyWeather.zfac`, which does not exist any more on main; replaced by a function
  defined in the test file. Unrelated to the scheme, but the GPU suite stopped there.

Tests with the scheme as default:

| Suite | Julia | Result |
|---|---|---|
| Main suite (`Pkg.test`, `--check-bounds=yes`), before the `reset_variable!` fix | 1.12.2 | 122067 pass, 2 fail (both in `dynamics/dispatch`) |
| `dynamics/dispatch` after the fix | 1.12.2 | 4 of 4 pass |
| **Main suite on the final code** | 1.12.2 | **122069 of 122069 pass** |
| **Main suite on the final code** (copy of the tree, set up as in CI) | 1.10.13 | **122069 of 122069 pass** |
| Enzyme, `full_diff_CI.jl`, final code | 1.10.13 and 1.11.7 | pass |
| Enzyme, extended (`test/differentiability/runtests.jl`), final code | 1.10.13 | all 23 test sets pass |
| Reactant suite (`test/reactant/runtests.jl`), CPU backend, before the `reset_variable!` fix | 1.12.2 | all 28 test sets pass |
| Reactant, CPU vs Reactant for the dry and wet model, committed code | 1.12.2 | 20 of 20 and 23 of 23 pass |
| GPU suite (A40, CUDA), before the `reset_variable!` fix | 1.12.2 | all pass up to `test/GPU/reactant.jl`, see below |
| Docs build | 1.12.2 | completes, no warning about the vertical diffusion page |

- `test/GPU/reactant.jl`: the BarotropicModel with `SpectralTransform` on Reactant's GPU backend does not
  compile ("cannot raise op to stablehlo"). The same test set fails in the same way on a clean checkout of main
  (`$D/tests/gpu_reactant_barotropic_main.log`), so it is not caused by this work. The suite stops there; the
  test sets after it in that file did not run.
- The Reactant suite has to run on the CPU backend as in CI (`CUDA_VISIBLE_DEVICES=` on login02); on the A40 it
  fails in `test_transform.jl` with a relative difference of 2·10⁻³, the reduced precision of GPU matrix products.
- Logs in `$D/tests/`. The Enzyme runs follow CI and `Pkg.develop` the packages into
  `SpeedyWeather/test/differentiability/Project.toml`; that file was restored with `git checkout` afterwards.

### 8.13 Default of `SpeedLimitDrag`

The drag is explicit in time: the wind is reduced by c·excess²·Δt per step, with Δt = 2Δt of the leapfrog. For
c·excess·Δt > 1 that is more than the excess over the speed limit, and for c·excess²·Δt > |V| the wind is reversed.
The time step at T32 (2Δt = 4800 s) is 16 times that of T512 (300 s), so a coefficient that T512 needs is too
strong at T32:

| Default | T32 | T512 | Result |
|---|---|---|---|
| `drag = 3e-6` (main) | overshoots above 150 m/s | 3·10⁻⁶ | T32 fine; T512L24 blows up in the top layer on day 39 (8.8) |
| `drag = 1e-5` at all resolutions (`89be0251`) | overshoots above 101 m/s, reverses above 133 m/s | 10⁻⁵, stable (8.8, 8.10) | **T32L24 with the scheme NaN on 2000-07-05** (stable for a year with 3e-6; job 3349009); `test/output/zarr_output.jl` fails |
| `drag = 3e-6`, `resolution_scaling = 0.5` (`b3c43ac4`) | as main | 1.2·10⁻⁵ | T32 identical to main by construction; **T512L24 with the scheme completes 60 days** (job 3349380): highest wind 109 m/s, coldest lowest-layer temperature −53 °C, strongest inversion 20 K, no point below −60 °C, 2.25 SYPD with 6-hourly output |

- The failing test runs T6 with four layers for a day, with a time step of 3 hours. That run blows up on main as
  well (winds of 3700 m/s after one day without the scheme and with `3e-6`); the test only checks that the output
  is finite, and with `1e-5` it no longer was.
- Side result of job 3349009 (`$D/validation/means_drag1e-5_fixed`): T256L8 **without** the scheme completes 300
  days with `drag = 1e-5`. With `3e-6` it went NaN between day 60 and 70 (8.9). So that NaN came from the weak drag,
  not from the missing vertical mixing. What the missing mixing does at T256L8 is still there without the scheme:
  5900–12100 superadiabatic layer pairs. With the scheme and `1e-5`: 5–46 pairs, 300 days completed.
- "Points below 170 K" is not a measure of the scheme: in all T256L8 runs that reach late winter, with and without
  the scheme and for every drag, the coldest points are in the top two layers on days 210–250 (163–170 K, 0–116
  points in the 10-daily snapshots).
- T256L8 with the committed scaled drag (8.5·10⁻⁶; job 3349399, `$D/validation/means`): both runs complete 300 days;
  superadiabatic layer pairs 6000–11100 without and 4–26 with the scheme, highest wind 96 m/s in both.
- Tests with the committed scaled drag (`b3c43ac4`): main suite 122074 of 122074 on Julia 1.12.2 and 1.10.13;
  CPU vs Reactant for the dry and wet model 20 of 20 and 23 of 23; Enzyme `full_diff_CI.jl` on 1.10.13 and 1.11.7
  and the extended Enzyme tests on 1.10.13 (23 test sets) pass; T32L24 with the scheme stable for the year again
  (job 3349399). Not repeated after the drag change: the GPU suite and the full Reactant suite (the drag kernel
  is unchanged, the coefficient is computed on the host).
- The three-year T512L24 runs use `drag = 1e-5` set explicitly. The new default gives 1.2·10⁻⁵ there, 20 % more.
  The 60-day run with it behaves like the one with 10⁻⁵ (8.8: 109 m/s, −53 °C). There is no multi-year run with
  exactly the new default.
- With the scaling the coefficient is 8.5·10⁻⁶ at T256, 6·10⁻⁶ at T128 and 4.2·10⁻⁶ at T64. The overshoot threshold
  rises with resolution (excess of 69 m/s at T32, 139 m/s at T128, 278 m/s at T512), since the time step shrinks
  faster than the coefficient grows.

#### 8.13.1 Damped drag instead of the resolution scaling

The user's objection to the scaling: it is a detour, since the time step and not the resolution makes the drag
unstable, and only for large c·excess·Δt; a wind just above the limit is never a problem. Instead, the drag is
time integrated as the surface fluxes of #1336, with a damping factor φ(z) of an `AbstractFluxIntegration`:

- The drag read the *current* leapfrog step (`which_prognostic_step(var, ::AbstractLeapfrog, ::AbstractDrag) = 2`),
  so it was a centred leapfrog step, not the forward step over 2Δt the table above assumes. A damping term −a·u
  in a centred leapfrog step has a computational mode with |λ| = aΔt + √(1 + (aΔt)²) > 1, which only the Robert
  filter (ν = 0.1) holds down, and only for small aΔt. `SpeedLimitDrag` therefore first read the previous step,
  like the parameterizations, a forward step over 2Δt (`41854301`). That made T32L24 with the scheme blow up
  (8.13.2), so it reads the **current** step again (`b7156e83`), as all drags. With φ the drag of one step is
  still bounded to about the excess (rate per step aΔt·φ < 0.5), but φ = 1/(1 + z) is no longer the exact solution.
- Without other tendencies the excess speed e decays as de/dt = −c·e², a relaxation of the speed towards the speed
  limit on the time scale τ = 1/(c·e). The explicit drag is multiplied by φ(z), z = 2Δt·c·e, of the new option
  `time_integration`. Its default `ImplicitFluxes()`, φ = 1/(1 + z), is the **exact** solution of the quadratic
  drag over the time step: e(Δt) = e/(1 + c·e·Δt) = e − Δt·c·e²·φ(z). `ExponentialFluxes()` integrates the
  relaxation with c·e frozen exactly, e·exp(−z), which decays faster than the quadratic drag does;
  `RationalFluxes()` approximates that and lets the speed undershoot the limit by up to 10 % of the excess. With
  each of them the wind is never decelerated beyond the speed limit, for any c, e or Δt.
- Hence one coefficient for all resolutions, `drag = 1e-5` (what T512L24 needs, 8.8). `resolution_scaling` and
  `speed_limit_drag_coefficient` are removed. At T512 (2Δt = 300 s) and 20 m/s excess z = 0.06, so the drag is
  6 % weaker than the undamped one of the three-year runs. At T32 (2Δt = 4800 s) z = 0.96 and φ = 0.51 for
  20 m/s excess.
- `AbstractFluxIntegration`, its four subtypes and `flux_damping` moved from `parameterizations/tendencies.jl` to the
  new `time_stepping/flux_integration.jl` (since renamed, see the revision log of 2026-11-21), which is included before the dynamics. `flux_damping(::ExponentialFluxes, z)`
  now guards numerator and denominator separately: the drag has z = 0 wherever it is inactive, where the former
  `ifelse(small, 1, -expm1(-z)/z)` evaluated 0/0 in the branch that is not taken (NaN derivatives).
- Tests: `test/dynamics/forcing_drag.jl` (antiparallel with `ExplicitFluxes()`, the step the drag reads, one step
  of every time integration against its closed form for z from 0.02 to 200) and
  `test/parameterizations/surface_flux_time_integration.jl` pass, with the previous and with the current step.
- With the previous-step version (`41854301`): `test/output/zarr_output.jl` passes (116 of 116; it failed with the
  undamped 1e-5), CPU vs Reactant 20 of 20 (dry) and 23 of 23 (wet) (`$D/tests/reactant_correctness_dragdamped.log`),
  Enzyme `full_diff_CI.jl` on 1.11.7 passes (run from a copy of the environment, the tracked
  `test/differentiability/Project.toml` untouched). Not repeated for the current step, which only changes the
  step that is read.
- Validation runs: 8.13.2.

#### 8.13.2 T32L24 blow-ups and the leapfrog step of the drag

First validation, with the damped drag on the **previous** step (`41854301`; T512L24 still running):

| Run | Without the scheme | With the scheme |
|---|---|---|
| T32L8, 360 days (GPU; with the scheme also CPU) | completes | completes |
| T32L24, 360 days, GPU (`$D/validation/means_drag_damped`, job 3371652) | NaN 2000-08-10 | NaN 2000-07-21 |
| T32L24, 360 days, CPU (`$D/validation/means_drag_damped_cpu`) | completes | NaN 2000-07-23 |
| T256L8, 300 days, GPU | completes, max wind 99 m/s, 6000–10800 superadiabatic layer pairs | completes, max wind 96 m/s, 2–30 pairs |
| T512L24, 60 days, GPU (`output/run_t512l24_vdiff_dragdamped_0001`) | – | **PENDING** (started 15:37, loaded the previous-step code) |

Where it blows up (`$D/validation/t32_blowup_probe.jl`, `logs/t32-probe-spinup140.log`): a single 200-day run of
T32L24 with the scheme went NaN on 2000-05-23. Restarted from its state on 2000-05-20 (day 140) and continued in one
`run!`, it goes NaN on 2000-05-22T15. Not in the jet aloft (top layers stay near 100 m/s) but in the **lowest
layers** (20–24 of 24) at the Antarctic coast, 72–80°S, 150–170°E (Ross Sea, Victoria Land), in southern winter:
a lowest-layer wind of about 115 m/s for a day, then 148 → 186 → 279 → 480 → 1285 → 8238 m/s in six steps (×1.6 per
40-minute step). The sign does not alternate between steps, so it is not the computational mode.

Restart test, 5 days from that state for every drag variant (`t32_blowup_probe.jl restart`,
`logs/t32-probe-restart140.log`):

| Drag | Step | Result | Peak wind in the lowest layers [m/s] |
|---|---|---|---|
| 3e-6, undamped (main) | current | completes | 109 |
| 1e-5, undamped (`89be0251`) | current | completes | 106 |
| 1e-5, damped (now default) | current | completes | 108 |
| 1e-5, damped (`41854301`) | previous | NaN | 115, then grows |
| 1e-5, `RationalFluxes()` | previous | NaN | 116, then grows |
| 1e-5, undamped | previous | NaN | 174 within 6 h, then grows |
| 3e-6, damped | previous | NaN | 113, then grows |
| none (`drag = 0`) | – | NaN | – |
| 1e-5, damped, without the scheme | previous | completes (the event does not occur) | – |

The drag on the previous step does act (the top-layer winds of the run with it are 100–110 m/s, 115–123 m/s without
drag), but every variant on the current step keeps the lowest-layer jet below 110 m/s and every variant on the
previous step lets it grow past 113 m/s and blow up, like no drag at all. Why is not understood. A guess (not
tested): the jet near the surface becomes unstable above roughly 110–115 m/s at this time step, and a drag on the
previous step lags one step behind a jet that grows within hours.

One-year ensemble with the scheme (CPU, `$D/validation/run_t32_ensemble.sh`, `logs/t32ens/`, code after
`8f1520f0`): four realizations per drag that differ only in the length of the `run!` calls (every call starts the
leapfrog again with an Euler step, which perturbs the trajectory):

| Drag | One call | 10-day calls | 30-day calls | 90-day calls |
|---|---|---|---|---|
| 3e-6, undamped, current step (main) | NaN 2000-06-10 | NaN 2000-07-15 | completes | completes |
| 1e-5, damped, current step (now default) | completes | completes | completes | completes |
| 1e-5, damped, previous step (`41854301`) | NaN 2000-05-23 | NaN 2000-07-23 | completes | completes |

- **T32L24 with the scheme is marginal with main's drag too**: 2 of 4 years blow up. The earlier stable year with
  3e-6 (8.6, 8.13) was one lucky realization. The 2000-07-05 NaN with the undamped 1e-5 (8.13) can no longer be
  taken as evidence against that coefficient.
- With the damped 1e-5 on the previous step 2 of 4 years blow up, as with main's drag (maximum winds of 103–106 m/s
  at the ends of the calls of the two that complete).
- With the damped 1e-5 on the current step all four years complete, with maximum winds of 89–95 m/s at the ends of
  the calls. Four realizations do not prove stability; the instability in the lowest layers at the Antarctic coast
  is still there and only less likely (0 of 4 against 2 of 4 for the other two drags).
- Open: the cause of that instability (scheme, thin lowest layers, orography at T32; without the scheme it did not
  occur in the restart test and one of two years blew up in the first validation).

### 8.14 Tibet at L24

The draft of #1283 reported T31L24 NaN after about 12 hours over Tibet and attributed it to the explicit coupling
of land and a thin lowest layer. Here T32L24 with the scheme is stable for a year (8.6, 8.13). Looked at directly
with `$D/validation/tibet_oscillation.jl`: lowest-layer and soil temperature at every time step for 10 days after a
20-day spin-up, at the 10 land points of 28–38°N, 78–100°E above 3000 m. A 2Δt oscillation would make the second
difference in time about twice the first difference.

| | Lowest-layer T: mean second / first difference | largest second difference [K] | Soil T: mean second / first difference |
|---|---|---|---|
| T32L8 without the scheme | 0.32 | 0.30 | 0.19 |
| T32L8 with | 0.33 | 0.29 | 0.19 |
| T32L24 without | 0.30 | 0.86 | 0.19 |
| T32L24 with | 0.37 | 0.93 | 0.19 |

No oscillation with or without the scheme, and no NaN. Ten points and ten days in January; not a proof for all
seasons.

### 8.15 Audit of the scheme's variables (2026-11-20, user request)

Question: are all six variables necessary, or can some be computed on the fly? Each 3D array is 58 MB at T512L24
in Float32. Outcome: the scheme's dedicated 3D arrays went from five to two.

| Variable (before) | Finding | Now |
|---|---|---|
| `boundary_layer_height` (2D) | diagnostic only, read nowhere | kept, it is 2D |
| `vertical_diffusivity_momentum`, `_heat` | needed as storage between the loop that computes them and the four solves | replaced by `vertical_conductance_momentum`, `_heat`, the product ρK/Δz [kg/m²/s] |
| `density_over_distance` | only ever read as a product with the diffusivity | **removed**, part of the conductances |
| `tridiagonal_upper`, `tridiagonal_rhs` | inherent to the Thomas algorithm, but no dedicated variable needed | **removed**, the model's general scratch arrays `vars.scratch.grid.a` and `.b` are used (write before read within the call, as convection and the radiation transmissivities use them) |

Decisions of the user along the way:

- Computing the diffusivity on the fly inside the solves: not done. With four separate solves it would be
  computed four times, and the virtual temperature would need an array of its own for the duration of the solves.
- One diffusivity for momentum and heat: **not wanted** (it fixes the turbulent Prandtl number at 1). The two
  stability functions stay.
- Storing conductances instead of diffusivities: done (`8f1520f0`). The diffusivity in m²/s is no longer stored;
  it is K = conductance · Δz/ρ. The unit tests recompute ρ/Δz from the state to check K. `vertical_diffusivity!` is
  now `vertical_conductance!`, and `implicit_vertical_diffusion!` takes the conductance.

Checks, for the scratch arrays (`6bc132e2`) and then for the conductances (`8f1520f0`), because the Reactant
workaround of 8.3 is sensitive to how the loops read their arrays:

| | General scratch arrays | Conductances |
|---|---|---|
| Unit tests of the scheme, `--check-bounds=yes` | 68 of 68 | 71 of 71 |
| Reactant, whole `vertical_diffusion!` as one kernel, raised against not raised (wet coupled, wet uncoupled, dry coupled) | relative differences 10⁻⁷ to 7·10⁻⁵ | 10⁻⁷ to 6·10⁻⁵ |
| Same check: boundary layer height and the four tendencies against the version before | same digits | same digits |
| Reactant, CPU vs Reactant over 100 steps, PrimitiveDryModel / PrimitiveWetModel | 20 of 20 / 23 of 23 | 20 of 20 / 23 of 23 |

- `Variables` of the dry and wet model: no scratch array of the scheme any more, two 3D parameterization
  variables (the conductances) and the 2D boundary layer height.
- Each solve reads one array per layer where it read two; the effect on the GPU timing is not measured.
- Not repeated with these two changes: main suite, CUDA, Enzyme, GPU suite. Scripts and logs:
  `$D/reactant_conductance_check.jl`, `$D/reactant_correctness_conductance.jl`, `$D/tests/*scratch_ab*`,
  `$D/tests/*conductance*`.

### 8.16 Boundary layer height: Ri_c = 10 and the surface excess terms (2026-11-20)

Typical observed boundary layer heights (radiosonde climatologies, e.g. Seidel et al. 2010, 2012): stable
nocturnal over land 50–300 m, polar winter 10–200 m, convective daytime over mid-latitude land 1–2 km in summer,
subtropical deserts 3–5 km, stratocumulus regions 0.5–1 km, trade-wind cumulus 1–2 km, mid-latitude oceans
0.5–1.5 km. With Ri_c = 1 the T32 runs of 8.6 had global means of 0.9–1.2 km.

With Ri_c = 10 (user's runs, not in this protocol) the tropical boundary layer is 1–2 km deep and the tropics
improve markedly (revision log); under the mid-latitude jet h exceeds 10 km. Explanation: the bulk Richardson
number of the height uses the full wind speed |V(z)|² in the denominator, without any buoyancy excess of the
parcel. In the tropics winds are weak and Ri reaches 1 just above the lowest layers; under the jet |V|² is large
throughout the troposphere and Ri stays below 10 up to the tropopause.

Tried (not committed, reverted by the user): the standard remedy of the literature, with options for each term,
- b u⋆² (b = 100) added to the squared wind,
- a thermal excess θ_T = b′(w′θᵥ′)₀/w_m (b′ = 8.5, at most 3 K) of the rising parcel for an upward surface
  buoyancy flux from the current sensible heat and humidity fluxes, w_m = (u⋆³ + 0.6w⋆³)^⅓ with w⋆ from the
  boundary layer height of the previous time step (to avoid the circularity h ↔ θ_T),
- the wind difference to the lowermost layer |V(z) − V_N|² instead of |V(z)|² (`relative_wind = true`).

Result (user's run): the tropics fell back to the state of Ri_c = 1, "a sluggish soup": little dynamics, low
winds, no cyclones, large areas of weak convective rain. Not investigated further. A likely cause, not tested:
the relative wind raises Ri where the lowermost layer itself is windy (the trades) and the shear above is weak,
so the tropical h became shallower, not deeper. The improvements of Ri_c = 10 thus seem to come from the
full-wind formulation with a large Ri_c, not from the physics that the excess terms represent. Whether the thermal
excess alone (with the full wind) would help was not tried.

Open: the effect of the deep mid-latitude h (profile diffusivity through the troposphere, momentum mixing across
the jet) on the jet strength and the T512L24 stability (8.8, 8.10), which was validated with Ri_c = 1.

### 8.17 Performance audit (2026-11-21, user request)

> State 2026-11-21 18:00: **audit complete**, no run is left. All measurements are done and the unit tests pass
> with the bounds-check patches, both with `@inbounds` active and with `--check-bounds=yes` (8.17.5). Nothing is
> implemented; candidates and recommendation in 8.17.6. Not done: Enzyme, and the Reactant and GPU suites with
> the patches.

Request: audit the vertical diffusion for further optimisations, for both CPU and GPU. Nothing in the source is
changed by the audit. All variants are prototypes outside the repository, in
`$A = /p/tmp/maxgelbr/code/Speedy2/claude_vdiff_perf_2026-11-21/`, each a function of its own that is compared
against the scheme of `007f04d9` (frozen copies `$A/tree` for CPU, `$A/tree_gpu` for GPU and Reactant). Two
sessions work on it in parallel, because the user's "continue" forked the session: `speedyweather-jl-a8` does the
CPU side (8.17.2), `speedyweather-jl-d6` GPU and Reactant (8.17.3, 8.17.4) and this protocol. After a second
"continue" at 13:15 `speedyweather-jl-6a` took over from d6, which stopped. Session a8 ended at 13:18 with three
of its runs unfinished; 6a restarted them (8.17.5).

#### 8.17.1 Where the time goes

CPU: EPYC 9354 (login02), one thread, Float32, default `PrimitiveWetModel`, T32 L8 (3168 grid points), after 3
days of spin-up, minimum of 30 calls (`$A/cpu_baseline.jl`, `cpu_inbounds.jl`, `cpu_profile.jl`).

- A time step takes 13.8 ms with the scheme and 10.9 ms without: **the scheme is 21 % of a T32 L8 time step on
  the CPU** (2.9 ms). On the H100 it is 1.4–2.8 % (8.17.3), because there the transforms and the dynamics
  dominate.
- It is the most expensive column parameterization:

  | Parameterization | Time [ms] | Share of the time step |
  |---|---|---|
  | `vertical_diffusion` | 3.07 | 22 % |
  | `radiation` | 2.57 | 19 % |
  | `boundary_layer` | 1.68 | 12 % |
  | `convection` | 0.65 | 5 % |
  | `large_scale_condensation` | 0.51 | 4 % |
  | `albedo`, the three surface fluxes | 0.62 | 4 % |
  | all column parameterizations | 9.1 | 66 % |

- Inside the scheme (loops with `@inbounds`, 2.65 ms): the conductances 1.62 ms (61 %), the four solves 0.23 to
  0.25 ms each (0.95 ms, 36 %). That is 121 ns per grid point and layer with the bounds checks of the package's
  loop, 104 ns without.
- Line profile of the conductances: per interface about 16 divisions, 4 square roots, 2 `exp` (mixing length decay,
  background diffusivity) and 1 `log` (stability correction, evaluated for every interface although it is only
  used in the surface layer, the lowest 10 % of the boundary layer). `exp` and `log` are about a fifth of the
  scheme's samples. Per column, `bulk_richardson_surface` (two `saturation_humidity`) is computed a second time:
  the drag of the boundary layer scheme has just computed the same surface Richardson number.
- The solves are limited by the latency of the divisions in the pivot recurrence (3 divisions per layer and
  variable), not by memory: four separate sweeps, and u and v (and T and q) recompute identical pivots.

GPU (H100, `007f04d9`, 8.17.3): the scheme adds 0.52 ms to the 0.79 ms fused column kernel at T256 L24 (2.8 % of
the 18.5 ms time step) and 1.83 ms to 2.74 ms at T512 L24 (2.3 % of 78.9 ms). In 8.4, before the conductance
and scratch changes, it was 0.62 ms and 2.25 ms.

#### 8.17.2 CPU (written by session `speedyweather-jl-a8`, pending numbers completed by `speedyweather-jl-6a`)

login02 (AMD EPYC 9354, AVX-512), one thread, Julia 1.12.2, Float32, default `PrimitiveWetModel` after 3 days, frozen
copy of `007f04d95` (`$A/tree`, `$A` as above). Times are the
minimum of 30–40 calls. Scripts: `$A/cpu_baseline.jl` (cost in the model), `cpu_inbounds.jl`, `cpu_variants.jl`
and `cpu_range.jl` (prototypes in `variants.jl` and `range_variants.jl`, which do not override the package).

**Cost in the model today** (`cpu_baseline.jl`; "scheme" is `column_parameterization_cpu!` of the scheme as the
model calls it, i.e. with bounds checks, see 8.17.5):

| | Time step [ms] | without the scheme [ms] | All column parameterizations [ms] | Scheme [ms] | Scheme, share of the step |
|---|---|---|---|---|---|
| T32 L8 | 13.87 | 10.90 | 9.13 | 3.14 | 22 % |
| T32 L24 | 37.85 | 29.84 | 20.92 | 8.65 | 23 % |
| T64 L8 | 52.63 | 42.68 | 31.45 | 10.74 | 20 % |
| T128 L8 | 226.0 | 187.5 | 116.4 | 38.3 | 17 % |

- On the CPU the scheme is the most expensive parameterization (T32 L8: 3.1 ms, radiation 2.7 ms, boundary layer
  1.7 ms, convection 0.65 ms), and the column parameterizations together are two thirds of a time step. This is
  the opposite of the GPU, where the scheme is 1.4–2.8 % of a step (8.17.3).
- It costs 100–108 ns per grid point and layer at all four configurations (without bounds checks), so it is
  limited by arithmetic, not by memory: T128 L8 is hardly slower per point than T32 L8.
- Parts at T32 L8 without bounds checks (2.65 ms): conductances 1.63 ms (61 %), the four solves 0.92 ms (35 %).

**Variants, whole scheme** (`cpu_variants.jl`, `cpu_range.jl`; all loops with `@inbounds`, so the first column is
already faster than the model today). Largest relative difference to the package (relative to the largest value
of the field) over all four configurations: tendencies 9·10⁻⁵, conductances 2·10⁻⁵, boundary layer height
9·10⁻⁷; the solver changes alone give 6·10⁻⁷.

| [ms] | Package | Column variant | Range form, Base `exp`/`log` | Range form, vectorised `exp`/`log` and height search | same, surface Ri recomputed |
|---|---|---|---|---|---|
| T32 L8 | 2.66 | 1.92 | 1.27 | 0.18 | 0.29 |
| T32 L24 | 7.57 | 5.41 | 3.95 | 0.48 | 0.58 |
| T64 L8 | 9.28 | 6.59 | 4.40 | 0.62 | 1.00 |
| T128 L8 | 34.8 | 25.9 | 16.3 | 2.34 | 3.74 |
| relative to the package | 1 | 0.71–0.74 | 0.47–0.52 | 0.06–0.07 | 0.08–0.11 |

- **Column variant** (same structure as now, one column at a time): the four changes of the conductances below,
  and the solves of u and v, and of dry static energy and humidity, in one sweep each with a shared pivot
  (`implicit_vertical_diffusion_shared!`: one upper diagonal and one division for the pivot per layer for both
  variables, right-hand side stored in the tendency for the coupled solve, so one scratch array is enough).
- **Range form** (`range_variants.jl`): the scheme as a function of a range of columns, layers in the outer
  loops and columns in the inner `@simd ivdep` loops. What is carried from layer to layer (pivot, right-hand
  side, conductance of the interface above) is read from the arrays at k−1, the per-column scalars (friction
  velocity, stability ratio, stability correction of the surface layer) are in three 2D scratch arrays. Called
  with `1:npoints` the inner loops vectorise across columns; called with `ij:ij` it is the per-column scheme of
  the GPU kernel (8.17.3). With Base's `exp` and `log` only the solves and the virtual temperature vectorise
  (LLVM does not vectorise a loop with these calls); with branch-free Float32 versions (`exp_branchfree`,
  `log_branchfree`, after Cephes) and the boundary layer height searched with the layers in the outer loop,
  everything but the surface layer does.
- The range form run column by column on the CPU (as the GPU kernel would run it) takes 2.40 / 6.76 / 8.21 /
  31.6 ms: 10 % faster than the package, 20–25 % slower than the column variant, because the carried values go
  through memory. Blocks of 64 to 1024 columns are not faster than all columns at once.

**Parts at T32 L8** [ms]:

| Conductances (column loop) | | Solves | |
|---|---|---|---|
| package | 1.63 | package, 4 solves | 0.92 |
| same code as a variant (control) | 1.56 | shared pivots, 2 sweeps | 0.57 |
| reciprocals instead of about half of the divisions per interface | 1.54 | one sweep for all four variables | 0.50 |
| boundary layer height without a division per layer | 1.52 | range form over all columns (vectorised) | 0.04 |
| `log` of the stability correction only in the surface layer (`if`, not `ifelse`) | 1.44 | | |
| surface Richardson number and friction velocity read, not recomputed | 1.39 | | |
| all four | 1.25 | | |
| all four, without the two `exp` (wrong result, bound) | 0.93 | | |
| range form over all columns, Base `exp`/`log` | 1.21 | | |
| range form, vectorised `exp`/`log` and height search | 0.14 | | |

- The solves are limited by the chain pivot → upper diagonal → next pivot with one division per layer. Two
  variables per sweep share that chain (−38 to −41 % over the four configurations); all four in one sweep run
  both chains at once (−44 to −47 %). Vectorised over the columns the chain no longer matters: 23 times faster
  at T32 L8, 21 times at the other three configurations.
- `bulk_richardson_surface` is evaluated twice per column and time step, in the drag of the boundary layer and
  again in `surface_layer` of the scheme (two saturation humidities each). Storing the surface Richardson number
  (and the friction velocity) in the boundary layer scheme saves 11 % of the conductances in the column loop at
  L8 and 4 % at L24 (it is a cost per column), and 17 % (L24) to 38 % (L8) of the whole scheme once everything
  else is vectorised.
- The logarithm is only needed for interfaces below the surface layer height (0.1 h), at L8 almost none. With
  `ifelse` it is evaluated for every interface (7–8 % of the conductances).
- Fewer divisions gain 1–4 % only: the divider is not the bottleneck in the column loop.
- The two exponentials (mixing length decay and background diffusivity) are 20 % of the conductances. Both are
  functions of the height of the interface only; they could only be removed by changing the formulation
  (profiles from reference heights, or one height scale for both), which was not tried.
- `Base.sqrt` checks its argument and throws, which blocks vectorisation. The range form uses `Base.sqrt_llvm`;
  `sqrt(abs(x))` has the same effect with public API.
- A solver written with `map` and closures over tuples of fields was 3–4 times slower than the package's
  (no allocations; why was not looked into). With recursion over the tuples instead (`forward_steps!`) it is
  faster than hand-written code for two variables (0.57 against 0.74–0.77 ms).

Not measured on the CPU: Float64 (the vectorised `exp`/`log` of the prototype are Float32 only and fall back to
Base otherwise), Enzyme through the range form, more than one thread.

#### 8.17.3 GPU (sessions `speedyweather-jl-d6`, `speedyweather-jl-6a`)

`$A/gpu_variants.jl`, `submit_gpu.sh`: every variant as a kernel of its own over all columns, compared with the
package's result on the same state, minimum and median of 30 calls; additionally the fused column kernel of the
model with the variant in place of the scheme.

- A40 on login02, T64 L8 and T128 L8 (`$A/gpu_variants_a40_*.log`): all variants compile and agree with the
  package (tendencies 2·10⁻⁶, conductances 6·10⁻⁶ to 4·10⁻⁵ relative, the shared-pivot solves alone 10⁻⁷),
  including the range form called per column, also with the layer-outer height search and the branch-free
  `exp`/`log`. No ranking at these sizes: a kernel that does nothing but takes the same arguments
  (`vars`, `diffusion`, `model`) already needs 0.15 ms per launch there, as much as every variant.
- H100, T256 L24, T512 L24 and T512 L8 (job 3387439, `$A/gpu_variants-3387439.log`), minimum of 30 calls. The
  median is 0.12 ms higher for every kernel; the differences between variants are the same in both.

  | | T256 L24 | T512 L24 | T512 L8 |
  |---|---|---|---|
  | Time step with / without the scheme [ms] | 18.52 / 18.00 (0.52, 2.8 %) | 78.91 / 77.09 (1.81, 2.3 %) | 40.71 / 40.12 (0.59, 1.4 %) |
  | Fused column kernel with / without the scheme [ms] | 1.311 / 0.787 | 4.572 / 2.744 | 1.734 / 1.153 |
  | Empty kernel with the same arguments (launch) [ms] | 0.114 | 0.106 | 0.117 |

  "Own" is the variant as a kernel of its own (launch included), "fused" the model's fused column kernel of all
  parameterizations with the variant in place of the scheme, both in ms:

  | Variant | T256 L24 own | T256 L24 fused | T512 L24 own | T512 L24 fused | T512 L8 own | T512 L8 fused |
  |---|---|---|---|---|---|---|
  | package | 0.552 | 1.314 | 1.755 | 4.575 | 0.479 | 1.740 |
  | package, `@propagate_inbounds` accessors (8.17.5) | 0.541 | 1.286 | 1.719 | 4.493 | 0.466 | 1.688 |
  | shared pivots (u with v, T with q) | 0.465 | 1.198 | 1.507 | 4.189 | 0.416 | 1.606 |
  | one sweep for all four | 0.511 | 1.259 | 1.707 | 4.402 | 0.427 | 1.612 |
  | reciprocals and division-free height loop, package solves | 0.467 | | 1.482 | | 0.413 | |
  | **reciprocals + shared pivots** | **0.387** | **1.080** | **1.257** | **3.749** | **0.350** | **1.508** |
  | reciprocals + one sweep | 0.434 | 1.136 | 1.455 | 3.955 | 0.358 | 1.505 |
  | reciprocals, `log` branch, stored surface Ri + shared pivots | 0.379 | | 1.231 | | 0.336 | |
  | range form per column: reciprocals, `log` branch | 0.397 | | 1.266 | | 0.379 | |
  | range form per column: + stored surface Ri | 0.386 | | 1.239 | | 0.356 | |
  | range form per column: + layer-outer height search | 0.431 | | 1.381 | | 0.399 | |
  | range form per column: + branch-free `exp`/`log` | 0.443 | | 1.411 | | 0.409 | |
  | solves only: package / shared pivots / one sweep | 0.346 / 0.165 / 0.239 | | 1.085 / 0.622 / 0.727 | | 0.321 / 0.163 / 0.211 | |
  | conductances only: package / all four options | 0.173 / 0.128 | | 0.560 / 0.383 | | 0.179 / 0.123 | |

  - Reciprocals with shared pivots take the scheme inside the fused kernel from 0.52 ms to 0.29 ms at T256 L24,
    from 1.83 ms to 1.00 ms at T512 L24 (−45 % at both) and from 0.58 ms to 0.35 ms at T512 L8 (−40 %), which
    is 1.3 %, 1.0 % and 0.6 % of a time step. About half of that comes from the solves and half from the
    conductances: the arithmetic matters on the GPU as well, the kernel is not limited by memory traffic alone
    as assumed in 8.4.
  - The single conductance options are not additive on the GPU: reciprocals alone are slower than the control as
    a kernel of their own (0.663 against 0.559 ms at T512 L24), while reciprocals together with the division-free
    height loop save 0.27 ms in the kernel of the whole scheme. Only the combinations in the table are reliable.
  - One sweep for all four variables is slower than two sweeps with shared pivots at L24 in every pairing and
    equal at L8 (on the CPU it is slightly faster, 8.17.2).
  - **The range form called per column is as fast as the column variant with the same options** (0.386 against
    0.379 ms, 1.239 against 1.231 ms, at L8 0.356 against 0.336 ms), so reading the carried values back from the
    arrays costs at most 0.02 ms on the H100. Its two CPU-oriented options (height search with the layers
    outside, branch-free `exp`/`log`) cost 0.05 ms at T256 L24 and T512 L8 and 0.15 ms at T512 L24 and belong on
    the CPU only.
  - The `log` branch and the stored surface Richardson number make no measurable difference on the GPU.
- Second job 3390072 (`$A/gpu_extra.jl`, `$A/gpu_extra-3390072.log`), cumulative stages of the package's scheme,
  minimum of 30 calls in ms. The package measured twice in a row differs by 0.003 ms at most. Stages 1 to 3 give
  results identical to the package.

  | Stage | T256 L24 scheme kernel | T256 L24 fused | T512 L24 scheme kernel | T512 L24 fused | Bounds checks, scheme / fused |
  |---|---|---|---|---|---|
  | 0 package | 0.558 | 1.313 | 1.758 | 4.595 | 26 / 32 |
  | 1 + `@propagate_inbounds` on `diffused_value`, `column_humidity` | 0.545 | 1.289 | 1.725 | 4.510 | 10 / 16 |
  | 2 + `@propagate_inbounds` on `get_step`, `get_prognostic_step`, `get_tendency_step` | 0.543 | 1.271 | 1.726 | 4.436 | 0 / 0 |
  | 3 + virtual temperature computed where it is read, no scratch loop | 0.529 | 1.246 | 1.683 | 4.350 | 0 / 0 |
  | 4 reciprocals + shared pivots in place of the scheme, with 1 and 2 | 0.392 | 1.070 | 1.255 | 3.708 | 8 / 8 (from the prototype) |

  - The two bounds-check fixes take 3.2 % (T256 L24) and 3.5 % (T512 L24) off the fused column kernel of all
    parameterizations, 0.04 and 0.16 ms, which is 0.2 % of a time step. Stage 2 does not change the scheme's own
    kernel; its gain in the fused kernel is in the other parameterizations.
  - Without the scratch loop for the virtual temperature (possible once stage 1 has removed the cause of the
    Reactant failure, 8.17.4) the fused kernel is another 0.025 and 0.086 ms faster.
  - Stage 4 does not include stage 3; together they would put the fused kernel at about 1.05 and 3.62 ms, i.e. the
    scheme at half of its present cost (not measured as one variant).

#### 8.17.4 Reactant (sessions `speedyweather-jl-d6`, `speedyweather-jl-6a`)

`$A/reactant_variants.jl`, as in 8.3: the variant as one kernel over all columns, compiled without and with
raising on the same state, the outputs compared (conductances, boundary layer height, the four tendencies).
Reactant 0.2.308, CPU backend, T32 L8, `PrimitiveWetModel`, coupled solve (`$A/reactant_*.log`, lines `RESULT`).
Two variants are added to those of 8.17.2: "merged", the virtual temperature computed where it is read in both
loops instead of the scratch loop (what made the raised kernel write nothing in 8.3), and "carried", the
package's solver with `x = x_below` carried through the loop (what could not be raised in 8.5). Each is run with
the package's `@inline` accessors `diffused_value` and `column_humidity` and with `@propagate_inbounds` ones
(8.17.5).

| Variant | Accessors `@inline` (package) | Accessors `@propagate_inbounds` |
|---|---|---|
| package | agrees, 1.1·10⁻⁴ | agrees, 1.1·10⁻⁴ |
| merged (no scratch loop for the virtual temperature) | **raised kernel writes nothing** | agrees, 1.1·10⁻⁴ |
| carried (`x = x_below`) | **does not compile** ("cannot raise op to stablehlo", `scf.while`) | agrees, 1.1·10⁻⁴ |
| shared pivots | agrees, 1.1·10⁻⁴ | agrees, 1.1·10⁻⁴ |
| one sweep | agrees, 1.1·10⁻⁴ | agrees, 1.1·10⁻⁴ |
| reciprocals + shared pivots | | agrees, 1.6·10⁻⁴ |
| + `log` branch (`if`) | | agrees, 1.6·10⁻⁴ |
| + stored surface Ri | | agrees, 1.3·10⁻⁴ |
| range form per column | | agrees, 1.6·10⁻⁴ |
| range form + layer-outer height search | | agrees, 1.6·10⁻⁴ |
| range form + branch-free `exp`/`log` | | agrees, 1.6·10⁻⁴ |

("agrees": largest difference between raised and not raised, relative to the largest value of the field.)

- **The two Reactant failures of 8.3 and 8.5 are caused by the bounds checks of the two accessors.** Both are
  reproduced with the package's accessors and gone with `@propagate_inbounds`. The mechanism matches what is
  known from the Reactant test suite, which cannot raise kernels under `--check-bounds=yes`. Consequences: the
  separate virtual temperature loop (the workaround of 8.3, with its `TODO` on the scratch memory) and the rule
  "no layer values carried through a loop" (8.5) are not needed once the accessors are fixed, and the committed
  code only works because of the workaround, not because the cause was removed.
- Every prototype is raised correctly, including the range form with its `@simd ivdep` loops, `Base.sqrt_llvm`
  and the bit operations of the branch-free `exp`/`log`.
- `PrimitiveDryModel` with `@propagate_inbounds` accessors: package, merged and carried agree, 1.1·10⁻⁵
  (`$A/reactant_propagate_dry.log`).
- Not checked: the uncoupled solve of the prototypes, Reactant on the GPU backend, CPU against Reactant over 100
  steps (`test_correctness.jl`), Enzyme.

#### 8.17.5 Bounds checks

Three independent gaps in the `@propagate_inbounds` chain, all of them also on main. None changes a result.

1. **CPU, all column parameterizations.** `column_parameterization_cpu!` and everything it calls is
   `@propagate_inbounds`, but nothing at the root of the chain says `@inbounds`
   (`column_parameterizations_cpu!` and `_column_parameterization_cpu_by_name!` are plain functions), so on the
   CPU every column parameterization runs with bounds checks, also before the `SPLIT_COLUMN_LOOPS` workaround.
   The GPU kernel has `inbounds = true`. Effect of `@inbounds` on the loop over the grid points, T32 L8
   (`$A/cpu_inbounds.jl`):

   | Parameterization | as in the model [ms] | `@inbounds` loop [ms] | |
   |---|---|---|---|
   | boundary layer | 1.68 | 1.63 | −3 % |
   | surface momentum / heat / humidity flux | 0.05 / 0.11 / 0.18 | 0.03 / 0.05 / 0.13 | −46 / −53 / −28 % |
   | vertical diffusion | 3.10 | 2.65 | −15 % |
   | large-scale condensation | 0.50 | 0.44 | −11 % |
   | convection | 0.65 | 0.57 | −12 % |
   | albedo | 0.29 | 0.29 | 0 % |
   | radiation | 2.68 | 2.36 | −12 % |
   | sum | 9.22 | 8.15 | −12 % |

2. **The scheme, all backends.** `diffused_value` and `column_humidity` are `@inline`, not
   `@propagate_inbounds`, so their `getindex` keeps its bounds check under any `@inbounds` (it only reaches one
   level of inlining): three checked reads per layer and variable in the solver and one in the virtual
   temperature loop. They are the only helpers of the parameterizations that index an array without propagating.
   Effect of `@propagate_inbounds` on the two alone: on the CPU nothing in the scalar column loop (four solves
   at T32 L8 0.92 ms with either), but their checks are what kept the layer-outer loops of the range form from
   vectorising; on the H100 the fused column kernel 1.314 → 1.286 ms at T256 L24, 4.575 → 4.493 ms at T512 L24
   and 1.740 → 1.688 ms at T512 L8 (the four solves as a kernel of their own −18 to −33 %); **on Reactant it
   removes the cause of the two raising failures** (8.17.4).
3. **Every column parameterization, all backends.** `get_step`, `get_prognostic_step` and `get_tendency_step`
   (`time_stepping/steps.jl`) are `@inline` as well, so the `view` they create keeps its `checkbounds`, once per
   column and call. They are called 98 times in `src/parameterizations`. All bounds checks left in the scheme's
   kernel after item 2 are these (`$A/gpu_boundschecks.jl`: 4 in `vertical_conductance!`, 4 in
   `vertical_diffusion!`, 2 in `bulk_richardson_surface`). Calls of `throw_boundserror` in the CUDA device code
   (A40, T64 L8, `$A/gpu_extra_a40_test.log`), with identical results in all three rows:

   | | Kernel of the scheme | Fused column kernel of all parameterizations |
   |---|---|---|
   | package | 26 | 32 |
   | + `@propagate_inbounds` on `diffused_value`, `column_humidity` | 10 | 16 |
   | + `@propagate_inbounds` on `get_step`, `get_prognostic_step`, `get_tendency_step` | 0 | 0 |

   The six checks of the fused kernel that are not in the scheme belong to other parameterizations and go away
   with the step functions as well. On the H100 the step functions take another 1.4–1.6 % off the fused column
   kernel (1.289 → 1.271 ms at T256 L24, 4.510 → 4.436 ms at T512 L24, 8.17.3).

Patches, each applied to a copy of `007f04d9`: `$A/inbounds_fix.patch` (items 1 and 2, copy `$A/tree_fix`) and
`$A/get_step_fix.patch` on top (item 3, 18 methods, copy `$A/tree_fix2`). End to end on the CPU
(`cpu_baseline.jl`, the copies measured one after the other):

| | `007f04d9` | items 1 and 2 | items 1 to 3 |
|---|---|---|---|
| T32 L8 time step [ms], four rounds (items 1 to 3: the last two) | 13.87 / 13.75 / 13.74 / 13.90 | 12.95 / 13.05 / 12.96 / 13.04 (−5.9 %) | 12.87 / 12.93 (−6.7 %) |
| T32 L8 all column parameterizations [ms] | 9.13 / 9.00 / 9.01 / 9.12 | 8.18 / 8.18 / 8.19 / 8.18 (−9.7 %) | 8.08 / 8.08 (−10.9 %) |
| T32 L8 scheme [ms] | 3.14 / 3.05 / 3.03 / 3.04 | 2.60 / 2.60 / 2.60 / 2.61 | 2.56 / 2.56 |
| T32 L24 time step / all column parameterizations [ms] | 37.85 / 20.92 | | 35.18 / 18.63 (−7.1 / −10.9 %) |
| T64 L8 time step / all column parameterizations [ms] | 52.63 / 31.45 | | 48.93 / 27.78 (−7.0 / −11.7 %) |
| T128 L8 time step / all column parameterizations [ms] | 226.0 / 116.4 | | 212.7 / 102.7 (−5.9 / −11.8 %) |

- The three fixes together make a CPU time step 6–7 % faster at all four configurations, all of it in the column
  parameterizations (−11 to −12 %). Item 3 contributes 0.09 ms at T32 L8 (0.7 % of the time step).
  Logs `$A/cpu_baseline_tree*`, `cpu_baseline_fix*` (chain `$A/run_chain_6a.sh` and a8's `run_chain.sh`).
- Unit tests on `tree_fix2` (items 1 to 3) with the `@inbounds` active, i.e. without `--check-bounds=yes`
  (`$A/run_tests.jl`, environment `SpeedyWeather/Project.toml`, `$A/tests_fix2_default.log`): **54066 passed,
  0 failed, 1 error** in 19 files of `test/parameterizations` (12069 tests, among them the 71 of the scheme) and
  13 files of `test/dynamics` (41997 tests). The error is `dynamics/copy_variables.jl` loading `JLArrays`, a
  test-only dependency that this environment does not have (its 27 tests before that line pass); it has nothing
  to do with the patch. Not included: the test files that need other test-only dependencies (`dispatch.jl`,
  `numericalradiation.jl`, `terrarium_coupling.jl`) and the output tests.
- With `--check-bounds=yes`, as in CI, the three fixes are without effect by construction. The same files give
  the same result there: **54066 passed, 0 failed, 1 error** (again the `JLArrays` load), 142 min,
  `$A/tests_fix2_checkbounds.log`. `parameterizations/surface_fluxes.jl` counts 0 tests in both runs: the file
  has two test sets and no `@test`, so it only checks that the prescribed-flux models run.
- A first run of session a8 on `tree_fix` was killed when that session ended.
- Not run: the Enzyme tests. Item 1 changes the loop that the `SPLIT_COLUMN_LOOPS` workaround wraps for
  llvm/llvm-project#219885, so the extended differentiability tests have to be run before it is merged. Also not
  run: the Reactant and GPU suites with the patches.

#### 8.17.6 Candidates and recommendation

"Scheme" percentages are relative to the scheme's own cost, which is 21 % of a T32 L8 time step on the CPU and
1.4–2.8 % of a T256/T512 time step on the H100.

| # | Candidate | CPU | H100 | Reactant | Result changes by | Size of the change |
|---|---|---|---|---|---|---|
| 1 | `@inbounds` on the grid point loop of `column_parameterization_cpu!` | with 2: all column parameterizations −10 to −12 %, time step −6 to −7 % (four configurations) | none (the kernel has it) | not affected | nothing | one line |
| 2 | `@propagate_inbounds` instead of `@inline` on `diffused_value`, `column_humidity` | nothing by itself, but required for 9 | scheme −5 % (solves −20 %) | **removes the cause of the failures of 8.3 and 8.5** | nothing | three lines |
| 3 | `@propagate_inbounds` on `get_step`, `get_prognostic_step`, `get_tendency_step` | time step −0.7 % at T32 L8 | fused column kernel of all parameterizations −1.5 %; with 2 no bounds check is left in it (32 before) | fewer checks to raise in every parameterization | nothing | 18 one-line methods in `steps.jl` |
| 4 | With 2: remove the virtual temperature scratch loop (the workaround of 8.3) | not measured | scheme −5 % | raised correctly, wet and dry | nothing (identical on the H100) | a few lines, removes a `TODO` |
| 5 | Solves with shared pivots: u with v, T with q, right-hand side in the tendency | solves −38 to −41 %, scheme −13 % | scheme −22 % | raised correctly | 10⁻⁷ | new solver for N variables; the uncoupled solve needs one more scratch array or the present path |
| 6 | Reciprocals in the conductances, boundary layer height without a division per layer | conductances −1 to −4 % | with 5: scheme −40 to −45 % | raised correctly | 10⁻⁵ (rounding) | rewrite of two loops, same formulas |
| 7 | `log` of the stability correction only in the surface layer (`if`) | conductances −8 % | none | raised correctly | nothing | a few lines |
| 8 | Surface Richardson number and friction velocity stored by the boundary layer scheme instead of recomputed | conductances −11 % at L8, −4 % at L24 | none | raised correctly | nothing | one or two new 2D variables, couples the two schemes |
| 9 | Range form: the scheme as a function of a range of columns, layers in the outer loops; all columns at once on the CPU (`@simd`), one column in the kernel. Includes 5–7 | scheme −50 % with Base `exp`/`log`; −89 to −94 % with branch-free `exp`/`log` and a layer-outer height search | same as 5 + 6 (−40 to −45 %); the two CPU options must stay off | raised correctly, also with the CPU options | 10⁻⁵ | rewrite of the scheme, three 2D scratch arrays, two compile-time switches for the CPU |
| – | One sweep for all four variables | solves −44 to −47 % | slower than 5 | raised correctly | 10⁻⁷ | not recommended: 9 makes it irrelevant on the CPU |
| – | The two `exp` per interface | bound: conductances −25 % | bound: small | | changes the physics | not tried: needs profiles from reference heights or one height scale |

What it adds up to, if all of 1–9 were done: on the CPU at T32 L8 the time step goes from 13.9 ms to about
10.6 ms (−23 %; −16 % with Base `exp`/`log`), of which the bounds checks of 1–3 are 0.9 ms for all
parameterizations and the scheme 2.3 ms. On the H100 the time step gets 0.6–1.3 % faster.

Recommendation (shared by both sessions; nothing of it is implemented, the decision is the user's):

1. **Bounds checks first, as a change of their own** (1–3, optionally 4): small, no change of results, a gain for
   every column parameterization on the CPU, and it removes the cause of the Reactant failures instead of working
   around it. Before merging: the Enzyme tests including the extended ones, because 1 changes the loop that the
   `SPLIT_COLUMN_LOOPS` workaround wraps; the Reactant and GPU suites.
2. **Then the scheme in the range form (9)** as a separate step with its own plan section: it is the only
   candidate that changes the CPU picture (the scheme from 21 % of a T32 L8 time step to 3 %, or 11 % with
   Base's `exp`/`log`), it costs nothing on the GPU, and it is one code path for CPU, GPU and Reactant. To
   decide there: whether the CPU gets the branch-free `exp`/`log` (a8's prototype is Float32 only and its
   accuracy is not characterised), how the uncoupled solve and the options `diffuse_*` are handled, and Enzyme
   through the `@simd` loops.
3. If the range form is not wanted: 5 to 8 in the present per-column structure give −28 % on the CPU (about
   −21 % without 8) and −40 to −45 % on the GPU.
4. 8 is independent of the rest and only matters once the scheme is cheap (then it is 17–38 % of what is left on
   the CPU).

Outside the scheme, seen on the way and not followed up: `boundary_layer` costs 1.7 ms per T32 L8 step on the
CPU for a 2D computation (12 % of the time step), and launching a kernel with `vars` and `model` as arguments
costs 0.11 ms on the H100 node and 0.15 ms on login02 before any work is done.

### 8.12 Next

0. Fill the item marked PENDING in `PR_bulk_richardson_diffusion.md` (multi-year runs) and set the PR number in
   `CHANGELOG.md`.
1. Commit the changelog (PR number) and this document.
2. Follow the four multi-year runs (8.10) to the end; continue preempted ones with `START_RUN_NUMBER`. Decide on
   the T512L24 setup from them.
3. **To-do: is `HyperDiffusion(time_scale = Hour(2))` really needed?** Evidence so far: it removes the Antarctic
   coastal cold event in two restart tests (8.10), the event occurred once each in two of three multi-year runs
   and never ended a run. To decide: compare the multi-year run with it (`drag1e-5_hyper2`) against the one
   without (`drag1e-5`) over the full three years (events, drift, highest winds), and look at what the stronger
   hyperdiffusion costs in resolved small scales (e.g. kinetic energy spectrum), which has not been looked at.
4. Tibet at L24: done at T32 (8.14).
5. With `critical_Richardson = 10` (8.16): repeat the stability checks that were made with Ri_c = 1, at least the
   T32L24 ensemble (8.13.2) and a T512L24 run, and look at the maximum winds and the jet with the deep
   mid-latitude boundary layer. Compare the `blh` output with the ERA5 boundary layer height climatology.
6. Performance (8.17.6), to be decided by the user: (a) the three bounds-check fixes as a change of their own
   (patches `inbounds_fix.patch` and `get_step_fix.patch` in the audit directory; the Enzyme tests including the
   extended ones and the Reactant and GPU suites have to be run with them first); (b) the scheme in the range
   form, with a plan section of its own.

Slurm: `gpushort` allows two running jobs per user and 24 hours; the multi-year runs other than the first use
`gpupreempt`, which has no such limit but can be preempted.
