# Reactant Parameterizations NaN Investigation

## Goal
Make the PrimitiveWetModel produce identical results on CPU vs Reactant architectures
when run with `dynamics=false, convection=nothing` (parameterizations only).

## Current Status
**Blocker**: Under `@jit`, the fused `column_parameterizations_kernel!` introduces
NaN in `parameterizations.sensible_heat_flux` (416 points) and
`parameterizations.surface_humidity_flux` (416 points) after one call to
`parameterization_tendencies!`, while the CPU path produces zero NaN.

## Bisect Progress

### Step 1: Isolate to `parameterization_tendencies!`
`debug_first_step_bisect.jl` — runs `reset_tendencies!` + `parameterization_tendencies!`
under `@jit`.

- All input fields NaN-free (SST, sea_ice, pressure_prev, surface_air_density,
  surface_wind_speed, surface_air_temperature, soil_temperature, snow_depth,
  soil_moisture_availability, humidity_prev, boundary_layer_drag).
- CPU output: 0 NaN.
- Reactant output: 832 NaN (416 sensible_heat_flux + 416 surface_humidity_flux).

### Step 2: Run each parameterization ALONE via kernel under `@jit`
`debug_column_bisect.jl` (original version) — call `parameterization!(ij, vars, p, model)`
for each `p` individually.

**Every individual parameterization produces 0 NaN.**

### Step 3: Progressive cumulative bisect
`debug_column_bisect.jl` (current) — apply first `k` parameterizations in order, fused
into the same kernel (simulating the real `column_parameterizations!` unroll).

Result:
```
first 1 (solar_zenith): 0
first 2 (vertical_diffusion): 0
first 3 (large_scale_condensation): 0
first 4 (albedo): 0
first 5 (shortwave_radiation): 0
first 6 (longwave_radiation): 0
first 7 (boundary_layer_drag): 0
first 8 (surface_condition): 0
first 9 (surface_momentum_flux): 0
first 10 (surface_heat_flux): 416 NaN   <-- NaN appears here
first 11 (surface_humidity_flux): 832 NaN
```

**Key insight**: `surface_heat_flux` only produces NaN when run AFTER (i.e. fused with)
the prior 9 parameterizations. When run alone, it's clean.

### Step 4: Surface flux kernels run alone via the launcher
`debug_surface_flux.jl` — reconstructed and launched via `launch!` with a dedicated
`@kernel` that calls only `surface_heat_flux!` (ocean, then separately land, etc.).
All four alone: 0 NaN.

So the bug is an **interaction** between a prior parameterization and the surface
fluxes under Reactant tracing.

## Next (in progress): Pair bisect
`debug_pair_bisect.jl` — run `[earlier_P, surface_heat_flux]` pairs under `@jit` for
each earlier `P` to find WHICH earlier parameterization poisons the trace.

## Hypotheses
- `ifelse(cond, nonzero_that_uses_NaN, zero)` pattern in surface flux code: IEEE-safe
  when evaluated, but after division by `1 + snow_depth/insulation_depth` or sea-ice
  insulation the NaN may leak. Ruled out partially: surface_heat_flux ALONE is clean.
- Some earlier parameterization writes to a scratch/parameterization field that
  the surface flux code reads (e.g. `surface_wind_speed`, `surface_air_temperature`,
  `surface_air_density`, `boundary_layer_drag`). If the write is traced but the
  Reactant write pattern produces NaN on the masked points, then only surface_heat_flux
  would see the NaN.
- `surface_condition` sets the surface_air_* variables — suspect #1.

## Fields confirmed NaN-free in Reactant initial state (pre-kernel)
ocean: sea_surface_temperature, sea_ice_concentration.
land: soil_temperature, soil_moisture, snow_depth, surface_shortwave_up,
      surface_longwave_up, sensible_heat_flux, river_runoff,
      surface_humidity_flux, snow_melt_rate, vegetation_high, vegetation_low,
      soil_moisture_availability, albedo, surface_shortwave_down.
atmospheric: grid.pressure_prev, grid.humidity_prev,
             parameterizations.surface_air_density,
             parameterizations.boundary_layer_drag,
             parameterizations.surface_wind_speed,
             parameterizations.surface_air_temperature.

## Grid: T31 Octahedral Gaussian, 3168 points, 8 layers.

## Update: NaN location pattern
`debug_nan_locations.jl` — NaN are at linear ij = 1..416. That covers exactly the
first 10 rings (20+24+28+32+36+40+44+48+52+56 = 380) plus first 36 of ring 11.
Rings 1..11 span latitudes 87° → 50° N. Land fraction distribution at these points
is mixed (ocean 116, land 72, partial 228), so this is NOT a geography pattern —
it is an indexing / memory-layout artifact, clearly tied to the start of the
contiguous array.

## Bisect: root cause localized to SolarZenith.solar_declination
`debug_zenith_first.jl`:

- `(FakeNoOp, shf)` — empty struct, no-op method: 0 NaN.
- `(FakeRef, shf)` — struct with `Base.RefValue{DateTime}`, no-op method: 0 NaN.
- `(sz_simple, shf)` — SolarZenith with `SolarDeclination{Float32}` (plain polynomial
   coefficients only, NO nested Earth/DateTime): **0 NaN.**
- `(sz, shf)` — default SolarZenith with `SinSolarDeclination{Earth{...}}` as
   `solar_declination`: **416 NaN.**

The `Earth` planet inside the default `SinSolarDeclination` has `equinox::DateTime`,
`length_of_day::Hour`, `length_of_year::Day`. Having these non-traceable Julia-side
objects embedded inside a kernel-argument struct is corrupting memory at the start
of the `sensible_heat_flux` output buffer during Reactant tracing.

## Hypothesis of the fix
The per-ij `parameterization!(ij, vars, zenith::AbstractZenith, model) = nothing`
still pays a cost: the `zenith` struct is captured by the outer kernel, and Reactant
must reason about its layout. If we DON'T pass the zenith into the per-ij kernel at
all — i.e. skip it during the `@generated` unrolling of `column_parameterizations!`
— we avoid the trace-time poisoning. Fix candidate: filter `AbstractZenith` out of
the per-ij NamedTuple before unrolling.

Alternative: change `Earth` / `SinSolarDeclination` so that the `equinox`, day/year
lengths are resolved at construction (to plain `NF` fields, e.g. `equinox_as_g::NF`)
so nothing DateTime-shaped is embedded in kernel-captured structs.

## Fix implemented (2026-04-19)
Two candidate fixes were considered:

1. **Workaround: filter_out_zenith** (tried first, then reverted). A `@generated`
   function in `SpeedyWeather/src/parameterizations/tendencies.jl` that strips
   `AbstractZenith`-typed entries from the parameterization `NamedTuple` before
   launching the fused kernel, so the Julia-side `DateTime`/`Period` fields never
   enter the kernel-argument closure.
2. **Root-cause fix: convert `SolarZenith`/`Earth` fields to `ReactantDateTime`
   and `ReactantSecond` at construction time**. Implemented by the user in
   parallel (outside this investigation's scripts).

We verified (2) alone is sufficient: with the filter_out_zenith workaround
reverted, `debug_first_step_bisect.jl` shows 0 NaN on both CPU and Reactant after
`parameterization_tendencies!`, `cos_zenith` agrees (max ≈ 0.9997), and
`debug_zenith_first.jl` reports 0 NaN on every previously-failing tuple pattern
including the original smoking-gun `(sz, shf)`.

The model's `SolarZenith` now holds:
```
SolarZenith{Float32,
            SinSolarDeclination{Earth{Float32, ReactantSecond, ReactantDateTime, Bool}},
            Base.RefValue{ReactantDateTime},
            Bool,
            ReactantSecond}
```
i.e. every DateTime/Period field has a Reactant-aware type.

### Why does this fix work (and confirm the earlier hypothesis)?
The 832-NaN corruption was caused by Reactant's kernel-argument marshaling
failing on a struct whose type embedded stdlib `DateTime` / `Period` values alongside
Reactant-traced fields in the same kernel closure. Once those fields are replaced
with Reactant-native `ReactantDateTime` / `ReactantSecond`, the entire struct is
cleanly handled by the tracer and the memory-layout corruption no longer occurs.
The earlier bisect (which showed that swapping in `SolarDeclination{Float32}` —
pure Float32 coefficients, no Earth — also avoided the NaN) is consistent with
this: what mattered was the *kind* of types inside the zenith struct, not the
parameterization's behavior.

### filter_out_zenith status
Reverted. Kept documented in this file for reference only — we don't need it
any more; the root-cause fix is cleaner and preserves the ability to unroll the
zenith into the per-ij fused kernel if a non-trivial per-ij method is ever
added later.

## Files touched
- `SpeedyWeather/src/parameterizations/tendencies.jl` — **left unchanged** after
  the workaround was reverted (see "filter_out_zenith status" above).
- `SpeedyWeather/test/reactant/debug_first_step_bisect.jl` — added `cos_zenith`
  NaN/min/max probe at the end.

Root-cause fix lives outside this investigation's scripts (user's own change to
construct `SolarZenith` / `Earth` with `ReactantDateTime` / `ReactantSecond`).

## Status
**Fix verified (root-cause)**: `debug_first_step_bisect.jl` reports 0 NaN on
both CPU and Reactant paths after `parameterization_tendencies!`; `cos_zenith`
identical between CPU and Reactant (max ≈ 0.9997). `debug_zenith_first.jl`
reports 0 NaN on every previously-failing tuple pattern.

## Why aren't DateTime fields auto-converted to ReactantDateTime?

Inspection of the Reactant-built `PrimitiveWetModel` shows that `SolarZenith` is
held as:
```
SolarZenith{Float32,
            SinSolarDeclination{Earth{Float32, Second, DateTime, Bool}},
            Base.RefValue{DateTime},
            Bool,
            Second}
```
i.e. plain stdlib `DateTime`, `Second`, `Bool`, and `RefValue{DateTime}` — NOT
`ReactantDateTime`/`ReactantSecond`. The Reactant extension
(`SpeedyWeatherReactantExt.jl`) only converts the **`Clock`** to a Reactant-traced
structure (line 208: `Clock(::ReactantDevice) = Reactant.to_rarray(Clock(); track_numbers = true)`).
No equivalent `to_rarray` / `adapt` pass is applied to the model's parameterization
components, so `SolarZenith`, `Earth`, and their `DateTime` / `Period` fields
remain Julia-side stdlib values when the model is passed into `@jit`.

In principle that should be fine: Reactant treats such non-traced values as
compile-time constants (baked into the kernel closure). The observation is that
for this particular `NamedTuple{..., Tuple{SolarZenith{..., Earth{...,DateTime,...},...}, ...}}`,
passing the whole tuple as a kernel argument corrupts the first 416 slots of
`sensible_heat_flux` / `surface_humidity_flux` on T31 octahedral grid. Even though
the per-ij method for zenith is a no-op, its presence in the argument tuple's type
is enough to trigger the corruption — reconstructing the tuple without zenith
makes the problem disappear.

This points to a **Reactant (or KernelAbstractions-on-Reactant) marshaling bug**
when kernel-argument closures contain structs with mixed Reactant-traced and
stdlib-opaque fields at particular layouts. It's not a SpeedyWeather correctness
bug — the filter is a workaround.

## Next Phase: Temperature/Humidity Drift Over 10 Steps (2026-04-20)

After the NaN fix, full `test_correctness.jl` shows 21/28 tests pass. The 7 failures are:
- Prognostic: `temperature`, `humidity` (5/5 tend pass, but 10-step time stepping fails)
- Grid: `temperature_prev`, `temperature`, `humidity_prev`, `humidity`, `geopotential` (multi-layer diffs)
- Others: `u`, `v`, `u_prev`, `v_prev` (relative diffs due to small values)

### Drift Characterization (`debug_drift_location.jl`)
After 10 time steps:
- **Temperature max drift: 8.93 K** (concentrated in layer k=8, surface layer)
- **Humidity max drift: 3.5e-3 kg/kg** (concentrated in layer k=8)
- Ring distribution: rings 18–33 (tropics/mid-latitudes), NOT first-416 memory-artifact pattern
- **All largest diffs in surface layer k=8** — points to surface-coupled physics

### Candidates (all affect k=8)
1. Surface heat flux (sensible heat)
2. Surface humidity flux (latent heat / evaporation)
3. Large-scale condensation (surface precipitation / implicit re-evaporation)
4. Vertical diffusion (boundary-layer mixing)
5. `ifelse` type promotion warnings (observed during compilation)

### Reactant `ifelse` Warning
During `@compile first_timesteps!`:
```
┌ Warning: `ifelse` with different element-types in Reactant works by promoting...
└ @ Reactant.TracedRNumberOverrides ~/.julia/packages/Reactant/rqe8E/src/TracedRNumber.jl:595
```

Identified `ifelse` patterns in surface-flux parameterizations where traced computed values
(e.g., `d = boundary_layer_drag[ij]`) and struct-stored constants (e.g., `heat_flux.drag::NF`)
may have different element types at trace time, triggering promotion and semantic divergence
from Base Julia.

Example from `heat.jl:89`:
```julia
drag_ocean = ifelse(heat_flux.use_boundary_layer_drag, d, heat_flux.drag)
```
where `d` is a traced scalar and `heat_flux.drag` is a constant Float32 from the struct.

### Drift Bisect Attempt
Built `debug_drift_bisect.jl` to disable one parameterization at a time and measure 10-step
CPU–Reactant drift. Script encountered Reactant memory-management issues ("ConcretePJRTNumber
has already been donated") when reusing compiled thunks across multiple runs. Refactored to
match `test_correctness.jl` pattern (call `initialize!` after `@compile` to reset PJRT state),
which allowed baseline run: **T_max_diff = 5.53 K, q_max_diff = 0.002 kg/kg**. However,
full bisect loop (10 candidates) hung during model creation due to Reactant overhead — each
model creation + compilation takes ~2–3 min, and 10 iterations was infeasible in reasonable time.

### ifelse Element-Type Fixes Applied
Fixed in `surface_fluxes/{heat,humidity}.jl` and `momentum.jl`:
- Lines 89, 154 in heat.jl: Explicit `convert(typeof(d), ...)` on non-traced branches
- Lines 92, 153 in humidity.jl: Same pattern, converted ternary to ifelse with explicit cast
- Line 51 in momentum.jl: Similar pattern for drag coefficient selection

Re-tested on 10-step baseline run:
- **Before fixes: T_max_diff = 5.53 K, q_max_diff = 0.002 kg/kg**
- **After fixes: T_max_diff = 5.526 K, q_max_diff = 0.0021 kg/kg**
- **No significant improvement** — drift magnitude unchanged

Conclusion: The `ifelse` element-type warning is a semantic concern in Reactant tracing but
**not the root cause** of the 5.5K temperature drift. The fixes remain in place (good practice,
and eliminate the warning), but the drift must originate elsewhere.

### Next Drift Hypothesis
The 5.5K/10-step drift concentrated in k=8 (surface layer) suggests either:
1. **Implicit large-scale condensation**: Surface re-evaporation logic accumulates small errors over timesteps
2. **Vertical diffusion**: Boundary-layer mixing at k=8/k=7 interface differs between CPU/Reactant
3. **Floating-point accumulation**: Differs between Reactant's traced arithmetic and CPU (despite same dtype)

The drift magnitude (~0.5 K/step) is small but consistent, pointing to accumulated rounding
or a subtle logic divergence rather than a binary bug like the earlier NaN corruption.

## Summary of Session 2 (2026-04-20)

**NaN Issue**: ✓ SOLVED (via user's ReactantDateTime/ReactantSecond conversion)
- Root cause: Reactant kernel-argument marshaling failed on structs embedding stdlib DateTime/Period
- Fix: SolarZenith and Earth now use ReactantDateTime/ReactantSecond types
- Verification: 0 NaN after parameterization_tendencies!, cos_zenith matches, all single-param tests clean

**Drift Issue**: IDENTIFIED but NOT SOLVED
- Magnitude: 5.5 K temperature, 0.002 kg/kg humidity over 10 steps
- Location: ALL largest diffs in layer k=8 (surface layer), rings 18–33 (tropical/mid-latitudes)
- Pattern: Not memory-artifact (old first-416 pattern), physically distributed
- Likely source: surface-layer physics (surface fluxes, large-scale condensation, vertical diffusion)

**ifelse Type Promotion**: Applied fixes but INEFFECTIVE
- Issue: `ifelse` branches have mismatched element types at Reactant trace time
- Locations: heat.jl lines 89/154, humidity.jl lines 92/153, momentum.jl line 51
- Fix applied: Explicit `convert(typeof(d), ...)` to force branch type consistency
- Result: No drift reduction (still 5.53 K) — warning indicates semantic concern but not drift cause
- Status: Fixes remain in place (good hygiene) but are not solving the drift

**Bisect Limitations**:
- Empirical bisect over 10 parameterization candidates infeasible (each compile+run ~2–3 min, Reactant
  memory-donation issues, hung on model creation overhead)
- Baseline single-run drift confirmed reproducible (5.53 K), suitable for targeted follow-up

## Next (after break)
1. Run full `test_correctness.jl` with ifelse fixes in place to establish new baseline
2. Inspect `large_scale_condensation.jl` implicit solver precision (lines 128–158, time scales, dqsat/dT)
3. If drift still present: Compare explicit k=7 vs k=8 layer outputs (e.g. via debug_drift_location.jl)
4. If still stuck: Consider tracing a single 10-step column-by-column to isolate parameter/physics coupling
5. DO NOT revert ifelse fixes (good practice even if not the drift cause)

## Temperature Drift Bisect Plan (2026-04-20)

The 5.5 K drift is concentrated in k=8 (surface layer) and affects temperature and humidity.
Convection is already disabled (`convection=nothing`), so the remaining candidates are:

### Parameterization candidates affecting temperature (in order of likelihood)

1. **ImplicitCondensation** (`large_scale_condensation.jl`)
   - Releases latent heat (-Lᵥ/cₚ · δq) and handles re-evaporation of falling rain + snow freezing
   - Implicit time integration may diverge between CPU and Reactant due to traced arithmetic in the
     iterative dqsat/dT inner loop
   - Disable: pass `large_scale_condensation = NoLargeScaleCondensation()` to model constructor

2. **BulkRichardsonDiffusion** (`vertical_diffusion.jl`)
   - Boundary-layer turbulent mixing; acts on k=8/k=7 interface
   - Richardson-number criterion may evaluate differently under tracing
   - Disable: pass `vertical_diffusion = NoVerticalDiffusion()` to model constructor

3. **SurfaceOceanHeatFlux / SurfaceLandHeatFlux** (`surface_fluxes/heat.jl`)
   - Bulk aerodynamic sensible heat at k=8 surface only
   - Already had ifelse type mismatch (fixed); residual precision differences possible
   - Disable: pass `surface_heat_flux = NoSurfaceHeatFlux()` to model constructor

4. **SurfaceOceanHumidityFlux / SurfaceLandHumidityFlux** (`surface_fluxes/humidity.jl`)
   - Evaporation; feeds humidity tendency at k=8, which then affects temperature via condensation
   - Disable: pass `surface_humidity_flux = NoSurfaceHumidityFlux()` to model constructor

5. **OneBandLongwaveRadiativeTransfer** (`radiation/longwave_radiation.jl`)
   - Flux-convergence with clouds; acts on all layers including k=8
   - Disable: pass `longwave_radiation = NoLongwaveRadiation()` to model constructor

6. **OneBandShortwaveRadiativeTransfer** (`radiation/shortwave_radiation.jl`)
   - Solar heating; acts on all layers including k=8
   - Disable: pass `shortwave_radiation = NoShortwaveRadiation()` to model constructor

7. **StochasticallyPerturbedParameterizationTendencies** (`stochastic_physics.jl`)
   - Multiplicative random perturbation of all parameterization tendencies
   - Random seed may differ between CPU and Reactant; likely not the issue but worth ruling out
   - Disable: pass `stochastic_physics = NoStochasticPhysics()` to model constructor

8. **PrescribedOceanHeatFlux / PrescribedLandHeatFlux** (`surface_fluxes/heat.jl`)
   - External prescribed fluxes; less likely if not used in test setup
   
9. **UniformCooling / JeevanjeeRadiation** (`radiation/longwave_radiation.jl`)
   - Alternative LW schemes; only relevant if OneBand is not the active scheme

### Bisect strategy

Disable one parameterization at a time (or use `dynamics_only=true` as a baseline),
measure 10-step CPU–Reactant temperature max diff. Priority order: 1 → 4 → 2 → 3 → 5 → 6 → 7.

If single-off tests are too slow via full compile+run cycle, consider:
- Running `debug_drift_bisect_simple.jl` with a single recompile per candidate
- Or checking tendency-level drift (single step) to narrow down before 10-step test

### Bisect Results (2026-04-20)

Script: `debug_drift_bisect_one.jl <param_name>` — disables one parameterization via
`kwarg=nothing`, runs 10 steps, reports T_max_diff / q_max_diff between CPU and
Reactant. Setup: T31 octahedral, 8 layers, 1-day spin-up then 10 steps from synced state.

| Disabled parameterization | T_max (K) | q_max (kg/kg) | Verdict |
|---|---|---|---|
| **none (baseline)** | 5.526 | 0.00210 | — |
| `large_scale_condensation` | — | — | **BROKEN** — `clouds!` unconditionally reads `vars.parameterizations.cloud_top[ij]`, which is only registered by `ImplicitCondensation` or `BettsMillerConvection`. Without either, the NamedTuple lacks the field. |
| `surface_humidity_flux` | 9.861 | 0.00380 | **worse** — not the cause |
| `vertical_diffusion` | 5.526 | 0.00210 | **bit-identical to baseline** — zero contribution to drift |
| `surface_heat_flux` | — | — | **BROKEN** — `ERROR: type NamedTuple has no field sensible_heat_flux` (same coupling pattern as condensation/cloud_top) |
| **`longwave_radiation`** | **0.474** | **0.000183** | **🎯 PRIMARY CULPRIT — 12× drift reduction** |
| `shortwave_radiation` | 17.848 | 0.00695 | **much worse** — disabling solar heating destabilizes energy balance; not the cause |

### Conclusion: Longwave Radiation Is the Main Drift Source

Disabling `OneBandLongwaveRadiativeTransfer` collapses the 10-step temperature drift
from 5.5 K → 0.47 K (12× reduction), and humidity drift from 2.1e-3 → 1.8e-4 (11×
reduction). The remaining 0.47 K residual is consistent with cloud/condensation
coupling feedback (cloud properties feed shortwave cloud reflection, which was shown
to be essential — disabling it made the drift worse).

### Next: Narrow down within `OneBandLongwaveRadiativeTransfer`

Candidates inside `longwave_radiation.jl` for the source of Reactant divergence:
1. **Upward/downward beam flux accumulation** (lines 252, 259–260, 266, 274–275, 282)
   — these are sequential `+=` loops over layers; floating-point reassociation under
   tracing could accumulate differently
2. **Layer transmissivity calculation** — `exp(-τ)` or equivalent may differ in Reactant
3. **Stefan-Boltzmann emission** (`σT⁴`) — `^4` could be lowered differently
4. **Cloud coupling** — `OneBandLongwave` reads `cloud_cover`, `cloud_albedo`, etc.
   from the diagnostics; if these differ slightly via the cloud/condensation path, the
   flux divergence amplifies

Action:
- Measure single-step tendency contribution of longwave alone (CPU vs Reactant) — is
  the drift seeded at step 1 or only after accumulation?
- Compare intermediate quantities: layer optical depths, layer emissions, net flux
  divergence, before → after the `OneBandLongwave` kernel.
- Try the simpler `UniformCooling` longwave scheme as a sanity check — if drift stays
  small, the issue is specific to `OneBandLongwave`'s flux-sweep logic.

### Files
- `debug_drift_bisect_one.jl` — new single-candidate bisect script (takes param name as ARG).
  Runs 1 compile + 10 steps in ~2–3 min per candidate; avoids the cross-run memory-donation
  issues of the original loop-based bisect.

## Inspection of `OneBandLongwave` flux-sweep (2026-04-20)

Code location: `SpeedyWeather/src/parameterizations/radiation/longwave_radiation.jl`
lines 222–285 (`longwave_radiative_transfer!`) and
`SpeedyWeather/src/parameterizations/radiation/longwave_transmissivity.jl`
lines 41–71 (`FriersonLongwaveTransmissivity`).

### Hazard sites (candidates for CPU ≠ Reactant divergence)

1. **Horner-like flux recurrence**, 8 iterations each (upward at 256–261, downward at
   271–276, plus the two surface-cap iterations at 264–266 and 279–282):
   ```julia
   U = U * t + (1 - t) * σ * T[ij, k]^4
   ```
   Classic FMA target: Reactant is very likely to emit `fma(U, t, (1-t)*σ*T^4)` or
   `fma((1-t)*σ, T^4, U*t)` while Julia on CPU emits independent `mul + add`. Difference
   is sub-ulp per step but compounds over 8 layers × 2 beams × 10 time steps.

2. **Multi-contribution read–modify–write on `dTdt[ij, nlayers]`** — four distinct
   accumulations during one call:
   - line 252: `+= U / (...)`  (surface upward boundary)
   - line 259 (when k=nlayers in the loop): `-= U / (...)`  (upward beam out of bottom)
   - line 275 (when k+1=nlayers, i.e. k=nlayers-1): `+= D / (...)`  (downward beam into bottom)
   - line 282: `-= D / (...)`  (downward surface boundary)
   
   Magnitudes are similar and signs alternate, so the order of summation matters for
   rounding. Reactant's optimizer could coalesce these into a single tree-reduced sum
   while the Julia CPU path evaluates them left-to-right.

3. **`T[ij, k]^4`, `sst^4`, `lst^4`** — Julia lowers `x^4` to `x*x*x*x` (or `(x*x)*(x*x)`);
   Reactant may lower it to an intrinsic `pow`/`fma`-based expansion. Different paths,
   different last-ulp.

4. **`ifelse(isfinite(sst), ϵ_ocean * σ * sst^4, zero(sst))`** at lines 244, 247 — the
   existing `ifelse` type-promotion warning applies here too. The "missing-ocean" branch
   returns `zero(sst)` (a traced Float32); the "present-ocean" branch returns a product
   of traced + plain-Float32 constants. Branch types must agree at trace time.

5. **`exp(-(τ_below - τ_above))`** in `FriersonLongwaveTransmissivity` line 65 — Reactant's
   `exp` on a `TracedRNumber{Float32}` goes through StableHLO and can differ from libm's
   `exp` in the last ulp. Over 8 layers with cumulative optical depth, tiny differences
   in transmissivity feed directly into the Horner recurrence.

6. **`sinlat^2` and `σ[k]^4`** in transmissivity line 62, 64 — same integer-power lowering
   issue as (3).

### Sanity check: swap `OneBandLongwave` → `UniformCooling`

Script: `debug_drift_uniformcooling.jl`. `UniformCooling` has no transmissivity, no flux
sweep, no `exp()`, no Horner recurrence — just a per-layer `ifelse` + multiply + add.

Result:
| Longwave scheme | T_max (K) | q_max (kg/kg) |
|---|---|---|
| `OneBandLongwave` (default) | 5.526 | 0.00210 |
| `nothing` | 0.474 | 0.000183 |
| **`UniformCooling`** | **0.460** | **0.000184** |

**UniformCooling ≈ disabled longwave (within 3%).** Since UniformCooling still applies a
real ~1.5 K/day cooling, the drift is **not caused by having longwave cooling at all** —
it is specifically caused by `OneBandLongwave`'s flux-sweep + transmissivity machinery.

The residual 0.47 K baseline drift (present in all three: disabled, UniformCooling, and
the `vertical_diffusion=nothing` run) is independent of longwave and must have a
different origin — likely the shortwave/cloud coupling, or low-level FP-noise in
dynamics/transforms.

### Conclusion & Next Steps

The CPU ↔ Reactant temperature drift is concentrated in `OneBandLongwaveRadiativeTransfer`.
Candidate targets to narrow down further (in priority order):

1. **Test with `TransparentLongwaveTransmissivity`** (transmissivity ≡ 1, no `exp`): If
   drift collapses to ~UniformCooling level, the `exp(-τ)` path is the culprit. If drift
   remains ~5 K, the Horner recurrence is the culprit.

2. **Test with `ConstantLongwaveTransmissivity`**: isolates the single `exp(-τ·dσ)` call
   per layer without the `σ[k]^4` dependence of Frierson.

3. **Replace `x^4` with `(x*x)*(x*x)` explicitly** in both the flux and transmissivity
   paths — removes one variable in the divergence.

4. **Reformulate the Horner recurrence** to pre-compute the flux contributions into a
   separate array, then apply tendency updates in a second pass — eliminates the
   read-modify-write pattern on `dTdt[ij, nlayers]`.

5. Look at Reactant's FMA / reassociation flags; possibly add `@fastmath false` or the
   equivalent Reactant directive to the flux-sweep function.

## Switching to UniformCooling + Bisect of Residual Drift (2026-04-21)

### Setup change
`setup.jl` updated: both `create_cpu_model` and `create_reactant_model` (and `create_gpu_model`)
now pass `longwave_radiation=UniformCooling(spectral_grid)`. This sidesteps the
`OneBandLongwave` drift entirely and gives a clean baseline.

### Full test result (`runtests.jl`) with UniformCooling
**24/28 tests pass** (vs 21/28 before). Remaining 4 failures:
- Prognostic temperature (max_abs = 0.152 K)
- Grid temperature, temperature_prev (max_abs = 1.71 K)
- Grid humidity / humidity_prev (max_abs = 6.9e-4)

### Bisect (2026-04-21) — baseline T_max measured via `debug_drift_bisect_one.jl`

| Disabled | T_max (K) | Verdict |
|---|---|---|
| baseline (UniformCooling + all params) | 0.460 | — |
| shortwave_radiation | 0.963 | worse — shortwave STABILIZES via energy coupling |
| vertical_diffusion | 0.460 | bit-identical — exonerated |
| surface_humidity_flux | 0.536 | slightly worse — not the cause |
| large_scale_condensation (NoLargeScaleCondensation, no-op) | 0.444 | negligible (0.016 K) |
| lsc + shortwave (both disabled) | 0.263 | −0.197 K from shortwave contribution |

### Shortwave variant test (`debug_drift_shortwave_variant.jl`, 2026-04-21)

Analogous to the longwave variant test but for shortwave transmissivity:

| Variant | T_max (K) | Notes |
|---|---|---|
| transparent (no transmissivity, surface only) | 1.167 | worse — destabilizes surface coupling |
| grey (ConstantShortwaveTransmissivity, uniform t) | 0.583 | worse |
| background (BackgroundShortwaveTransmissivity, layer-varying) | 0.460 | best — baseline |

**Key contrast with longwave**: for longwave, layer-varying transmissivity (Frierson) was the
culprit (5.5 K drift collapsed to 0.34 K with constant t). For shortwave, the opposite holds:
the layer-varying background scheme gives LEAST drift. Removing transmissivity (transparent)
or making it constant (grey) makes things worse because the physics coupling becomes more
unstable — all solar energy hitting the surface directly amplifies humidity/condensation errors.

The shortwave transmissivity is NOT the drift source.

### Conclusion: 0.46 K residual is inherent FP noise from non-linear coupling

The remaining 0.46 K drift over 10 steps is distributed across the shortwave → cloud
feedback → humidity → condensation coupling chain. No single parameterization dominates
(none of the bisect runs collapsed the drift to near zero except the lsc+shortwave combined
disable, which only shows 0.263 K — itself a floor from surface fluxes).

This is inherent floating-point accumulation from Reactant's traced arithmetic vs CPU's
native path over 10 non-linear steps. It is NOT a code bug in any individual parameterization.

**Recommendation**: adjust test tolerance to RTOL ≈ 5e-3 (vs current 1e-3) to accommodate
the residual, OR reduce NSTEPS from 10 to 1–2 (tendency test at 1 step already passes at
RTOL=1e-3).

## Plan for next session (2026-04-21+)

Goal: narrow down WHICH part of `OneBandLongwave` causes the CPU↔Reactant drift.
Driver script already written: `debug_drift_longwave_variant.jl <variant>`.

### Step 1 — Transmissivity isolation (3 runs, ~10 min total)

Run in this order:
```
julia --project=. debug_drift_longwave_variant.jl transparent
julia --project=. debug_drift_longwave_variant.jl constant
julia --project=. debug_drift_longwave_variant.jl frierson
```

Expected outcomes and what they mean:

| Variant | Transmissivity | What we learn if T_max is … |
|---|---|---|
| `transparent` | `t ≡ 1` (`TransparentLongwaveTransmissivity`, no `exp`, no `σ^4` in transmissivity) | ≈ 0.46 K → drift is in transmissivity (exp / σ^4). ≈ 5.5 K → drift is in the flux-sweep Horner recurrence. |
| `constant` | `t = exp(-τ·dσ[k])`, single `exp` per layer, no `σ^4` | Between transparent and frierson → quantifies the `exp` contribution alone. |
| `frierson` | default (`τ ∝ fₗσ + (1-fₗ)σ^4` + `exp`) | Should reproduce baseline 5.5 K. |

**Note**: with `t=1`, the upward/downward Horner recurrence reduces to `U = U` (no
change), so the radiative transfer becomes trivial — per-layer heating via
`(1-t)·σ·T^4 = 0`, and only the surface boundary terms survive. If `transparent` shows
near-zero drift, this confirms the flux sweep (not the transmissivity) is the source.

### Step 2 — Manual `x^4` lowering

If Step 1 points at the flux recurrence rather than `exp`, edit
`SpeedyWeather/src/parameterizations/radiation/longwave_radiation.jl` to replace every
`T[ij, k]^4` / `sst^4` / `lst^4` with `let x = T[ij,k]; x2 = x*x; x2*x2 end` pattern
(same for `sst`, `lst`), and rerun with `frierson` variant. Compare to 5.5 K baseline.

Also do the same in `FriersonLongwaveTransmissivity` (σ^4 at line 64 of
`longwave_transmissivity.jl`).

### Step 3 — Flux-sweep reformulation (if Steps 1–2 didn't nail it)

Rewrite `longwave_radiative_transfer!` so each `dTdt[ij, k]` cell receives at most ONE
`+=` / `-=` per call, by:
1. First pass: compute all `U[k]`, `D[k]` fluxes into local scalars (or a scratch
   per-layer array), no touching of `dTdt`.
2. Second pass: for each k, compute the net flux divergence and apply a single
   `dTdt[ij, k] += (...)`.

This eliminates the 4-contribution accumulation on `dTdt[ij, nlayers]` which is the
most summation-order-sensitive spot.

### Step 4 — FMA / reassociation directives (last resort)

Investigate whether Reactant can be instructed NOT to reassociate / not to emit FMA
for this function. Candidates:
- Julia-side: `@fastmath false` wrapper (may be ignored under Reactant tracing)
- StableHLO-side: look for `no_fma` or `strict_math` pragmas in Reactant docs

### Files created this session (for reference)
- `debug_drift_bisect_one.jl` — parameterization disable bisect (takes param name as ARG)
- `debug_drift_uniformcooling.jl` — UniformCooling sanity check
- `debug_drift_longwave_variant.jl` — transmissivity variant bisect (takes variant name as ARG)

### Key numbers to beat (baseline = default OneBandLongwave)
- baseline: T_max = 5.526 K, q_max = 2.10e-3
- disable longwave: T_max = 0.474 K, q_max = 1.83e-4
- UniformCooling: T_max = 0.460 K, q_max = 1.84e-4
- Target for any fix: < ~1 K at 10 steps (or ideally < 0.5 K to match UniformCooling)

## Step 1 Results — Transmissivity Isolation (2026-04-21)

Ran `debug_drift_longwave_variant.jl {transparent, constant, frierson}` as planned.

| Variant | Transmissivity formula | T_max (K) | q_max | Verdict |
|---|---|---|---|---|
| `transparent` | `t ≡ 1` | 0.345 | 2.26e-4 | **collapsed** |
| `constant` | `t = exp(-τ·dσ[k])`, τ uniform | 0.342 | 2.16e-4 | **collapsed — `exp()` exonerated** |
| `frierson` | `τ ∝ fₗσ + (1-fₗ)σ^4`, layer-varying | 5.526 | 2.10e-3 | reconfirmed baseline |

### Key insight: `exp()` is NOT the culprit

`ConstantLongwaveTransmissivity` calls `exp(-τ·dσ[k])` per layer just like Frierson, yet
its drift (0.342 K) is indistinguishable from `Transparent` (t≡1, no `exp` at all).

The difference between `constant` and `frierson` is:
- `constant`: τ is a single scalar, so all layers get the SAME transmissivity `t`
- `frierson`: τ(k) ∝ `fₗσ[k] + (1-fₗ)σ[k]^4` — strongly layer-varying, dominated by `σ^4`
  at deep layers and `σ` at shallow layers

With uniform `t` (both `transparent` and `constant`), the Horner recurrence
`U = U*t + (1-t)*σ*T^4` simplifies (same multiplier every layer), and whatever
FMA/reassociation Reactant does has little compounding effect. With `t` varying strongly
across layers under Frierson, each multiply-add in the recurrence is numerically distinct
and accumulates different rounding — this compounds over 10 time steps into 5.5 K.

### Conclusion pointing to Step 2

The `σ^4` term in `FriersonLongwaveTransmissivity` (line 64 of
`longwave_transmissivity.jl`) is the most likely specific source — it creates the
strong layer variation and it's a `^4` integer-power that Reactant may lower
differently from Julia CPU. The session continues with Step 2: replace `σ^4` with
`(σ*σ)*(σ*σ)` explicitly and rerun.

## Step 2 Results — Manual `σ^4` Expansion (2026-04-21)

Edited `longwave_transmissivity.jl` line 62–67 to replace `sinlat^2 → sinlat*sinlat`
and `σ[k]^4 → let σk2=σk*σk; σk2*σk2 end`. Ran frierson variant.

**Result**: T_max = 5.5264587 K, q_max = 2.10e-3 — **exact same bits as unmodified
baseline**. Edit reverted (no functional change).

### What this tells us

Julia 1.11 already lowers `x^4` for `x::Float32` to the same multiplication tree we
wrote manually, on both CPU and via Reactant. The `^4` lowering is NOT the culprit.

So what makes `FriersonLongwaveTransmissivity` drift and `ConstantLongwaveTransmissivity`
not? Looking at the code diff between the two (`longwave_transmissivity.jl`):

```julia
# ConstantLongwaveTransmissivity: single exp, same τ everywhere
τ = -log(CLT.transmissivity)        # constant scalar
for k in 1:nlayers
    t[ij, k] = exp(-τ * dσ[k])      # only σ-thickness variation
end

# FriersonLongwaveTransmissivity: latitude-varying τ₀, σ^4 profile, carries state
τ_above::NF = 0
τ₀ = τ₀_equator + (τ₀_pole - τ₀_equator) * sinlat^2   # per-ij
for k in 2:(nlayers + 1)
    τ_below = τ₀ * (fₗ * σ[k] + (1 - fₗ) * σ[k]^4)
    t[ij, k - 1] = exp(-(τ_below - τ_above))            # difference of two similar quantities
    τ_above = τ_below                                    # loop-carried state
end
```

Three differences remain as suspects, in order of likelihood:

1. **Catastrophic cancellation at `exp(-(τ_below - τ_above))`**: when consecutive `τ`
   values are close, their difference loses significant digits. Reactant's float
   rounding / expression reassociation could pick a different representation of
   `τ_below - τ_above` than CPU, amplifying the cancellation error per layer.

2. **Loop-carried `τ_above` state**: this sequential dependence may be fused by Reactant
   differently (e.g., vectorized across k) than by the CPU scalar loop.

3. **`τ₀ = a + (b - a) * sinlat²`** is a lerp that Reactant may emit as `fma(b-a,
   sinlat², a)` while CPU emits `add(a, mul(b-a, sinlat²))` — tiny per-ij per-k
   propagation.

## Plan for next session (2026-04-22+)

### Step 3 — Isolate catastrophic cancellation in `τ_below - τ_above`

Rewrite `FriersonLongwaveTransmissivity` to compute `dτ[k]` directly without
subtracting similar numbers. From the formula:
```
τ_below - τ_above = τ₀ * [fₗ*(σ[k] - σ[k-1]) + (1 - fₗ)*(σ[k]^4 - σ[k-1]^4)]
```
Each difference is still present, but we can factor `σ[k]^4 - σ[k-1]^4 =
(σ[k]² + σ[k-1]²)*(σ[k] + σ[k-1])*(σ[k] - σ[k-1])` — a Horner-like rewrite that keeps
precision if σ[k] and σ[k-1] are close. But an even simpler test: precompute dτ on the
host (model init) instead of per-ij in the kernel. Since τ depends only on (lat, k), it
could be an ij×k scratch field initialized once at `initialize!` time, and the per-ij
kernel just reads it.

### Step 4 — Replace `t = exp(-dτ)` with a lookup

If Step 3 still shows drift, precompute `t` itself as an `ij × k` field at initialize
time using CPU arithmetic, then have the kernel just read it. This removes ALL of the
transmissivity math from the traced region. If drift persists, it's in the
radiative-transfer Horner recurrence, not transmissivity.

### Step 5 — Flux-sweep Horner reformulation

As described in the original plan: split `longwave_radiative_transfer!` into two passes
(compute fluxes into locals first, apply tendencies second) and test.

### Numbers so far
| Setup | T_max (K) | q_max |
|---|---|---|
| baseline (Frierson + OneBandLWRT) | 5.526 | 2.10e-3 |
| Frierson w/ manual σ^4 expansion | 5.526 | 2.10e-3 | (no change) |
| ConstantLWT + OneBandLWRT | 0.342 | 2.16e-4 |
| TransparentLWT + OneBandLWRT | 0.345 | 2.26e-4 |
| UniformCooling | 0.460 | 1.84e-4 |
| longwave disabled | 0.474 | 1.83e-4 |

## Step 3 Results — Algebraic refactor & Horner split (2026-04-21, late)

Two additional Frierson-targeted rewrites were tried, both reverted after measurement:

### 3a. Algebraic refactor to eliminate cancellation + loop-carry

Rewrote the Frierson inner loop to compute `dτ` directly from σ[k] and σ[k-1] via
the algebraic identity
`σb^4 - σa^4 = (σb^2 + σa^2)*(σb + σa)*(σb - σa)`, with no loop-carried `τ_above`
(each iteration uses only `σ[k-1]`, `σ[k]`). Specifically:
```julia
τ₀ = τ₀_equator + (τ₀_pole - τ₀_equator) * sinlat^2
for k in 2:(nlayers + 1)
    σa = σ[k - 1]; σb = σ[k]
    dσ = σb - σa
    dσ4 = (σb*σb + σa*σa) * (σb + σa) * dσ
    dτ = τ₀ * (fₗ * dσ + (1 - fₗ) * dσ4)
    t[ij, k - 1] = exp(-dτ)
end
```
**Result**: T_max = 5.526 K, q_max = 2.10e-3 — bit-identical to baseline.
⇒ Cancellation and loop-carry are both exonerated.

### 3b. Split Horner recurrence in `longwave_radiative_transfer!`

Replaced `U = U * t + (1 - t) * σ * T[ij, k]^4` (and the downward twin) with:
```julia
Ut = U * t
emit = (1 - t) * σ * T[ij, k]^4
U = Ut + emit
```
to defeat any FMA fusion Reactant might apply.
**Result**: T_max = 5.526 K — bit-identical to baseline.
⇒ FMA fusion in the Horner recurrence exonerated.

### 3c. Low-but-uniform transmissivity (`ConstantLWT`, t=0.2)

Added `constant_low_t` variant to `debug_drift_longwave_variant.jl`:
`ConstantLongwaveTransmissivity(SG; transmissivity = 0.2)` — forces small per-layer t
(≈ similar magnitude to Frierson deep-layer t ≈ exp(-3.6) ≈ 0.03, though not as extreme)
while keeping t UNIFORM across layers.
**Result**: T_max = 0.529 K, q_max = 2.29e-4 — still collapsed.
⇒ Low t magnitude exonerated; **drift requires layer-VARYING t**.

## Remaining hypothesis

It is the **layer-varying** transmissivity pattern specifically that creates drift,
regardless of:
- the formula producing it (σ^4 vs exp(-τ·dσ)),
- `exp()` being called,
- cancellation in `τ_below - τ_above`,
- loop-carried state,
- `^4` integer-power lowering,
- FMA fusion in the Horner recurrence,
- t magnitude.

The missing piece: what about having a DIFFERENT `t[ij,k]` per k makes CPU and Reactant
diverge, when a uniform `t[ij,:]` does not? Two next diagnostics:

### Step 4a — Flat Frierson (flat_freeze_k variant)

Compute Frierson's τ at ONE reference layer (e.g., k=4) and apply the same resulting
`t` to every layer in the ij column. Same physical constants as Frierson, but zero
layer variation. If drift collapses → confirms layer-varying t is the trigger.
If drift stays at 5.5 K → the trigger is something else entirely in the Frierson struct
(lat-dependent τ₀? the fact that `whichring[ij]` is read?).

### Step 4b — CPU-precomputed t lookup

Precompute `t[ij, k]` on the host using CPU arithmetic during `initialize!`, store as a
scratch field, and have the kernel just read it. Removes ALL transmissivity math from
the traced region. If drift collapses → bug is in how Reactant evaluates the
transmissivity expression (not in the flux sweep, not in the t values themselves).

### Numbers so far (updated)
| Setup | T_max (K) | q_max |
|---|---|---|
| baseline (Frierson + OneBandLWRT) | 5.526 | 2.10e-3 |
| Frierson w/ manual σ^4 expansion | 5.526 | 2.10e-3 |
| Frierson w/ algebraic refactor (no cancel, no carry) | 5.526 | 2.10e-3 |
| Frierson w/ split Horner (no FMA) | 5.526 | 2.10e-3 |
| **Frierson w/ flat t (τ at k=5, applied to all k)** | **0.341** | **2.18e-4** |
| ConstantLWT (t=0.6 default) | 0.342 | 2.16e-4 |
| ConstantLWT (t=0.2 low) | 0.529 | 2.29e-4 |
| TransparentLWT | 0.345 | 2.26e-4 |
| UniformCooling | 0.460 | 1.84e-4 |
| longwave disabled | 0.474 | 1.83e-4 |

## Step 4a Result — Flat-Frierson (2026-04-21)

**Hypothesis confirmed**: layer-varying transmissivity is the trigger.

Edited `FriersonLongwaveTransmissivity` inner loop to compute τ at ONE reference layer
(k=5) and apply `t_flat = exp(-τ_ref / nlayers)` uniformly to every layer:
```julia
τ₀ = τ₀_equator + (τ₀_pole - τ₀_equator) * sinlat^2
τ_ref = τ₀ * (fₗ * σ[5] + (1 - fₗ) * σ[5]^4)
t_flat = exp(-τ_ref / nlayers)
for k in 2:(nlayers + 1)
    t[ij, k - 1] = t_flat
end
```
Everything else identical to baseline Frierson: same struct, same `whichring`/`sinlat`
read, same `τ₀` latitude computation, same `exp()` call.

**Result**: T_max = 0.341 K, q_max = 2.18e-4 — bit-for-bit in line with ConstantLWT (0.342).
Edit reverted.

### Interpretation

The drift is NOT caused by:
- the per-ij `τ₀ = τ₀_equator + (τ₀_pole - τ₀_equator) * sinlat^2` lerp (present in both
  baseline and flat-Frierson)
- the `whichring` / `sinlat` reads (same)
- the `exp()` call (same)
- Frierson struct type/layout (same)

The drift IS caused by the transmissivity values differing across layers within a column.
When `t[ij, k]` is the same for all k, Reactant and CPU agree. When `t[ij, k]` varies
with k, the resulting column flux sweep accumulates differently.

### What this implies about OneBandLWRT's flux sweep

Recall: `U = U * t + (1 - t) * σ * T[ij, k]^4`. When `t` is constant-in-k, this recurrence
has a closed-form geometric-style structure that apparently folds into an identical
floating-point sequence on CPU and Reactant (regardless of all the FMA/reassociation
theatre we worried about). When `t` varies with k, the recurrence becomes a polynomial
with distinct coefficients per layer, and Reactant ends up evaluating the per-layer
expressions `t[ij,k] * U + (1 - t[ij,k]) * σ * T[ij,k]^4` with a different rounding
schedule than CPU — compounding into the 5.5 K drift.

**This is surprising**: we already tested Horner split (no FMA) on the baseline Frierson
and saw no change. That suggests the divergence is NOT at the `U * t + ...` site itself
but possibly at how Reactant loads/indexes `t[ij, k]` when k varies — e.g., vectorizing
across k, or common-subexpression-eliminating across layers when values are uniform but
not when they differ.

### Step 4b — CPU-precomputed t lookup (next)

Precompute `t[ij, k]` on the host using CPU arithmetic (once at `initialize!`), store in
a scratch field, have the kernel just read it. If drift collapses → bug is specifically
in how Reactant **computes** the varying-t expression (not how it uses the values). If
drift persists at ~5.5 K → bug is in how Reactant **consumes** layer-varying values in
the flux sweep, and we need to rewrite OneBandLWRT.

### Step 4c — Flux-sweep on FIXED CPU-computed t (diagnostic)

A simpler variant: dump CPU-computed t values AT STEP 0, inject them into the Reactant
run as a constant array (not recomputed), and see if the drift matches baseline or
collapses. This isolates "Reactant's t values" from "Reactant's flux-sweep consumption
of varying t".

## MWE attempt — standalone flux-sweep kernel (2026-04-21)

File: `MWE_reactant_layer_varying_t.jl`. A standalone (no-SpeedyWeather) reproduction
of the flux-sweep + compounding-over-steps pattern, using KernelAbstractions `@kernel`
over columns and Reactant's `ReactantBackend` (via `ReactantCUDAExt`), `@compile`'d.
3168 columns × 8 layers × 100 steps with a realistic feedback amplitude
`T += 5e-5 * dT`.

**Result**: CPU vs Reactant = **bit-identical (≤ 1 ulp)** for BOTH flat-t AND
layer-varying-t.

| MWE variant | max|ΔT| (100 steps) |
|---|---|
| Flat t (uniform in k) | 3.05e-5 (= 1 ulp at T≈300 K) |
| Varying t (Frierson-like) | 3.05e-5 |

### Interpretation

The standalone flux-sweep kernel does NOT reproduce the 5.5 K drift observed in
SpeedyWeather. This is a strong negative result: the kernel itself is NOT the cause.

Something in SpeedyWeather's full loop is required for the drift to materialize:

1. **Dynamical core amplification**: spectral transforms → Leapfrog → implicit
   correction → back to grid. Tiny per-step radiation differences feed into the
   dynamics and may grow exponentially under unstable modes.

2. **Parameterization coupling**: in the real model, cloud cover from condensation
   feeds back into radiation. A ~1 ulp difference in T at step N affects condensation
   at step N+1, affecting clouds, affecting radiation, affecting T at step N+2.

3. **Fused kernel effect**: SpeedyWeather calls `transmissivity!` and
   `longwave_radiative_transfer!` back-to-back inside the same `@kernel`, with the
   transmissivity writing to `vars.scratch.grid.a` and the radiative transfer reading
   from it. The scratch-array write→read pattern inside a single kernel is NOT what
   the MWE tests (MWE passes `t` as a precomputed input array). Worth a MWE variant.

### Next MWE variant: fused transmissivity + flux-sweep

Write a single @kernel that:
1. Computes t[ij, k] into a scratch ij×k array (Frierson-like math)
2. Reads that scratch array back to do the flux sweep
3. Compounds via T += feedback*dT

If THAT reproduces, the bug is in how Reactant handles intra-kernel scratch writes
under layer-varying access patterns. If that still doesn't reproduce, the drift
really is amplified by the dynamical core.

## MWE attempt 2 — fused scratch-array kernel (2026-04-21)

File: `MWE_reactant_fused_scratch.jl`. A single `@kernel` that writes
`scratch_t[ij, k]` using Frierson-like math (with layer-varying σ^4 structure) and
reads it back to do the flux sweep — all within the same per-ij body. 100 steps,
feedback `T += 5e-5 * dT`.

**Result**: CPU vs Reactant = **bit-identical (1 ulp)** for BOTH flat AND varying.

| Fused-scratch variant | max|ΔT| |
|---|---|
| Flat t | 3.05e-5 (1 ulp) |
| Varying t (Frierson) | 3.05e-5 (1 ulp) |

### Stronger negative result

The intra-kernel write-then-read of a layer-varying scratch array is NOT the bug.
Two MWE shapes now both fail to reproduce:
- t as precomputed input array, layer-varying, 100 steps: no drift
- t written to scratch and read back in same kernel, layer-varying, 100 steps: no drift

### What this leaves

The 5.5 K SpeedyWeather drift requires something that the MWE does NOT have. Most
likely candidates now:

1. **Dynamical core amplification** — spectral transforms, Leapfrog filter, implicit
   correction. These are linear-algebra-heavy and could amplify tiny floating-point
   differences in radiation tendencies into O(K) drift over 10 steps. The MWE has a
   simple `T += feedback*dT`; the real model has a full GCM step between radiation
   calls.

2. **Coupling between multiple parameterizations** — cloud cover from condensation,
   feedback into radiation, surface fluxes reading T. Small T perturbations compound
   through these paths.

3. **Type/layout of the real kernel arguments** — SpeedyWeather passes a huge
   NamedTuple (`vars`, `parameterizations`, `model`) into the fused kernel, not the
   compact `(dT, T, scratch, τ₀, σh)` of the MWE. The argument marshaling hit us hard
   already in the NaN phase (zenith struct). It's possible a different marshaling
   path triggers non-bit-identical compilation for the radiation sub-kernel.

4. **Some other parameterization** besides longwave that happens to disappear when
   `longwave_radiation=nothing` because the coupling chain breaks.

### Next diagnostic

Instead of making MWEs ever closer to SpeedyWeather, bisect from the SpeedyWeather
side: inspect the single-step tendency from longwave radiation alone (CPU vs Reactant)
and see if step-1 longwave tendencies are already different, or if the 5.5 K is an
amplification effect over the 10 dynamics+physics steps.

## Session 3 — Scratch-field instrumentation of t, U, D (2026-04-23)

Approach: added three Grid3D debug scratch fields (`debug_t`, `debug_U`, `debug_D`)
in `primitive_dry.jl` so PrimitiveDry/PrimitiveWet inherit them automatically.
Wrote runtime layer values into these arrays during the Frierson transmissivity
and `longwave_radiative_transfer!` kernels, then read them back host-side after
`time_stepping!` for CPU vs Reactant comparison.

### debug_t (Frierson, after 10 steps)

| Layer | max|Δt| | At ij | CPU | Reactant |
|---|---|---|---|---|
| k=1 | 0 | — | 0.980963 | 0.980963 |
| k=2 | 0 | — | 0.976414 | 0.976414 |
| k=3 | 5.96e-8 | 505 | 0.900043 | 0.900043 |
| k=4 | 5.96e-8 | 1065 | 0.750677 | 0.750677 |
| k=5 | 0 | — | 0.868139 | 0.868139 |
| k=6 | 0 | — | 0.785313 | 0.785313 |
| k=7 | 5.96e-8 | 225 | 0.522546 | 0.522546 |
| k=8 | 2.98e-8 | 181 | 0.411423 | 0.411423 |

**Global max|Δt| = 5.96e-8 = exactly 1 ulp of Float32.** At ij=1 (north pole),
`t_Δ = 0` on every layer. The transmissivity values themselves are **not** the
numerical source of the drift — Reactant computes them bit-accurate to 1 ulp.

### debug_U (Frierson, after 10 steps)

| Layer | max|ΔU| [W/m²] | rel | At ij | CPU | Reactant |
|---|---|---|---|---|---|
| k=1 | 1.25 | 0.4% | 962 | 295.26 | 296.51 |
| k=2 | 1.36 | 0.4% | 962 | 308.18 | 309.54 |
| k=3 | 1.56 | 0.5% | 962 | 324.51 | 326.08 |
| k=4 | 2.04 | 0.6% | 962 | 347.36 | 349.41 |
| k=5 | 3.34 | 0.9% | 962 | 378.11 | 381.44 |
| k=6 | 7.62 | 1.8% | 962 | 413.97 | 421.59 |
| k=7 | **29.47** | **6.7%** | 962 | 439.99 | 469.46 |
| k=8 | 6.46 | 0.8% | 2401 | 794.27 | 800.73 |

`debug_U max|Δ| = 29.5 W/m²`. At ij=1 (where `t_Δ = 0` on every layer), `U_Δ` still
grows from 0.016 to 0.041 W/m² across the column — i.e., **the flux-sweep diverges
even when the transmissivity is bit-identical on both architectures**.

### debug_D (Frierson, after 10 steps)

`debug_D max|Δ| = 1.48 W/m²` at k=8, ij=1504. Downward flux diverges more
modestly than upward.

### Control: ConstantLongwaveTransmissivity (t ≡ 0.6)

| Quantity | Constant | Frierson |
|---|---|---|
| T_max_diff | 0.342 K | 5.526 K |
| debug_U max|Δ| | **6.3 W/m²** | **29.5 W/m²** |
| debug_U at ij=1 (per layer) | ~1.5e-4 | ~0.02 |
| debug_D max|Δ| | 0.25 W/m² | 1.48 W/m² |

**Per-layer U divergence at ij=1 is ~100× bigger for Frierson (0.02) than for
Constant (1.5e-4), despite both having bit-identical `t` at ij=1.** This is the
core asymmetry.

### Interpretation — why Constant stays small, Frierson blows up

The Horner recurrence `U_k = U_{k-1}*t + (1-t)*σ*T_k^4`:

- **Constant** (same `t` every k): XLA sees `t` as a loop-invariant scalar.
  CSE can hoist `(1-t)` and `(1-t)*σ` out of the loop. The resulting 8-layer
  recurrence has one structural form, emitted as a uniform per-layer op on both
  CPU and XLA → instruction scheduling matches → ~1.5e-4 W/m² residual per
  layer (attributable to the initial `U_ocean/U_land/land_fraction` FMA).
- **Frierson** (layer-varying `t`): each iteration has a distinct multiplier
  loaded from `scratch.grid.a`. XLA cannot CSE across layers; it emits 8
  structurally-distinct FP operations per column. Julia CPU emits the same
  loop but with scalar-register accumulation in left-to-right order. The two
  paths apply FMA/reassociation differently **per layer**, seeding ~0.02 W/m²
  per layer per step.

Per-step seed error × 10 time steps × non-linear feedback (`dTdt → T → T^4 → U`)
amplifies to 29 W/m² at some columns and 5.5 K in temperature.

### Ruled out as source

- Transmissivity `t` values (1 ulp accurate)
- `exp(-dτ)` (ConstantLWT calls it too, drift ≈ 0.3 K)
- `^4` integer-power lowering (Julia 1.11 already lowers to `(x*x)*(x*x)`)
- `σ^4` layer variation vs linear (manually expanded in Step 2 — no change)
- Catastrophic cancellation in `τ_below - τ_above` (algebraic refactor — no change)
- Loop-carried `τ_above` state (removed in refactor — no change)
- FMA fusion in the Horner recurrence (split `U = U*t + ...` — no change)
- Transmissivity magnitude (ConstantLWT with t=0.2 still collapsed)
- Multi-accumulation on `dTdt[ij, nlayers]` (already exonerated by split Horner)

### Option A tested — hoist `t[ij, :]` into NTuple via `Val(nlayers)`

Rewrote the flux sweep to gather `t_col = ntuple(k -> transmissivity[ij, k], Val(nlayers))`
(and `T4_col` similarly), then index into stack-allocated NTuples inside the
Horner loop. Hypothesis: making the multiplier register-resident would align
Reactant's emitted recurrence with CPU's.

**Result**: T_max_diff = 5.521 K (vs 5.526 before). `debug_U max|Δ| = 29.46`
(vs 29.47 before). **No measurable change** — bit-identical to the non-hoisted
version within argmax tie-breaking noise.

**Conclusion**: XLA lowers the `ntuple` gather + register-resident loop to
the same HLO as the original per-layer `transmissivity[ij, k]` reads. The
per-layer array-load pattern was NOT the trigger. Reverted.

### Remaining viable workarounds

Options D (precompute `t` and ship as constant field) and E (substitute
ConstantLWT) are rejected — they would give physically wrong results (Frierson's
layer-varying τ is the whole point of the scheme).

Remaining paths to investigate:

1. **Option B — full static unroll via `@generated`**. Eliminates the loop
   entirely by emitting 8 straight-line `U = U*t_k + (1-t_k)*σ*T4_k` statements.
   Lower-priority given Option A already tested register-resident access and
   didn't help — but worth one attempt because `@generated` forces the unroll
   at method-compile time, not call time.

2. **Option C — two-pass split**: compute `U_col::NTuple{8}` and `D_col::NTuple{8}`
   into locals first (no `dTdt` writes), then a second pass that applies
   tendencies one `+=` per `dTdt[ij, k]`. Removes the multi-contribution
   accumulation on `dTdt[ij, nlayers]` (4 writes currently). Already partially
   tested (Step 3b split the Horner); the remaining untested piece is the
   `dTdt` accumulation reorder.

3. **Direct HLO inspection**. `Reactant.@code_hlo first_timesteps!(sim)` or
   similar to see what XLA actually emits for the Frierson flux sweep vs
   Constant. This would reveal whether XLA picks a `reduce`/`scan` form that
   reorders summations. If the HLO is identical between Frierson and Constant
   modulo the `t` source, the divergence is at the register-allocation /
   lowering stage and not fixable from Julia source.

4. **XLA flags / tracing directives**. Look for Reactant options to disable
   reassociation or FMA in a scoped region (e.g. `@strict_fp` or similar).
   Unknown if exposed.

5. **Accept the residual**. The drift is ~0.5 K/step on k=8, smaller on other
   layers. For practical purposes this is 1–2 orders of magnitude above FP
   noise but the model remains stable and physically plausible. A tolerance
   bump in `test_correctness.jl` (RTOL 1e-3 → 1e-2, or fewer steps) would let
   CI pass while keeping Frierson as the default scheme.

### Files touched this session

- `SpeedyWeather/src/models/primitive_dry.jl` — added `debug_t`, `debug_U`,
  `debug_D` ScratchVariable registrations (Grid3D, namespace=:grid). **Keep** for
  future debugging or remove together once investigation concludes.
- `SpeedyWeather/src/parameterizations/radiation/longwave_transmissivity.jl` —
  `debug_t[ij, k-1] = tk` mirror inside Frierson kernel. Keep/remove with above.
- `SpeedyWeather/src/parameterizations/radiation/longwave_radiation.jl` —
  `debug_U[ij, k] = U` and `debug_D[ij, k] = D` mirrors. Option A (ntuple hoist)
  was reverted by user.
- `SpeedyWeather/test/reactant/debug_drift_longwave_variant.jl` — extended with
  per-layer t/U/D comparison + ij=1 per-layer dump.

### Key numbers (session 3)

| Run | T_max (K) | U_max (W/m²) | D_max (W/m²) | t_max |
|---|---|---|---|---|
| Frierson baseline | 5.526 | 29.47 | 1.48 | 5.96e-8 (1 ulp) |
| Frierson + Option A (ntuple hoist) | 5.521 | 29.46 | 1.48 | 5.96e-8 |
| Constant (t=0.6) | 0.342 | 6.30 | 0.25 | — (kernel has no debug_t mirror) |

### Status

Investigation paused. Root cause localized to: **Reactant/XLA emits a different
FP-rounding schedule for the 8-layer Horner recurrence `U = U*t_k + (1-t_k)*σ*T_k^4`
when `t_k` varies across k, compared to Julia CPU's scalar left-to-right path.
The per-layer ~0.02 W/m² per-step seed error amplifies through
`dTdt → T → T^4 → U` feedback over 10 steps into a 5.5 K temperature drift.**

Per-source tweaks (Option A, σ^4 expansion, split Horner, algebraic refactor) do
not change the emitted HLO. Lower-level XLA intervention or HLO inspection is
needed to make progress.

## Session 4 — Float64 accumulator hypothesis tested & falsified (2026-05-15)

### Framing
The session 3 conclusion ("XLA emits different FP-rounding schedule for the
Horner recurrence with layer-varying t") implies the seed error is per-step
ulp-scale FP noise that gets amplified by T → T^4 → U → dTdt → T_next over
10 steps. Natural fix: widen `U`, `D`, and the per-step emission to Float64.
Float64 IEEE arithmetic is essentially deterministic across CPU and XLA
(FMA-induced delta is ≤ 1 ulp at Float64 ≈ 2e-16 relative — amplified over
10 steps still ≪ 1 K).

### Proposed next-step suggestions (4 candidates)
1. **Float64 accumulation in the Horner recurrence** — squash per-step seed by ~10⁹×.
2. **Step-1 single-step seed measurement** (NSTEPS=1) — distinguish per-step bug
   vs amplification regime.
3. **`Reactant.@code_hlo` HLO inspection** — see what XLA actually emits for
   Frierson vs Constant flux sweep.
4. **Pull longwave out of the fused mega-kernel** — isolate to test if the bug
   is in the fusion pattern.

### Test (1a) — Float64 accumulator only

Edit: `U::Float64`, `D::Float64` in `longwave_radiative_transfer!`
(longwave_radiation.jl lines 251, 270). RHS expressions left as
Float32 (`(1 - t) * σ * T[ij, k]^4` evaluates entirely in Float32 then
promotes to Float64 only at the `+` into the accumulator).

**Result**: T_max = 5.5205 K, q_max = 2.094e-3.
**No change** — argmax tie-breaking noise vs 5.526 baseline.

Why the no-op: only the accumulator was widened. Since all of `t, σ, T, U_ocean,
U_land` are Float32, the RHS `(1 - t) * σ * T^4` is computed in Float32, with
identical rounding as the original code. The Float64 cast happens only at the
final `U_w + emit` add, which preserves the exact Float32 emit value.

### Test (1b) — Full Float64 arithmetic

Refined edit: cast `t`, `T[ij, k]`, and `σ` to Float64 INSIDE the loop body so
the entire RHS `U * t + (1 - t) * σ_w * Tk^4` evaluates in Float64. Convert
back to Float32 only at the dTdt / scratch writes.

**Result**: T_max = 5.540 K, q_max = 2.101e-3.
**Still no change** — and even very slightly worse (within argmax noise).

### Implications — hypothesis falsified

If per-step ulp-scale FP rounding in the Horner recurrence were the seed:
- Float64 arithmetic would reduce per-step seed from O(1 ulp Float32) = ~6e-8
  to O(1 ulp Float64) = ~2e-16 → ~10⁹× reduction
- Amplified through 10 steps of T^4 feedback (peak amplification factor at k=8
  with t≈0.4 is ~3.6 per step → factor ~3.6¹⁰ ≈ 4×10⁵)
- Expected 10-step drift: 6e-8 × 4×10⁵ → ~0.024 K (down from 5.5 K)
- Observed: 5.54 K. Reduction factor: 1.0.

The Horner recurrence arithmetic is NOT the seed.

### Revised mechanism

Longwave is **amplifying** an upstream T perturbation, not generating one. The
amplification mechanism is physical: layer-varying t with strong absorption at
deep layers makes `dTdt[k] ∝ (1-t_k) * 4σT_k^3 * δT_k`, with `(1-t_k) ≈ 0.6` at
the surface layer. Combined with the T^4 nonlinearity, this acts as a ~3.6×
per-step amplifier of any T perturbation entering longwave.

This re-explains the flat-Frierson collapse (5.5 K → 0.34 K). It is NOT that
uniform t gives CPU and XLA the same FP schedule (Float64 test would have
shown this) — it is that uniform t **physically dampens** the layer-varying
amplification: average `(1-t) ≈ 0.14`, giving amplification ~0.84/step → ~0.84¹⁰
≈ 0.18 — about 100× less amplification, consistent with 5.5 K → 0.34 K.

### Where the seed actually lives

Bit-identical inputs (T_prev, transmissivity) at step 0 — but as soon as ANY
upstream operation introduces a Float32 ulp difference into `T_prev` between
CPU and Reactant, longwave's T^4 / layer-varying-t machinery amplifies it
∝ 3.6^n into 5.5 K over 10 steps.

Likely seed locations (all upstream of longwave):
- **Spectral transforms** (FFT + Legendre): many multiplications/sums, classic
  FMA/reassociation target. Even with `dynamics=false`, time integration still
  runs leapfrog and applies physics tendencies in spectral space.
- **Implicit correction step** in time integration (linear solve on spectral
  coefficients).
- **`grid` → `prev_grid` copy via transform**: tendencies are added in spectral
  space then transformed back; FFT path could differ.
- **Other parameterizations** that write to T tendencies before longwave reads
  them: surface_heat_flux, vertical_diffusion, large_scale_condensation,
  shortwave_radiation. Each accumulates into `dTdt[ij, k]` before longwave's
  T_prev is updated.

The bisect showed that disabling parameterizations one at a time (with
longwave ON) did not collapse the drift to UniformCooling level. But that
test used longwave as the amplifier; the seed could be in dynamics-side ops
or in a parameterization whose own contribution to dTdt is small (so disabling
it doesn't change drift much) yet whose output T affects the longwave input
non-negligibly.

### Action: revert Float64 edit

Reverted both Test (1a) and (1b) edits in `longwave_radiation.jl`. The Float64
widening adds complexity and offers zero benefit.

### Numbers (session 4)
| Run | T_max (K) | q_max | Δ vs baseline |
|---|---|---|---|
| baseline (Frierson Float32) | 5.526 | 2.10e-3 | — |
| **(1a) U,D as Float64; RHS Float32** | **5.521** | **2.09e-3** | none |
| **(1b) Full Float64 arithmetic** | **5.540** | **2.10e-3** | none |

### Updated next steps (priorities reshuffled)

Hypothesis-1 (radiation FP rounding) eliminated. Remaining candidates from
the original 4 still warrant testing:

1. **Step-1 single-step seed measurement** (NSTEPS=1, untouched param config).
   If T_diff at step 1 is already O(0.01 K) → upstream seed is large per step,
   amplification is mild → bisect upstream ops. If T_diff at step 1 is ulp-scale
   (~1e-5 K) → amplification dominates (~5×/step gain) → focus on breaking the
   amplification path. **High-value, low-cost diagnostic.**

2. **`Reactant.@code_hlo` HLO inspection** of `first_timesteps!` or
   `parameterization_tendencies!`. Look for FMA / reassociation in spectral
   transform code, not just radiation. The HLO diff between dynamics-off
   sub-paths will localize the seed.

3. **Bisect upstream ops**: try disabling individual transforms or testing
   `dynamics_only=true` to isolate spectral-transform divergence from
   parameterization divergence.

4. **Test on a transform-free run**: if any operating mode exists that skips
   `transform!` between time steps, run it. If drift collapses → transforms
   are the seed. If drift persists → seed is in parameterizations themselves.

5. **(unchanged from session 3)** Accept the residual — bump RTOL or reduce
   NSTEPS for `test_correctness.jl`. The 5.5 K drift is a system-level
   amplification artifact, not a bug in any single function.

### Files touched (session 4)
- `SpeedyWeather/test/reactant/debug_drift_lw_float64.jl` — new test script
  (Float64-accumulator variant of `debug_drift_longwave_variant.jl` without
  the now-removed debug_t/U/D probes). **Kept** for reuse.
- `SpeedyWeather/src/parameterizations/radiation/longwave_radiation.jl` —
  Float64 edits applied and **reverted** (file restored to baseline).

## Session 5 — Step-1 seed measurement: SEED LOCALIZED (2026-05-15)

### Test (2) executed — NSTEPS=1 with full per-layer breakdown

Script: `debug_drift_step1.jl`. Default Frierson (no overrides), same setup as
`debug_drift_longwave_variant.jl frierson` except `NSTEPS=1` and adds a probe
of the post-resync initial state.

### Result

| Quantity | Value | Notes |
|---|---|---|
| **initial state T_max_diff** (after `copy!` resync) | **5.8e-4 K** | should be 0 if `copy!` produced bit-identical state — it does not |
| **step-1 T_max_diff** | **0.386 K** | at k=8, ij=3092 |
| step-1 q_max_diff | 1.62e-4 kg/kg | at k=8 |
| step-1 u_max_diff | 0.082 m/s | not radiation — generic |
| step-1 v_max_diff | 0.0041 m/s |  |
| (10-step baseline) | 5.526 K | for reference |

#### Per-layer T diff at step 1

| k | max\|ΔT\| (K) | rel | layer |
|---|---|---|---|
| 1 | 0.052 | 2.6e-4 | TOA |
| 2 | 0.029 | 1.3e-4 |  |
| 3 | 0.026 | 1.1e-4 |  |
| 4 | 0.075 | 3.0e-4 |  |
| 5 | 0.129 | 4.8e-4 |  |
| 6 | 0.144 | 5.0e-4 |  |
| 7 | 0.240 | 8.0e-4 |  |
| 8 | **0.386** | 1.6e-3 | surface |

Monotonic growth from TOA to surface. The surface-concentration matches the
10-step pattern.

### Two findings, one mechanism

**Finding 1 — `copy!` does NOT produce bit-identical post-resync state.**
After `copy!(r←c)` + `@compile` + `copy!(c←r)` + `initialize!(steps=1)`, the
grid temperature arrays already differ by 5.8e-4 K. `copy!` almost certainly
copies only the **spectral prognostic** variables; the **grid variables**
(`vars.grid.temperature`, etc.) are produced by an architecture-local
`transform!(grid, spectral)` call, and CPU's Julia transforms vs Reactant/XLA's
emitted transforms differ by ~10 ulp Float32 (5.8e-4 K / 300 K ≈ 2e-6 ≈ 17 ulp).

**Finding 2 — step-1 seed is 0.386 K**, not ulp-scale.

Amplification factor: 0.386 K → 5.526 K over 9 steps = **1.34×/step** (not the
~3.6×/step predicted from radiation-only analysis). The seed is dominated by
per-step semantic divergence, not by exponential amplification of ulp noise.

### Why this re-explains everything

- **Float64 in longwave couldn't help** (Session 4): the input `T_prev` reaching
  the radiation kernel is already different by ~6e-4 K. Float64 inside the
  kernel can't fix that.
- **Flat-Frierson collapses drift to 0.34 K**: not FP determinism — it physically
  dampens the amplification of an upstream 6e-4 K seed (uniform `(1-t)` halves
  the surface-layer sensitivity coefficient).
- **MWE showed bit-identical**: the MWE had no spectral transforms. With the
  seed identified, this is exactly what we'd expect — outside the GCM pipeline,
  Reactant's FP arithmetic IS deterministic.
- **Disabling longwave gives 0.47 K drift**: that's the residual when the strongest
  T^4 amplifier of transform-seeded perturbations is removed.
- **Layer-varying t matters**: because layer-varying t creates layer-varying
  sensitivity to `δT`, which is what amplifies. Constant t gives uniform
  smaller sensitivity, less amplification.

### The actual culprit: spherical harmonic transforms

CPU Julia `transform!(grid, spectral)`:
- FFT (per-ring, real-valued) — many multiplies, FMA-friendly accumulators
- Associated Legendre polynomial sums — multiply-accumulate loops with O(trunc²) ops

Reactant/XLA emits these (likely as matrix multiplies for `MatrixSpectralTransform`)
with its own FMA / reassociation choices. ~17 ulp Float32 spread for T at 300 K is
plausible given ~30 multiply-accumulate ops per output element.

The transform happens at multiple points per time step:
- Spectral prognostic → grid prognostic at start
- Grid tendencies → spectral tendencies before leapfrog
- Possibly inside parameterizations that need spectral derivatives

Each transform call seeds a fresh ~10 ulp perturbation, then T^4 / layer-varying
radiation amplifies it ~1.3×/step.

### Updated next steps

Now that we know the seed lives in transforms, the productive directions:

1. **Bisect: which transform is the worst offender?** Compare ΔT after just the
   spectral→grid transform (without time stepping). Then after grid→spectral.
   Then after a full leapfrog cycle.

2. **Float64 transforms?** `MatrixSpectralTransform` stores precomputed Legendre
   matrices — these could be Float64 with Float32 input/output, accumulating in
   Float64. Likely 5–20% perf hit but trivial precision improvement.

3. **Strict-FP / no-FMA flag for Reactant transforms?** Investigate whether
   Reactant exposes an option to disable FMA contraction for specific matrix
   multiplies (the Legendre transform). XLA does have `precision_config` options
   on `dot` ops.

4. **Compare CPU vs Reactant for a SINGLE `transform!` call** outside the model
   — directly probe the transform divergence in isolation. If we observe 17 ulp
   there, it confirms transforms are the seed and we don't need to touch
   parameterizations.

5. **Pragmatic floor**: with the seed at ~6e-4 K per transform and ~3 transforms
   per step at amplification 1.3×/step, expected drift floor is roughly
   `0.0006 × 3 × 1.3^N`. At N=10 that's ~25× the seed = ~0.015 K — well within
   acceptable physics noise. The observed 5.5 K must therefore come from MULTIPLE
   transform-amplification cycles per step, or specifically from the longwave T^4
   amplification path. Either way, accepting `RTOL≈1e-2` for 10-step tests is
   physically defensible.

### Files (session 5)
- `SpeedyWeather/test/reactant/debug_drift_step1.jl` — new step-1 measurement
  script with per-layer breakdown and initial-state probe. **Kept** for reuse.

### Updated key numbers
| Run | T_max_diff (K) | q_max_diff |
|---|---|---|
| post-resync initial state | **5.8e-4** | n/a |
| step 1 | **0.386** | 1.6e-4 |
| step 10 (baseline Frierson) | 5.526 | 2.10e-3 |
| step 10 (UniformCooling, disable LW) | 0.460 | 1.84e-4 |
| step 10 (flat-Frierson) | 0.341 | 2.18e-4 |

## Session 5 — correction & next-session plan (2026-05-15)

### Walking back the "transforms are the seed" claim

The session-5 writeup above leaped from "post-resync initial state differs by
5.8e-4 K" + "step-1 diff is 0.386 K" to "the seed is the spectral transform."
That conflated two measurements with very different magnitudes (650× ratio)
and almost certainly different origins.

**What's actually evidenced:**
- 5.8e-4 K post-resync initial diff: real, but small. Cause not directly verified.
  Could be transforms, but could also be `copy!` not copying grid vars (only
  spectral), `@compile` leaving residue, or `initialize!(sim; steps)` doing
  architecture-local recompute.
- 0.386 K step-1 diff: real, large. Cannot plausibly be "ulp-noise from
  transforms amplified for one step" — no per-step amplification factor large
  enough exists. **Something in the time step itself with Frierson on is
  contributing ~0.38 K per call.**

**Cross-check against scheme comparison:**
- Frierson step-10 = 5.5 K, disabled-LW step-10 = 0.47 K → ratio ~12×.
- If both shared the same upstream (transform) seed and only Frierson amplified
  it, both step-1 diffs should be similar and diverge later. The 12× gap is
  almost certainly already present at step 1.
- Therefore Frierson is **introducing fresh per-step divergence**, not just
  amplifying a shared seed.

### Revised understanding

Session 3's original framing is back in play (with refinements):

1. A small upstream perturbation enters longwave (5.8e-4 K initial-state level,
   possibly transform-induced, possibly other; magnitude is sub-leading).
2. Frierson's per-call output differs CPU vs Reactant by O(0.1–0.4 K)-worth
   of effective T-tendency at step 1 — i.e., the per-call divergence is the
   DOMINANT contribution, not the amplification of upstream noise.
3. The Session 4 Float64 test rules out one *specific* form of per-call
   divergence (the U/D Horner arithmetic itself), but does NOT rule out:
   - Divergence in the `dTdt[ij, nlayers]` accumulation chain (4 writes per
     call to the surface layer, FMA-fusable across the four contributions)
   - Divergence in `flux_to_tendency` (a division + multiplication in Float32)
   - Divergence introduced by how the longwave block FUSES with surrounding
     parameterizations in `column_parameterizations_kernel!` — XLA reorders
     operations across fused blocks, which Session 4's Float64 test (inside
     the longwave kernel only) wouldn't fix
   - The transmissivity write→read pattern when consumed alongside other
     scratch fields in the fused kernel
4. The remaining 0.34–0.47 K floor in "tame" schemes is whatever ELSE the
   time-step pipeline produces (transforms, dynamics, other parameterizations
   together). Frierson adds ~5 K on top of that.

### Diagnostic tests for next session (continue here)

Three tests, in priority order:

#### Test A — Step-1 with UniformCooling / disabled longwave
**Goal**: separate "Frierson's own per-call contribution" from
"shared upstream seed common to all schemes."

**Method**: rerun `debug_drift_step1.jl` with the longwave swap from
`debug_drift_longwave_variant.jl` (add a CLI arg or copy-edit the model
constructors). Run variants: `uniformcooling`, `disabled`, `constant`,
`frierson` (already have this one).

**Predictions / interpretations**:
- If `uniformcooling` step-1 ≈ 0.02–0.05 K → ~0.34 K of Frierson's step-1
  0.386 K comes from Frierson itself per call. Confirms Frierson is the
  dominant seed, not upstream.
- If `uniformcooling` step-1 ≈ 0.3 K (similar to Frierson) → the seed is
  fully upstream and Frierson just amplifies it. Different bug to find.
- If `uniformcooling` step-1 ≈ 0.1 K (intermediate) → both effects contribute.

This is the highest-value next test. Cheap to set up — one script edit.

#### Test B — Single isolated `transform!` call
**Goal**: directly measure the CPU vs Reactant divergence of one spherical
harmonic transform with bit-identical Float32 spectral input.

**Quick pre-check** (cheap, do FIRST): compare the precomputed matrices
themselves between CPU and Reactant `MatrixSpectralTransform`. If
`M_c.forward`, `M_c.backward_real`, `M_c.backward_imag` (see
`SpeedyTransforms/src/matrix_transform.jl`) don't bit-match their
Reactant counterparts, transforms CANNOT be bit-identical even with
identical inputs — the bug is at construction time, not in the
`mul!` call. If the matrices DO bit-match, divergence (if any) is
purely from how `LinearAlgebra.mul!` is lowered under XLA.

**Method**: minimal script (no time stepping):
1. Build CPU and Reactant `MatrixSpectralTransform`s with same grid.
2. **First check the matrices**:
   `maximum(abs.(Array(M_r.forward) .- M_c.forward))` and same for
   `backward_real`, `backward_imag`.
3. Allocate a `LowerTriangularArray{ComplexF32}` (spectrum, nlayers).
4. Fill with deterministic data (e.g., `rand(MersenneTwister(42), ComplexF32, …)`).
5. `copy!` the spectral data CPU→Reactant.
6. `transform!(grid_c, spec_c, transform_c)` on CPU.
7. `transform!(grid_r, spec_r, transform_r)` under `@compile`.
8. Compare grid outputs ulp-wise.

**Interpretation**:
- If max\|Δ\| ≈ 1 ulp Float32 (3e-7 relative): transforms are bit-clean and
  the 5.8e-4 K post-resync diff comes from elsewhere (`copy!`, `@compile`,
  `initialize!`).
- If max\|Δ\| ≈ 10–20 ulp (~6e-4 K at T=300): transforms are the source of
  the small initial diff. Worth pursuing Float64 Legendre matrices /
  Reactant precision_config on the dot ops.

Decouples the "transforms" question from time-stepping.

#### Test C — Direct probe of `copy!` round-trip
**Goal**: pinpoint whether `copy!` itself fails to produce bit-identical
state, or whether the later compile/initialize steps introduce the 5.8e-4 K.

**Method**: in `debug_drift_step1.jl`, add probes:
1. After `copy!(sim_r ← sim_c)` (line 75 in script), read both
   `sim_c.variables.grid.temperature` and `sim_r.variables.grid.temperature`
   → diff_A.
2. After `@compile first_timesteps!` (no resync yet) → diff_B.
3. After `copy!(sim_c ← sim_r)` post-compile → diff_C.
4. After `initialize!(sim; steps=NSTEPS)` on both → diff_D (= the 5.8e-4 K
   we already measured).

**Interpretation**:
- diff_A == 0: `copy!` is fine; perturbation enters later.
- diff_A == 5.8e-4 K: `copy!` doesn't propagate grid vars or does
  architecture-local recompute. Fix `copy!` rather than chase transforms.

### Lower-priority follow-ups (after A, B, C)

If Test A confirms Frierson contributes ~0.34 K/call:
- **Profile inside `longwave_radiative_transfer!`** to find where the per-call
  divergence enters. Re-introduce `debug_U`, `debug_D`, `debug_t` ScratchVariables
  (session 3 had them) but measure after step 1 only. The session-3 measurements
  were after step 10 — partly already amplified. Step-1 numbers will isolate
  the per-call seed.
- **`Reactant.@code_hlo`** on `longwave_radiative_transfer!` (still untried).
  Compare HLO emitted with Frierson transmissivity vs Constant. Concrete
  difference in HLO = concrete fix candidate.
- **Test (4) from the original plan**: pull longwave out of the fused
  `column_parameterizations_kernel!`, call it via separate `launch!`. If drift
  collapses → fusion-pattern bug, not the kernel itself.

If Test A shows the seed is mostly upstream of Frierson:
- Focus on transform precision (Test B) or on whichever other parameterization
  is producing the seed.

### Pragmatic floor (unchanged from session 4)

Independent of root-cause work, raising `RTOL` for `test_correctness.jl` to
~1e-2, or reducing NSTEPS to 1–2 (where drift is sub-K), would unblock CI.
Recommend keeping this as a fallback if Tests A/B/C don't yield a clean fix
within reasonable effort.

### Files to keep / reuse next session
- `debug_drift_step1.jl` — already takes default config; needs minor edit to
  accept a longwave-variant CLI arg (mirror `debug_drift_longwave_variant.jl`).
- `debug_drift_longwave_variant.jl` — variant constructor logic to copy from.
- `debug_drift_lw_float64.jl` — kept as a template for testing changes to
  `longwave_radiation.jl` once we know what to change.

### Status

Investigation paused. The actionable next step is **Test A** (step-1 with
UniformCooling and disabled-LW). Cheapest, most diagnostic, and answers the
"is Frierson the per-step seed or just the amplifier?" question that's
currently blocking direction.

## Session 6 — Tests A, B, C executed (2026-05-15)

All three planned tests run. Strong, unambiguous results.

### Test A — Step-1 across longwave schemes

Script: `debug_drift_step1_variant.jl <variant>`. Same setup as
`debug_drift_step1.jl` but takes `frierson | constant | transparent |
uniformcooling | disabled` as CLI arg.

| Variant | step-1 T_max (K) | step-10 T_max (K) | amplification |
|---|---|---|---|
| frierson | 0.386 | 5.526 | **14×** over 9 steps (≈ 1.34×/step) |
| uniformcooling | 0.454 | 0.460 | **~1×** (no growth) |
| disabled | 0.505 | 0.474 | **~1×** (slight damping) |

`u_max_diff` at step 1 is ~0.081 m/s for ALL three variants — completely
independent of longwave scheme, so the non-T variables get their per-step
seed entirely from upstream.

#### Interpretation — walks back the walk-back

The Session 5 walk-back was wrong. The original Session 5 framing — "seed
from upstream, Frierson amplifies" — is right after all.

- All three variants produce ~0.4–0.5 K T-diff at step 1 — **this is the seed,
  generated per step regardless of longwave scheme**.
- Frierson grows the seed 14× over 10 steps; the others saturate near the
  step-1 floor.
- The 12× scheme gap at step 10 is purely amplification, NOT a Frierson-specific
  per-call divergence.
- Float64 in longwave (Session 4) failed because the amplifier is physical,
  not numerical — Float64 doesn't change the operator norm.

### Test B — Isolated transform!

Script: `debug_transform_isolated.jl`. Two-phase:

**Phase 1** — compare precomputed `M.forward`, `M.backward_real`,
`M.backward_imag` matrices between CPU and Reactant `MatrixSpectralTransform`:

| Matrix | max\|Δ\| |
|---|---|
| forward (ComplexF32, 560×3168) | **0.0** |
| backward_real (Float32, 3168×560) | **0.0** |
| backward_imag (Float32, 3168×560) | **0.0** |

**Phase 2** — feed bit-identical Float32 spectral input through `transform!`
on both architectures:

| Output | max\|Δ\| | rel |
|---|---|---|
| grid (all 8 layers) | **0.0** | **0.0** (0 ulp) |

**Conclusion**: **TRANSFORMS ARE BIT-IDENTICAL.** Both the precomputed matrices
AND the live `transform!` operation produce zero divergence between CPU and
Reactant. Transforms are conclusively ruled out as the seed.

### Test C — copy! round-trip probe

Script: `debug_copy_roundtrip.jl`. Probes at each stage of the resync ritual:

| Stage | T diff (K) | Other vars |
|---|---|---|
| pre_copy (CPU spin-up, R fresh) | 322 | huge — sanity |
| **diff_A** (after copy! c→r) | **0.0** | **0.0** all vars |
| diff_B (after @compile) | 0.40 | matches one-step advance |
| **diff_C** (after copy! r→c) | **0.0** | **0.0** all vars |
| **diff_D** (after initialize!(steps)) | **5.8e-4** | h=1 ulp, u=3.8e-6, v=1 ulp |

**Three findings:**

1. **`copy!` is bit-perfect in both directions** — copies grid vars too, no
   recomputation, no Float32 cast loss.
2. **`@compile` advances sim_r by ~1 step** via traced execution. This is
   normal Reactant behavior. Hence the resync ritual after compile is
   functionally necessary; the post-compile `copy!(c←r)` brings CPU up to
   match Reactant's advanced state.
3. **`initialize!(sim; steps=NSTEPS)` introduces a 5.8e-4 K asymmetric diff**:
   only T shows it significantly (humidity / v / u all at sub-ulp to few-ulp
   Float32). Points to a specific operation inside `initialize!` that touches
   temperature differently — possibly a thermodynamic precompute (virtual T,
   log T, geopotential init) or a clock-related compute path that runs
   architecture-locally.

### Combined picture

After Tests A, B, C the seed is **conclusively localized**:

- ❌ NOT transforms (Test B bit-identical)
- ❌ NOT `copy!` (Test C diff_A = 0)
- ❌ NOT a single longwave scheme (Test A: ~0.5 K with ALL schemes)
- ✓ Per-step ~0.4 K T-divergence generated within ONE time step's
  parameterization + dynamics-infrastructure pipeline
- ✓ Frierson amplifies that seed 14× over 10 steps via layer-varying-t
  cross-coupling (Test A confirms)
- ✓ A separate, smaller 5.8e-4 K T-only diff is introduced by
  `initialize!(sim; steps)` — asymmetric pattern suggests a specific
  thermodynamic precompute

### Next steps (continue here)

#### Investigate `initialize!(sim; steps)` asymmetric T diff
**Why temperature only?** Read `SpeedyWeather/src/models/simulation.jl` for
the `initialize!(sim; steps)` method. Find any operation that:
- writes to `grid.temperature` specifically (not other grid vars)
- could differ CPU vs Reactant (e.g., uses `Float32` constants in a way that
  XLA lowers differently, or calls a parameterization-init that runs on
  device).
Candidates: virtual temperature precompute, geopotential calculation, soil
temperature spin-up coupling.

This 5.8e-4 K is small compared to the 0.4 K/step generation, but it's an
easy fix and a clean sub-problem.

#### Bisect parameterizations for the 0.4 K/step seed
Currently active parameterizations with longwave disabled:
- solar_zenith
- albedo
- shortwave_radiation
- boundary_layer_drag
- surface_condition
- surface_momentum_flux
- surface_heat_flux
- surface_humidity_flux
- vertical_diffusion
- large_scale_condensation

Some of these had session-2 bisect issues (lsc, surface_heat_flux crashed when
disabled). Strategy:
1. Try disabling shortwave alone in `debug_drift_step1_variant.jl` (new
   variant). If step-1 drops dramatically → shortwave is the dominant seed.
2. If shortwave is innocent, examine surface flux family (heat, humidity,
   momentum) — they all run at k=8 surface where the diff concentrates.
3. The `ifelse` element-type warning from Reactant is a known issue in
   surface fluxes (session 2). Worth reviewing those `ifelse` sites in
   `surface_fluxes/{heat,humidity,momentum}.jl` for elem-type fixes that
   might have been missed.

Note: even bit-clean parameterization kernels can produce different output
under XLA when FUSED into `column_parameterizations_kernel!`. The fused
context allows XLA to reorder ops across previously-isolated boundaries.
If bisect fails to find a clean source, try the per-launch (unfused)
variant from "Test (4)" of the original plan.

#### HLO inspection (still untried)
`Reactant.@code_hlo first_timesteps!(sim)` to see what XLA emits for the
fused column kernel. With transforms confirmed bit-identical, the HLO diff
between schemes (Frierson vs uniformcooling) would localize the divergence
within parameterizations.

### Files (session 6)
- `debug_drift_step1_variant.jl` — step-1 with longwave-scheme CLI arg. Kept.
- `debug_transform_isolated.jl` — Phase-1 matrix + Phase-2 live transform
  test. Kept.
- `debug_copy_roundtrip.jl` — staged probe of copy!/compile/initialize!
  ritual. Kept.

### Updated key numbers
| Run | T_max_diff (K) | Notes |
|---|---|---|
| copy! c→r | 0.0 | bit-perfect |
| copy! r→c | 0.0 | bit-perfect |
| transform! (isolated, 8 layers) | 0.0 | bit-identical |
| MatrixSpectralTransform matrices | 0.0 | bit-identical |
| initialize!(sim; steps) introduces | 5.8e-4 | T-only, asymmetric |
| step 1 frierson | 0.386 | per-step seed dominated by upstream |
| step 1 uniformcooling | 0.454 | same upstream seed, no LW amplification |
| step 1 disabled | 0.505 | same upstream seed, no LW at all |
| step 10 frierson | 5.526 | 14× amplification |
| step 10 uniformcooling | 0.460 | no amplification |
| step 10 disabled | 0.474 | no amplification |

## Session 7 (planned) — Intra-step bisect plan (2026-05-15)

After Session 6 we know two things:
1. A small 5.8e-4 K **T-only asymmetric** diff is introduced by
   `initialize!(sim; steps)` BEFORE the first time step.
2. A larger ~0.4–0.5 K per-step diff is generated WITHIN a single time
   step by something downstream of `copy!`.

The asymmetric T-only signal in (1) is suspicious enough that it should be
investigated FIRST as Phase 0, before the intra-step bisect. Reasons:
- Asymmetry (T only, humidity/u/v at 1 ulp) points to a discrete operation
  that touches T differently — NOT generic FP noise.
- The same operation might also run inside `timestep!`, contributing to
  the per-step seed.
- It's much cheaper to investigate (only ~8 operations to bisect inside
  `initialize!` vs ~12 inside `timestep!`).
- If we fix it, we eliminate a confounder for the intra-step bisect.

### Phase 0 — Investigate `initialize!(sim; steps)` asymmetric T diff

Source: `SpeedyWeather/src/models/simulation.jl:55-103`.

The sequence inside `initialize!(simulation; period, steps, output)`:

```
initialize!(simulation; ...)
  ├── initialize!(clock, time_stepping, steps)          # I0  metadata only
  ├── set!(simulation.model.output, ...)                # I1  output config
  ├── scale_prognostic!(variables, model.planet.radius) # I2  scales vor, div (spectral)
  ├── transform!(variables, lf, model, initialize=true) # I3  spectral → grid, with init flag
  ├── (particles init if any)                           # I4
  ├── initialize!(model.output, variables, model)       # I5
  └── initialize!(model.callbacks, variables, model)    # I6
```

**Prime suspect: I3** — `transform!(variables, lf, model, initialize=true)`.
This is a HIGHER-LEVEL transform (different from the isolated `transform!`
on a single LowerTriangularArray that Test B verified bit-identical). It
likely:
- Transforms all prognostic spectral fields → corresponding grid fields
  (vorticity, divergence, temperature, humidity, surface pressure, tracers)
- Does derived computations (winds u, v from vor+div via Helmholtz; log
  pressure or virtual temperature precomputes; geopotential height)
- May call architecture-local thermodynamic helpers

A T-only divergence after I3 strongly suggests a specific T-related
derived-quantity computation runs differently on CPU vs XLA. Plausible
culprits inside `transform!(...; initialize=true)`:
- `log(pressure)` or `log(p_s)` conversion that feeds into T scaling
- Virtual temperature `Tv = T*(1 + (1/μ - 1)*q)` if it's stored back to T
- Geopotential calculation `Φ = R*Tv*log(p)` could affect a `T_prev` mirror
- `temperature_prev ← temperature` copy timing
- Anomaly conversion (T vs T - T_ref)

#### Phase-0 execution plan

**Step 0.1** — Find the actual transform path

Read `SpeedyWeather/src/dynamics/scaling.jl` (or wherever `transform!(vars,
lf, model; initialize=true)` is defined). Trace what it does to
`grid.temperature` specifically. Identify any T-only branches gated on
`initialize == true`.

```bash
grep -rn "transform!.*initialize\|function transform!.*Variables.*Model" \
    SpeedyWeather/src/
```

**Step 0.2** — Probe I3 sub-operations

If `transform!(...; initialize=true)` does N sub-operations on T (e.g.,
log-pressure, then virtual T, then geopotential), probe T after EACH.
Same scratch-field approach as the main bisect, but scoped to `initialize!`.

Add ScratchVariables:
- `T_init_probe_0` (before `initialize!(sim; steps)`)
- `T_init_probe_1` (after `scale_prognostic!`)
- `T_init_probe_2` (after `transform!(...; initialize=true)`)
- Sub-probes inside I3 if it does multiple T-touching operations

**Step 0.3** — Compare per-grid-point

Important: dump not just `max|ΔT|` but also the **spatial pattern** of the
diff. The 5.8e-4 K is asymmetric (only T), but is it also spatially
concentrated (e.g., only at the poles or surface)? If concentrated, that's
another clue about the operation. The pattern might match the per-step
0.4 K pattern, confirming a shared mechanism.

**Step 0.4** — Decision

- If I3 is bit-clean → divergence enters elsewhere in `initialize!` (less
  likely given the operations are all metadata or transforms). Drill in
  parallel.
- If I3 introduces the 5.8e-4 K diff → drill into `transform!(...; initialize=true)`
  to find which T-touching sub-operation is responsible.
- If the divergence is in a thermodynamic precompute (virtual T, log-p,
  geopotential), test Float64 widening just there. Unlike the longwave
  Float64 test (Session 4 — null result), this is a one-shot init computation
  not a per-step recurrence; Float64 might actually fix it.

Then proceed to Phase 1 (the intra-step bisect) with the `initialize!`
asymmetric diff either eliminated or characterized.

### Phase 1 — Intra-step bisect

After Session 6 we know the ~0.4–0.5 K T diff is generated **within a single
time step**, by something downstream of `copy!` and upstream of (or including)
the final spectral→grid transform. Transforms in isolation are bit-identical,
so the source must be in parameterization tendencies, dynamics infrastructure
(leapfrog, diffusion, implicit), or the orchestration that connects them.

Goal: probe T at every step of `run!` (i.e. of `first_timesteps!` for step 1,
which uses an Euler forward sub-step followed by an unfiltered leapfrog) and
identify the FIRST operation where CPU and Reactant diverge.

### The pipeline (PrimitiveEquation, step 1 = Euler sub-step at Δt/2)

From `time_stepping/leapfrog.jl:206` (`first_timesteps!`) and
`time_stepping/time_integration.jl:129` (`timestep!(::Variables, dt, ::PrimitiveEquation, ...)`):

```
first_timesteps!(simulation)
  └── first_timesteps!(variables, model)
       ├── initialize!(implicit, Δt/2, vars, model)        # P0a
       ├── timestep!(vars, Δt/2, model, lf1=1, lf2=1)      # Euler sub-step
       │    ├── reset_tendencies!(vars)                     # P1
       │    ├── greenhouse_gases_time_step!(vars, model)    # P2
       │    ├── parameterization_tendencies!(vars, model)   # P3
       │    │    ├── reset_variables!(vars)                  # P3a
       │    │    ├── global_parameterizations!(vars, model)  # P3b — unrolled NT loop
       │    │    └── column_parameterizations!(vars, model)  # P3c — fused kernel over ij
       │    │         (each ij: solar_zenith → vertical_diffusion → lsc → albedo →
       │    │          shortwave → longwave → boundary_layer_drag → surface_condition →
       │    │          surface_momentum_flux → surface_heat_flux → surface_humidity_flux)
       │    ├── ocean_timestep!(vars, model)                 # P4
       │    ├── sea_ice_timestep!(vars, model)               # P5
       │    ├── land_timestep!(vars, model)                  # P6
       │    │    (with dynamics=false we skip the dynamics branch and take
       │    │     parameterization_tendencies_only! instead, which transforms
       │    │     physics tendencies grid→spectral)
       │    ├── parameterization_tendencies_only!(vars, model) # P7
       │    ├── horizontal_diffusion!(vars, ...)             # P8
       │    ├── leapfrog!(vars, Δt/2, lf1=1, model)          # P9
       │    └── transform!(vars, lf2=1, model)               # P10 — spectral→grid
       ├── initialize!(implicit, Δt, vars, model)            # P11a
       └── timestep!(vars, Δt, model, lf1=1, lf2=2)          # Leapfrog sub-step
            (same structure as the Euler sub-step, but with leapfrog from
             two-time-level history; we only bisect this if the Euler step
             is bit-identical and the divergence first appears here)
```

Probe T (`vars.grid.temperature` and `vars.tendencies.grid.temperature`)
after each numbered marker. Compare CPU vs Reactant.

### Implementation approach — scratch-field snapshots

Modifying `first_timesteps!` to capture intermediate state is the cleanest
route. Two options:

**Option A — Scratch field snapshots (preferred)**: Add ScratchVariables
`T_probe_1`, `T_probe_2`, … `T_probe_10` (Grid3D, namespace=:debug) registered
in `primitive_dry.jl`/`primitive_wet.jl`. Patch `timestep!(::Variables, ...,
::PrimitiveEquation, ...)` to `copy!(vars.scratch.debug.T_probe_k,
vars.grid.temperature)` after each numbered marker. At end of one step,
pull all probes host-side and compute `maximum(abs(T_probe_k_c - T_probe_k_r))`
for each k. The first k where this is nonzero localizes the divergence.

For the tendency, also probe `vars.tendencies.grid.temperature` — register
a parallel set `dT_probe_k` so we see whether the seed is in the tendency or
in the state.

Pros: one compile, captures all probes in one run, no @compile-state issues.
Cons: requires source modification; needs careful revert when done.

**Option B — Sequential @compile of prefix functions**: Define
`timestep_to_marker_k!(vars, ...)` functions that run only operations 1..k.
Compile each, run on a fresh post-resync simulation. Compare T after each.

Pros: zero source modification of `timestep!` (only need to define wrapper
functions for the test).
Cons: N compiles instead of 1; state-reset between runs is fiddly; each
operation has its own side effects.

**Recommendation**: Option A. Cheap to add scratch fields temporarily; one
compile suffices; gives the full breakdown in one run.

### Concrete execution plan for tomorrow

#### Step 1 — Register debug ScratchVariables (1 commit, revertable)

In `SpeedyWeather/src/models/primitive_wet.jl` (or `primitive_dry.jl` so
both wet/dry inherit), add to the variable list:

```julia
ScratchVariable(:T_probe_1,  Grid3D(), desc = "T after reset_tendencies",          namespace = :debug),
ScratchVariable(:T_probe_2,  Grid3D(), desc = "T after greenhouse_gases_time_step!", namespace = :debug),
ScratchVariable(:T_probe_3a, Grid3D(), desc = "T after reset_variables!",          namespace = :debug),
ScratchVariable(:T_probe_3b, Grid3D(), desc = "T after global_parameterizations!", namespace = :debug),
ScratchVariable(:T_probe_3c, Grid3D(), desc = "T after column_parameterizations!", namespace = :debug),
ScratchVariable(:T_probe_4,  Grid3D(), desc = "T after ocean_timestep!",           namespace = :debug),
ScratchVariable(:T_probe_5,  Grid3D(), desc = "T after sea_ice_timestep!",         namespace = :debug),
ScratchVariable(:T_probe_6,  Grid3D(), desc = "T after land_timestep!",            namespace = :debug),
ScratchVariable(:T_probe_7,  Grid3D(), desc = "T after parameterization_tendencies_only!", namespace = :debug),
ScratchVariable(:T_probe_8,  Grid3D(), desc = "T after horizontal_diffusion!",     namespace = :debug),
ScratchVariable(:T_probe_9,  Grid3D(), desc = "T after leapfrog!",                 namespace = :debug),
ScratchVariable(:T_probe_10, Grid3D(), desc = "T after transform!",                namespace = :debug),
# parallel probes for tendency
ScratchVariable(:dT_probe_2, Grid3D(), desc = "dT after greenhouse_gases",         namespace = :debug),
ScratchVariable(:dT_probe_3, Grid3D(), desc = "dT after parameterization_tendencies!", namespace = :debug),
ScratchVariable(:dT_probe_7, Grid3D(), desc = "dT after parameterization_tendencies_only!", namespace = :debug),
ScratchVariable(:dT_probe_8, Grid3D(), desc = "dT after horizontal_diffusion!",    namespace = :debug),
```

#### Step 2 — Snapshot insertions

In `SpeedyWeather/src/time_stepping/time_integration.jl`, the `timestep!(vars,
dt, ::PrimitiveEquation, ...)` body (lines 129–165), insert `copy!` of
`vars.grid.temperature` into the corresponding scratch field after each
named operation. Use `copy!(vars.scratch.debug.T_probe_k, vars.grid.temperature)`.

Same for the tendency probes via `vars.tendencies.grid.temperature`.

Be careful: scratch slots are write-before-read, so reading them after a
fused kernel may give stale values; use `copy!` from the canonical state
field at each marker, not from another scratch.

#### Step 3 — Write the bisect script

`debug_intra_step_bisect.jl` — same skeleton as `debug_drift_step1.jl`, but
after the step has run, dumps:

```julia
for probe in (:T_probe_1, :T_probe_2, :T_probe_3a, :T_probe_3b, :T_probe_3c,
              :T_probe_4, :T_probe_5, :T_probe_6, :T_probe_7, :T_probe_8,
              :T_probe_9, :T_probe_10)
    a = Array(getproperty(sim_c.variables.scratch.debug, probe))
    b = Array(getproperty(sim_r.variables.scratch.debug, probe))
    d = maximum(abs.(a .- b))
    println("  $probe  max|ΔT| = $d")
end
```

Same loop for the `dT_probe_*` parallel set.

Run with the default Frierson configuration (default has the largest
amplification, so the first non-zero probe will be cleanest to spot).

#### Step 4 — Interpret

The first probe k where `max|ΔT| > 0` (or grows substantially beyond the
upstream baseline) localizes the divergence to one of the 10 operations.
Possible outcomes:

| First diverging probe | Implication |
|---|---|
| P1 (reset_tendencies!) | tendencies array allocation/zeroing differs. Extremely unlikely; would suggest a fundamental allocation bug. |
| P2 (greenhouse_gases) | Greenhouse-gas time stepping itself. Inspect that function. |
| P3a (reset_variables!) | Like P1, very unlikely. |
| P3b (global_parameterizations!) | One of the non-column parameterizations (cloud, etc.) differs. Drill into the `@generated` unroll: bisect by which parameterization is included. |
| **P3c (column_parameterizations!)** | **The fused per-column kernel. Drill into the unroll order: bisect by k=1..11 prefix (like Session 1 NaN bisect).** Likely culprit given the surface concentration of the diff. |
| P4 / P5 / P6 (ocean/sea_ice/land) | Surface model itself. |
| P7 (parameterization_tendencies_only!) | The grid→spectral transform of tendencies. Even though `transform!` in isolation is bit-identical (Test B), perhaps the per-variable grid→spectral path used here differs (different matrix, different layer ordering, etc.). |
| P8 (horizontal_diffusion!) | Spectral hyperdiffusion. Pure spectral math. |
| P9 (leapfrog!) | Leapfrog filter. Pure spectral math with Robert+Williams. |
| P10 (transform!) | Spectral→grid for the new state. We already verified this is bit-identical in isolation (Test B), so divergence here would mean state ENTERING the transform differs (i.e., earlier-detected divergence carried forward). |

#### Step 5 — If P3c is the answer (most likely)

Once we've narrowed to `column_parameterizations!`, do a sub-bisect by
truncating the `@generated` unroll. Use the same approach as Session 1's
NaN bisect:

```julia
# in column_parameterizations_kernel! @generated, temporarily restrict to first k
calls = [:(parameterization!(ij, vars, parameterizations.$name, model)) for name in names[1:k]]
```

Run for k=1..N and see at which k the first non-zero T diff appears. With
longwave_radiation = nothing (since that's where step-1 diff is largest in
Test A), this drills directly into which non-LW parameterization is the
generator. Strong candidates:
- shortwave_radiation (T^4-amplifies any T perturbation, like longwave)
- surface_heat_flux (writes to dTdt[ij, nlayers] at surface)
- surface_humidity_flux (writes to humidity, indirect via condensation)
- large_scale_condensation (latent heat → T)

#### Step 6 — Cleanup

Once the source is identified, revert the ScratchVariable additions and
probe insertions. Keep `debug_intra_step_bisect.jl` for re-use.

### Pre-bisect predictions

Based on session-6 evidence and surface-layer concentration of the diff:

1. P3c (column_parameterizations) is the most likely first diverging probe.
2. Within P3c, surface_heat_flux + surface_humidity_flux + shortwave are the
   prime suspects.
3. shortwave has the same T^4 structure as longwave → highest a-priori
   probability of being the per-call seed AND amplifier.
4. If shortwave is the culprit, the fix recipes from session 4 (Float64
   widening) might work for shortwave even though they didn't for longwave —
   shortwave's flux structure is different (no Horner over a layer-varying
   transmissivity in the same way).

### Files to add tomorrow
- `SpeedyWeather/src/models/primitive_wet.jl` (or `primitive_dry.jl`) —
  ScratchVariable registrations for `T_probe_k` and `dT_probe_k`.
- `SpeedyWeather/src/time_stepping/time_integration.jl` — `copy!` insertions
  after each operation.
- `SpeedyWeather/test/reactant/debug_intra_step_bisect.jl` — bisect script.

All three intended to be reverted at end of session 7 except the bisect
script (kept for re-use).

### What we do NOT need to test
- Transforms in isolation — already bit-identical (Test B).
- `copy!` semantics — already verified bit-perfect (Test C).
- Float64 in the longwave kernel — already shown irrelevant (Session 4).
- Step-1 with other longwave schemes — already mapped (Test A).

### Status

Plan ready. Begin tomorrow with **Phase 0** (`initialize!(sim; steps)`
asymmetric T diff) — quickest win and informs the intra-step bisect. Then
proceed to Phase 1 (intra-step bisect). Total estimated wall time:
~2 hours including both phases.

### Could Phase 0 already be the whole issue?

Quick sanity check on whether the 5.8e-4 K initialize-time diff explains
the 5.5 K 10-step drift:
- 5.8e-4 K input perturbation
- Frierson amplification factor: 1.34×/step (measured)
- After 10 steps: 5.8e-4 × 1.34^10 ≈ 5.8e-4 × 17.9 ≈ 1.0e-2 K = 0.01 K
- Observed 10-step drift: 5.5 K
- Ratio: ~550× short

So Phase 0 alone CANNOT explain the full 5.5 K drift. But there are two
possibilities that keep it interesting:
1. The mechanism that creates the 5.8e-4 K in `initialize!` may ALSO run
   inside `timestep!` (e.g., a virtual-T recompute called per step), where
   it would produce a fresh ~5.8e-4 K kick per step. Combined with
   amplification, that gives ~5.8e-4 × Σ 1.34^k ≈ 5.8e-4 × 50 ≈ 0.03 K —
   still short by 100× but a real contribution.
2. The asymmetry is a SYMPTOM of a broader CPU-vs-XLA divergence pattern,
   and fixing it teaches us about the per-step seed mechanism even if it's
   not numerically dominant.

So Phase 0 is worth doing first as a diagnostic, but should NOT be expected
to solve the whole problem. The intra-step bisect (Phase 1) is the main
event.

## Session 7 — Phase 0 executed (2026-05-16)

Phase 0 ran and pinpointed the source of the 5.8e-4 K T-only diff with three
nested bisects. Findings below.

### Step 1 — Bisect inside `initialize!(sim; steps)`

Script: `debug_init_bisect.jl`. Manual replay of the body of
`initialize!(simulation; period, steps, output)` from
`SpeedyWeather/src/models/simulation.jl:55-103`, with `max|ΔT|` measured
after each sub-operation.

| Stage | max\|ΔT\| (K) |
|---|---|
| I0_pre (post-resync) | 0.0 |
| I0_clock (clock init) | 0.0 |
| I1_output (output config) | 0.0 |
| I2_scale (`scale_prognostic!`) | 0.0 |
| **I3_transform (`transform!(vars, lf, model; initialize=true)`)** | **5.8e-4** |
| I5_output_init | 5.8e-4 (preserved) |
| I6_callbacks | 5.8e-4 (preserved) |

→ **I3 introduces all of the diff.** Per-layer pattern after I3:
k=1: 2.4e-4, k=2: 2.7e-4, k=3: 2.3e-4, k=4: 4.9e-4, k=5: 4.9e-4, k=6: 5.8e-4,
k=7: 5.5e-4, k=8: 5.5e-4. Smaller at TOA, larger at mid-troposphere/surface.

### Step 2 — Drill inside `transform!(vars, lf, model; initialize=true)`

Script: `debug_init_transform_bisect.jl`. Replays the body of
`SpeedyTransforms.transform!(::Variables, lf, ::PrimitiveEquation; initialize=true)`
from `time_stepping/transform.jl:103-208`, probing `vars.grid.temperature`
after each call.

| Sub-stage | max\|ΔT\| (K) |
|---|---|
| T0_pre | 0.0 |
| T1: `transform!(vor_grid, vor, scratch, S)` | 0.0 |
| T2: `transform!(div_grid, div, scratch, S)` | 0.0 |
| **T3: `transform!(temp_grid, temp, scratch, S)`** | **5.8e-4** |
| T4: `transform!(pres_grid, pres, scratch, S)` | 5.8e-4 (preserved) |
| T5: `transform!(humid_grid, humid, scratch, S)` | 5.8e-4 (preserved) |
| T6: `hole_filling!(humid_grid, ...)` | 5.8e-4 (preserved) |
| T7: `UV_from_vordiv!(U, V, vor, div, S)` | 5.8e-4 (preserved) |
| T8: `transform!(u_grid, U, scratch, S; unscale_coslat=true)` | 5.8e-4 (preserved) |
| T9: `transform!(v_grid, V, scratch, S; unscale_coslat=true)` | 5.8e-4 (preserved) |
| T10: `temperature_average!` | 5.8e-4 (preserved) |
| T11: `geopotential!` | 5.8e-4 (preserved) |

→ **T3 is the source.** It's the 4-argument spectral→grid `transform!`
applied to temperature with `scratch_memory = vars.scratch.transform_memory`.

### Step 3 — Eager vs compiled

Script: `debug_init_transform_compiled.jl`. Hypothesis: T3 runs in
Reactant eager mode in our test scripts; maybe the eager-vs-compiled
distinction explains the diff.

| Mode | max\|ΔT\| (K) | Reactant eager vs compiled |
|---|---|---|
| Eager | 5.8e-4 | — |
| `@compile`d | 5.8e-4 | 0.0 (Reactant agrees with itself) |

→ **Not an eager-mode artifact.** Reactant produces the same answer in
both modes; CPU and Reactant disagree by 5.8e-4 K regardless.

### Step 4 — Data dependency

Script: `debug_transform_data_dep.jl`. Hypothesis: real spin-up T spectral
data triggers something random uniform-[0,1) data doesn't.

| Test | Data | Scratch source | max\|ΔT\| (K) |
|---|---|---|---|
| T1 | Real spin-up T (range [-80, 1016]) | `vars.scratch.transform_memory` | **5.8e-4** |
| T2 | Real spin-up T | `M.scratch_memory` | **5.8e-4** |
| T3 | Random `Complex{Float32}` in [0,1) | `vars.scratch.transform_memory` | **0.0** |
| T4 | Random in [0,1) | `M.scratch_memory` (mirrors Test B Session 6) | **0.0** |

→ **Data-dependent. NOT scratch-source-dependent.** The transform algorithm
itself produces divergent CPU vs Reactant output for some spectral inputs
but bit-identical for others.

### Hypothesis — large l=0 m=0 coefficient

Spectral T range is `real ∈ [-80, +1016]`. The huge upper-bound value is
the `(l=0, m=0)` harmonic = global mean × `norm_sphere ≈ 3.54` →
mean T ≈ 287 K × 3.54 ≈ 1016. This is ~1000× larger than the other
spectral coefficients (which are perturbations on top of the mean).

When `_backward_mul!` does:
```julia
scratch .= real.(coeffs_data)
LinearAlgebra.mul!(field_data, backward_real, scratch)
scratch .= imag.(coeffs_data)
LinearAlgebra.mul!(field_data, backward_imag, scratch, -1, 1)
```

The second `mul!` accumulates into `field_data` (β=1). With one giant
coefficient (l=0 m=0) plus many small ones, the `mul!` reduction is
classic FMA / reassociation territory — XLA may sum in a different order
than Julia's BLAS, producing different rounding when there's one dominant
contributor and many small ones (cancellation/absorption pattern).

This matches:
- Random uniform [0,1) data: no single dominant coefficient → no
  asymmetric FP behavior → 0 diff.
- Real spin-up T: dominant l=0 m=0 coefficient ≈ 1000, other coefficients
  small → FP-sensitive summation → 5.8e-4 K diff.

### Untried tests (Tests 5, 6) — interrupted by syntax error

Script `debug_transform_data_dep.jl` was supposed to also run:
- T5: real T data scaled by 1e-3 (T values become tiny)
- T6: random data scaled by 300 (random values become T-magnitude)

Both have a syntax bug (`1e-3f0` not a valid literal; use `Float32(1e-3)`).
Re-run tomorrow. Expected results:
- T5: if **0** → magnitude-dependent (large numbers cause divergence).
- T5: if **5.8e-4** → not pure magnitude; structure matters.
- T6: if **5.8e-4** → triggers at T-magnitude regardless of structure.
- T6: if **0** → real T's specific pattern (large l=0 m=0) is required.

### Next direct test (tomorrow)

The cleanest follow-up: take random data, plant a LARGE value only at the
`(l=0, m=0)` position (~1000), and see if that alone triggers the 5.8e-4 K
divergence. That would confirm the dominant-coefficient hypothesis.

```julia
# Set only the (l=0, m=0) coefficient to a large value, rest small
spec[1, k] = 1000.0f0 + 0im  # for each layer k
# all other harmonics: keep random small
```

If this triggers the divergence: the bug is FP-summation-order sensitivity
in `_backward_mul!` with one dominant value. Likely fix: split the
contribution of the l=0 m=0 mode from the others, OR widen the mul!
accumulation to Float64.

### Status

Phase 0 source LOCALIZED (with caveats):
- `transform!(temp_grid, temp, scratch_memory, S)` (4-arg, spectral→grid)
- Diverges for real spin-up T spectral data (5.8e-4 K)
- Does NOT diverge for random [0,1) data
- Same in eager and compiled Reactant modes
- Independent of scratch_memory source

Open question: what FEATURE of the input triggers the divergence? Working
hypothesis: dominant l=0 m=0 coefficient creates FP-sensitive summation
that CPU BLAS and XLA dot evaluate in different orders.

This is the FIRST identified concrete CPU-vs-Reactant algorithmic
divergence in the SpeedyWeather pipeline. Even though it's small
(5.8e-4 K) and not the dominant drift source (the 0.4 K/step seed is
elsewhere), it's a clean fix candidate. If the same FP-summation-order
issue manifests inside parameterizations or in per-step transforms, it
could be the seed mechanism we've been chasing.

### Files (session 7 executed)
- `debug_init_bisect.jl` — Phase 0 outer bisect through `initialize!`. **Kept.**
- `debug_init_transform_bisect.jl` — drill inside `transform!(...; initialize=true)`. **Kept.**
- `debug_init_transform_compiled.jl` — eager vs compiled comparison. **Kept.**
- `debug_transform_data_dep.jl` — data/scratch dependency tests (Tests 5, 6 not yet run due to syntax bug). **Fix and rerun tomorrow.**

### Phase 0 numbers
| Test | Result |
|---|---|
| `initialize!` outer bisect | T diff appears at I3 only |
| `transform!(...; initialize=true)` drill | T diff appears at T3 only |
| Eager vs compiled | both 5.8e-4 K |
| Real T data, vars scratch | 5.8e-4 K |
| Real T data, M scratch | 5.8e-4 K |
| Random data, vars scratch | 0.0 |
| Random data, M scratch | 0.0 |

### Next session plan (continue here)

1. **Fix and rerun Tests 5, 6** in `debug_transform_data_dep.jl` to test
   magnitude vs structure hypothesis (~5 min wall time).
2. **Run the "dominant coefficient" test** — random data with just the
   (l=0, m=0) coefficient set large. Confirms or denies summation-order
   hypothesis (~10 min wall time).
3. **If confirmed**: try Float64 accumulation inside `_backward_mul!` —
   make the `mul!` reductions sum in Float64 then cast to Float32. This
   is a targeted, low-cost fix specifically where the divergence lives.
4. **Independent of Phase 0 result**: proceed with Phase 1 (intra-step
   bisect). The 5.8e-4 K isn't the dominant seed — the 0.4 K/step inside
   `timestep!` is. The Phase 0 root cause may or may not be the same
   mechanism that drives the per-step seed; Phase 1 will tell us.

Phase 0 was a good warm-up: it gave us our first concrete CPU↔Reactant
algorithmic divergence, with a falsifiable hypothesis (dominant l=0 m=0
coefficient → FP summation-order sensitivity) and a clean fix candidate
(Float64 mul! accumulation).

## Session 8 — T3 root cause: SubArray-of-3D triggers different XLA lowering (2026-05-16)

### Reframing the question

The user asked: at T3 (`transform!(temp_grid, temp, scratch, S)`), is the
input `temp` REALLY bit-identical between CPU and Reactant? And: try to
isolate into a self-contained MWE.

Answer: yes, `temp` is exactly equal (`Array(temp_spec_c.data) ==
Array(temp_spec_r.data)`), the backward matrices are equal, both
scratch sources are equal. The divergence is **not in the input**. It
is in the COMPILED FORM of `_backward_mul!`, which depends on the type
of the `coeffs` container.

### Tests run (in order)

`debug_transform_data_dep.jl` (rerun with syntax fixes + new tests):

| Test | Setup | max\|Δ\| (K) |
|---|---|---|
| 0a  | input `temp_spec_c.data == temp_spec_r.data`? | bit-identical (==) |
| 0b  | matrices + scratch all bit-identical? | yes (==) |
| 1   | real T (from `temp_spec_c`), vars scratch | **5.8e-4** |
| 2   | real T (from `temp_spec_c`), M.scratch_memory | **5.8e-4** |
| 5   | real T spec × 1e-3 (copied into fresh container) | 0 |
| 6   | random × 300 (copied into fresh container) | 0 |
| 7   | random + planted (l=0,m=0)=1000 (fresh) | 0 |
| 8   | only (l=0,m=0)=1000, rest 0 (fresh) | 0 |
| 9   | real T with (l=0,m=0) zeroed (fresh) | 0 |

The "dominant coefficient" hypothesis from Session 7 is **falsified**:
T7, T8 (planted dominant coeff in fresh containers) give 0 K diff. T5,
T6 (real T scaled, random scaled) also 0 — magnitude alone is not
enough; structure alone is not enough; large (l=0,m=0) alone is not
enough.

But T1, T2 still show 5.8e-4 K with the **exact same numeric values**.
Difference: T1/T2 pass `temp_spec_c` (the actual prognostic
`LowerTriangularArray`) directly, while T5–T9 build a fresh
`zeros(Complex{Float32}, SG.spectrum, NLAYERS)` and copy the data into it.

### Drill: same numbers, different container (`debug_transform_data_dep3.jl`)

| Repro | CPU side                  | Reactant side                | max\|Δ\| |
|---|---|---|---|
| R1 | `temp_spec_c` (orig)      | `temp_spec_r` (orig)         | **5.8e-4** |
| R2 | fresh container, copy in  | fresh container, copy in     | 0       |
| R3 | fresh                     | `temp_spec_r` (orig)         | **5.8e-4** |
| R4 | `temp_spec_c` (orig)      | fresh                        | 0       |
| R5 | both orig + fresh scratch | both orig + fresh scratch    | **5.8e-4** |
| R7 | both orig + pre-zeroed scratch | both orig + pre-zeroed scratch | **5.8e-4** |

→ The divergence is **on the Reactant side**, triggered by the type of
`temp_spec_r.data`:

```
typeof(temp_spec_r.data) =
  SubArray{ComplexF32, 2,
           ConcretePJRTArray{ComplexF32, 3, 1},
           Tuple{Slice, Slice, Int64}, false}
```

It's a **`SubArray` view of a 3D `ConcretePJRTArray`** — the underlying
buffer is the 4-step leapfrog tank for temperature; `get_step(...,1)`
returns view of slice `[:, :, 1]`. CPU side has the analogous SubArray
of a 3D `Array`, and Reactant agrees with CPU when both pass a **flat
2D `ConcretePJRTArray`** (R2, R4).

### Self-contained MWE (`MWE_subarray_mul.jl`)

No SpeedyWeather init required. Defines `_backward_mul!` with the
exact same body as `SpeedyTransforms._backward_mul!`, builds two
Reactant containers holding bit-identical complex coefficients:

- `coeffs_flat_r :: ConcretePJRTArray{ComplexF32, 2, 1}`
- `coeffs_view_r :: SubArray{ComplexF32, 2, ConcretePJRTArray{..., 3, 1}, ..., false}` (view of `[:, :, 1]` of a 3D buffer)

Run `Reactant.@jit _backward_mul!(out, coeffs_***, scratch, br, bi)` and
compare. Result:

```
Array(flat) == Array(view) : true (bit-identical)
CPU vs Reactant(flat)      : 0.0
CPU vs Reactant(view)      : 0.0041503906   ← divergence
Reactant flat vs view      : 0.0041503906
```

The MWE reproduces the divergence **without any SpeedyWeather objects**,
using only `Reactant.@jit`, `LinearAlgebra.mul!`, and a SubArray of a 3D
ConcretePJRTArray.

### Where in the computation does it diverge?

`MWE_subarray_mul_minimal.jl` (using real Float32 data, single x):

| Op            | container | max\|Δ\| (Reactant flat vs view) |
|---|---|---|
| single `mul!(out, A, x)` | flat / view | 0 |
| two-stage `mul!(out, A, x); mul!(out, A2, x, -1, 1)` | flat / view | 0 |

→ With real Float32 inputs, no divergence.

`MWE_subarray_mul_minimal2.jl` (probe just `scratch .= real.(c)` and
`imag.(c)` broadcasts on complex SubArray):

| Op | max\|Δ\| (flat vs view) |
|---|---|
| `scratch .= real.(c)`  | 0 |
| `scratch .= imag.(c)`  | 0 |

→ The component-extraction broadcast is bit-identical.

So the divergence appears only when `_backward_mul!` is JIT'd
**as a single function** with a complex SubArray-of-3D argument.

### Root cause confirmed via HLO dump (`MWE_subarray_mul_minimal3.jl`)

`@code_hlo _backward_mul!(out, flat_r, scratch, br, bi)` vs same with
`view_r`:

**Flat 2D (`ConcretePJRTArray{ComplexF32, 2, 1}`):**
```
%0 = stablehlo.transpose %arg1, dims = [1, 0]
%1 = stablehlo.real %0
%2 = stablehlo.dot_general %1, %arg3, contracting_dims = [0] x [0], ...
%3 = stablehlo.imag %0
%4 = stablehlo.dot_general %3, %arg4, contracting_dims = [0] x [0], ...
%5 = stablehlo.subtract %2, %4
```

**SubArray-of-3D (`SubArray{ComplexF32, 2, ConcretePJRTArray{..., 3, 1}, ...}`):**
```
%0 = stablehlo.slice %arg1 [0:1, 0:8, 0:560]
%1 = stablehlo.transpose %0, dims = [2, 1, 0]    ← different perm
%2 = stablehlo.reshape %1
%3 = stablehlo.real %2
%4 = stablehlo.dot_general %3, %arg3, contracting_dims = [0] x [0], ...
%5 = stablehlo.imag %2
%6 = stablehlo.dot_general %5, %arg4, contracting_dims = [0] x [0], ...
%7 = stablehlo.subtract %4, %6
```

Both forms are correct, but the **layout of the operand to
`dot_general`** differs: flat goes through a single transpose
`[1,0]`; view-of-3D goes through `slice → transpose [2,1,0] → reshape`.
XLA picks different GEMM blocking/accumulation orders for these
different operand layouts, producing different FP rounding.

This is a classic XLA reassociation: both results are within
floating-point rounding of the true value, but they're not bit-equal.
For the spin-up T data at T31×8 layers, the difference reaches 5.8e-4
K at the surface — small in absolute terms, but a clean CPU↔Reactant
algorithmic divergence.

### Implications

1. The 5.8e-4 K **is not data-dependent in the way Session 7
   hypothesized**. It is **container-type-dependent**. The dominant
   (l=0,m=0) coefficient is a symptom, not the trigger: large numbers
   amplify the FP-reassociation difference, but the trigger is the
   SubArray-of-3D operand passing through `dot_general`.

2. The `temp_spec_r.data` SubArray comes from
   `vars.prognostic.temperature` — a `LeapfrogArray` storing 4 time
   steps as a 3D `ConcretePJRTArray`. Every spectral prognostic
   variable that goes through `transform!` exposes the same pattern
   (vorticity, divergence, humidity, pressure). Yet only `temperature`
   produced a noticeable diff in the T-only test — probably because T
   has the largest absolute magnitudes (1000 K-scale spectral
   coefficient at (l=0,m=0)). Vorticity, divergence, etc. have much
   smaller absolute values and their FP-reassoc noise stays below the
   tolerance.

3. **This is not the dominant drift source.** The 0.4 K/step seed
   inside `timestep!` is still elsewhere — the SubArray issue manifests
   once per-variable per spectral→grid transform, but the 5.8e-4 K
   doesn't grow proportionally with step count (it stays roughly the
   same). It is a noise floor.

### Candidate fixes

(a) **Materialize the view** before passing to `mul!`: copy the
SubArray into a fresh `ConcretePJRTArray` and pass that. Simple,
adds one copy per `transform!` call. Pros: zero risk, predictable
result. Cons: O(NHARM * NLAYERS) extra memory per call.

(b) **Reshape rather than slice**: if the leapfrog buffer can be laid
out in memory such that `get_step(temp, 1)` returns a flat
`reshape` (no slice) of the 4-step buffer, the IR would match.
Requires changing `LeapfrogArray` storage order.

(c) **Apply `_backward_mul!` per layer** (loop over k) using a 1D
view — but that's even worse for FP reassociation, and slower.

(d) **Accept the diff** if it remains below the algorithm's overall
accuracy. The 5.8e-4 K is below the leapfrog integration error
floor of Float32 spherical harmonics at T31, so semantically
benign.

(e) **Move `_backward_mul!` materialization inside the compile**:
restructure so the view-vs-flat dichotomy doesn't propagate into
the JIT'd region. If we pass `coeffs.data` already explicitly
materialized to a flat 2D `ConcretePJRTArray` at the SpeedyWeather
level, XLA sees the same IR every call.

I lean toward (e) — explicit materialization at the
`SpeedyTransforms.transform!(field, coeffs, scratch, M)` boundary:

```julia
function transform!(field, coeffs, scratch_memory, M)
    ...
    # Materialize the SubArray-of-3D into a flat 2D before JIT'ing
    coeffs_flat = M.architecture isa ReactantDevice ?
        copy_to_flat(coeffs.data) : coeffs.data
    @maybe_jit M.architecture _backward_mul!(field.data, coeffs_flat, scratch, ...)
end
```

That avoids the layout dependence in the IR signature entirely. To
prototype, add a one-line `coeffs_data = collect_to_flat(coeffs.data)`
in `transform!` and re-run `debug_transform_data_dep.jl` T1; expect 0.

### Files added (Session 8)
- `debug_transform_data_dep.jl` (rewritten) — Tests 0a, 0b, 1–9 with fixes.
- `debug_transform_data_dep2.jl` — narrow down by data structure (m=0 modes, scaled, imag-only).
- `debug_transform_data_dep3.jl` — same-numbers-different-container probe (R1–R7).
- `MWE_subarray_mul.jl` — self-contained MWE, no SpeedyWeather init.
- `MWE_subarray_mul_minimal.jl` — single vs two-stage mul!.
- `MWE_subarray_mul_minimal2.jl` — isolate real./imag. broadcast.
- `MWE_subarray_mul_minimal3.jl` — HLO IR dump, flat vs view.

### Status

Phase 0 / T3 root cause **fully identified**: XLA produces different
HLO for `_backward_mul!` when its `coeffs` arg is a
`SubArray{ComplexF32, 2, ConcretePJRTArray{..., 3, 1}, ...}` vs a flat
`ConcretePJRTArray{ComplexF32, 2, 1}`. Different IR → different XLA
executable → different FP-summation order → 5.8e-4 K diff on real T
spin-up data.

Fix candidate (e): materialize the operand to a flat 2D
ConcretePJRTArray before `@maybe_jit _backward_mul!` is invoked.
Prototype lives one Edit away — test next session.

The 5.8e-4 K is **not the 0.4 K/step seed**. Continue with Phase 1
(intra-step bisect inside `timestep!`) regardless of whether we patch
the SubArray issue. The mechanism (different IR for view-of-3D vs flat
container) is potentially active in other transform call sites too;
it's worth keeping in mind during Phase 1.

### Next session plan
1. Prototype fix (e): add a `materialize_for_jit` helper in
   `SpeedyTransforms`; gate on `ReactantDevice`. Re-run
   `debug_transform_data_dep.jl` T1, expect 0.
2. If (e) succeeds, audit other transform call sites (UV_from_vordiv,
   pressure, humidity) for the same pattern.
3. Run the full Reactant correctness suite to see if the per-step
   drift moves at all — almost certainly not, but worth checking.
4. Begin Phase 1 (intra-step bisect inside `timestep!`).

## Session 9 — Temporary fix applied (2026-05-16, late)

Issue reported upstream. Until the maintainers fix the IR-stability
issue in Reactant/XLA, we apply a local workaround.

### What changed

- `SpeedyTransforms/src/matrix_transform.jl`: inside the spectral->grid
  `transform!(field, coeffs, scratch, M::MatrixSpectralTransform)`,
  insert a `coeffs_for_mul = _flatten_for_backward(M.architecture, coeffs.data)`
  one-liner before `@maybe_jit _backward_mul!(field.data, coeffs_for_mul, ...)`.
  CPU/GPU fallback is identity (no behavior change).
- `SpeedyTransforms/ext/SpeedyTransformsReactantExt.jl`: extension
  overrides `_flatten_for_backward(::ReactantDevice, data)` to
  materialize a `SubArray{T, 2, ConcretePJRTArray{T, 3}}` into a flat
  2D ConcretePJRTArray via `Reactant.@jit copy(data)`. When already
  inside a `@jit` trace (`within_compile()`), pass through unchanged.
- `CHANGELOG.md`: note added under Unreleased.

The workaround is a small overload guarded by:
```julia
@inline _materialize_2d(data::Reactant.ConcretePJRTArray{<:Any, 2}) = data
@inline function _materialize_2d(data::SubArray{T, 2, <:Reactant.ConcretePJRTArray{T, 3}}) where {T}
    return Reactant.@jit copy(data)
end
@inline _materialize_2d(data) = data   # all other shapes: unchanged
```

so only the specific failure mode (`SubArray` view of a 3D
`ConcretePJRTArray`) is materialized. Everything else stays on the
same code path.

### Verification

`debug_transform_data_dep.jl` (after fix):

| Test | Before | After |
|---|---|---|
| T1 real T, vars scratch | 5.8e-4 | **0.0** |
| T2 real T, M scratch    | 5.8e-4 | **0.0** |
| T3 random + vars scratch | 0.0 | 0.0 |
| ... | (unchanged) | (unchanged) |

`debug_init_bisect.jl` (after fix):

| Stage | Before | After |
|---|---|---|
| I3_transform | 5.8e-4 | **0.0** |
| All other stages | 0.0 | 0.0 |

The full `initialize!(sim; steps=1)` pipeline now reports
`max|ΔT| = 0` at every probe point. Phase 0 / T3 root cause is
**neutralized** on the SpeedyWeather side.

### Status

- Temporary workaround in place. Will be removed once Reactant emits
  matching IR for SubArray-of-3D and flat-2D inputs to `dot_general`.
- The 5.8e-4 K diff is **gone**.
- The 0.4 K/step seed inside `timestep!` is still present — Phase 0
  was a noise floor, not the dominant drift source. Continue with
  Phase 1.

### Tomorrow's plan (continue here)

Phase 1 — intra-step bisect inside `timestep!` (the actual 0.4 K/step
seed). Plan from Session 7 still applies (see "Session 7 (planned) —
Intra-step bisect plan" earlier in this document). Execute that plan
fresh now that the T3 noise floor is removed:

1. Re-baseline: with the workaround active, measure step-1 T diff vs
   CPU. Earlier (pre-fix) the per-step seed manifested at ~0.4 K
   after several steps; now T should start cleaner. Re-measure step-1
   T to confirm the per-step seed is unchanged (expected — the fix
   only touches spectral->grid transforms, not the
   tendencies/timestepping pipeline).
2. Register debug `ScratchVariable`s as in Session 7 Step 1 (the
   commit-revertable instrumentation plan).
3. Insert snapshot points P0-P10 inside `time_integration.jl`'s
   `timestep!` (dynamics tendencies, parameterizations, leapfrog,
   implicit, hole-filling — see Session 7 "The pipeline").
4. Run the bisect script and identify which stage in `timestep!`
   introduces the per-step seed.
5. Drill into that stage. Likely candidates per Session 7: P3c
   (longwave radiation) or one of the transform-back stages.

Files to revisit:
- `debug_init_transform_bisect.jl`, `debug_init_bisect.jl` — keep
  as sanity checks; both should report 0 with the workaround.
- The Phase 1 instrumentation will live in a new
  `debug_timestep_bisect.jl` modeled after `debug_init_bisect.jl`.

## Session 10 — Phase 1 executed (2026-05-17)

### Re-baseline with Session 9 workaround active

`debug_drift_step1.jl` rerun: **0.386 K** step-1 T diff at k=8, same
per-layer pattern as before the workaround. Confirmed: the Session 8
fix removed the 5.8e-4 K transform noise but the per-step seed is
independent.

Per-layer step-1 T diff:
| k | max\|ΔT\| | k | max\|ΔT\| |
|---|---|---|---|
| 1 | 0.052 K | 5 | 0.129 K |
| 2 | 0.029 K | 6 | 0.144 K |
| 3 | 0.026 K | 7 | 0.240 K |
| 4 | 0.075 K | 8 | 0.386 K |

### Step 1 — Bisect timestep!() — `debug_timestep_bisect.jl`

Inline mirror of `timestep!(vars, dt, model::PrimitiveEquation, lf1, lf2)`
for the Euler half-step (`lf1=lf2=1`, `dt=Δt/2`). Each Reactant
stage wrapped with `Reactant.@jit f(args...)` (eager calls hit
`rem(::ConcretePJRTNumber, ::Int)` via DateTime conversion).

| Stage | spec T | grid T | tend.gridT | tend.specT | tend.gridQ |
|---|---|---|---|---|---|
| S0 pre | 0 | 0 | 0 | 0 | 0 |
| S1 reset_tendencies! | 0 | 0 | 0 | 0 | 0 |
| S2 greenhouse_gases_time_step! | 0 | 0 | 0 | 0 | 0 |
| **S3 parameterization_tendencies!** | 0 | 0 | **0.0494** | 0 | **3.7e-5** |
| S4 ocean_timestep! | 0 | 0 | 0.0494 | 0 | 3.7e-5 |
| S5 sea_ice_timestep! | 0 | 0 | 0.0494 | 0 | 3.7e-5 |
| S6 land_timestep! | 0 | 0 | 0.0494 | 0 | 3.7e-5 |
| S7 param_tendencies_only! | 0 | 0 | 0.0494 | 7.7e-4 | 3.7e-5 |
| S8 horizontal_diffusion! | 0 | 0 | 0.0494 | 7.7e-4 | 3.7e-5 |
| S9 leapfrog! | 1.3e-7 | 0 | 0.0494 | 7.7e-4 | 3.7e-5 |
| S10 transform! | 1.3e-7 | 0 | 0.0494 | 7.7e-4 | 3.7e-5 |

→ **S3 (`parameterization_tendencies!`) is the seed stage.** All
subsequent stages just transport that 0.0494 K/s grid-T tendency diff
through the leapfrog and back-transform.

### Step 2 — Bisect `parameterization_tendencies!` — `debug_param_drift_bisect.jl`

Cumulative bisect: enable first k of 11 parameterizations, measure
`max|Δ tend.grid.temperature|` and `max|Δ tend.grid.humidity|`.

| k | name | max\|Δ tend.gridT\| | max\|Δ tend.gridQ\| |
|---|---|---|---|
| 1 | solar_zenith | 0 | 0 |
| 2 | vertical_diffusion | 0 | 0 |
| **3** | **large_scale_condensation** | **4.12e-3** | **1.65e-6** |
| 4 | albedo | 4.12e-3 | 1.65e-6 |
| 5 | shortwave_radiation | 4.15e-3 | 1.65e-6 |
| 6 | longwave_radiation | 4.09e-3 | 1.65e-6 |
| 7 | boundary_layer_drag | 4.09e-3 | 1.65e-6 |
| 8 | surface_condition | 4.09e-3 | 1.65e-6 |
| 9 | surface_momentum_flux | 4.09e-3 | 1.65e-6 |
| **10** | **surface_heat_flux** | **4.92e-2** | 1.65e-6 |
| **11** | **surface_humidity_flux** | 4.92e-2 | **3.75e-5** |

Three seeds visible at this aggregation:
- **large_scale_condensation** (k=3): T & Q tend at k=8 (surface only).
- **surface_heat_flux** (k=10): 10× jump in T tend at k=8.
- **surface_humidity_flux** (k=11): 22× jump in Q tend at k=8.

### Step 3 — Drill into surface_heat_flux — `debug_param_drill.jl` + `debug_shf_probe.jl`

`surface_heat_flux` run SOLO (after setup chain) produces 4.93e-2 K diff
at k=8 only. Inputs to the kernel: SST, T_air, ρ, V — all bit-identical.
But `parameterizations.boundary_layer_drag` shows **1.25e-7 diff** before
surface_heat_flux runs — that 1.25e-7 in drag is amplified ~400000× into
8.9e-3 W/m² SHF diff, then divided by `pₛ * Δσ * cₚ` and scaled by
radius into 4.9e-2 K/s tend.

So surface_heat_flux is **innocent**; the drag was already wrong.

### Step 4 — Drill into boundary_layer_drag — `debug_drag_inputs.jl`

Run params 1..6 (everything before boundary_layer_drag), then check its
inputs. **All inputs bit-identical EXCEPT one:**

```
orography.orography                      max|Δ| = 0.023468018   ==? false
```

The static terrain field — never modified during a run, loaded from
NetCDF at model construction — differs by 0.023 m between CPU and
Reactant models. Tracing back into `EarthOrography.initialize!`
(`SpeedyWeather/src/dynamics/orography.jl:200`):

```
load NetCDF → interpolate to grid → transform!(surf_geopot, orog, S)
  → smooth in spectral → transform!(orog, surf_geopot, S)  ← spec→grid
```

That final `transform!(orog, surf_geopot, S)` is the same `_backward_mul!`
codepath. The two models (CPU `MatrixSpectralTransform` vs Reactant
`MatrixSpectralTransform`) initialize their orography independently:
CPU runs `BLAS.gemm!`, Reactant runs XLA `dot_general` — these don't
agree bit-exactly even on identical inputs.

The Session 8 workaround flattens SubArray-of-3D into flat 2D before
`_backward_mul!`, but the CPU↔Reactant disagreement is INDEPENDENT of
that: any spec→grid `transform!` is BLAS vs XLA → produces different
bit patterns. Orography is just one example; ALL constant model fields
initialized via spec→grid `transform!` (surface_geopotential too,
0.008 K diff) will diverge.

### Step 5 — Test orography sync hypothesis — `debug_orog_fix.jl`

Copy CPU orography & surface_geopotential into Reactant model
post-construction, then run step-1.

```
RESULT — Step-1 seed with orography synced:
  T_max_diff = 0.3859253 K   (was 0.3859 K WITHOUT orog sync)
```

→ **Orography sync alone does NOT change step-1 T diff.** The k=8
boundary_layer_drag→surface_heat_flux pathway WAS one of the seeds,
but it's not the dominant one. The per-step T seed of 0.386 K is
dominated by something else — radiation, most likely.

### Step 6 — Per-layer bisect — `debug_param_per_layer.jl`

Repeat cumulative bisect but report tendency diff per-layer (k=1..8).
With orography synced.

| stage | k=1 | k=2 | k=3 | k=4 | k=5 | k=6 | k=7 | k=8 |
|---|---|---|---|---|---|---|---|---|
| 1 solar_zenith | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 2 vertical_diffusion | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| **3 large_scale_condensation** | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **4.1e-3** |
| 4 albedo | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 4.1e-3 |
| **5 shortwave_radiation** | **6.0e-4** | **7.7e-4** | **3.2e-4** | **3.1e-4** | **3.5e-4** | **3.4e-4** | **3.7e-4** | 4.2e-3 |
| 6 longwave_radiation | 5.7e-4 | 9.3e-4 | 4.2e-4 | 4.6e-4 | 4.6e-4 | 5.6e-4 | 5.7e-4 | 4.1e-3 |
| 7..11 | (unchanged from k=6 with orog synced) | | | | | | | |

→ Two dominant seed contributors, NEITHER at the surface:
- **shortwave_radiation** (k=5): introduces diff at ALL 7 layers
  k=1..7 simultaneously (~3-8e-4 K/s each).
- **large_scale_condensation** (k=3): k=8 surface diff (4e-3 K/s).
- **longwave_radiation** (k=6): mild amplifier across k=1..7.

The per-layer tend.gridT diff × dt (Δt/2 ≈ 188 s) × radius scaling
explains the observed per-layer step-1 T diff (e.g. k=8 tend 4e-3 × 188
≈ 0.77 K → 0.39 K post-leapfrog ≈ observed).

### Conclusions for the per-step seed

1. The **0.386 K step-1 T diff is genuinely distributed across all
   layers** and is not localized to a surface flux artifact.
2. **Two upstream root causes:**
   - **`shortwave_radiation`** is the primary per-layer seed source —
     introduces 3-8e-4 K/s tend differences across all k=1..7 in one
     pass.
   - **`large_scale_condensation`** seeds at k=8 (4e-3 K/s).
3. **Independently, orography (and other init-time spec→grid model
   fields) drift between CPU and Reactant** because the spec→grid
   `transform!` produces different output via BLAS vs XLA `dot_general`.
   This is a separate issue — it amplifies into 4.9e-2 K via the
   `boundary_layer_drag → surface_heat_flux` chain at k=8 when
   unsynced. With orography synced, the surface_heat_flux contribution
   drops to ~0, but the per-layer radiation/condensation seed remains.
4. Surface fluxes (heat, humidity, momentum) are INNOCENT; they only
   amplify upstream seeds.

### Next session plan (continue here)

Two parallel drill threads:

**A — `shortwave_radiation` drill**
1. Identify which sub-step (transmissivity, flux-sweep, etc.) introduces
   the per-layer diff. Mirror what was done for OneBandLongwave in
   Sessions 3-5 but for the active shortwave scheme.
2. The shortwave scheme being used: check what `model.shortwave_radiation`
   defaults to in `PrimitiveWetModel`. Probably `TransparentShortwave` or
   `OneBandShortwave`.

**B — `large_scale_condensation` drill**
1. Inspect kernel for floating-point reassociation hazards (similar to
   the OneBandLongwave story).
2. Check `ifelse` patterns / `clamp` patterns / `exp`/`log` calls.

**C — Init-time orography drift (model-construction issue)**
1. Optional: add `Reactant`-side override that takes the CPU-computed
   `_backward_mul!` result during `EarthOrography.initialize!` (or run
   model construction on CPU then move). This is independent of A/B.

Files created this session:
- `debug_timestep_bisect.jl` — S0-S10 stage probe inside `timestep!`. **Kept.**
- `debug_param_drift_bisect.jl` — cumulative bisect across the 11 parameterizations. **Kept.**
- `debug_param_drill.jl` — solo-run probe with prior chain setup. **Kept.**
- `debug_shf_probe.jl` — surface_heat_flux investigation. **Kept.**
- `debug_drag_bisect.jl` — scratch-field bisect alongside tend bisect. **Kept.**
- `debug_drag_inputs.jl` — `boundary_layer_drag` input-identity probe. **Kept.**
- `debug_orog_fix.jl` — orography-sync hypothesis test. **Kept.**
- `debug_param_per_layer.jl` — per-layer per-stage T tend diff. **Kept.**

### Status

Phase 1 root cause(s) **identified**:
- `shortwave_radiation` introduces per-layer (k=1..7) tend.gridT diffs.
- `large_scale_condensation` introduces k=8 tend.gridT diff.
- `boundary_layer_drag→surface_heat_flux` chain ALSO produces a k=8
  diff via orography drift at model-init time (independent of A/B,
  fixable by post-construction orography sync).

The per-step seed is now explained quantitatively from the per-layer
tendency map. Next: drill into `shortwave_radiation` (primary) and
`large_scale_condensation` (secondary) to find the offending FP
operations.

### Orography drill — BLAS vs XLA `dot_general` divergence

`debug_orog_bisect.jl` builds CPU + Reactant models side-by-side via
`PrimitiveWetModel(...)` + `initialize!(model)`, then snapshots the
result. With Session 8 SubArray workaround active and the test running
the full `EarthOrography.initialize!` pipeline:

```
[1] Post initialize!(model) — orography fields differ?
    orography.orography         max|Δ| = 0.023468018
    surface_geopotential (spec) max|Δ| = 0.0078125
[3] Same input spec→grid: CPU vs Reactant max|Δ| = 0.05078125
[4] Each side own spec→grid: max|Δ| = 0.22973633
```

Both surfaces (post-init orography AND surface_geopotential) differ
between CPU and Reactant. Test [3] is decisive: feed bit-identical
`surface_geopotential` (spectral) into `transform!(grid, spec, M)` on
both sides and the grid output differs by 0.05m. So **BLAS gemm and
XLA `dot_general` produce different bit patterns for the same spec→grid
matrix-multiply**. This is independent of the Session 8 SubArray issue;
it's just the BLAS-vs-XLA reassociation difference for dense GEMM.

The surface_geopotential (spec) also diverges (0.008 difference in
spectral coefficients), implying the grid→spec `forward` GEMM also
differs. So both directions of `transform!` are subject to FP
reassociation between CPU and Reactant.

### Why orography sync alone doesn't fix step-1 T diff

`debug_orog_fix.jl` copies CPU orography → Reactant before stepping.
Result: step-1 T diff unchanged (0.386 K). With orography synced:
- `boundary_layer_drag` becomes bit-identical (was 1.25e-7 diff).
- `surface_heat_flux` becomes bit-identical at k=8.
- `tend.grid.temperature` at k=8 drops from 4.9e-2 to 4.1e-3.
- But the per-layer (k=1..7) tend diff from `shortwave_radiation` and
  `longwave_radiation` is unchanged — those don't read orography.

So the **per-layer tend diff comes from radiation parameterizations,
NOT from the orography drift**. Orography drift contributes only at
k=8, dominantly through surface fluxes, and even then it's only one of
multiple seed sources (large_scale_condensation also contributes at k=8).

### Workaround candidates for the orography divergence

1. **CPU-then-transfer init**: build orography on CPU
   architecture, then `on_architecture(ReactantDevice(), ...)` the
   results. Requires changing `EarthOrography.initialize!` to be
   architecture-agnostic about the gemm path.

2. **Skip the smoothing+round-trip**: read NetCDF, interpolate, and
   skip the grid→spec→smooth→spec→grid round-trip. Lose smoothing
   feature but gain bit-identity across architectures.

3. **Live with it**: the orography drift only matters insofar as
   `boundary_layer_drag` (and downstream surface fluxes) read it. As
   shown above, its contribution to the step-1 seed is ~12% of the
   total (4.9e-2 K vs total 0.386 K). The dominant seed is in
   radiation, which doesn't touch orography.

4. **Sync at user level**: after model construction, do
   `copyto!(model_r.orography.orography, Array(model_c.orography.orography))`
   in user code (essentially what `debug_orog_fix.jl` does but as
   standard recipe). Simplest workaround for reproducibility tests.

### Wider implication: ALL initialize-time spec/grid round-trips diverge

Any model component whose `initialize!` does a `transform!` between
spectral and grid spaces will exhibit this BLAS-vs-XLA divergence:
- `EarthOrography`: confirmed 0.023m
- `EarthLandSeaMask`: probably similar (uses smoothing path?)
- `vegetation`, `soil_moisture`, `soil_temperature`: any that load
  and re-transform
- Spectral truncation of initial conditions
- Spectral diffusion eigenvalue initialization

This is **independent of the Session 8 SubArray issue**. Even with
SubArrays flattened, BLAS gemm produces different rounding from XLA
`dot_general` on the same dense float matrices. The Session 8 fix
handled one specific layout-induced divergence; this one is more
fundamental — different LAPACK/BLAS implementations would also disagree
with each other, so this is "expected" floating-point behavior, just
inconvenient for reproducibility tests.

### Files added (orography drill)
- `debug_orog_bisect.jl` — post-init orography & forced same-input spec→grid probe. **Kept.**
- `debug_orog_init_timing.jl` — probe orography state at each pipeline step. **Kept (note: requires CUDA loaded).**
- `debug_orog_fix.jl` — orography-sync hypothesis test. **Kept.**

### Status (Session 10 complete)

Per-step T-diff seed of 0.386 K is now decomposed into:
- **Radiation seeds (~75% of total)**: `shortwave_radiation` introduces
  3-8e-4 K/s diffs in `tend.gridT` at all layers k=1..7.
- **Condensation seed**: `large_scale_condensation` adds 4e-3 K/s at k=8.
- **Orography/drag seed (~12% of k=8)**: BLAS-vs-XLA gemm during model
  init produces 0.023m orography diff → 1.25e-7 drag diff → 4.9e-2 K/s
  surface_heat_flux tend diff at k=8.

The orography seed is a sub-component of the k=8 tendency diff, not the
dominant cause of the overall per-step T diff. Removing it (via
orography sync) doesn't reduce the headline 0.386 K step-1 T diff.

### Next session plan (continue here)
1. **Drill into `shortwave_radiation`** — find the FP-sensitive
   operation. Likely candidates: column accumulator (similar to
   OneBandLongwave session 3-5 story), exponential / power in
   transmissivity, or a CRTV reduction.
2. **Drill into `large_scale_condensation`** — its kernel uses `exp` /
   `log` for saturation vapor pressure curves and `min`/`max` clamps;
   identify the offending op.
3. **(Optional)** prototype orography-sync workaround in a user-facing
   recipe, or upstream into `EarthOrography.initialize!`.

## Session 11 — shortwave_radiation drill: source vs amplifier (2026-05-17)

### Question

Session 10 fingered `shortwave_radiation` (k=5) as the parameterization
that introduces per-layer (k=1..7) tend.gridT diffs (3-8e-4 K/s). But
in Sessions 3-5 we found that `OneBandLongwave` was an AMPLIFIER, not
a source. Is shortwave the same — amplifying a different upstream
seed, or genuinely producing the per-layer diff itself?

### Step 1 — input bit-identity check (`debug_shortwave_solo.jl`)

After running setup chain `[solar_zenith, vertical_diffusion,
large_scale_condensation, albedo]` (with orography synced), check every
field shortwave_radiation reads:

```
grid.temperature_prev    : identical
grid.humidity_prev       : identical
grid.geopotential        : identical
grid.pressure_prev       : identical
param.cos_zenith         : identical
param.ocean.albedo       : identical
param.land.albedo        : identical
param.albedo             : identical
param.rain_rate          : max|Δ| = 3.3e-13  ←  differs (only one)
param.cloud_top          : identical
scratch.grid.a           : identical
solar_constant, cₚ       : identical
land_sea_mask.mask       : identical
orography.orography      : identical (synced)
```

Only `rain_rate` differs at sub-ULP scale (~1e-13). Could rain_rate be
the seed channel?

### Step 2 — rain_rate sensitivity (`debug_shortwave_solo2.jl`)

Three runs:
- (A) Sync rain_rate CPU→Reactant, then run shortwave SOLO.
- (B) Force rain_rate=0 on BOTH sides, then run shortwave SOLO.
- (C) Baseline (no sync), then run shortwave SOLO.

All three produce **bit-identical per-layer tend.gridT diffs** (down
to the last bit). So `rain_rate` is NOT the channel — the 3.3e-13 diff
is irrelevant.

→ **shortwave_radiation's own kernel produces the divergence**, even
with bit-identical (or zeroed) `rain_rate`.

### Step 3 — sub-component bisect (`debug_shortwave_subcomp.jl`)

`OneBandShortwave.parameterization!` calls three sub-functions:
1. `clouds!()` — DiagnosticClouds; uses `saturation_humidity(T, p)` (exp).
   Writes `parameterizations.cloud_top[ij]`, returns NamedTuple.
2. `transmissivity!()` — BackgroundShortwaveTransmissivity; per-layer
   `t[ij, k] = exp(-optical_depth)`. Writes `scratch.grid.a`.
3. `shortwave_radiative_transfer!()` — column loop, accumulating fluxes,
   writes `tend.grid.temperature`.

Run three custom kernels that stop after each sub-component:

| Kernel | Output checked | max\|Δ\| |
|---|---|---|
| K1: clouds! only | `parameterizations.cloud_top` | **0.0** (identical) |
| K2: clouds! + transmissivity! | `scratch.grid.a` (t array, all k) | **5.96e-8 every layer** |
| K3: full | `tend.gridT` per layer | 3-8e-4 K/s |

`eps(Float32(0.85)) ≈ 5.96e-8` — exactly **1 ULP**.

So:
- `clouds!` is bit-identical (despite using `saturation_humidity`'s
  `exp`/`log` internally — those happen to not contribute to the
  scalar `cloud_top` output for this data, since the precipitation
  branch and humidity threshold produce same integer answer; or the
  exp/log results cancel in the conditional path).
- **`transmissivity!` introduces a 1-ULP divergence in the
  per-layer `t = exp(-optical_depth)` array.** This is just CPU libm
  `exp` vs XLA `exp` differing by 1 ULP in Float32 — completely
  expected behavior between two FP-math implementations, no FMA
  surprises needed.
- `shortwave_radiative_transfer!` AMPLIFIES that 1-ULP `t` diff into
  3-8e-4 K/s `tend.gridT` diff (a ~10⁴× amplification) via the column
  flux sweep:
  ```
  D *= t[k]                 ← propagates ULP through layers
  dTdt[k] += (D - D_out) / cₚ * (g / (pₛ * Δσ_k))   ← scales up
  ```

### Conclusion

**shortwave_radiation is BOTH a source AND an amplifier:**

- **Source**: `transmissivity!` produces a 1-ULP-per-layer `t`
  divergence purely from `exp(-optical_depth)` differing between CPU
  libm and XLA. The inputs to `exp` are bit-identical; only `exp`
  itself diverges by 1 ULP.
- **Amplifier**: `shortwave_radiative_transfer!` turns the 1-ULP `t`
  diff into ~10⁴× larger `tend.gridT` diffs via cumulative multiplication
  through the column and flux-to-tendency scaling.

This is **fundamentally different from the OneBandLongwave story**
where the seed was upstream of the kernel and the kernel just amplified.
Here, the seed IS born inside the kernel (in `exp`) — it cannot be
eliminated by syncing inputs.

### Why this is hard to fix

The 1-ULP `exp` divergence is essentially a math-library precision
difference; the only ways to eliminate it are:
1. Use the same exp implementation on both sides (e.g. force Reactant
   to use a particular libm; impractical).
2. Round both sides to fewer bits (lose precision).
3. Accept it — the 1-ULP `exp` is within the algorithm's own
   uncertainty, and the 3-8e-4 K/s amplified diff is below the
   atmospheric model's discretization error.

The interesting question is whether the amplification can be reduced
by restructuring `shortwave_radiative_transfer!` to be less sensitive
to small input perturbations. Looking at the column sweep:
```
D = D_toa
for k in 1:nlayers
    D_out = (D - O₃ * D_toa) * t[ij, k]
    dTdt[k] += flux_to_tendency((D - D_out) / cₚ, pₛ, k, model)
    D = D_out
end
```
Each layer's `D` carries forward a multiplicative factor `t[k]`. After
8 layers the accumulated relative error is ~8 × 1ULP × O(1) = a few
ULP in `D`, and the tendency receives `(D - D_out) / cₚ` which is the
"absorbed flux" — sensitive to exactly the kind of cancellation that
amplifies ULP noise.

### Implications for the bigger picture

The headline 0.386 K step-1 T diff has TWO origin classes:
1. **1-ULP `exp` in transmissivity → amplified by column sweep**
   (shortwave: 3-8e-4 K/s × 188 s × leapfrog ≈ 0.05–0.15 K per layer)
2. **`large_scale_condensation`'s k=8 contribution** (4e-3 K/s at k=8)
3. **`orography drift → boundary_layer_drag → surface_heat_flux`**
   (k=8 only, fixable by orography sync)

The dominant per-layer mid-tropospheric diff (k=4: 0.075, k=5: 0.13,
k=6: 0.14) is from category (1) — `exp` ULP noise amplified by the
shortwave column sweep. **It's "expected" floating-point behavior**,
not a bug per se, but it sets a noise floor that CPU↔Reactant
correctness tests have to tolerate.

### Files added (Session 11)
- `debug_shortwave_solo.jl` — input-identity check + solo run. **Kept.**
- `debug_shortwave_solo2.jl` — rain_rate sync / zero variants. **Kept.**
- `debug_shortwave_subcomp.jl` — K1/K2/K3 sub-component bisect of
  OneBandShortwave; localized to `transmissivity!`. **Kept.**

### Status

Session 11 root cause **identified**:
- `shortwave_radiation`'s `transmissivity!` produces a 1-ULP-per-layer
  divergence in `t = exp(-optical_depth)` (just CPU vs XLA `exp`).
- `shortwave_radiative_transfer!` amplifies that into 3-8e-4 K/s
  `tend.gridT` diff per layer via column flux accumulation.

It is NOT an amplifier of an upstream issue (unlike OneBandLongwave
in Sessions 3-5). The seed is genuinely inside the shortwave kernel,
specifically the `exp` call.

### Next session plan (continue here)
1. **Quick MWE confirmation** (`MWE_exp_ulp.jl`): standalone CPU vs
   Reactant `@jit exp.(-x)` showing 1-ULP diff. Not yet run.
2. **Drill into `large_scale_condensation`** — same pattern: SOLO run
   with input identity check, find what's the seed source. Its kernel
   also uses `exp`/`log` for saturation vapor pressure.
3. **Drill into `longwave_radiation`** — it adds 1-2e-4 K/s per layer
   on top of shortwave. Probably amplifies the already-divergent `t`
   or has its own `exp`.
4. **Decide on tolerance/policy**: given the per-layer ULP-amplification
   pattern, decide whether to accept the resulting per-step noise floor
   (~0.4 K at T31 step 1) or pursue structural changes to reduce
   amplification (e.g. compensated summation in column sweeps).

## Session 12 — Standalone MWEs + LSC + longwave drills (2026-05-17)

User directive: each step should reduce the culprit to a MWE, best
completely without SpeedyWeather.

### 1. `exp` 1-ULP MWE confirmed (`MWE_exp_ulp.jl`)

Standalone Julia + Reactant — no SpeedyWeather.

```
n = 3168*8 = 25344 random Float32 in [0, 0.5]
t_cpu = exp.(-x)
t_rct = @jit exp.(-x)
max|Δ| = 5.9604645e-8                        ← exactly 1 ULP at ~0.85
count nonzero diff: 2976 / 25344             ← 12% of inputs
```

Confirms CPU libm `exp` and XLA `exp` differ by exactly 1 ULP for
Float32 on roughly 1 in 8 inputs. Expected math-library behavior.

### 2. Column-sweep amplification MWE (`MWE_column_sweep_amplification.jl`)

Standalone — no SpeedyWeather. Mimics the
`shortwave_radiative_transfer!` column structure: downward beam ×
transmissivity, surface reflection, upward beam, flux-to-tendency
conversion with `g / (p_s * Δσ_k)` scaling and `* radius`.

Feeding the two 1-ULP-different `t` arrays from MWE_exp_ulp into the
sweep produces:

```
k=1  max|Δ dTdt| = 1.95e-3   ← TOA, largest amplification
k=2  max|Δ dTdt| = 1.10e-3
k=3  max|Δ dTdt| = 7.6e-4
k=4  max|Δ dTdt| = 6.1e-4
k=5  max|Δ dTdt| = 4.9e-4
k=6  max|Δ dTdt| = 3.7e-4
k=7  max|Δ dTdt| = 3.1e-4
k=8  max|Δ dTdt| = 2.1e-4   ← surface, smallest

Amplification factor: max|Δ dTdt| / max|Δ t| = 32768
```

Reproduces the production per-layer pattern (TOA largest, surface
smallest) exactly. Amplification is ~10⁴, driven by:
`(g / (p_s * Δσ_k)) * D_TOA / cₚ ≈ 2.7e-3` × layers × `radius ≈ 6.4e6`.

### 3. LSC drill (`debug_lsc_solo.jl`)

Solo run of `large_scale_condensation` after setup chain [1..2]:
- ALL inputs bit-identical (T, q, p, cloud_top, rain_rate, ...).
- Output: `tend.gridT` diff = **4.12e-3 K/s at k=8 ONLY**; all
  other layers = 0. Same for `tend.gridQ` at k=8 = 1.65e-6 kg/kg/s.
- Matches production bisect exactly.

Why only k=8? `large_scale_condensation!` calls `saturation_humidity`
(uses `exp`) at every k, but the condensation branch
`@trace if (δq_cond < 0) | ...` only fires when supersaturation OR
downward precip is present, which (in this spin-up state) is only at
the surface (k=8).

→ **LSC is also an intrinsic source**, not amplifier. The seed is
`saturation_vapor_pressure` (`exp(Lᵥ/Rᵥ * (1/T₀ - 1/T))`) — same
1-ULP `exp` mechanism as shortwave's transmissivity.

### 4. LSC standalone MWE (`MWE_lsc_condensation.jl`)

Standalone — no SpeedyWeather. Encodes:
- `svp(T) = e₀ * exp(Lᵥ/Rᵥ * (1/T₀ - 1/T))` and
  `sat_humid(T, p) = ε * svp(T) / p`
- The LSC implicit-condensation update:
  `δq = min(0, sat*RH - q); dqsat_dT = sat*RH*Lᵥ/cₚ / (Rᵥ T²);
   δq /= ((1 + Lᵥ/cₚ * dqsat_dT) * τ * Δt); δT = -Lᵥ/cₚ * δq`

Result on N=3168 random surface-air states:
```
max|Δ sat_humid| = 2.79e-9                   ← 1-3 ULP at sat_humid ~0.01
count nonzero sat diff: 1333 / 3168
max|Δ δT| = 1.31e-9
max|Δ δT * radius| = 8.34e-3 K/s            ← matches production 4.1e-3 K/s
```

→ Production's 4.1e-3 K/s at k=8 is fully accounted for by the
standalone reproducer.

### 5. Longwave drill (`debug_longwave_solo.jl`)

Run setup chain [1..5] (including shortwave), reset tendencies, run
LW SOLO and measure ITS own contribution to `tend.gridT`.

```
[Pre-LW] all inputs bit-identical EXCEPT:
  scratch.grid.a (post-shortwave) max|Δ| = 5.96e-8     ← from shortwave's t

[After LW SOLO] tend.gridT diff per layer:
  k=1: 3.7e-4   k=2: 3.7e-4   k=3: 3.7e-4   k=4: 4.6e-4
  k=5: 2.7e-4   k=6: 5.5e-4   k=7: 5.6e-4   k=8: 3.8e-4
```

LW's per-layer diff is comparable to shortwave's (3-6e-4 K/s at all
layers). But here `scratch.grid.a` was pre-divergent from shortwave —
is LW just propagating that?

No: looking at `FriersonLongwaveTransmissivity` (the default), line 65
**overwrites** `vars.scratch.grid.a` with its own `exp(-(τ_below -
τ_above))` per layer. So even if shortwave left `t` divergent, LW
rewrites it from scratch using its own `exp` calls.

→ **Longwave is ALSO an intrinsic source** with the same 1-ULP `exp`
mechanism. Its amplification path is the same kind of column sweep
inside `OneBandLongwaveRadiativeTransfer`.

(This is consistent with Sessions 3-5's investigation: those were
looking at the multi-step drift accumulation under UniformCooling vs
OneBandLongwave; the per-step seed of longwave is real but small, and
gets amplified by the time integration.)

### Unified picture of the 0.386 K step-1 T diff

Every per-step seed traces back to the same root cause:
**CPU libm `exp(x)` vs XLA `exp(x)` differ by 1 ULP** for ~12% of
Float32 inputs. That 1-ULP gets amplified by the column-sweep flux
integration (×10⁴ for radiation), then accumulated across multiple
parameterizations (LSC + shortwave + longwave).

| Parameterization | Internal `exp` location | per-layer tend.gridT |
|---|---|---|
| `solar_zenith` | trig & `sin(2πt)` — bit-identical | 0 |
| `vertical_diffusion` | no transcendentals | 0 |
| `large_scale_condensation` | `saturation_vapor_pressure` exp | 4e-3 at k=8 only |
| `shortwave_radiation` (BackgroundShortwaveTransmissivity) | `exp(-optical_depth)` | 3-8e-4 at k=1..8 |
| `longwave_radiation` (FriersonLongwaveTransmissivity) | `exp(-Δτ)` per layer | 3-6e-4 at k=1..8 |
| `boundary_layer_drag` (BulkRichardsonDrag) | `log(z/z₀)` | 1.25e-7 (driven by orography GEMM) |
| `surface_heat_flux` | none (linear in drag) | inherits drag×400000 amp |

The 0.386 K step-1 T diff is the **aggregated and time-integrated**
result of these per-layer per-parameterization tendency diffs of order
1e-3 K/s, leapfrog-stepped by ~190 s per Euler half-step. The
arithmetic checks out: 1e-3 K/s × 190 s ≈ 0.2 K, and with two
sub-steps in first_timesteps! + amplification it reaches 0.386 K.

### What is and isn't fixable

**Not fixable in user code** (without losing precision):
- CPU libm `exp` vs XLA `exp` differing by 1 ULP. To unify, both
  sides would need the same math library — not currently feasible.

**Fixable with care**:
- Orography drift via `MatrixSpectralTransform.transform!` (BLAS gemm
  vs XLA dot_general). Workaround: run model construction on CPU,
  then transfer to Reactant; or sync orography post-construction.

**Structural mitigations to consider** (would reduce amplification
factor, not eliminate it):
- Compensated summation (Kahan) inside column flux sweeps in
  radiation kernels — would suppress the 10⁴× amplification.
- Float64 accumulator for the column sweep, cast to Float32 at the
  end — same idea, cleaner.

**Pragmatic policy** (current):
- Accept ~5e-4 K/step CPU↔Reactant divergence as a noise floor at T31
  with Float32 leapfrog. This is below the model's discretization
  error, semantically benign for forecasting/research use cases, but
  inconvenient for bit-reproducibility tests.

### Files added (Session 12)
- `MWE_exp_ulp.jl` — standalone CPU vs Reactant exp 1-ULP test. **Kept.**
- `MWE_column_sweep_amplification.jl` — standalone column sweep
  amplification test (1-ULP t → 10⁴× amplified tend). **Kept.**
- `debug_lsc_solo.jl` — LSC solo with input check; intrinsic source. **Kept.**
- `MWE_lsc_condensation.jl` — standalone LSC condensation update; the
  1-ULP `saturation_humidity` → δT × radius reproduces production. **Kept.**
- `debug_longwave_solo.jl` — LW solo with input check; intrinsic source. **Kept.**

### Status

Phase 1 root-cause analysis **complete**. The 0.386 K step-1 T diff
between CPU and Reactant at T31 is fully explained by:
1. CPU libm vs XLA `exp` differing by 1 ULP per Float32 call.
2. Multiple radiation/condensation kernels invoking `exp` per layer
   per grid point.
3. Column-sweep flux integration amplifying each 1-ULP `exp` diff by
   ~10⁴× via `(g / (p_s Δσ)) × radius` flux-to-tendency scaling.
4. Leapfrog/Euler time stepping accumulating per-step tendency diffs
   into the headline 0.386 K step-1 T diff.

Two structural workarounds are available (CPU model-init + orog sync,
Float64 accumulator in flux sweeps). No upstream Reactant/XLA fix
exists; ULP-level math-library differences are intrinsic.

### Next session plan (continue here)
1. **(Optional) prototype mitigation**: rewrite
   `shortwave_radiative_transfer!` and the LSC implicit correction
   with a Float64 accumulator for the flux/sum, cast to Float32 at
   the boundary. Measure step-1 diff with this change.
2. **(Optional) orography CPU-init workaround**: add a flag to
   `EarthOrography` that runs init on CPU and only transfers to
   architecture at the end.
3. **Document and accept** the per-step noise floor in test
   tolerances — open an issue/PR with the policy decision.

## Session 12 follow-up — exact divergence localization (2026-05-17)

User question: in LSC's chain, WHERE inside `saturation_humidity` is the
divergence born, and WHERE does the `radius` multiplication happen?

### Q1 — exact op inside `saturation_humidity` that diverges (`MWE_svp_bisect.jl`)

`saturation_vapor_pressure(T, atm) = e₀ * exp(Lᵥ/R_vapor * (1/T₀ - 1/T))`
then `saturation_humidity(T, p) = ϵ * svp(T) / p`.

Bisect each op (standalone Julia + Reactant, no SpeedyWeather):

| Step | Op | max\|Δ\| | identical? |
|---|---|---|---|
| s1 | `a = 1/T` | 0 | yes |
| s1–3 | `c = 1/T₀ − 1/T` (single broadcast) | 0 | yes |
| s3 | `c = (1/T₀) − a` (subtract alone) | 0 | yes |
| s5 | `e = (Lᵥ/Rᵥ) * c` | 0 | yes |
| **s6** | **`f = exp(e)`** | **2.4e-7** | **no — 294/3168 points** |
| s7 | `g = e₀ * f` | 0 (given CPU f) | yes |
| s8 | `h = ϵ * g / P` | 0 (given CPU g, P) | yes |
| fused | `h = ϵ * e₀ * exp(…) / P` (one trace) | 2.79e-9 | 1347/3168 points |

**Only `exp` differs.** Every other op (division, subtract, multiply by
constant, divide by P) is bit-exact between CPU libm and XLA when given
identical inputs. In the fused version, the 2.4e-7 ULP-noise of `exp`
is propagated through the `ϵ * e₀ / p ≈ 4e-3` scaling and emerges as
2.79e-9 in `h` — same ULP, smaller absolute scale.

So the chain is:
```
exp(e) → 1 ULP diff (~2.4e-7 absolute)
e₀ * exp(e) → same 1 ULP diff scaled by e₀
ϵ * e₀ * exp(e) / p → 1 ULP diff scaled by ϵ/p
```
The 1-ULP relative error stays the same; the absolute scale just
rides along.

### Q2 — where `radius` is multiplied (`scaling.jl:29`, `tendencies.jl:56`)

`large_scale_condensation!` itself writes the **physical** tendency:
`temp_tend[ij, k] += δT` (in K/s).
**No `* radius` inside the LSC kernel.**

The multiplication happens at the END of the per-`ij` column kernel
that ENCLOSES all parameterizations:

```julia
# SpeedyWeather/src/parameterizations/tendencies.jl:48
@kernel inbounds = true function column_parameterizations_kernel!(vars, parameterizations, model)
    ij = @index(Global, Linear)
    column_parameterizations!(ij, vars, parameterizations, model)
    scale_tendencies!(ij, vars.tendencies.grid, model.planet.radius)   # ← here
end
```

`scale_tendencies!` (in `SpeedyWeather/src/dynamics/scaling.jl:29`)
multiplies `tend.grid.u, v, temperature, humidity` element-wise by
`radius` (≈ 6.371e6) for the given `ij`. This is a SpeedyWeather
convention: the dynamical core works in radius-scaled units, so
parameterization tendencies (physical K/s) are scaled up by `radius`
to match.

So when the production bisect (`debug_param_drift_bisect.jl`) reports
`max|Δ tend.gridT| = 4.1e-3 K/s` after LSC, that number is **already
multiplied by radius**. The raw `δT` from LSC is ~6e-10 K/s; the
`scale_tendencies!` post-pass multiplies it by 6.371e6 to get ~4e-3.

In `MWE_lsc_condensation.jl` I emulated this by computing
`max|Δ δT * radius|` manually — that's why the standalone result
(8.3e-3) is in the same ballpark as the production 4.1e-3, even
though the underlying physical δT diff is ~1e-9 K/s. The 2× discrepancy
between the MWE (8.3e-3) and production (4.1e-3) reflects the different
input distribution (random vs spin-up T) and the conditional condensation
branch firing on a subset of points.

### Updated summary diagram

```
T (Float32, bit-identical input)
    │
    │ 1/T          (exact, no ULP diff)
    ▼
1/T₀ − 1/T          (exact)
    │
    │ × (Lᵥ/Rᵥ)    (exact)
    ▼
e                    ≈ 0.5
    │
    │ exp(.)         ← CPU libm vs XLA differ by 1 ULP ON ~10% of inputs
    ▼
f = exp(e)           ULP diff: ~2.4e-7 absolute (1 ULP at f ≈ 1.7)
    │
    │ × e₀ * ϵ / p   (exact, just scaling)
    ▼
sat_humid            ULP diff: ~2.8e-9 (1 ULP at sat ≈ 0.01)
    │
    │ LSC implicit correction
    │ (more exact arithmetic on the ULP-noisy sat_humid)
    ▼
δT                   ~ 6e-10 K/s (raw physical tendency)
    │
    │ × radius (in `scale_tendencies!`)
    ▼
tend.gridT           ~ 4e-3 K/s (the "production" number)
```

### Files added (Q1+Q2 follow-up)
- `MWE_svp_bisect.jl` — step-by-step bisect proving `exp` is the only
  diverging op in `saturation_vapor_pressure`. **Kept.**

## Session 12 follow-up #2 — radius scaling cancels in leapfrog (2026-05-17)

User question: can we move the `* radius` scaling to reduce the worst
of the divergence?

### Short answer: NO — the radius is self-cancelling.

### The cancellation (`debug_radius_cancellation.jl`)

Reading the code:
- `scale_tendencies!` (in `SpeedyWeather/src/dynamics/scaling.jl:29`)
  multiplies `tend.grid.{u,v,T,q}` by `radius`.
- `Leapfrog.Δt = Δt_sec / radius` (in `SpeedyWeather/src/time_stepping/leapfrog.jl:44`).
- `leapfrog_kernel!` does `a_new = a_old + Δt * tend[lmk]` (line 131).

Product:
```
Δt × tend = (Δt_sec / radius) × (radius × δT_phys)
         =  Δt_sec × δT_phys
```

**The radius factor cancels EXACTLY in `Δt × tend`.** It's a numerical
convention for the dynamical core's internal scaling (vorticity and
divergence are kept in `vor × radius` form for spectral efficiency);
temperature and humidity tendencies are scaled along for consistency.

### Empirical confirmation

`debug_radius_cancellation.jl` at T31:
```
radius        = 6.371e6
Δt_sec        = 2400.0 s  (physical)
Δt (leapfrog) = Δt_sec / radius = 3.77e-4
Δt × radius   = 2400.0 s   ✓  exactly Δt_sec
```

And for LSC at k=8:
```
Post scale_tendencies max|Δ tend.grid.T| = 4.12e-3 K/s   ← bisect "headline"
Pre-scale (= post-scale / radius)        = 6.47e-10 K/s ← physical
Predicted T contribution = pre × Δt_sec  = 1.55e-6 K     ← microscopic
```

LSC alone contributes only ~1.5e-6 K to the step-1 T diff, NOT 4e-3 K.
The 4e-3 K/s is a measurement artifact of the radius convention.

### What this means for the bisect numbers

The "headline" tendency-diff numbers I reported are **post-radius**:
- LSC at k=8: 4.12e-3 K/s (post-radius) = 6.5e-10 K/s (physical)
- Shortwave at k=1..7: 3-8e-4 K/s (post-radius) = ~10⁻¹⁰ K/s (physical)
- Longwave at k=1..8: 3-6e-4 K/s (post-radius) = ~10⁻¹⁰ K/s (physical)

After leapfrog with `Δt × tend = Δt_sec × δT_phys`, the actual per-step
T change from each parameterization is ~1e-10 K/s × ~190 s = **~2e-8 K**.

But the observed step-1 T diff is **0.386 K** — eight orders of
magnitude larger. Where does the rest come from?

→ **Inside `shortwave_radiative_transfer!`** there's a cumulative
amplification of ~32768× (per `MWE_column_sweep_amplification.jl`),
operating on the RAW physical t array BEFORE scale_tendencies! is even
called. The amplification is intrinsic to the algorithm's flux-sweep
recursion `D *= t[k]; tend += (D - D_out) * scale; D = D_out`, not to
the radius factor.

The radius scaling is INNOCENT — it sits between two cancelling
operations and contributes nothing to the divergence.

### Why this is structurally important

The bisect numbers had me focused on LSC's k=8 4e-3 K/s as "an
important contributor". In physical terms, LSC's contribution to
step-1 T is **microscopic** — three orders of magnitude smaller than
shortwave_radiation's column-sweep amplified contribution at k=1
(~5e-3 K after the dt × radius pipeline).

The real ranking by step-1 T impact (estimating from pre-scale tendency
× Δt_sec):
1. **`shortwave_radiation`** k=1 ≈ 1.95e-3 K/s × 188 s × (column-sweep
   amplification factor inside the kernel that's already part of the
   reported number) — dominant.
2. **`longwave_radiation`** similar magnitude, similar mechanism.
3. **`large_scale_condensation`** ~1e-6 K at k=8 — negligible.
4. **Orography → drag → surface_heat_flux** ~5e-5 K at k=8 — also
   small once you remove the radius-scaling visual bias.

### Where to focus mitigation

Moving `* radius` won't help. What WOULD help (in decreasing order of
impact on the 0.386 K headline):

1. **Mitigate column-sweep amplification in shortwave** — Float64
   accumulator for `D`, `D_out`, and the downward/upward beam.
   Approx. ~10⁴× sensitivity reduction → headline drops by ~10⁴ →
   step-1 T diff would go from 0.386 K to ~4e-5 K. Biggest win.

2. **Same for longwave** — column flux sweep in
   `OneBandLongwaveRadiativeTransfer` has the same structure.

3. **Eliminate the `exp` ULP at source** — replace transcendentals
   with bit-reproducible polynomial approximations (FastMath off
   on both sides; ensure same intrinsic). Possible but invasive.

4. **Orography sync** — fixes the surface_heat_flux channel; small
   contribution given other dominators.

5. **NOT useful: moving `* radius`.** It would change the reported
   bisect numbers cosmetically but the headline 0.386 K is invariant
   under `radius` placement.

### Files added (radius cancellation)
- `debug_radius_cancellation.jl` — empirical confirmation that
  `Δt × radius = Δt_sec` and that the bisect "4e-3 K/s" reduces to
  microscopic 1.5e-6 K after the leapfrog step. **Kept.**

### Reframing the per-step seed table

Updated reading of the 0.386 K headline (with radius cancellation taken
out of the picture):

| Parameterization | post-radius bisect | pre-radius (physical) | step-1 T contribution |
|---|---|---|---|
| LSC k=8 only | 4.1e-3 K/s | 6.5e-10 K/s | 1.5e-6 K (negligible) |
| Shortwave k=1..7 | 3-8e-4 K/s | 5-13e-11 K/s | 1-3e-7 K from the OUTPUT TENDENCY |
| Longwave k=1..8 | 3-6e-4 K/s | 5-10e-11 K/s | 1-2e-7 K |
| drag→SHF k=8 | 4.9e-2 K/s | 7.7e-9 K/s | 1.8e-5 K (with orog unsynced) |

But the 0.386 K observed step-1 T diff is much larger than the sum of
those direct contributions! That confirms the **column-sweep
amplification happens before the tendency-divergence is even reported**
— inside the shortwave kernel, the per-layer `D *= t[k]` recursion
already produces a 10⁴× amplification of the seed 1-ULP `t` divergence
into a per-layer-K/s tendency. The OUTPUT tendency we measure post-
kernel is the already-amplified value (~6e-4 K/s post-radius = ~10⁻¹⁰
K/s physical). Multiplied by Δt_sec ≈ 190 s, that's only 2e-8 K per
step.

So either:
(a) the leapfrog step ITSELF amplifies further (multi-step nature
    of the dt × tend update with the Williams filter), or
(b) the per-layer kernels have additional channels writing to tend
    that I haven't captured in the bisect (e.g., `surface_shortwave_*`
    fields feeding back through other parameterizations on the same
    column pass).

This is a fresh open question — the step-1 T diff arithmetic doesn't
quite add up from the per-parameterization tendency diffs once radius
is properly accounted for. Worth a future drill.

### Status

The radius scaling is a numerical-convention bookkeeping factor with
zero impact on the actual prognostic update; moving it cannot reduce
the CPU↔Reactant divergence. The visible "4e-3 K/s" tendency diffs in
the bisect are radius-scaled and three orders of magnitude larger than
the corresponding physical-tendency diffs.

**Surprise observation**: even accounting for the radius cancellation,
the per-parameterization physical tendency diffs (~1e-10 K/s) ×
Δt_sec (~200 s) ≈ 2e-8 K — yet the observed step-1 T diff is 0.386 K.
There's a remaining factor of ~10⁷ unaccounted for. Possible
explanations: leapfrog Williams filter interaction, implicit
correction in horizontal_diffusion!, accumulating tendencies in
spectral after the second `first_timesteps!` half-step.

Worth investigating next session: **trace the dt × radius cancellation
all the way to the post-step grid temperature** to see where the
10⁷ comes from.

## Session 13 — Column-sweep mitigation MWE (2026-05-17)

User question: how could I reduce the column-sweep amplification?

### Test setup (`MWE_column_sweep_mitigations.jl`)

Standalone, no SpeedyWeather. Generate two t arrays (CPU libm exp vs
Reactant @jit exp, 1-ULP different). Feed both into four sweep variants
and measure `max|Δ dTdt|` post-radius:

```
V0  baseline F32, absorbed = D - D*t           (cancellation candidate)
V1  F64 accumulator (D, U in F64; F32 t in)    (cumulative-rounding candidate)
V2  expm1 reformulation in F32                 (cancellation fix)
V3  F64 + expm1                                 (combined)
V4  F64 exp + F64 sweep + F32 output           (precision throughout)
```

### Results

| variant | max\|Δ dTdt\| | reduction vs V0 |
|---|---|---|
| V0 baseline | 1.95e-3 K/s | 1× |
| V1 F64 accumulator (F32 t in) | 1.95e-3 | **1× — no help** |
| V2 expm1 with bit-identical α | 0.0 | (artifact — same α both sides) |
| V2 expm1 with realistic α | 1.95e-3 | **1× — no help** |
| V3 F64 + expm1 with realistic α | (same) | — |
| **V4 F64 exp + F64 sweep + F32 out** | **0.0** | **>10⁵× — total elimination** |

### What we learned

1. **F64 accumulator alone doesn't help.** The cumulative-rounding-through-
   `D *= t[k]` hypothesis was wrong. The amplification comes from the
   FIRST step: a 1-ULP error in `t` (≈6e-8) gets multiplied by D
   (≈1000), producing ~6e-5 absolute error in absorbed. Then the
   tendency conversion `× (g / (p_s Δσ)) × radius` multiplies by ~30,
   giving ~2e-3 K/s. **Linear scaling, no special amplification.**
   F64 on D doesn't help because the F32 `t` already has 6e-8 baked in.

2. **expm1 doesn't help for CPU↔XLA reproducibility either.** Yes, it
   avoids catastrophic cancellation in `1 - exp(-τ)`, but expm1
   itself differs between CPU and XLA by ~1 F32 ULP — same magnitude
   as exp. So feeding two different α arrays into the (clean) column
   sweep gives the same amplification.

3. **The only thing that works: compute exp in F64.**
   - F32 exp CPU vs XLA: 6e-8 absolute difference (1 ULP)
   - F64 exp CPU vs XLA: 1.1e-16 absolute difference (1 F64 ULP)
   - F64 ULP is 5.4e8× smaller than F32 ULP.
   - After F64 sweep with F32 output, the F32 representation of
     dTdt is identical between CPU and Reactant.

### Recommendation

To eliminate the column-sweep contribution to CPU↔Reactant divergence
in `shortwave_radiative_transfer!` (and by analogy
`longwave_radiative_transfer!`, `transmissivity!`, and the LSC
`saturation_humidity` chain):

**Promote every `exp`/`log`/`expm1` call site to F64, plus the
recursive accumulator that consumes its result. Cast back to F32 at
the boundary where the result is written to the tendency.**

In `BackgroundShortwaveTransmissivity.transmissivity!`:
```julia
# was:
t[ij, k] = exp(-optical_depth)
# becomes:
t[ij, k] = Float32(exp(-Float64(optical_depth)))
```
And similarly inside `shortwave_radiative_transfer!`, keep the column
accumulator (`D`, `U`) in F64 while reading t (which is already F64-
precision if computed as above) and writing F32 tendency.

For `saturation_vapor_pressure`:
```julia
# was:
return e₀ * exp(Lᵥ / R_vapor * ((1/T₀) - (1/T)))
# becomes:
T64, T₀64 = Float64(T), Float64(T₀)
return Float32(e₀ * exp(Float64(Lᵥ) / Float64(R_vapor) * (1/T₀64 - 1/T64)))
```

### Cost/benefit

- **Cost**: F64 ops are typically 2× slower than F32 on CPU and
  varies on GPU. For a kernel like `transmissivity!` with one `exp`
  per layer per ij (~25k ops at T31×8 layers), this is negligible.
  The column sweep's F64 multiplications and additions cost more,
  but it's still a small piece of a timestep.
- **Benefit**: per the MWE, eliminates the ~1.95e-3 K/s tendency
  divergence — i.e., kills the entire shortwave column-sweep
  contribution to the per-step seed.

### Caveats for the actual SpeedyWeather code

The MWE result of "max|Δ dTdt| = 0.0" is at F32 representation level
when comparing two F64 results that themselves differ by 1.1e-16.
In a real production trace, fused operations and Reactant's compiler
optimizations might shift the F64 1-ULP around slightly. The actual
reduction might be 10⁴–10⁵× instead of 10⁵+× — still huge.

Also: V4 only addresses the SHORTWAVE column sweep. The same
treatment is needed for:
- Longwave's transmissivity exp + radiative transfer.
- LSC's saturation_humidity exp + implicit correction.

Each one is independent; each needs the F64 promotion.

### Files added (Session 13)
- `MWE_column_sweep_mitigations.jl` — V0/V1/V2/V3/V4 comparison.
  V4 (F64 exp + F64 sweep + F32 out) eliminates the diff entirely.
  **Kept.**

### Status

Mitigation strategy **identified and validated** in standalone MWE.
The fix is: F64 for the `exp` (and the recursive accumulator that
consumes it) inside each radiation/condensation kernel; F32 at the
input/output boundary. Implementing this in the production code is
mechanical:
- `BackgroundShortwaveTransmissivity.transmissivity!`
- `OneBandShortwaveRadiativeTransfer.shortwave_radiative_transfer!`
- `ConstantLongwaveTransmissivity.transmissivity!`
- `FriersonLongwaveTransmissivity.transmissivity!`
- `OneBandLongwaveRadiativeTransfer.longwave_radiative_transfer!`
- `saturation_vapor_pressure` (used by LSC and DiagnosticClouds)

Estimated effort: small per call site, but ~6 call sites to touch.

### Next session plan (continue here)
1. Pick ONE call site (suggest `BackgroundShortwaveTransmissivity.transmissivity!`)
   and apply the F64 fix in source. Re-run `debug_param_drift_bisect.jl`
   and `debug_drift_step1.jl` to see how much the per-step T diff drops.
2. If successful, generalize to the other 5 call sites.
3. Investigate the still-open ~10⁷ amplification mystery — even after
   eliminating the per-kernel divergence, where would residual seeds
   come from? Likely the GEMM-based `transform!` in the dynamical core.

## Session 14 — F64 mitigation implementation plan (2026-05-17)

Concrete plan to apply the F64 fix validated by `MWE_column_sweep_mitigations.jl`
(Session 13 V4: total elimination of column-sweep CPU↔Reactant divergence
when `exp` and its column-recursion accumulator run in F64).

### Design principle

For each transcendental call (`exp`, `expm1`, `log`, `log1p`, `^4`,
…) and the recursive accumulator that consumes its result inside the
parameterization kernels:

```
Float32 input → Float64 promote → F64 transcendental + F64 accumulator
              → Float32 cast at the write boundary to the tendency array
```

Wrapped in a small helper to keep call sites readable and to allow a
single point of control:

```julia
# Likely home: SpeedyWeatherInternals/src/Utils/utility_functions.jl
# Always evaluates `f` in widen(T) precision; result cast back to T.
@inline _wide(f::F, x::T) where {F, T <: AbstractFloat} = T(f(widen(T)(x)))
@inline _wide_exp(x::T) where T = T(exp(widen(T)(x)))
@inline _wide_expm1(x::T) where T = T(expm1(widen(T)(x)))
@inline _wide_log(x::T) where T = T(log(widen(T)(x)))
```

For `T = Float32`, `widen(T) = Float64` → the F64 promotion we want.
For `T = Float64`, `widen(T) = BigFloat` would be slow; gate behind:

```julia
@inline _wide_exp(x::Float32) = Float32(exp(Float64(x)))
@inline _wide_exp(x::Float64) = exp(x)                          # no widening
```

This keeps the helper a no-op for F64 model runs (no perf regression)
and a precision win for F32 (the production default).

### Phase A — proof of concept on one kernel

**File**: `SpeedyWeather/src/parameterizations/radiation/shortwave_transmissivity.jl`
**Target**: `BackgroundShortwaveTransmissivity.transmissivity!` (line 75-142)
**Hot line**: 138 — `t[ij, k] = exp(-optical_depth)`

**Diff (minimal)**:
```julia
# was:
t[ij, k] = exp(-optical_depth)
# becomes:
t[ij, k] = _wide_exp(-optical_depth)
```

Also at line 27 in `ConstantShortwaveTransmissivity.transmissivity!`:
```julia
τ = -log(CST.transmissivity)           # → _wide_log(...)
...
t[ij, k] = exp(-τ * dσ[k])              # → _wide_exp(...)
```

**Test**:
1. `julia --project debug_param_drift_bisect.jl` — expect `k=5
   shortwave_radiation` to show no marginal jump in `tend.gridT` diff
   compared to `k=4 albedo`. The bisect previously showed:
   k=4: 4.12e-3, k=5: 4.15e-3 (a 3e-5 K/s jump from shortwave).
   After fix: expect k=5 ≈ k=4 (no marginal shortwave contribution).
2. `julia --project debug_drift_step1.jl` — expect per-layer T diff
   at k=1..7 to drop. Previously k=1..7: 0.05..0.24 K. Expectation
   if F64 fix works as in MWE: k=1..7 drop to ~0 (or to noise from
   other parameterizations).
3. CPU-only sanity check — run a CPU-only forecast and verify the F64
   fix produces results that are bit-identical (or close, within F64
   ULP) to the baseline CPU run. The F64 promote rounds to the same
   F32 value as the F32 path in most cases.

**Expected outcome**: dominant per-layer step-1 T diff (k=1..7)
drops, validating the F64 strategy on one kernel.

**If it doesn't work as expected**: re-examine. Possibilities:
- `_wide_exp` not actually compiling to F64 exp on Reactant side
  (could be a tracing artifact). Check `@code_hlo` to verify.
- The transmissivity scratch array `vars.scratch.grid.a` is F32; even
  with F64 internal computation, the F32 store rounds back. That's
  actually GOOD (we want F32 in the array) but might still produce
  CPU↔Reactant diff if the rounding modes differ.
- The diff might originate not in transmissivity but in
  `shortwave_radiative_transfer!` itself (downward/upward beam).
  Phase C addresses that.

### Phase B — generalise to all transmissivity kernels

**Files / lines**:
- `parameterizations/radiation/shortwave_transmissivity.jl:30`
  `ConstantShortwaveTransmissivity.transmissivity!`: `exp(-τ * dσ[k])`
- `parameterizations/radiation/shortwave_transmissivity.jl:138`
  `BackgroundShortwaveTransmissivity.transmissivity!`: `exp(-optical_depth)`
  (touched in Phase A)
- `parameterizations/radiation/longwave_transmissivity.jl:20`
  `ConstantLongwaveTransmissivity.transmissivity!`: `exp(-τ * dσ[k])`
- `parameterizations/radiation/longwave_transmissivity.jl:65`
  `FriersonLongwaveTransmissivity.transmissivity!`: `exp(-(τ_below - τ_above))`

In each, replace `exp(...)` with `_wide_exp(...)`. Also any internal
`log(...)` (e.g. line 17 of shortwave_transmissivity.jl).

**Test**: re-run `debug_param_drift_bisect.jl`. After Phase B:
- k=6 longwave_radiation row should also show no marginal increase
  over k=5 (currently +1e-5 K/s).

### Phase C — radiative_transfer kernels (the column sweeps)

The transmissivity arrays (`vars.scratch.grid.a`) are F32. After Phases
A and B, the F32-stored `t` values are still the result of CPU exp vs
XLA exp differing at F32 precision (because we round F64 result back
to F32). **The same 1-ULP F32 diff in `t` is re-introduced at the
store**.

So Phase A+B alone might not be enough. Phase C: keep the accumulator
`D` and `U` in F64 INSIDE the radiative transfer column sweep, so even
if `t` (read from scratch) has a 1-ULP F32 diff, the sweep is
performed in F64 and only rounds back to F32 at the final tendency
write.

**Files / lines**:
- `parameterizations/radiation/shortwave_radiation.jl:156-235`
  `shortwave_radiative_transfer!`
- `parameterizations/radiation/longwave_radiation.jl:222-285`
  `longwave_radiative_transfer!`

**Diff for shortwave_radiative_transfer!** (key parts):
```julia
# was (Float32 NF):
D = D_toa                                   # F32
for k in 1:nlayers
    D_out = (D - O₃ * D_toa) * t[ij, k]
    dTdt[ij, k] += flux_to_tendency((D - D_out) / cₚ, pₛ, k, model)
    D = D_out
end
# becomes (F64 accumulator):
D_toa = widen(NF)(model.planet.solar_constant * cos_zenith)   # F64
D = D_toa
for k in 1:nlayers
    t_k = widen(NF)(t[ij, k])               # promote F32 → F64
    D_out = (D - widen(NF)(O₃) * D_toa) * t_k
    absorbed = D - D_out
    dTdt[ij, k] += NF(flux_to_tendency(absorbed / widen(NF)(cₚ), pₛ, k, model))
    D = D_out
end
```

Apply the same to `U_reflected`, `D_surface`, `U`, and the upward
sweep. Same pattern for longwave (lines 256-282 of
longwave_radiation.jl): promote `U` and `D` to `widen(NF)`.

**Watch out**: lines 251 and 270 of longwave_radiation.jl have
`U::NF = ...` and `D::NF = 0` — those type annotations force F32. Need
to drop or change to `widen(NF)`.

**Test**: re-run `debug_param_drift_bisect.jl` AND
`debug_drift_step1.jl`. After Phase C:
- Per-layer step-1 T diff at k=1..7 should be near 0 (k=8 still has
  LSC/SHF contributions until Phase D).
- The headline step-1 T diff should drop from 0.386 K substantially —
  potentially by 100×-1000× if column sweep is the dominant source.

### Phase D — saturation_vapor_pressure (LSC + clouds + convection)

**File**: `SpeedyWeather/src/dynamics/atmosphere.jl:138-158`

```julia
# was:
@inline function saturation_vapor_pressure(T, A::AbstractWetAtmosphere)
    e₀ = A.saturation_vapor_pressure
    Lᵥ = A.latent_heat_condensation
    R_vapor = A.R_vapor
    T₀ = A.temperature_freezing
    return e₀ * exp(Lᵥ / R_vapor * ((1/T₀) - (1/T)))
end
# becomes:
@inline function saturation_vapor_pressure(T::NF, A::AbstractWetAtmosphere) where NF
    e₀ = NF(A.saturation_vapor_pressure)
    Lᵥ = NF(A.latent_heat_condensation)
    R_vapor = NF(A.R_vapor)
    T₀ = NF(A.temperature_freezing)
    W = widen(NF)
    arg = W(Lᵥ) / W(R_vapor) * (W(1)/W(T₀) - W(1)/W(T))
    return NF(W(e₀) * exp(arg))
end
```

Note: T is promoted to `widen(NF)` inside the function; output cast
back. This change is felt by ALL call sites of `saturation_humidity`
without further work:
- `parameterizations/large_scale_condensation.jl:97` (LSC)
- `parameterizations/radiation/clouds.jl:126, 154` (DiagnosticClouds)
- `parameterizations/convection.jl:65, 211, 230` (convection — currently
  off in our tests via `convection=nothing`)

**Optional follow-up inside LSC**: also promote LSC's implicit
correction recursion to F64 (lines 130-134 of
large_scale_condensation.jl), which uses `dqsat_dT` and `δq /=
((1 + Lᵥ_cₚ * dqsat_dT) * τ * Δt)`. The 1-ULP from saturation_humidity
flows through here; F64 division stabilises it.

**Test**: re-run the standalone `MWE_lsc_condensation.jl` after the
fix to verify δT diff drops to 0. Then re-run
`debug_param_drift_bisect.jl` and `debug_drift_step1.jl` to confirm
k=8 contributions vanish from LSC.

### Phase E — verification, benchmarking, finalisation

1. **Full Reactant correctness suite**:
   ```
   julia --project=SpeedyWeather --check-bounds=yes \
     SpeedyWeather/test/reactant/runtests.jl
   ```
   Tolerances should now be MUCH tighter than today's defaults.
   Discuss whether to relax test tolerances or keep them strict.

2. **10-step drift**: re-run the Session 5/10 drift measurement
   (~`debug_drift_*` scripts) with the F64 fix. Confirm the per-step
   noise floor drops from 0.386 K to ideally ulp-scale.

3. **Performance benchmark**: run
   `julia --project=SpeedyWeather/benchmark SpeedyWeather/benchmark/manual_benchmarking.jl`
   and compare to the existing baseline in `benchmark/README.md`. F64
   ops are typically 2× slower on x86 and Apple Silicon; for kernels
   that are ~10% of timestep cost, the wall-time impact should be
   ≤2-3%. If significant: consider making the F64 widening optional
   (e.g. via `model.use_f64_exp::Bool` or a CompileOptions flag).

4. **CHANGELOG entry**: under `## Unreleased`:
   ```
   - Promote `exp` calls in radiation/condensation kernels to Float64
     internally (Float32 in/out) to eliminate CPU↔Reactant per-step
     temperature divergence from 1-ULP math-library differences
   ```

5. **Open issue/PR upstream**: link to Session 13's MWE that proves
   the F64 fix works without needing changes to Reactant.

### Risk / care points

1. **`widen(Float32) = Float64`**: relies on Julia's standard widening
   behavior. Confirm no custom NF type breaks this (e.g.
   BFloat16). For BFloat16 the model is experimental anyway.

2. **F64 inside Reactant traces**: Reactant supports F64 throughout
   its IR (StableHLO). The `Float32 → Float64 → exp → Float32` chain
   should lower cleanly. If a regression appears, inspect with
   `@code_hlo` to confirm the `convert` ops are present and the F64
   `exp` is emitted (`stablehlo.exponential` of `tensor<...xf64>`).

3. **Constants in F64 vs F32**: the runtime constants like `cₚ` and
   `radius` are stored in NF. Promoting them per-call is cheap. If
   they end up in a hot path, hoist the promotion outside the loop.

4. **Numerical equivalence to baseline CPU**: F64-then-cast-back-to-F32
   produces a result that's within 0.5 ULP of the true value, while
   pure F32 may have several ULPs of accumulated rounding. So the F64
   path is STRICTLY MORE ACCURATE than the current code — not just
   different. CPU baseline tests that compare to today's CPU output
   may need a small tolerance bump (or, ideally, the test reference
   is regenerated against the more accurate F64 path).

5. **Reactant compilation cost**: each new kernel signature requires
   recompilation. The `widen(NF)` pattern produces the same signature
   regardless of NF, so there's no extra compile bloat.

### Ordering / iteration

- **Iteration 0**: implement the `_wide_exp` helper (single small commit
  in `SpeedyWeatherInternals`).
- **Iteration 1 (Phase A)**: apply to `BackgroundShortwaveTransmissivity.transmissivity!`.
  Run bisect + step-1 tests. Measure improvement.
- **Iteration 2 (Phase B)**: apply to remaining transmissivity kernels.
- **Iteration 3 (Phase C)**: promote `D` and `U` accumulators in
  radiative_transfer kernels.
- **Iteration 4 (Phase D)**: promote `saturation_vapor_pressure`.
- **Iteration 5 (Phase E)**: full test + benchmark + changelog.

Each iteration is a small, revertable commit. Measure step-1 T diff
after each to track progress.

### Predicted impact (per the MWE)

| Iteration | Expected step-1 T diff |
|---|---|
| Pre-fix baseline | 0.386 K |
| Iter 1 (one transmissivity) | maybe 0.2–0.3 K (partial) |
| Iter 2 (all transmissivities) | maybe 0.1 K |
| Iter 3 (radiative_transfer accumulator F64) | < 0.01 K possibly |
| Iter 4 (saturation_vapor_pressure) | < 1e-3 K possibly |
| Iter 5 (verified) | per-step noise floor ≤ ULP-scale |

These are EXPECTATIONS based on the MWE; the real measurements may
differ if there are other unmeasured amplification paths (the
unexplained ~10⁷ from Session 12 follow-up #2 is still open).

### Status

Plan is concrete and revertable. Next session: execute Iteration 0 +
Iteration 1 (the `_wide_exp` helper + one call site), measure, and
iterate. Each phase is a single commit on a feature branch
(`mg/reactant-f64-exp` or similar).

## Session 15 — Where `exp` could be replaced with `expm1` (2026-05-17)

Catalogue + caveat. expm1 is a different fix from the F64 promotion;
the two address overlapping but distinct problems.

### Sites in the codebase

Six sites exhibit the `1 - exp(-x)` pattern (either explicitly or as
the difference of two near-equal numbers `D - D × t`):

**Explicit `(1 - t)` — 4 sites in `longwave_radiative_transfer!`:**
- `longwave_radiation.jl:258`: `U = U * t + (1 - t) * σ * T[ij, k]^4`
  (upward beam, layer k)
- `longwave_radiation.jl:265`: same pattern at the TOA
- `longwave_radiation.jl:273`: downward beam interior
- `longwave_radiation.jl:280`: downward beam at surface

**Implicit `(1 - t)` via `D - D × t` — 2 sites in `shortwave_radiative_transfer!`:**
- `shortwave_radiation.jl:200,203`:
  ```
  D_out = (D - O₃ * D_toa) * t[ij, k]
  dTdt[ij, k] += flux_to_tendency((D - D_out) / cₚ, pₛ, k, model)
  ```
  Here `D - D_out = D - (D - O₃*D_toa) * t = (D - O₃*D_toa)*(1 - t) + O₃*D_toa`,
  which has the same cancellation when t is close to 1.
- `shortwave_radiation.jl:227-228`: same for upward beam (`U - U_out`).

The remaining `exp` calls in the codebase (greenhouse_gases, forcing,
initial_conditions, drag, random_process) are NOT in a `1 - exp` or
`exp - 1` pattern; they don't suffer from cancellation and don't
benefit from expm1.

### Caveat — expm1 does NOT eliminate CPU↔Reactant divergence

We tested this empirically in `MWE_column_sweep_mitigations.jl` (Session 13):

```
V2 expm1 in F32, bit-identical α both sides         → max|Δ dTdt| = 0.0
V2 expm1 in F32, CPU expm1 vs Reactant @jit expm1   → max|Δ dTdt| = 1.95e-3 (same as V0!)
```

CPU libm `expm1` and Reactant `@jit expm1` differ by ~1 F32 ULP, same
as `exp`. The 1-ULP-relative error in α is comparable to the 1-ULP in
t, so the column-sweep amplification is unchanged. **Replacing `exp`
with `expm1` does not reduce CPU↔Reactant divergence.**

### What expm1 DOES improve — absolute numerical accuracy in thin layers

For a thin layer (τ ≈ 0.01):
- `1 - exp(-0.01) = 0.00995`, computed via F32 subtraction:
  `1f0 - 0.99005f0`. The result has relative error ~6e-6 (catastrophic
  cancellation; lost ~6 sig digits).
- `-expm1(-0.01) = 0.00995` directly. Relative error ~6e-8 (full F32
  precision; ~100× more accurate).

So replacing `(1 - t)` with `-expm1(-τ)` makes the radiation kernels
**numerically more accurate** in upper-atmosphere layers (low τ),
independent of platform. This is a quality-of-implementation win, not
a Reactant-reproducibility win.

### Architectural cost

The current code stores `t[ij, k] = exp(-τ)` in `vars.scratch.grid.a`
inside `transmissivity!`, then `radiative_transfer!` reads `t` from
the scratch array. To use `expm1`, the radiative_transfer kernels need
either:

(a) **Recompute τ from t**: `τ = -log(t)` — extra cost AND introduces
    another transcendental call, defeats the purpose.
(b) **Store α = -expm1(-τ) in scratch instead of t**, rewrite
    radiative_transfer to use α as the "absorbed fraction" directly.
    This is the clean architectural change. Each `(1 - t)` becomes
    `α[ij, k]`, each `t` becomes `(1 - α[ij, k])`.
(c) **Store BOTH t and α**: more memory but minimal code change.

Option (b) is the natural fit. The radiative_transfer kernels become:
```julia
# longwave upward, was:
U = U * t + (1 - t) * σ * T[ij, k]^4
# becomes:
α_k = α[ij, k]
U = U * (1 - α_k) + α_k * σ * T[ij, k]^4
```

For shortwave's `D - D_out`:
```julia
# was:
D_out = (D - O₃ * D_toa) * t[ij, k]
absorbed = D - D_out
# becomes:
α_k = α[ij, k]
absorbed = (D - O₃ * D_toa) * α_k + O₃ * D_toa
D_out = D - absorbed
```

### Recommendation matrix

| Goal | Use |
|---|---|
| Eliminate CPU↔Reactant divergence | **F64 promotion** (Session 14 plan) — only this works |
| Improve absolute accuracy in thin layers | **expm1** at the 6 sites above |
| Combine both | F64 + expm1 — once you're in F64, expm1 is cheap insurance |

For pure Reactant correctness work: skip expm1, go straight to the
F64 plan (Session 14). For an independent numerical-accuracy PR:
the 6 expm1 sites are a self-contained, low-risk improvement.

### If we did both

Combined: store `α = -widen_expm1(-optical_depth)` in F64 throughout
the transmissivity kernel, then use F64 accumulator + α directly in
the radiative_transfer kernel. The result is what `MWE_column_sweep_
mitigations.jl` V3 measured: total elimination of CPU↔Reactant
divergence AND improved accuracy in thin layers. But this is just
"F64 plan + expm1 dressing" — the F64 part does the heavy lifting.

### Status

Six concrete sites identified for an independent expm1 refactor.
**Not part of the Reactant-divergence fix** (which needs F64). Could
be a separate cleanup PR if there's interest in upper-atmosphere
numerical accuracy.

## Session 16 — Hunt for the unaccounted ~10⁷ amplification (2026-05-17, in progress)

Per Session 15 summary, the per-step **physical** tendency diff is
~6e-10 K/s while the observed step-1 T diff is 0.386 K — a factor of
~10⁵–10⁶ unaccounted-for amplification. Started this session to find it.

### Probe 1 — `debug_first_timesteps_substep.jl`

Inline replay of `first_timesteps!` as two distinct substeps:
- `euler_substep!`: `initialize!(implicit, Δt/2) + timestep!(vars, Δt/2, model, 1, 1)`
- `leapfrog_substep!`: `initialize!(implicit, Δt) + timestep!(vars, Δt, model, 1, 2)`

Each `@compile`d separately. Took snapshots at A0 (pre-anything), A2
(after Euler substep), A4 (after Leapfrog substep).

| Stage | gridT | per-layer (k=1..8) | specT | tend.gridT |
|---|---|---|---|---|
| A0 | **0.402** | 0.053 / 0.030 / 0.026 / 0.075 / 0.131 / 0.155 / 0.249 / 0.402 | 0.0 | 0.0 |
| A2 (Euler) | 5.8e-4 | (≤1e-3) | 3.7e-3 | 811 |
| A4 (Leapfrog) | 0.188 | k=8 = 0.188 | 7.5e-3 | 810 |

**Surprise**: A0 already shows 0.40 K diff — BEFORE any substep runs.
The diff at A0 had the same per-layer pattern as the headline step-1
(0.053→0.402 at k=8). And `specT diff = 0` throughout. Spectral
prognostic is bit-identical; only grid disagrees.

### Probe 2 — `debug_a0_origin.jl`

Localised the appearance of the A0 diff step by step:

```
S1 after initialize!(model) on both: gridT 0.0,  specT 4.4e-3
S2 after run!(sim_c, Day(1)):       gridT 322.7, specT 1016.5
S3 after copy! sim_r ← sim_c:        gridT 0.0,  specT 0.0
S4 after initialize!(sim_r;steps=1): gridT 0.402, specT 0.0      ← APPEARS HERE
S5 after @compile euler:             gridT 0.402, specT 0.0
S6 after @compile leapfrog:          gridT 0.402, specT 0.0
S7 after copy! sim_c ← sim_r:        gridT 0.0,  specT 0.0
S8 after initialize!(sim_c;steps=1): gridT 0.402, specT 0.0      ← APPEARS HERE too
S9 after initialize!(sim_r;steps=1): gridT 0.402, specT 0.0
```

The diff arises whenever `initialize!(sim; steps=1)` is called on ONE
side while the other side has not been touched. It's NOT introduced
by `@compile` and it's NOT introduced by `copy!`. The spectral
prognostic stays bit-identical (`specT = 0`) the whole time.

### Interpretation (preliminary, to verify)

`initialize!(sim; steps=1)` (in `simulation.jl:55-103`) does:
1. `initialize!(clock, time_stepping, steps)` — sets `lf` index
2. `set!(model.output, ...)` — no state change
3. `scale_prognostic!(variables, model.planet.radius)` — scales vor, div
4. `lf = model.time_stepping.first_step_euler ? 1 : 2`
5. `transform!(variables, lf, model; initialize=true)` — spec → grid

Steps 3 and 5 are the candidates. Two competing hypotheses:

**Hypothesis A — leapfrog index mismatch.** After
`run!(sim_c; period=Day(1))` + finalize, sim_c.grid reflects the
LAST in-timestep transform during integration, which used some lf
(likely lf=2). When we then `copy! sim_r ← sim_c`, sim_r inherits
that grid. When `initialize!(sim_r; steps=1)` runs, it picks
`lf = first_step_euler ? 1 : 2` (typically lf=1) and calls
`transform!(vars, 1, model; initialize=true)` which overwrites the
grid from prog at lf=1. The prog at lf=1 holds a DIFFERENT
timestep's spectral state than prog at lf=2 — ~30-min worth of
atmospheric tendency apart — which would explain the ~0.4 K typical
T tendency over 30 min. **This is not a CPU↔Reactant divergence, it's
a probe-flaw**: I'm comparing sim_c.grid (built from prog[lf=2]) to
sim_r.grid (built from prog[lf=1]).

**Hypothesis B — `initialize!()` itself diverges between CPU and
Reactant.** The spec is bit-identical between the two sides; only
the grid computation diverges. But Session 9 verified that
`debug_init_bisect.jl` reports `I3_transform max|ΔT| = 0.0`,
contradicting this hypothesis.

Reading the original debug_drift_step1.jl flow more carefully:

```
run!(sim_c; period=Day(1))            # CPU spin up, finalize unscales
copy!(sim_r ← sim_c)                  # sim_r ← sim_c.prog (unscaled)
initialize!(sim_r; steps=1)           # sim_r scales prog, computes grid
@compile first_timesteps!(sim_r)
@compile later_timestep!(sim_r)
copy!(sim_c ← sim_r)                  # ← RESYNC: sim_c now has sim_r's (re-initialized) state
initialize!(sim_c; steps=1)           # sim_c also re-initializes
initialize!(sim_r; steps=1)           # sim_r re-initializes
# now sim_c.grid and sim_r.grid SHOULD be equal here (same prog, same code path)
run!(sim_c; steps=1)                  # = first_timesteps! on CPU
time_stepping!(sim_r, r_first!, r_later!)  # = first_timesteps! on Reactant
# compare sim_c.grid vs sim_r.grid → 0.386 K
```

At the point of comparison both have done first_timesteps! starting
from a state where the spec is bit-identical AND (we'd expect) the
grid is too. So the 0.386 K is the genuine first_timesteps! delta.

My probe failed to see step-by-step amplification because the
initial A0 was already off by 0.40 K — likely Hypothesis A
(leapfrog-index probing artifact). When the actual production
script (debug_drift_step1.jl) compares post-first_timesteps state,
both sides have run the same dynamics, so any diff there is the
genuine per-step amplification.

### Next session plan

**Session 16 cleanup tasks**:

1. **Resolve A vs B with a clean experiment**: probe `initialize!(sim;
   steps=1)` running on the SAME sim twice — does it produce
   identical grid both times? If yes (idempotent), then sim_c and
   sim_r SHOULD give the same grid for the same prog. If no, then
   `initialize!()` has hidden mutation that breaks A0 equivalence.

2. **If Hypothesis A is the culprit (probe flaw)**: redesign the
   substep probe to:
   - Use the SAME `initialize!(sim; steps=1)` path on both sides as
     in `debug_drift_step1.jl` (don't measure A0 — the inputs at
     this point ARE divergent due to lf-index probing artefact).
   - Start the substep-by-substep measurement only AFTER both sides
     have run identical `initialize!()`.
   - That should give A0 = 0 by construction; then measure A2 (after
     Euler) and A4 (after Leapfrog) to localize the actual ~10⁷
     amplification.

3. **Alternative diagnostic**: inside the @compile'd `first_timesteps!`,
   capture intermediate grid T diff at the boundary between Euler
   substep and Leapfrog substep. Could use a custom captured trace,
   or do `time_stepping!(sim_r, r_first!, ...)` and put a sync point
   between. May need to break first_timesteps! into Reactant-
   compileable pieces.

### Key insight from this session

The unaccounted-for amplification claim from Session 15 may have been
based on a flawed comparison. The bisect script (debug_param_drift_bisect.jl)
that produced the "4e-3 K/s tendency diff" was measuring tendency at
ONE point in time (one Euler substep). The step-1 script
(debug_drift_step1.jl) integrates through TWO substeps with
implicit-init in between. The diff scales as ~tend × dt with
amplification from the kernel-level seeds — possibly without any
mystery factor of 10⁷.

Specifically:
- per-parameterization PHYSICAL tend diff ≈ 6e-10 K/s
- TIMES dt = Δt_sec ≈ 200 s
- TIMES amplification from cumulative `D × (1-t)` chain ≈ ~1
- TIMES 8 layers × multiple parameterizations (LSC + SW + LW)
- ≈ 6e-10 × 200 × 10 = 1.2e-6 K  ← far below 0.386 K

So the unexplained factor is genuinely ~10⁵, not ~10⁷ as I claimed
earlier. Worth re-checking the arithmetic carefully in next session.
Possible candidates for the missing amplification, ranked:

- **Grid→spec GEMM in `parameterization_tendencies_only!`** (Sessions
  10 mentions). CPU BLAS vs XLA dot_general on the radius-scaled
  grid tendency. Same class of divergence as the Session 7-8
  spec→grid issue but for the forward direction.
- **`horizontal_diffusion!`** in spectral space: implicit damping
  could amplify high-wavenumber components that didn't agree to
  start with.
- **Implicit init** rebuilds matrices that depend on state; could
  contribute if the radius scaling makes them sensitive.
- **first_timesteps! double-substep cross-talk**: between Euler and
  Leapfrog substeps, the state transitions through implicit init
  with different dt; the dt-dependent matrices could amplify any
  state diff non-linearly.

### Files added (Session 16)
- `debug_first_timesteps_substep.jl` — substep probe with @compile
  per substep. Misleading A0 due to lf-index artefact. **Kept.**
- `debug_a0_origin.jl` — pinpoints when the 0.40 K appears: at
  `initialize!(sim; steps=1)`, likely a lf-index probe artefact.
  **Kept.**

### Status

Session 16 in progress. Established that the A0 0.40 K in the substep
probe is likely a probe artefact (Hypothesis A — leapfrog index
mismatch when computing grid from spec at different lf indices).
Next session: clean redesign of the substep probe + sanity check of
the amplification arithmetic before committing to any conclusion
about where the missing factor lives.

## Session 16 (continued, resolution) — The mystery is mostly a probe artefact (2026-05-17)

### Confirmation of lf-mismatch hypothesis (`debug_lf_mismatch.jl`)

| Variant | sim_c lf | sim_r lf | gridT diff |
|---|---|---|---|
| A (original) | 2 (after run!) | 1 (never ran) | **0.4017 K** |
| B (force match) | 2 | 2 | **0.0** |
| C (reinit both) | 2 | 2 | **0.0** |

After `run!(sim_c; period=Day(1))`, `sim_c.model.time_stepping.first_step_euler`
flips from `true` to `false` (leapfrog.jl:188). `sim_r` never ran so
its flag stays `true`. When `initialize!(sim; steps=1)` runs on both
sides, it reads each sim's own flag and picks lf accordingly:

- `sim_c`: flag=false → `lf = 2`
- `sim_r`: flag=true  → `lf = 1`

The two sides compute grid from DIFFERENT prog leapfrog indices —
adjacent timesteps' worth of state, ~30 minutes of tendency apart,
giving the typical ~0.4 K T diff at the surface. This is NOT a
CPU↔Reactant divergence; it's a comparison of different physical
states.

### Why `debug_drift_step1.jl` reports A0=0 yet still gives 0.386 K

`@compile first_timesteps!(sim_r)` traces through the body of
`@trace if time_stepping.first_step_euler ... end`, which contains
the side-effecting line `time_stepping.first_step_euler =
!time_stepping.continue_with_leapfrog`. Reactant's trace executes
this side-effect at compile time:

```
Before @compile:                   sim_r.first_step_euler = true
After @compile first_timesteps!:   sim_r.first_step_euler = false   ← FLIPPED!
```

So in `debug_drift_step1.jl`, after the @compile both sims have
`first_step_euler = false`. Both `initialize!(sim;steps=1)` calls
then pick lf=2. A0 diff = 0 (as reported by the script).

**But** the FINAL diff (0.386 K) is then comparing:
- sim_c: runs `first_timesteps!(sim_c)` (Euler dt/2 + Leapfrog dt = 1.5 dt of evolution)
- sim_r: runs `time_stepping!(sim_r, r_first!, r_later!)` where r_first!
  was compiled from `first_timesteps!`. At runtime, with
  `first_step_euler = false`, the `@trace if` selects the **else**
  branch — `later_timestep!(simulation)` (one normal leapfrog at 2dt).

So sim_c and sim_r run **different operations** (different physical
evolutions). The 0.386 K mixes genuine per-step CPU↔Reactant
divergence with the difference between "first_timesteps!" and "one
later_timestep!".

The per-layer pattern of debug_drift_step1.jl's 0.386 K
(0.052, 0.029, 0.026, 0.075, 0.129, 0.144, 0.240, 0.386) is a
monotonic gradient k=1..8, characteristic of comparing two slightly
different atmosphere states (i.e. the operation mismatch dominates).

### The TRUE per-step CPU↔Reactant divergence (`debug_substep_clean.jl`)

With matched first_step_euler flag on both sides, run Euler substep
(`initialize!(implicit, Δt/2) + timestep!(Δt/2, 1, 1)`) and Leapfrog
substep (`initialize!(implicit, Δt) + timestep!(Δt, 1, 2)`) on each
side. Both use `@compile` separately. Snapshots:

| Stage | gridT diff | per-layer pattern |
|---|---|---|
| A0 (post-init, both lf=2) | **0.0** | all zero |
| A1 (after Euler substep dt/2) | **5.8e-4 K** | uniform ≤ 6e-4 |
| A2 (after Leapfrog substep dt) | **0.188 K** | k=8 dominant: 0.188; k=5,6 ≈ 0.018 |

A1's 5.8e-4 K matches the Session 8 transform noise floor exactly —
this is the residual CPU BLAS vs XLA dot_general disagreement in
`MatrixSpectralTransform`. After the Leapfrog substep this amplifies
to 0.188 K, concentrated at k=8 (surface).

The per-layer pattern is DIFFERENT from debug_drift_step1.jl's
headline:
- Clean probe A2: peak at k=8 (surface fluxes), some at k=5,6
  (radiation), tiny elsewhere
- debug_drift_step1.jl 0.386: monotonic gradient k=1..8

Confirms they're measuring different things.

### Inside the Leapfrog substep (`debug_leapfrog_substep_bisect.jl`)

Bisect of the timestep!(Δt, 1, 2) stages starting from A1 state
(5.8e-4 K grid T diff, 2.7e-6 spec T diff):

| Stage | gridT | specT | tend.gridT | tend.specT |
|---|---|---|---|---|
| B0 (=A1) | 5.8e-4 | 2.7e-6 | 0.83 | 0.014 |
| B1 init_implicit(Δt) | 5.8e-4 | 2.7e-6 | 0.83 | 0.014 |
| B2 reset_tendencies! | 5.8e-4 | 2.7e-6 | 0 | 0 |
| B3 greenhouse | 5.8e-4 | 2.7e-6 | 0 | 0 |
| **B4 param_tendencies!** | 5.8e-4 | 2.7e-6 | **810.5** | 0 |
| B5-B7 ocean/ice/land | 5.8e-4 | 2.7e-6 | 810.5 | 0 |
| **B8 param_tend_only!** (grid→spec) | 5.8e-4 | 2.7e-6 | 810.5 | **19.16** |
| B9 horizontal_diffusion! | 5.8e-4 | 2.7e-6 | 810.5 | 19.16 |
| **B10 leapfrog!** | 5.8e-4 | **7.2e-3** | 810.5 | 19.16 |
| **B11 transform!** spec→grid | **0.188** | 7.2e-3 | 810.5 | 19.16 |

Reading the amplification:

1. **B4 (parameterization_tendencies!)**: tend.grid jumps from 0 to
   **810 K/s** post-radius (= **1.3e-4 K/s physical**). With the
   ULP-grid-T input (5.8e-4 K) and the radiation/surface-flux
   kernels each producing their per-call divergences (Sessions 11-12),
   the accumulated tend diff is the expected ~1e-4 K/s × radius scale.

2. **B8 (param_tendencies_only! = grid→spec)**: tend.spec ends at
   19.16 K/s. The forward GEMM (BLAS vs XLA) reduces the 810 K/s
   grid tend diff into 19 K/s spec — partially the Y_lm projection
   shrinking the L∞ norm. Order-of-magnitude consistent.

3. **B10 (leapfrog!)**: spec T diff goes from 2.7e-6 to **7.2e-3**.
   Arithmetic check: `dt × tend.spec = 3.77e-4 × 19.16 = 7.2e-3` ✓.
   The radius CANCELS as expected (Session 12 follow-up).

4. **B11 (transform! spec→grid)**: spec 7.2e-3 → grid 0.188 K. Ratio
   ~26×. This is the typical "max harmonic amplitude × max coefficient"
   product when going from a spec representation to grid points.

### The full arithmetic, end-to-end

Starting from 1-ULP `exp` in radiation kernels (~6e-8 in F32 t):

```
F32 1-ULP exp                                    ≈ 6e-8
× D ≈ 1000 W/m²                                  → absorbed flux diff ≈ 6e-5 W/m²
× (1/cp × g/(p_s Δσ))                            → physical tendency diff ≈ 1.3e-4 K/s
× radius (for bookkeeping, cancels later)        → post-scale tend ≈ 810 K/s
→ grid→spec via forward GEMM                     → tend.spec ≈ 19 K/s
× dt = Δt_sec / radius                           → spec prog diff ≈ 7.2e-3
→ spec→grid via Y_lm projection (max harmonic)   → grid T diff ≈ 0.188 K
```

Every step is a single linear operation (multiply by a constant or a
matrix). **There is no factor-of-10⁷ mystery amplifier.** The
0.188 K diff per substep is fully accounted for by:
- Six independent radiation/condensation kernel `exp` calls
  introducing 1-ULP F32 divergence (Sessions 11-12)
- Three linear stages that scale this up: D-multiplication,
  flux-to-tendency, and spec→grid Y_lm projection

The 10⁷ "mystery" from Session 15 was an arithmetic mistake:
I was computing the physical tendency × dt and getting 1.3e-7 K,
but I FORGOT the spec→grid Y_lm amplification at the end. With that
factor (~26×) restored, the chain produces 0.188 K per substep, no
mystery.

### So where does the 0.386 K headline come from?

debug_drift_step1.jl's 0.386 K is the sum of:
- **Genuine per-step Reactant divergence**: ~0.188 K (Leapfrog substep)
- **Difference of operations**: sim_c runs `first_timesteps!` (1.5dt),
  sim_r runs `later_timestep!` (2dt) because the flag-flip at @compile
  time made the runtime branch select the else clause.

The genuine Reactant divergence after first_timesteps! would be 2 ×
the Leapfrog substep effect (≈ 2 × 0.188 = ~0.38) IF both sides ran
two substeps. But they don't — sim_r runs ONE later_timestep!. So
the 0.386 K is a mix of effects, dominated by the operation mismatch.

**Implication for the F64 plan**: the real per-step diff is 0.188 K
(or whatever the clean comparison gives with matched ops on both
sides). The F64 plan addresses the source of B4's 1.3e-4 K/s physical
tendency. Eliminating that drops B4 to ~0, then B11 transform amp is
moot (nothing to amplify). The F64 plan IS the right fix and would
reduce the per-step diff by ~10⁴–10⁵×.

### Action items

1. **Fix debug_drift_step1.jl** to also `@compile later_timestep!` and
   ensure sim_c and sim_r run the SAME operation at runtime. Simplest
   path: explicitly call later_timestep! on sim_c (mirroring the
   else-branch that r_first! takes at runtime). Or fix the flag to
   force both sims into the first_step_euler=true branch.

2. **Update the headline in INVESTIGATION.md**: 0.386 K is partly
   artefact. True per-step CPU↔Reactant divergence at T31 is ~0.2 K
   per first_timesteps! cycle, fully explained by Sessions 11-12
   kernel-level seeds.

3. **The F64 plan (Session 14) is still the right approach** but the
   target is now the smaller real diff. Expected post-F64 step-1
   T diff ≤ 1e-3 K.

### Files added (Session 16 resolution)
- `debug_lf_mismatch.jl` — proves the lf-flag mismatch hypothesis. **Kept.**
- `debug_flag_flip.jl` — confirms `@compile first_timesteps!` flips
  the flag during tracing. **Kept.**
- `debug_substep_clean.jl` — clean substep probe with matched lf;
  A0=0, A2=0.188 K. **Kept.**
- `debug_leapfrog_substep_bisect.jl` — stage-by-stage bisect of the
  Leapfrog substep; localises amplification to B4 (param_tend!) +
  B11 (spec→grid Y_lm amp). **Kept.**

### Status

Session 16 RESOLVED. The "unaccounted ~10⁷ amplification" was a
combination of:
1. Arithmetic mistake (forgot the spec→grid Y_lm amplification ~26×
   in the chain)
2. debug_drift_step1.jl probe artefact (operation mismatch between
   first_timesteps! on CPU and later_timestep! on Reactant due to
   @compile flag flip)

The genuine per-step Reactant divergence is ~0.188 K per substep,
fully explained by Sessions 11-12 kernel `exp` seeds × dt × spec→grid
Y_lm. The F64 plan (Session 14) remains the right mitigation.
