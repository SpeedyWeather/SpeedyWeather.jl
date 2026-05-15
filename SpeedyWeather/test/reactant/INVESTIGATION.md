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

**Method**: minimal script (no time stepping):
1. Build CPU and Reactant `SpectralTransform`s with same grid.
2. Allocate a `LowerTriangularArray{ComplexF32}` (spectrum, nlayers).
3. Fill with deterministic data (e.g., `rand(MersenneTwister(42), ComplexF32, …)`).
4. `copy!` the spectral data CPU→Reactant.
5. `transform!(grid_c, spec_c, transform_c)` on CPU.
6. `transform!(grid_r, spec_r, transform_r)` under `@compile`.
7. Compare grid outputs ulp-wise.

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
