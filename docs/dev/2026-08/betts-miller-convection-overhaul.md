# Betts-Miller convection overhaul: branchless vertical loops, entrainment, convective snow

> Status: **in progress**. Plan approved; docs audit and code restructuring not yet started.

Date of initial draft: 2026-08-22

Base revision: `557b38d5` (`mc/convection`, off `main`)

## Originating prompt

> Switch to main, start a new branch mc/convection in that I want you to overhaul the Betts
> Miller convection scheme we have. Check that the code matches the equations in the docs, and
> make changes towards GPU performance, which you can benchmark locally but only benchmark
> convection in isolation and together with other parameterizations don't bother doing a full
> model benchmark. At the moment the convection scheme writes NaNs above the level of zero
> buoyancy but maybe it's more efficient to just write the environmental profile here and avoid
> the branching of doing the relaxation only up to a specific layer. But have a look at what else
> can be improved. Then have a look at PR 976 and include the entrainment profile also into this
> PR to do both together. And make a change to convection such that if the lowermost layer
> temperature is below freezing then turn the convective rain into snow. Plan with Opus, execute
> with Sonnet.

## Revision log

- **2026-08-22, initial draft.** Plan produced by an Opus sub-agent (given full context: current
  `convection.jl`, `docs/src/convection.md`, the PR #976 diff fetched from GitHub, the
  `large_scale_condensation.jl` snow pattern, GPU kernel-launch conventions in `tendencies.jl`,
  and CLAUDE.md conventions), reviewed by Sonnet. One open design question — whether entrainment
  should default to `NoEntrainment` (bit-identical to today) or PR #976's `LinearEntrainment(0.5)`
  (matches the PR, but changes model behavior for existing users) — was put to the user, who chose
  **`NoEntrainment`**. Plan approved as written; execution proceeds with Sonnet.

## Problem description

The Betts-Miller convection scheme (`SpeedyWeather/src/parameterizations/convection.jl`) runs as
one GPU thread per horizontal grid point (see `column_parameterizations_kernel!` in
`SpeedyWeather/src/parameterizations/tendencies.jl`), each thread serially walking the vertical
column. The reference-profile routines (`pseudo_adiabat!`, `dry_adiabat!`) mark "above the level
of zero buoyancy" (LZB) by writing `NaN`, and `convection!` then re-loops over
`level_zero_buoyancy:nlayers` — a per-column dynamic loop bound — five separate times in the wet
scheme (three in the dry scheme). Dynamic-bound loops per GPU thread are a known real hazard here:
the CHANGELOG already documents `#1193`, an AMDGPU/LLVM miscompile of exactly this pattern
elsewhere in the model.

Four things are needed, bundled into one PR at the user's request:
1. Verify the code against the equations documented in `docs/src/convection.md`.
2. Restructure the scheme to use uniform full-range loops (write the *environmental* profile
   above LZB instead of NaN, so downstream accumulations contribute exact zero there without a
   branch), for GPU performance, and look for any other GPU-unfriendly patterns.
3. Port the entrainment-profile feature from stale PR #976 (by @nviebig; its base predates the
   `physics/` → `parameterizations/` file rename) onto the restructured code.
4. Add: convective rain becomes snow when the lowermost layer's temperature is below freezing.

## Background

See the full plan (as approved) below — reproduced here as the working reference for execution.

<details>
<summary>Approved plan (click to expand)</summary>

### 1. Docs correctness audit (commit 1, docs-only)

Fix `docs/src/convection.md` — code is correct, docs have five gaps:
1. `\tau_{SBM} = 2h` stated as default → actually `Hour(4)` in code. Fix to 4h.
2. `q_{ref} = RH_{SBM}T_{ref}` → should be `q_{ref} = RH_{SBM}\,q^\star(T_{ref})`.
3. Integration-limit sign convention in the `P_q`, `P_T` and final precipitation integrals reads
   backwards for `p_{LZB} < p_0`; either flip the limits or state layer weights are positive.
4. Document the undocumented `max(δq·Δσ, 0)` "no reevaporation" clamp in the precip integral.
5. Document that the `deep_convection` multiplicative flag on precipitation is load-bearing (not
   redundant with the qref correction) because of the `max(·,0)` clamp.

No physics bugs found — the dry/pseudo adiabat, LZB detection, Pq/PT criteria, deep-convection eq
(5)-(6), and shallow qref eq (11)-(15) all match the code exactly.

### 2. Branchless GPU restructuring (commit 2)

Key invariant: `temp`/`humid` passed into `pseudo_adiabat!` *are* `temp_environment`/
`humid_environment` (same array, `convection.jl:46-47,65`), so writing
`temp_ref_profile[ij,k] = temp_environment[ij,k]` above LZB makes every downstream
`temp[ij,k] - temp_ref_profile[ij,k]` exactly `+0.0`, not approximately zero — all accumulations
(Pq, PT, Qref, tendencies, rain) contribute exact zero there with no branch needed.

**`pseudo_adiabat!` (wet):** prefill `temp_ref_profile`/`humid_ref_profile` with the environment
for `k in 1:nlayers` (replaces the NaN fill). Fold the humidity reference into the ascent loop
itself: in both the saturated and dry branches, `saturation_humidity` is already evaluated at
exactly `(temp_ref_profile[ij,k], σ[k]*pres)`, so `RH·sat_humid` there is bitwise identical to
what the old separate fill loop produced — this removes that loop and the `saturation_humidity`
calls it made above LZB. At loop exit, replace the `NaN32` write with
`ifelse(buoyant, ..., temp_environment[ij,k])` / same for humidity. Signature gains
`humid_ref_profile`, `relative_humidity`, and `entrainment` (§3) as arguments.

**Wet `convection!`:** the five dynamic-bound loops collapse to two full-range loops:
- **Loop A** (fuses old Pq/PT loop + old Qref loop): `for k in 1:nlayers`, accumulate
  `Pq += (humid - humid_ref)*Δσ`, `PT -= (temp - temp_ref)*Δσ` (exact zero above LZB, no mask
  needed), and `Qref -= ifelse(k >= lzb, humid_ref, 0)*Δσ` (mask **required** here — Qref must
  only sum the in-column part).
- Compute `deep_convection`/`shallow_convection`, early-return on `no_convection` (kept — cheap
  win, see limitations), then `Δσ_lzb`, and branchless scalar `ΔT = ifelse(deep, ..., ...)`,
  `fq = ifelse(deep, 1, 1 - Pq/Qref)`.
- **Loop B** (fuses old deep-correction + shallow-correction + tendency + precip loops):
  `for k in 1:nlayers`, compute `T_ref = temp_ref[k] - ifelse(k>=lzb, ΔT, 0)`,
  `q_ref = humid_ref[k] * ifelse(k>=lzb, fq, 1)` **without writing them back** to the scratch
  arrays, then apply tendencies and accumulate `rain_convection` exactly as today (all exact zero
  above LZB).
- Do **not** naively make the old *correction-writing* loops full-range — that would leave a
  spurious `ΔT`/`fq` baked into the reference above LZB and inject fake tendencies there. The fix
  is to never write the corrected profile back; only apply it inline in loop B.

**`dry_adiabat!` + dry `convection!`:** same treatment — prefill with environment, one full-range
PT-accumulation loop, one full-range tendency loop with inline `ifelse`-masked correction, no
write-back. Also hoist the per-level `/DBM.time_scale.value` division to a precomputed `τ⁻¹`
(matching the wet scheme, minor fix).

Also fix while touching this file: drop the stray `NaN32`/`NF(NaN)` type mismatch (subsumed by
the restructuring), drop the now-redundant `local` declarations at loop scope.

**Known, accepted remaining limitation:** the ascent `while buoyant && k > 1` loop keeps a
per-column dynamic trip count — that's physically the LZB search and can't be avoided without a
different algorithm. Document this as an explicit out-of-scope limitation, not a bug.

### 3. Entrainment profiles, adapted from PR #976 (commit 3)

Port PR #976 (by @nviebig, `nv/entrainment-profiles`, stale — its base predates the
`physics/` → `parameterizations/` rename) onto the restructured code:

- New `SpeedyWeather/src/parameterizations/entrainment.jl`, included in `SpeedyWeather.jl` right
  before `convection.jl`. Defines `AbstractEntrainment <: AbstractParameterization`,
  `NoEntrainment` (singleton, `(::NoEntrainment)(σ) = 0`, **default**), `LinearEntrainment{NF}`
  (ramps from `surface_entrainment` at σ=1 to 0 at `σ_entrainment`), `ConstantEntrainment{NF}`
  (flat rate) — matching PR #976's math exactly. Add an `entraining(::AbstractEntrainment)` trait
  (`false` for `NoEntrainment`, `true` otherwise) so the mixing block compiles away entirely for
  the default case.
- Add `entrainment::Entrainment <: AbstractEntrainment` field (`@component`, default
  `NoEntrainment()`) to both `BettsMillerConvection` and `BettsMillerDryConvection`, threaded
  through to `pseudo_adiabat!`/`dry_adiabat!`. Promote `BettsMillerDryConvection` from plain
  `@kwdef` to `@parameterized @kwdef` (needed for `@component`), replacing its hand-written
  `adapt_structure` with `Adapt.@adapt_structure`.
- Wet ascent: inside the `if saturated` branch, after computing `humid_parcel`, guarded by
  `if entraining(entrainment)`: mix `temp_parcel`/`humid_parcel` toward
  `temp_environment[ij,k]`/`humid_environment[ij,k]` by `ε = entrainment(σ[k])`, then **recompute**
  `sat_humid = saturation_humidity(temp_parcel, σ[k]*pres, atmosphere)` so the folded humidity
  reference (§2) stays exactly correct after mixing. Dry ascent: same idea, temperature only, in
  `dry_adiabat!`.
- New "Entrainment" section in `docs/src/convection.md` (PR #976 added none): physical motivation,
  the mixing equations, the three profiles, the resolution-dependence caveat (ε is per-model-level,
  not per unit height), and that `NoEntrainment` is the default so behavior is unchanged unless a
  user opts in.
- Credit PR #976 and @nviebig in the CHANGELOG entry and PR description. **Do not** comment on or
  close #976 — that's left to the user's judgement, not automated here.

### 4. Convective rain → snow below freezing (commit 4)

Mirror the existing pattern in `large_scale_condensation.jl` (`snow::Bool`, `@param
freezing_threshold`), but simpler — a single check, no falling-flux melt cascade (convective
precip is deposited immediately, not fluxed through layers):

- Add to `BettsMillerConvection`: `snow::Bool = true`,
  `@param freezing_threshold::NF = 273.15` — a **new** field, not a reuse of
  `atmosphere.temperature_freezing` (that constant is the Clausius-Clapeyron reference
  temperature; overloading it as a precip-phase switch would couple two unrelated knobs). Default
  273.15 K matches "below freezing" literally (LSC's 263 K default is for a different, cascading
  check and isn't the right anchor here).
- At the end of `convection!` (wet only — dry scheme has no precipitation), after computing
  `rain_convection`, phase-swap by `temp[ij, nlayers]` (verified: layer `nlayers` is the surface
  layer, per `temp_environment[ij, nlayers]` as the ascent's surface start point) against
  `freezing_threshold`, same swap idiom as LSC: `rain, snow = cond ? (snow, rain) : (rain, snow)`.
  Accumulate into new `snow_convection`/`snow_rate_convection`, and into the shared
  `vars.parameterizations.snow_rate[ij]` instead of `rain_rate[ij]` when snowing.
- `variables(::BettsMillerConvection)`: add `ParameterizationVariable` entries for
  `snow_convection`, `snow_rate_convection`, `snow_rate`.
- `SpeedyWeather/src/output/variables/precipitation.jl`: add `ConvectiveSnowOutput` /
  `ConvectiveSnowRateOutput`, copying the existing `LargeScaleSnowOutput`/
  `LargeScaleSnowRateOutput` pattern; register in `PrecipitationOutput()`.
- Fix latent bug found during the audit: `reset_variables!` in
  `SpeedyWeather/src/parameterizations/tendencies.jl` resets `rain_rate`/`snow_rate` but not
  `rain_rate_convection`, which is only assigned on the non-early-return path — a column that
  stops convecting keeps a stale rate. Add `rain_rate_convection` and (new) `snow_rate_convection`
  to `reset_variables!`.
- New docs subsection under "Convective precipitation": the check, mass conservation (rain+snow
  unchanged), no melt cascade (contrast with large-scale), `snow=false` disables it.

### 5. Benchmarking — local only, not committed to the repo

Not a `manual_benchmarking.jl`/README.md entry — a throwaway script (e.g. `/tmp/bench_convection.jl`,
run via `julia --project=SpeedyWeather/benchmark`). Uses
`SpeedyWeather._column_parameterizations_cpu!(vars, parameterizations_namedtuple, model)`, which
takes the parameterization set as an argument, so isolation needs no model surgery:

- **A. Isolation**: `_column_parameterizations_cpu!(vars, (convection = model.convection,), model)`.
- **B. Combined**: full `column_parameterizations!(vars, model)` with real physics.
- **C. Combined, no convection**: same with `convection = nothing`, for a B−C delta.
- Spin the model forward first (e.g. 100 steps) before timing so convection is actually active in
  some columns; record `count(cloud_top .< nlayers+1)` alongside timings so before/after runs are
  comparable. Sweep `nlayers ∈ (8, 16, 32)` at `truncation=31`, plus one `truncation=63, nlayers=8`
  point (this axis is what exposes the "full-range loops do more work for shallow-LZB columns"
  risk). Assert zero allocations. Run identically on `main` and on this branch.
- If GPU hardware is available locally, repeat A/B via `launch!(arch, LinearWorkOrder, (npoints,),
  column_parameterizations_kernel!, ...)` + `synchronize`; otherwise treat CPU timings as the
  primary signal and say so explicitly when reporting results.
- Separately, the real acceptance gate: run 200 steps on `main` vs. this branch (default
  `NoEntrainment`, so no simulated-behavior change) and diff `vars.prognostic` fields with `==`
  (exact, not `≈` — every changed term is an exact zero or a bitwise-identical reuse). Do this
  right after §2 (before entrainment/snow are added) and again after §3 to confirm entrainment's
  default path is still bit-identical.

### 6. Repo conventions

- Plan doc for this work goes to `docs/dev/2026-08/betts-miller-convection-overhaul.md` following
  the CLAUDE.md template, status updated as work lands.
- CHANGELOG.md, three bullets under `## Unreleased` (PR number filled in once opened): branchless
  loops, entrainment (crediting #976/@nviebig), convective snow.
- `SpeedyWeather/Project.toml`: bump `0.22.1` → `0.23.0-DEV` (new exported public types, new type
  parameter on `BettsMillerConvection`/`BettsMillerDryConvection`, new output variables — a
  breaking change per CLAUDE.md's versioning rule, even though keyword constructors stay
  source-compatible).
- Run Runic on every changed file before committing.

### 7. Testing (`--check-bounds=yes`)

Extend `SpeedyWeather/test/parameterizations/convection.jl` with sibling testsets (keep the
existing `Convection × Model` sweep as-is):
1. Entrainment functor unit tests (`NoEntrainment`/`LinearEntrainment`/`ConstantEntrainment`
   values, monotonicity, `@inferred`).
2. Short `run!` per entrainment type for wet+dry models (6 runs, not a full cross product).
3. Entrainment has a measurable effect (tendencies differ vs. `NoEntrainment`, cloud top no
   higher).
4. Reference-profile invariants directly on `pseudo_adiabat!`/`dry_adiabat!`: no NaNs, exact
   equality to environment above LZB, valid LZB range — this is the regression guard for the
   restructuring itself.
5. Snow: force a cold/warm surface layer, assert rain/snow route correctly and
   `rain_rate_convection + snow_rate_convection` is invariant to the `snow` toggle (mass
   conservation); assert new variables are registered.
6. Run `SpeedyWeather/test/parameterizations/all_parametrizations.jl` (NaN-freedom check, picks up
   new snow fields automatically) and `SpeedyWeather/test/output/netcdf_output.jl` (already
   exercises `PrecipitationOutput()`).
7. Manually run `test/parameters.jl` and `test/type_stability_test.jl` (excluded from `Pkg.test`
   by `runtests.jl`, but touched by the new `@param`s and the `BettsMillerDryConvection`
   `@parameterized` promotion).
8. The bit-identity gate from §5 is the final correctness sign-off, run and recorded in the plan
   doc, not a unit test.

### 8. Ordered steps

1. Write plan doc to `docs/dev/2026-08/`, status **planned**.
2. Commit: docs audit fixes only (§1).
3. Baseline benchmark + capture 200-step reference state on `main` (§5).
4. Commit: branchless restructuring (§2). Run bit-identity gate + convection/all_parametrizations
   tests.
5. Post-restructure benchmark; report before/after. Stop and report if CPU regresses >20%.
6. Commit: entrainment (§3). Re-run bit-identity gate (must still be exact with default
   `NoEntrainment`).
7. Commit: convective snow (§4).
8. Commit: CHANGELOG + version bump + Runic + plan doc status → **completed** with recorded
   benchmark numbers.
9. Full `Pkg.test("SpeedyWeather")` with `--check-bounds=yes`.

</details>

## Summary of changes

_To be filled in as work lands._

## Testing and verification

_To be filled in as work lands._

## Documentation changes

_To be filled in as work lands._

## Known limitations

- The ascent `while buoyant && k > 1` loop in `pseudo_adiabat!`/`dry_adiabat!` keeps a per-column
  dynamic trip count — inherent to the algorithm (the LZB search), not addressed by this PR.
- Entrainment mixing is applied in the saturated branch only (wet scheme), inherited from PR #976.
- Convective snow uses a single lowest-layer temperature check, no falling-flux melt cascade
  (unlike large-scale precipitation).

## Future work

- If GPU profiling later shows the ascent while-loop is still a divergence hotspot, a follow-up
  plan could explore always-ascend-to-top-and-mask strategies.
- Extending entrainment mixing to the dry-ascent branch of the wet scheme's pseudo-adiabat.
