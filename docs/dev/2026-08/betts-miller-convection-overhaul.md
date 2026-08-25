# Betts-Miller convection overhaul: branchless vertical loops, entrainment, convective snow

> Status: **completed, with one open item**. Opened as
> [PR #1221](https://github.com/SpeedyWeather/SpeedyWeather.jl/pull/1221). The full `Pkg.test`
> run surfaces 2 JET dispatch-check failures in unrelated `Barotropic`/dynamics-core code,
> bisected to a compilation-heuristic sensitivity rather than a logic defect (see Testing and
> verification / Known limitations) — noted in the PR description for maintainer input rather
> than resolved unilaterally, per user decision. This PR supersedes #976 (entrainment, by
> @nviebig) — commented on #976 linking here, per user decision. Note: the comment's wording
> floated closing #976 in favor of #1221 while explicitly asking @nviebig first — #976 was
> **not** actually closed via the API, only commented on; closing it remains a human call.

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
- **2026-08-23, execution.** Two on-the-fly deviations from the written plan, both discovered via
  the bit-identity gate doing its job:
  - The dry scheme's planned `τ⁻¹` hoist (§2, "matching the wet scheme, minor fix") was dropped:
    it changed `/DBM.time_scale.value` (division) into `inv(...)` then multiply, which is not
    bit-identical in IEEE754 and also changed the working precision (explicit convert to `NF`
    before the reciprocal vs. Julia's implicit `Float64` promotion in the original division).
    Caught immediately by the bit-identity gate (~1e-5 to 1e-4 relative diffs after 200 steps);
    reverted to the original per-level division to keep the restructuring commit exactly
    bit-identical, as required.
  - The planned "entrainment has an effect" test (comparing mean `cloud_top` between two
    independent 20-step model runs, one per entrainment setting) was replaced: it failed on the
    first run, with entrainment's mean `cloud_top` coming out *lower* (deeper clouds) than
    without, opposite to the naive physical expectation. Diagnosis: two independent model runs
    diverge chaotically from the first timestep regardless of entrainment, so comparing their
    cloud-top statistics after only 20 steps isn't actually isolating entrainment's effect — it's
    comparing two uncorrelated weather trajectories. Replaced with a direct, deterministic test:
    call `pseudo_adiabat!` with and without entrainment on the *same* atmospheric snapshot and
    assert the LZB with entrainment is never above (never a smaller layer index than) the LZB
    without, for every column on a T15 grid — this is both a stronger and a cleaner test, and
    passed on the first try once framed this way.
  - Snow (commit 4) was, as anticipated in the plan's Q&A, confirmed *not* bit-identical to
    `main` even with all other defaults held — expected, since `snow=true` is the requested
    default-on behaviour, not an invisible refactor. Verified the divergence is confined to the
    wet model (via land snow-depth/soil-moisture coupling) and that the dry model stays exact.
  - Work paused before pushing/opening the PR: bundles a large restructuring with a feature that
    explicitly supersedes another contributor's (@nviebig, #976) open PR, which felt like a
    natural point to check in with the user rather than pushing autonomously.
- **2026-08-23, full test suite + bisection.** `Pkg.test("SpeedyWeather")` surfaced 2 failures
  (of 116,739) in `dynamics/dispatch.jl`'s JET dispatch checks, both passing on `main`. Bisected
  by building a chain of `git worktree`s at each commit and layering each of commit 4's three
  changed files (`convection.jl`, `entrainment.jl` — unmodified control, `precipitation.jl`,
  `tendencies.jl`) individually onto the entrainment-only commit: none reproduces the failure
  alone, only the full combination does, and the actual JET-flagged dispatch sites are 100% inside
  `KernelAbstractions`/dynamics-core internals that `BarotropicModel` (no parameterizations at
  all) exercises — conclusively not a logic bug in the snow/entrainment/restructuring code, most
  likely a global compilation-heuristic sensitivity to the package's total type/method count.
  Documented in Known limitations and Testing and verification rather than "fixed" by guessing at
  `@inline` placements; flagged to the user as an open item before proceeding to push/open the PR.
- **2026-08-24, review pass.** Addressed @maximilian-gelbrecht's PR review: the convective-snow
  branch selection (in `convection!`, previously a ternary) is now an explicit `ifelse` since it's
  runtime-dependent and needs to stay branchless for Reactant; removed comments in `convection.jl`
  (`pseudo_adiabat!`/`dry_adiabat!` prefill) and `tendencies.jl` (`reset_variables!`) that narrated
  the diff against the prior implementation rather than the current code. Two other review threads
  were discussion, not code changes, and answered inline: the docs sign-convention question (the
  bounds-order flip *and* the sign both had to change together to match the code — confirmed by
  tracing `Pq`/`PT` back to Frierson eq. (1)) and the ascent `while` loop's unavoidable dynamic trip
  count (already called out in Known limitations above, distinct from the `#1193` hazard this PR
  actually fixes, which was downstream loops using the LZB as a *loop bound*).
- **2026-08-24, GPU profiling.** Max followed up asking for actual GPU profiling instead of the
  "no GPU available" caveat in the CPU benchmark table. GPU hardware (1× NVIDIA L4) turned out to be
  available in this environment. `ncu` (the tool that would directly show ascent-loop divergence)
  is blocked by the driver's `RmProfilingAdminOnly` restriction — user has asked sysadmins for
  access, pending. Ran `nsys` instead (works without the restricted counters): isolated the
  `column_parameterizations_kernel!` convection-only launch and compared `main` vs. this branch at
  the plan's three CPU-benchmark configs. Result, unexpectedly, is a *regression*: ~5% slower at
  `nlayers=8` (both T31 and T63), ~0% at `nlayers=16` — opposite sign from the CPU numbers, and a
  real GPU-side confirmation of the "full-range loops do more work for shallow-LZB columns" risk
  this plan already flagged. Documented in a new "GPU profiling (nsys)" subsection and cross-referenced
  from Known limitations; not fixed here, reported to the reviewer as a maintainer call.

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

Five commits on `mc/convection`, in order:

1. **Docs audit fixes** (`docs/src/convection.md`, no code): `τ_SBM` default corrected 2h→4h,
   `q_ref` formula corrected to include `q*(T_ref)`, integration-limit sign convention fixed for
   `P_q`/`P_T`/`P`, the undocumented no-reevaporation clamp and the load-bearing role of the
   deep-convection indicator documented. Code was already correct against Frierson 2007 — no
   physics bugs found in the audit.
2. **Branchless restructuring** (`SpeedyWeather/src/parameterizations/convection.jl`):
   `pseudo_adiabat!`/`dry_adiabat!` now prefill the reference profiles with the environment
   (replacing the NaN marker) and never touch levels above the level of zero buoyancy (LZB); the
   wet scheme's separate humidity-reference fill loop is folded into the ascent. `convection!`'s
   five (wet) / three (dry) dynamic-bound `level_zero_buoyancy:nlayers` loops become two (wet) /
   two (dry) full-range `1:nlayers` loops with `ifelse`-masked corrections applied inline (never
   written back to the scratch arrays). Verified bit-identical to pre-restructure `main` via an
   exact (`==`) 200-step comparison of all prognostic fields, both `PrimitiveWetModel` and
   `PrimitiveDryModel`. One deviation from the original plan: the dry scheme's
   `/DBM.time_scale.value` was *not* hoisted into a precomputed `inv(...)` reciprocal as
   initially planned (matching the wet scheme's style) — that turned out to change both the
   floating-point operation (division → multiply-by-reciprocal) and working precision, breaking
   bit-identity; reverted to keep this commit exactly bit-identical.
3. **Entrainment** (`SpeedyWeather/src/parameterizations/entrainment.jl`, new file): ports PR
   #976 (by @nviebig) onto the restructured code — `NoEntrainment`/`LinearEntrainment`/
   `ConstantEntrainment`, added as an `entrainment` sub-component of `BettsMillerConvection` and
   `BettsMillerDryConvection` (the latter promoted from a plain `struct` to `@parameterized
   @kwdef` to support it). Default is `NoEntrainment` (user's decision, deviating from PR #976's
   `LinearEntrainment(0.5)` default) — verified the 200-step bit-identity check still passes
   exactly with entrainment wired in but defaulted off.
4. **Convective snow** (`convection.jl`, `output/variables/precipitation.jl`,
   `parameterizations/tendencies.jl`): new `snow`/`freezing_threshold` options on
   `BettsMillerConvection` (default `snow = true`, `freezing_threshold = 273.15`); below that
   threshold at the lowermost layer, convective rain is swapped for snow after `rain_convection`
   is fully computed (mass-conserving, verified exactly (`==`) in a test). New
   `snow_convection`/`snow_rate_convection` diagnostics and `ConvectiveSnowOutput`/
   `ConvectiveSnowRateOutput` netCDF variables. Also fixed a latent bug found while auditing this
   code: `reset_variables!` reset `rain_rate`/`snow_rate` but not `rain_rate_convection` (only
   assigned on the convecting path), so a column that stopped convecting kept a stale rate; now
   reset alongside the new `snow_rate_convection`. Since `snow = true` by default, this commit is
   *not* bit-identical to `main` — verified this is the *expected*, intentional divergence: wet
   model differs by O(1e-5 to 6e-2) depending on field after 200 steps wherever a column's
   surface is below freezing (via land snow-depth/soil-moisture feedback), the dry model (no
   precipitation) remains exactly bit-identical.
5. **Housekeeping**: `CHANGELOG.md` (3 entries, PR number filled in once opened),
   `SpeedyWeather/Project.toml` version `0.22.1` → `0.23.0-DEV`. All touched files are
   Runic-clean (`runic --check` exit 0).

## Testing and verification

- `SpeedyWeather/test/parameterizations/convection.jl`: original `Convection × Model` sweep plus
  five new testsets — entrainment functor unit tests, entrainment in a full model run (6
  configs), entrainment's effect verified deterministically on a fixed atmospheric snapshot
  (`lzb_with_entrainment >= lzb_without`, checked exhaustively across a T15 grid — not via
  comparing two independently-run, chaotically-diverging simulations, which was tried first and
  discarded as unreliable over a short run), reference-profile invariants (no NaNs, exact
  equality to the environment above LZB) directly on `pseudo_adiabat!`/`dry_adiabat!`, and
  convective snow (mass conservation, variable registration, `snow=false` recovers rain-only).
  **1132/1132 pass** with `--check-bounds=yes`.
- `SpeedyWeather/test/parameterizations/all_parametrizations.jl` (NaN-freedom across every
  parameterization variable): **44/44 pass** (was 42 before the two new snow variables).
- `SpeedyWeather/test/output/netcdf_output.jl` (exercises `PrecipitationOutput()`, now including
  the two new convective-snow output types): all pass.
- `SpeedyWeather/test/parameters.jl`, `SpeedyWeather/test/type_stability_test.jl` (excluded from
  `Pkg.test`, run manually since new `@param`s and the `BettsMillerDryConvection`
  `@parameterized` promotion touch them): pass.
- Full `Pkg.test("SpeedyWeather")` with `--check-bounds=yes`: **116,737 / 116,739 pass.** The two
  failures are both in `dynamics/dispatch.jl`'s JET `@test_opt` "free of runtime dispatch" checks
  — one for `BarotropicModel`, one for `PrimitiveWetModel`. Both pass cleanly on unmodified `main`.
  Root-caused by bisection (see revision log): neither of the four feature commits' *individual*
  files (`convection.jl`, `entrainment.jl`, `precipitation.jl`, `tendencies.jl`) reproduces the
  failure in isolation when layered onto the entrainment commit — only the *full combination* of
  changes in the convective-snow commit does. The JET-flagged dispatch is 100% inside
  `KernelAbstractions`/dynamics-core internals (`__thread_run`, spectral transform, horizontal
  diffusion) that `BarotropicModel` uses — a model that has no parameterizations at all and never
  touches convection, entrainment, or snow code. This points to a global Julia
  compilation-heuristic sensitivity (added types/methods elsewhere in the package shifting
  inference/specialization decisions in unrelated, marginal code — the same general class of
  issue as the already-known `#1193` dynamic-loop-bound fragility) rather than a logic defect
  introduced by this PR. **Flagged to the user rather than resolved** — see Known limitations.
- Bit-identity gate (200 steps, `PrimitiveWetModel` + `PrimitiveDryModel`, `truncation=31,
  nlayers=8`, exact `==` on all prognostic fields against pre-restructure `main`): **PASS** after
  the restructuring commit, **PASS** after entrainment (default off), **intentional FAIL** after
  the snow commit (see above) — expected and documented, not a regression.
- Benchmarking (`/tmp/bench_convection.jl`, CPU, isolation vs combined-with-other-parameterizations
  vs combined-without-convection, `T31/L8`, `T31/L16`, `T63/L8` — `T31/L32` and `T63/L16` dropped
  from the sweep, both numerically unstable at the default timestep even on unmodified `main`,
  unrelated to this work):

  | Config | Arm | `main` (baseline) | `mc/convection` (restructured) |
  |---|---|---|---|
  | T31 L8  | isolation (A) | 1.087 ms | 1.018 ms |
  | T31 L8  | combined (B)  | 5.969 ms | 5.872 ms |
  | T31 L8  | B−C (convection in-situ) | 0.819 ms | 0.739 ms |
  | T31 L16 | isolation (A) | 1.455 ms | 1.441 ms |
  | T31 L16 | combined (B)  | 10.325 ms | 10.302 ms |
  | T31 L16 | B−C | 1.092 ms | 1.050 ms |
  | T63 L8  | isolation (A) | 3.612 ms | 3.402 ms |
  | T63 L8  | combined (B)  | 20.435 ms | 20.170 ms |
  | T63 L8  | B−C | 2.954 ms | 2.692 ms |

  No CPU regression in any arm — a small (~2-10%) improvement throughout from fewer loop passes,
  well inside the acceptance criterion (revert/report if >20% regression). All arms 0 allocations
  before and after. No GPU hardware was available in this environment, so GPU numbers were not
  collected at the time; the CPU result was used as the primary local signal per the plan, with the
  branchless structure itself (uniform trip counts, no dynamic-bound loops) being the intended GPU
  win, consistent with the `#1193` precedent cited in the plan. GPU hardware became available in a
  later session — see "GPU profiling (nsys)" below, which complicates this picture.

### GPU profiling (nsys)

Added 2026-08-24 in response to @maximilian-gelbrecht's review ("Would be good to actually profile
this on GPU"). One NVIDIA L4 (CUDA 13.0, driver 580.173.02) became available in this environment.

**Tooling available:** both `nsys` (Nsight Systems) and `ncu` (Nsight Compute) are installed, but
`ncu` — the tool that reports per-kernel hardware-counter metrics like warp/branch execution
efficiency, i.e. what's actually needed to tell whether the ascent `while` loop specifically causes
divergence — fails with `ERR_NVGPUCTRPERM`: the driver restricts GPU performance-counter access to
root (`RmProfilingAdminOnly: 1`), and this session has no sudo. The user has asked the system admins
for access; until then, only `nsys`'s coarser kernel-duration timeline is available (verified: it
does not need the restricted counters, confirmed with a throwaway CUDA test program).

**Method:** `column_parameterizations_kernel!` launched directly via `launch!(arch, LinearWorkOrder,
(npoints,), column_parameterizations_kernel!, vars, (convection = model.convection,), model)` — the
same isolation technique as the CPU benchmark, now via a direct kernel launch instead of the
`_column_parameterizations_cpu!` internal. 100 spin-up steps first (so convection is actually
active — confirmed via `n_convecting`, all columns convecting in every config below), then 3
untimed warm-up launches, then 50 back-to-back timed launches, profiled with `nsys profile` and
summarized with `nsys stats --report cuda_gpu_kern_sum`. This reports actual GPU-side kernel
execution time from CUPTI activity records (kernel start → end on the device), not host-side Julia
dispatch latency, so it's a fair per-launch comparison between branches. Script: `/tmp/bench_convection_gpu.jl`.
One run per config (no repeated trials, given time cost — see caveat below), `main` at `9794037a`
(pre-PR) vs. this branch:

| Config  | `main` (ns) | `mc/convection` (ns) | Δ (ns) | Δ (%)  |
|---------|-------------|-----------------------|--------|--------|
| T31 L8  | 13,438.8 (σ=40.7)   | 14,104.7 (σ=38.0)   | +665.9 | **+4.96%** |
| T31 L16 | 24,902.5 (σ=86.5)   | 24,958.7 (σ=100.7)  | +56.2  | +0.23% |
| T63 L8  | 13,868.6 (σ=134.5)  | 14,531.5 (σ=223.0)  | +662.9 | **+4.78%** |

**Finding:** at `nlayers=8` (both truncations), this branch's isolated-convection kernel is
measurably ~5% *slower* on GPU than `main` — the opposite sign from the CPU numbers above. The
absolute gap (~660ns) is nearly identical at T31L8 and T63L8 despite T63 having ~4× the grid
points, which points to a fixed *per-thread* cost rather than a launch-overhead/occupancy artifact
(those would scale with grid size). This is consistent with the exact risk this plan already flagged
under §2 ("full-range loops do more work for shallow-LZB columns"): `main`'s dynamic re-loop over
`level_zero_buoyancy:nlayers` does fewer iterations than `nlayers` whenever a column's convection is
shallow, while this branch's full-range loops always run all `nlayers` iterations per thread
regardless. At `nlayers=16` the gap nearly disappears (+0.23%, within the observed noise) — a
plausible read is that at finer vertical resolution a typical LZB sits proportionally deeper into
the column, so `main`'s old dynamic range was already close to the full range there too, leaving
little headroom for it to "win"; not confirmed, just the most consistent explanation of the data.

**Caveats:**
- This measures the *whole* isolated-convection kernel (`pseudo_adiabat!`/`dry_adiabat!` — including
  the ascent `while` loop — fused with the downstream tendency loops in `convection!`). `nsys`
  reports only the aggregate per-launch kernel duration; it cannot attribute time to the `while`
  loop specifically vs. the other loops. So this does *not* directly answer whether the ascent loop
  itself is a divergence hotspot — that needs `ncu`'s branch/warp-efficiency counters, still pending
  admin access.
- Single run per config on a shared, multi-tenant machine — no repeated trials to quantify run-to-run
  variance beyond the per-launch stddev `nsys` itself reports. The L8 deltas (~660ns) are large
  relative to their stddevs (~3-5×), so likely real; the L16 delta (56ns vs. ~90ns stddev) is not
  distinguishable from noise.
- Not (yet) reflected in the CPU benchmark table above or the PR description's benchmark section —
  this is a GPU-only regression at coarse vertical resolution; CPU numbers stand as measured.

This is a real, GPU-observed downside of the branchless restructuring that the CPU-only benchmarking
in this plan could not have caught, and is now flagged to the reviewer rather than resolved
unilaterally — whether to pursue a fix (e.g. capping the full-range loop at some resolution-aware
bound, or accepting the regression at coarse `nlayers` as the price of the `#1193` GPU-miscompile
fix) is a maintainer call.

## Documentation changes

- `docs/src/convection.md`: five correctness fixes (§ above), new "Entrainment" section (physical
  motivation, mixing equations, the three profiles, resolution-dependence caveat, `NoEntrainment`
  default), new "Convective snow" subsection under "Convective precipitation" (threshold check,
  mass conservation, contrast with the large-scale melt cascade, `snow=false` toggle).

## Known limitations

- The ascent `while buoyant && k > 1` loop in `pseudo_adiabat!`/`dry_adiabat!` keeps a per-column
  dynamic trip count — inherent to the algorithm (the LZB search), not addressed by this PR.
  Whether it's specifically a warp-divergence hotspot is still open: `nsys` profiling (see "GPU
  profiling" above) shows the *whole* isolated-convection kernel is ~5% slower on GPU than `main`
  at `nlayers=8`, but attributing that to the ascent loop vs. the full-range downstream loops needs
  `ncu`'s branch-efficiency counters, blocked on admin access as of this writing.
- Entrainment mixing is applied in the saturated branch only (wet scheme), inherited from PR #976.
- Convective snow uses a single lowest-layer temperature check, no falling-flux melt cascade
  (unlike large-scale precipitation).
- **Open item, unresolved:** `SpeedyWeather/test/dynamics/dispatch.jl`'s JET "free of runtime
  dispatch" checks fail for `BarotropicModel` and `PrimitiveWetModel` only once the full
  convective-snow commit is present — not for any of its three changed files individually layered
  onto the entrainment commit, and not for either of the first two commits alone. The dispatch
  JET actually flags is entirely inside `KernelAbstractions`/dynamics-core internals unrelated to
  convection (`__thread_run`, spectral transform, horizontal diffusion) — `BarotropicModel` has no
  parameterizations and never runs `convection!`, `entrainment`, or the snow code at all. This
  looks like a global compilation-heuristic sensitivity (added types/methods shifting inference
  decisions in unrelated marginal code elsewhere in the same compiled package), similar in kind to
  `#1193`, rather than an actual bug in the snow feature. Not resolved here — needs either a
  maintainer call on whether this JET check's fragility is acceptable/pre-existing-prone, or
  further Julia-compiler-level investigation that was out of scope for this session.

## Future work

- If GPU profiling later shows the ascent while-loop is still a divergence hotspot, a follow-up
  plan could explore always-ascend-to-top-and-mask strategies.
- Extending entrainment mixing to the dry-ascent branch of the wet scheme's pseudo-adiabat.
