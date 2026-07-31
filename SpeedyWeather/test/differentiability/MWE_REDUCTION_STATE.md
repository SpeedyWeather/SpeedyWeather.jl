# L30 (2026-07-07) — full_diff_CI.jl HANGS in Enzyme compile on 1.12 (full time_step!), passes on 1.10

`full_diff_CI.jl` = parameter AD (`Duplicated(model)` + `Duplicated(model.time_stepping)`) of the
WHOLE `time_step!` on a full `PrimitiveWetModel` (trunc5 nlayers1, after a 6h run).
- **1.10: PASSES** (1/1, 12.7 min compile, nonzero gradient).
- **1.12: HANGS in the Enzyme compile.** Confirmed a hang, NOT slow compile: RSS drops 194→113→36 MB
  while CPU stays ~80–99% for 60+ CPU-min with no output past ">>> entering autodiff" (a real
  compile grows to GBs). `maxtypeoffset!(4096)` does NOT help (identical churn). Even **Const(model)
  state AD** of the full `time_step!` hangs the same way → not specific to parameter AD.

**Not caused by the transform! rule per se:** `_mwe5` (transform!(::Variables), Const model) compiles
+ runs on 1.12; and full_diff_CI compiles fine on 1.10 with the identical rule. It is a 1.12-specific
Enzyme compile pathology on the *combined* full `time_step!`.

**Key difference from the passing pwet dyn-core tests:** those use `dynamics_only = true` (NO physics
parameterizations) and test components individually; full_diff_CI runs the FULL `time_step!` WITH
physics + implicit + leapfrog chained. Hypothesis under test (`_ts_dynonly.jl`): the hang is in the
physics-parameterization reverse (never exercised on 1.12) or the full chain, not the dyn-core alone.

⭐ **CI ERROR OBTAINED (2026-07-07) — Enzyme codegen ASSERTION, not a hang:**
```
Enzyme/DiffeGradientUtils.cpp:433 addToDiffe(...): Assertion `0 && "unhandled accumulate with partial sizes"'
  addToDiffe → visitCommonStore (AdjointGenerator.h:1219) → visitStoreInst → CreatePrimalAndGradient
  VT: i64 ... start=0 size=4 storeSize=8 ... extractvalue [3 x i64] %2299, 1
signal 6 (Aborted), exit 134
```
So the CI (assertions-enabled Enzyme_jll) ABORTS while generating the reverse of a STORE — Enzyme
can't accumulate a gradient into a value with "partial sizes" (a `[3 x i64]` aggregate; storing a
4-byte piece into an 8-byte slot). My local **release** Enzyme (0.13.173, assertions OFF) hits the
SAME bug but HANGS/UBs instead of asserting — hence "hang locally, fast abort in CI".

Diagnosis:
- **Upstream Enzyme codegen bug**, 1.12-specific (1.10 passes in 12.7 min). Not a type-analysis
  issue (maxtypeoffset irrelevant).
- In **natively-differentiated dyn-core** code (dynamics_only variant also hangs → NOT physics; and
  it's NOT the `transform!` rule, which is an AD boundary and passes on 1.10 with identical code).
- The `[3 x i64]` size aggregate + "partial sizes" store strongly suggests the **fused view-backed
  prognostic** (`LTA{SubArray}` of the 3-D fused parent) being differentiated through the leapfrog/
  `get_step` path — i.e. the #1099 fused-variables layout, not this PR's Enzyme rules. Consistent
  with the separate `Duplicated(view-LTA, Array-LTA)` make_zero error seen in barotropic
  `update_prognostic!`.
- Local reproduction is a HANG (release Enzyme), which makes SpeedyWeather-side bisection painful
  (can't see the clean assertion). CI (debug Enzyme) is the clean signal.

**Next options (need user steer):** (a) build a reduced MWE — differentiate `update_prognostic!`/
leapfrog on fused view-backed prognostic on 1.12 — to confirm the trigger and report upstream to
Enzyme with the assertion; (b) restructure the leapfrog/`get_step` on fused vars to avoid the
partial-view store (SpeedyWeather-side workaround, if localizable); (c) keep `full_diff_CI` gated
off 1.12 until Enzyme fixes it. maxtypeoffset restore is NOT the fix (offset didn't help).

### ⭐ ROOT CAUSE + FIX (2026-07-07) — view-preserving `make_zero(::Variables)`
User chose (b) restructure. Decisive test (`_viewshadow_test.jl`): building the shadow as a
**view-preserving** deepcopy (zeroed via the fuse parents, so the prognostic shadow stays
`LTA{SubArray}` matching the primal) makes `update_prognostic!` (leapfrog `get_step` over the fused
prognostic) **compile + autodiff-return on 1.12** — no hang, no assertion. Confirmed the root cause:
Enzyme's DEFAULT `make_zero` materialises the fused view-backed prognostic shadow to a plain `Array`,
so the shadow type mismatches the view-backed primal → the `get_step` store reverse over the
mismatched view/Array types is what hit `Duplicated(view,Array)` (barotropic) and "unhandled
accumulate with partial sizes" (full time_step! / full_diff_CI on 1.12).

**Fix (SpeedyWeatherEnzymeExt.jl):** override `Enzyme.make_zero(::Variables)` = `deepcopy` (preserves
view→fuse-parent aliasing + types) then zero the differentiable data via the fuse parents / non-view
leaves (mirrors the `to_vec(::Variables)` and `copy!(::Variables)` treatment of fused views).
Committed (`634df776`). Fixes the barotropic `Duplicated(view,Array)` `update_prognostic!` error.

### ⚠️ make_zero fix is NOT sufficient for full_diff_CI — SECOND trigger + platform insight (2026-07-07)
- CI (x86) ran commit `634df776` (WITH the make_zero fix) and STILL aborts with the identical
  assertion. My update_prognostic! success used **`Const(model)` (state AD)**; full_diff_CI is
  **`Duplicated(model)` (parameter AD)** of the FULL `time_step!` → a SECOND instance of the Enzyme
  codegen bug in the parameter-AD reverse, not the leapfrog make_zero one.
- ⭐ **PLATFORM (user): CI is x86, local dev is ARM.** x86 (assertions-enabled Enzyme_jll) ABORTS
  fast; ARM (release Enzyme_jll) HANGS on the same bug (no assertion). ⇒ this must be
  debugged/verified on an **x86 Julia 1.12 node** (fast aborts → iterable); ARM cannot observe it.
- Handoff written: [`HANDOFF_full_diff_CI_x86.md`](HANDOFF_full_diff_CI_x86.md) — full diagnosis,
  x86 reproduction recipe, localization plan (bisect Const-vs-Duplicated model, dynamics_only,
  per-component parameter AD), candidate restructures, and the list of committed fixes to keep.
- maxtypeoffset does NOT help (verified). The remaining store is `size=4→storeSize=8` (looks like a
  real→ComplexF32 partial store) — a complex-valued store in the parameter-AD reverse.

---

# L29 (2026-07-06) — FiniteDifferences `to_vec(::Variables)` fix (fused view-backed vars)

The pre-existing FD harness bug (`Cannot convert Array→SubArray` in `NamedTuple_from_vec`,
24 errors in the barotropic diff test) is fixed. Root cause: the fused `Variables` store some
prognostic/grid leaves as `SubArray`-backed views into a shared parent under `vars.fused.<sym>`
(a `FusedParent` struct with `.parent` + `.slot_map`). FiniteDifferences' generic struct `to_vec`
(a) double-counted view + parent and (b) reconstructed the view field as an `Array` then
`typeof(x)(vals)` forced `convert(::SubArray, ::Array)` → error.

**Fix** (`SpeedyWeather/ext/SpeedyWeatherFiniteDifferencesExt.jl`): custom `to_vec(::Variables)`
mirroring `Base.copy!(::Variables,::Variables)` — vectorize only the non-view float/complex array
leaves (recursing `FusedParent.parent`; skip `is_view_entry` leaves since their data lives in the
fuse parent), and reconstruct **in place on a `deepcopy` template** whose view→parent aliasing
`deepcopy` preserves, so filling the parents updates the views. Validated (`_tovec_variables_check.jl`,
8/8): round-trip `to_vec(from_vec(v))≈v`; a 2× perturbation propagates through the fused views
(`vorticity.data ≈ 2×`); `j′vp` runs with no Array→SubArray error. Note the to_vec length grew
1762→3914 once `FusedParent.parent` was included — the fused prognostic parent had been silently
MISSING from the vector before (wrong FD gradients, not just an error). Full barotropic FD-side
validation running (`_driver_baro_default_offset.jl`; expect 0 `to_vec` errors).
Scratch: `_tovec_repro.jl`, `_tovec_variables_check.jl`, `_alias_diag.jl`.

## L29b (2026-07-06) — diff tests at scale after the to_vec fix

**Barotropic** (`_driver_baro_default_offset.jl`, 1.12, default maxtypeoffset): **0 `to_vec`
Array→SubArray errors** (was 24). Summary 12 pass / 6 fail / 2 err / 1 broken. The 2 errors are a
SEPARATE pre-existing issue: `update_prognostic!` (barotropic.jl:84, seeds `:prognostic`) hits
`MethodError: no method Duplicated(::LTA{SubArray}, ::LTA{Array})` — `make_zero(vars)` materializes
the fused view-backed prognostic to Array, so the shadow type mismatches the view-backed primal
when Enzyme differentiates the leapfrog `get_step` views. This is an Enzyme-`make_zero` + fused-vars
interaction (NOT the FD to_vec bug, NOT my rules). The 6 fails are FD-tolerance.

**PrimitiveWet** (`_driver_pwet.jl`, both 1.10 and 1.12): the 3 dyn-core testsets now RUN cleanly —
**0 to_vec errors, 0 GC-corruption/SIGABRT, 0 EnzymeNoTypeErrors** on BOTH versions (vs L23 which
had all three). `dynamics_tendencies!` fails only the `@test all(isfinite)` (primitivewet.jl:19) —
the pre-existing 5-day-spinup NaN instability (issue #1), not an AD/FD bug; the FD comparison is
`@test_broken` anyway. (pwet FD `j′vp` is slow — ~50k primal evals/testset — so full run is long,
but no errors/crashes.)

**Net:** the FD `to_vec` fix unblocked the FD side of the diff tests on both Julia versions.
Remaining diff-test failures are all pre-existing and separate: (a) 5-day spinup NaN → non-finite
gradients; (b) FD tolerance (already `@test_broken`/known-marginal); (c) `make_zero`/`Duplicated`
view/Array mismatch for the fused prognostic in `update_prognostic!`. None are Enzyme type-analysis
errors or transform-reverse crashes anymore.

---

# MWE reduction state — EnzymeNoTypeError (Julia ≥ 1.11)

Goal: SpeedyWeather-free MWE of the bug that `Enzyme.API.looseTypeAnalysis!(true)` works around
(`update_prognostic!` / leapfrog kernel, Julia 1.12, Enzyme 0.13.173).

Method: downward reduction from the failing real call (all scripts force `looseTypeAnalysis!(false)`
to undo the ext; run with `julia +1.12 --project=SpeedyWeather --check-bounds=yes <file>`).

NOTE: the intermediate `_mwe_reduce*`/`_mwe_d*` scratch scripts referenced in the tables below
were deleted after the campaign finished (2026-07-02); the kept artifact is the final
[`MWE_enzyme_runtime_ntuple.jl`](MWE_enzyme_runtime_ntuple.jl).

## Results so far

| Step | File | Content | Result |
|------|------|---------|--------|
| R0 | `_mwe_reduce0.jl` | real `update_prognostic!`, flag off (sanity) | ❌ FAILED (expected) |
| R1 | `_mwe_reduce1.jl` | verbatim **local copy** of `update_prognostic!` + kernel, real SpeedyWeather types | ❌ FAILED → bug follows code **structure**, not the precompiled MethodInstance |
| R2 | `_mwe_reduce2.jl` | single delta on old non-reproducing MWE: `var_old, var_new = get_steps(var)` (**runtime-length ntuple**) | ❌ FAILED → **load-bearing ingredient found** |
| R3 | `_mwe_reduce3.jl` | single delta: the two `@boundscheck` lines | ✅ success (not needed) |
| R4 | `_mwe_reduce4.jl` | single delta: Clock/Leapfrog `ifelse` scalar plumbing (Δt, lf, w1, w2) | ✅ success (not needed) |
| R2b | `_mwe_reduce2b.jl` | R2 with **literal** `var_lf` (no runtime index anywhere) | ❌ FAILED → runtime-offset view **not needed**; ntuple alone triggers |
| R5 | `_mwe_reduce5.jl` | **SpeedyWeather-free** kitchen sink (Enzyme+KA+Dates only, mimics all ingredients) | ❌ FAILED on 1.12, ✅ SUCCESS on 1.10 (correct gradient) → **SpeedyWeather-free reproducer achieved** |
| R6 | `_mwe_reduce6.jl` | strip R5: 5 variants | V1/V2 (with constructor check) ✅, **V3 (no check) ❌**, V4 (raw views) ✅, V5 (no spectrum field) ✅ |
| R7a | `_mwe_reduce7a.jl` | V3 standalone (fresh process) | ❌ FAILED → confirmed, no ordering artifact |
| R7b | `_mwe_reduce7b.jl` | R7a with **literal** step offsets everywhere | ❌ FAILED → **no runtime index needed at all** |
| R7c | `_mwe_reduce7c.jl` | data-only wrapper (no spectrum field) | ✅ success → 2nd field required |
| R8a | `_mwe_reduce8a.jl` | 2nd field = bare `Vector{Int}` | ✅ success → nested struct required |
| R8b | `_mwe_reduce8b.jl` | CONTROL: literal-length `ntuple(f, 2)` | ✅ success → runtime length is THE trigger |
| R8c | `_mwe_reduce8c.jl` | plain Julia loop instead of KA kernel | ✅ success → KA kernel required |

## ✅ FINAL MWE: [`MWE_enzyme_runtime_ntuple.jl`](MWE_enzyme_runtime_ntuple.jl)

~90 lines, **only Enzyme + KernelAbstractions**. Verified: Julia 1.12.6 → `EnzymeNoTypeError`;
Julia 1.10.11 → SUCCESS with the exactly-correct gradient (`Δt(1+w1−w2) = 0.000995`).

Three load-bearing ingredients (each ablation-verified necessary):
1. **Runtime-length `ntuple`** of wrapped views: `ntuple(s -> view3(var, s), size(var, 3))` —
   the runtime-n branch ladder returns a small **Union of tuple types** of large view aggregates.
   Literal length → fine.
2. The views must go into a differentiated **KernelAbstractions CPU kernel** (plain loop → fine).
3. The wrapper must carry a **second, inactive nested-struct field** (Ints + `Vector{Int}`s; a
   data-only wrapper or a bare `Vector{Int}` field → fine).

NOT needed (all ablated): runtime step index / runtime-offset view (!), boundschecks, mutable
Const structs (Clock/Leapfrog), big Const model args, the LTA inner-constructor check
(in the *minimal* context the check even rescues it, R6-V1 vs V3 — presumably by forcing the
tuple through a typed heap allocation; the real SpeedyWeather code has the check and still fails
in its bigger context).

## Revised root cause

`get_steps(var) = ntuple(step -> get_step(var, step), size(var, 3))` in
`time_stepping/steps.jl` — not the runtime-offset `get_step` view previously blamed. The
small-union tuple of big view aggregates, passed into a differentiated KA kernel, exceeds what
Enzyme's type analysis can prove on Julia ≥ 1.11 (Memory-backed arrays make the view aggregates
larger/more nested; 1.10 layouts stay under the limit). Consistent with the error's "very large
sized registers" wording, with why `maxtypeoffset!` didn't help, and with why the old
`_literal_step` branch-ladder fix worked (it made the ntuple-lambda's views literal-offset,
collapsing the union).

## Source-fix campaign (follow-up, 2026-07-02)

Upstream issue filed: https://github.com/EnzymeAD/Enzyme.jl/issues/3275

Attempt 1 — compile-time-length tuple `get_steps(var, Val(2))` in the differentiated consumers:
**NOT sufficient.** The real function still failed (R1b, and R2c isolated it): with the *real*
`LowerTriangularArray` views (fat aggregates, full `Spectrum` inside) even a **compile-time**
`Tuple{V, V}` breaks Enzyme's type analysis — the runtime-length union was only the extra push
needed for the slim MWE wrapper. It is fundamentally an aggregate-size limit.

| Step | Content | Result |
|------|---------|--------|
| R1b | local copy + `get_steps(var, Val(2))` | ❌ still fails |
| D1/D3/D4/D5 | R1b minus boundscheck / runtime `var_lf` / clock scalars / tendency dispatch | ❌ all still fail |
| **R2c** | minimal failing base + `Val(2)` tuple | ❌ **compile-time tuple of real LTA views fails too** |
| **R2d** | minimal failing base, **no tuple**: `var_old = get_step(var, 1); var_new = get_step(var, 2)` | ✅ SUCCESS |
| R0 re-run | real `update_prognostic!` with the no-tuple fix, flag OFF, 1.12 | ✅ **SUCCESS, exact gradient 0.0022670224** |

Attempt 2 — **no tuple at all** (APPLIED): bind steps individually via `get_step(var, 1/2)` in
the Enzyme-differentiated consumers: `update_prognostic!` (steppers/leapfrog.jl),
`implicit_primitive_equations.jl`, `implicit_shallow_water.jl`. The generic `get_steps` (runtime
and `Val` methods) remain for host-side copies (initialize!, move_prognostic_grid_variables_back!).

## Flag-off verification (final, 2026-07-02)

Full driver (speedy_transforms + barotropic) with `looseTypeAnalysis!` DISABLED + no-tuple fix:

| | Julia 1.10 | Julia 1.12 |
|--|-----------|-----------|
| barotropic | 17 pass, 3 FD-fail, **0 error** (unregressed) | 8 pass, **3 errored** (was 4 pre-fix), 1 broken |
| update_prognostic! testset | ✅ | ✅ **fixed at source** (was error #1) |
| time_step!/parameters/NCycleLorenz `Duplicated(model)` | ✅ | ❌ error #2 persists — **not tuple-related** |

**Note on the final-config confirmation run (source fix + flag ON, 1.12):** 16 pass / 4 FD-fail /
**0 Enzyme errors**. The 4th FD-fail (`transform!(::Variables)`, barotropic.jl:163) is *marginal
FD flakiness*, not gradient corruption: across seeds the Enzyme gradient is deterministic
(e.g. 14.295202 twice) while the Float32 `central_fdm(11,1)` side wiggles (14.2938–14.2969),
right at the `atol=rtol=1e-3` boundary (verified with `_chk_transform.jl`, 3 seeds: pass/fail
flips seed-to-seed). Consider loosening that tolerance slightly.

**Conclusion / shipped configuration:**
1. No-tuple source fix (state AD differentiates natively on 1.12): `update_prognostic!`,
   `implicit_primitive_equations.jl`, `implicit_shallow_water.jl` bind step views individually.
2. `Enzyme.API.looseTypeAnalysis!(true)` (Julia ≥ 1.11) RESTORED in the ext — still required for
   parameter-sensitivity AD (`Duplicated(model)`, error #2 in `vorticity_flux_curldiv!`), which is
   the same class of type-analysis failure but not caused by `get_steps` tuples.
3. Unit test added: `test/dynamics/steps.jl` (get_step/get_steps runtime + Val methods, 20/20).
4. Upstream: https://github.com/EnzymeAD/Enzyme.jl/issues/3275 with `MWE_enzyme_runtime_ntuple.jl`.

---

# Error #2 reduction campaign (started 2026-07-02)

Target: `EnzymeNoTypeError` when differentiating w.r.t. the model (`Duplicated(model)`), Julia ≥ 1.11.
Probes force `looseTypeAnalysis!(false)`; scripts prefixed `_mwe2_*`.

| Step | File | Content | Result |
|------|------|---------|--------|
| E0 | `_mwe_err2.jl` | `dynamics_tendencies!`: `Const(model)` vs `Duplicated(model)` | Const ✅ / Duplicated ❌ (confirmed) |
| E1 | `_mwe2_e1.jl` | **direct** `vorticity_flux_curldiv!(div=false, add=true)` with `Duplicated(model)` | ❌ FAILED → smaller entry point |
| — | verbose log | `Cannot deduce type of copy` — 24-byte memcpy (Julia ≥1.11 `Vector` aggregate) at `getproperty`, **tendencies.jl:844 = `S = model.spectral_transform`** | root frame |
| E2c | `_mwe2_e2.jl` | first half only (kernel launch with model's `f`, `coslat⁻¹`, `whichring`) | ✅ not the trigger |
| E2b | `_mwe2_e2.jl` | second half only (`S = model.spectral_transform`; 2× `transform!`; `curl!`) | ❌ FAILED |
| E2a | `_mwe2_e2.jl` | verbatim local copy | ❌ FAILED → structure-borne |
| E3-d3 | `_mwe2_e3.jl` | bare `transform!` with `Duplicated(S)` + `Duplicated(scratch)` (S passed directly) | ✅ Duplicated S itself is FINE |
| E3-d1/d2 | `_mwe2_e3.jl` | same with `Const(scratch)` | ❌ `MethodError: no augmented_primal` — **rule-coverage gap**: SpeedyTransforms' `_fourier!` EnzymeRule has no method for Const scratch memory (side finding) |
| E4 | `_mwe2_e4.jl` | is the `S = model.spectral_transform` getproperty the trigger, and does a 1-field wrapper suffice? | ⏳ running |

Working hypothesis: the *getproperty load of the SpectralTransform out of the Duplicated model*
is the unprovable 24-byte Vector-field copy — same aggregate-size type-analysis class as error #1,
but on the model side. If a small wrapper suffices (E4b), a SpeedyWeather-free MWE mimicking a
many-Vector-field struct is within reach; if not, aggregate size of the surrounding model matters.

## Error #2, rounds E5–E9 — ⭐ trigger found: MUTABLE Duplicated parent × transform!

| Step | Content | Result |
|------|---------|--------|
| E5 | wrappers: S + 10/40 junk Vectors; S + geometry/output/feedback | ✅ all pass |
| E7 | wrappers: S + model-field halves A and B; plain KA-kernel/broadcast consuming `S.gradients.eigenvalues` under FULL model | ✅ all pass → transform! machinery required; no field subset triggers |
| E8 | wrappers with A+B (15 fields) and with **ALL 18 model fields** + S | ✅ both pass (!) — identical content to the real model, different type |
| **E9** | **`mutable struct MutModel; spectral_transform; end`** + transform! | ❌ **FAILED** — mutability is the trigger |
| E9 controls | mutable wrapper + plain broadcast; immutable wrapper + transform! | ✅ both pass |

**Error #2 root cause (shape):** loading a `SpectralTransform` out of a **`mutable struct`**
`Duplicated` argument and passing it through `transform!` breaks Enzyme's type analysis on
Julia ≥ 1.11 ("cannot deduce type of copy", a 24-byte `Vector` aggregate inside S). The same load
from an *immutable* wrapper — even one carrying all 18 fields of the real (mutable)
`BarotropicModel` — is fine, as is any plain (kernel/broadcast) use of S from the mutable parent.
`BarotropicModel` is `@kwdef mutable struct` — hence `Duplicated(model)` parameter AD always hits
this; `Const(model)` never does (no shadow loads).

Reproducer is now ~30 lines of SpeedyWeather-based code (`_mwe2_e9.jl`, case m1).
E10 (running): which transform! half needs it — `_fourier!` (custom EnzymeRules/FFTW plans) or
`_legendre!` (KA kernels) — decides the SpeedyWeather-free synthesis path.

Side finding (E3-d1/d2): SpeedyTransforms' `_fourier!` EnzymeRule has **no `augmented_primal`
method for `Const` scratch memory** → raw `MethodError` instead of a graceful path; transforms
currently require `Duplicated` scratch under AD.

## Error #2, rounds E10–E11 + current endgame

| Step | Content | Result |
|------|---------|--------|
| E10 | mutable wrapper + only `_fourier!` (custom rule/FFTW) / only `_legendre!` (KA kernels) | ❌ **both fail independently** → custom rules NOT required |
| E11 | SpeedyWeather-free synthesis attempt (mutable wrapper + 19-field transform-like struct + legendre-like KA kernel) | ✅ did not reproduce (first shot; needs downward reduction from `_legendre!` instead) |

**Kept artifacts:** [`MWE_enzyme_mutable_transform.jl`](MWE_enzyme_mutable_transform.jl) —
~60-line SpeedyWeather-based reproducer with mutability contrast controls (m1 ❌ / m1b ✅ / i1 ✅);
[`_mwe2_e10.jl`](_mwe2_e10.jl) — the `_fourier!`/`_legendre!` split, starting point for the next
descent. All other `_mwe2_*`/probe scripts deleted (results recorded above).

**Next steps for a SpeedyWeather-free error #2 MWE** (not yet done): downward-reduce `_legendre!`
with a local verbatim copy under the mutable wrapper (the method that cracked error #1), then
swap SpeedyTransforms types for local ones. Note E7/E9-m1b: a *simple* kernel using *one* S array
passes — the many-field kernel use of the fat S is part of the trigger.

**Upstream:** worth appending the mutability finding to
https://github.com/EnzymeAD/Enzyme.jl/issues/3275 (same type-analysis class, second shape:
`getfield` of a large immutable struct from a **mutable** Duplicated parent, consumed by kernels).

## Error #2, L-series (2026-07-02, last result before pausing)

`_legendre!` (forward) is a **plain CPU loop — no KernelAbstractions, no custom rules!**
Yet E10 showed the real `_legendre!` fails under the mutable wrapper. L-series (`_mwe2_L1.jl`):

| Step | Content | Result |
|------|---------|--------|
| L1 | **verbatim local copy** of `_legendre!` + inner accumulate, S from mutable wrapper | ✅ SUCCESS (!) |
| L2 | L1 without the two `@boundscheck` lines | ✅ SUCCESS |
| L3 | L2 with trivial inner loop | ✅ SUCCESS |

⚠️ **Unlike error #1** (where the local copy R1 reproduced), the local `_legendre!` copy does NOT
reproduce. The failure is tied to the *real* `_legendre!` — candidate deltas between L1 and the
failing E10 call `_legendre!(spec, f_north, f_south, scratch.column, S)`:
  (a) the real signature's **type annotations** (`specs::LowerTriangularArray`,
      `f_north::AbstractArray{<:Complex, 3}`, `S::SpectralTransform`) — L1's copy is unannotated;
  (b) in E10 the column scratch was extracted `scratch.column` from a Duplicated ScratchMemory
      INSIDE the function — L1 passed the column as a separate Duplicated argument;
  (c) the real `_legendre!` MethodInstance was inferred during **precompilation** (pkgimage);
      the local copy is freshly inferred.

## Error #2, L-series round 2 (2026-07-03) — delta pinned, real types confirmed necessary

Naming note: the steps below are separate scripts (`_mwe2_L2.jl`, `_mwe2_L3.jl`), not further
ablations of `_mwe2_L1.jl` (whose *internal* steps were also labelled L1/L2/L3 above — do not
conflate the two).

`_mwe2_L2.jl` matches E10's exact harness: `S = w.spectral_transform` read from the **mutable**
model wrapper, and `scratch.north`/`.south`/`.column` all extracted from **one shared Duplicated
`ScratchMemory` parent inside the closure** (the aliasing ingredient L1 lacked); real
`SpectralTransform`/`LowerTriangularArray`/`ScratchMemory` types kept throughout, only the
`_legendre!` implementation varies:

| Step | Content | Result |
|------|---------|--------|
| L1d | control: REAL `_legendre!`, E10-exact scratch/model harness | ❌ FAILED (expected) |
| L2a | LOCAL copy, real type annotations restored | ❌ FAILED |
| L2b | LOCAL copy, no annotations (plain args) | ❌ FAILED |

→ Both local copies now fail. This **rules out (a) type annotations and (c) precompilation** as
the L1 delta. Confirmed: **(b) the shared-parent getfield extraction** — pulling multiple fields
(`.north`/`.south`/`.column`) out of one Duplicated struct inside the differentiated closure,
alongside the mutable-wrapper `getfield` from E9 — is sufficient to make even a fresh,
unannotated, freshly-inferred local function reproduce the failure.

`_mwe2_L3.jl` takes L2's confirmed pattern (shared-parent scratch getfield + mut/imm model
wrapper) and swaps every SpeedyWeather type for a local mimic of matching shape: `MyTransform`
(26-ish fields incl. `Ptr`-carrying `FakePlan`s standing in for FFTW plan objects, so type
analysis sees the same "opaque" fields), `MyScratch`/`MyColumnScratch`, `MyLTA` with a fat
`MyLTASpectrum` 2nd field (mirroring error #1's LTA trigger). Real array *data* is copied in;
only the *types* are local.

| Step | Content | Result |
|------|---------|--------|
| mut | local mimic types, mutable model wrapper | ✅ SUCCESS (did not reproduce) |
| imm | local mimic types, immutable model wrapper | ✅ SUCCESS (did not reproduce) |

→ Matching the getfield pattern and aggregate *shape* is **not sufficient** once the real
SpeedyWeather types are gone — this contradicts the naive extrapolation from error #1 (where a
shape-matched, SpeedyWeather-free synthesis reproduced immediately). Real type identity of
`SpectralTransform`/`ScratchMemory`/`LowerTriangularArray` (or something in Enzyme's
type-analysis/rule dispatch specific to those types) is still load-bearing. Note also: the
mut-vs-imm contrast from E9 **disappears** here (both pass), so whatever the mimics are missing
also explains why mutability alone doesn't flip the outcome without the real types present.

**Implication for next steps:** bisect by swapping real types back in one at a time against the
L3 all-mimic baseline (real `ScratchMemory` + mimic transform/LTA; real `SpectralTransform` +
mimic scratch/LTA; real `LowerTriangularArray` output + mimic scratch/transform) to find which
single real type flips L3 from SUCCESS to FAILED — same bisection style as error #1's R7c/R8a/R8b.

## Error #2, L4 (2026-07-03) — ⭐ trigger localized: the `gradients` NamedTuple inside S

`_mwe2_L4.jl`: single-type bisection + field-level probes, all under the confirmed harness
(shared-parent scratch getfield + S from MUTABLE wrapper; parametric mimic transform so probe
fields can be swapped one at a time). Candidate suspects read off the real `SpectralTransform`
definition (spectral_transform.jl:12-84): abstract-eltype `Vector{AbstractFFTs.Plan}` fields (×4),
the Type-valued `LegendreShortcut` field, the nested `scratch_memory::ScratchMemory`, fat nested
spectrum/grid.

| Step | Content | Result |
|------|---------|--------|
| C0 | all-mimic control (L3 re-run) | ✅ (baseline reproduced) |
| C1 | all-real control (L2-equiv) | ❌ (baseline reproduced) |
| S1 | real `SpectralTransform` only (mimic scratch + spec) | ❌ **S is the trigger** |
| S2 | real `ScratchMemory` only | ✅ |
| S3 | real spec `LowerTriangularArray` only | ✅ |
| P1 | mimic S + fat spectrum as `S.spectrum` | ✅ |
| P2 | mimic S + Type-valued field | ✅ |
| P3/P3m | mimic S + 4 abstract-eltype plan vectors (immutable/MUTABLE elements) | ✅ both |
| P4 | mimic S + nested scratch inside S | ✅ |
| P5 | kitchen sink P1+P2+P3m+P4 | ✅ (!) |
| R1–R5 | real grid / real `Spectrum` / real legendre LTA / real plan vectors / real nested `ScratchMemory` swapped into mimic S | ✅ all five |
| **R6** | **real `gradients` NamedTuple** swapped into mimic S | ❌ **FAILED — the trigger** |

**Finding:** the real `gradients` NamedTuple (`gradient_arrays`, gradient_arrays.jl:162 — 8
`LowerTriangularArray`s + 2 `Vector`s, each LTA carrying data + a full nested `Spectrum` with 4
more Vector fields) is the single field whose presence in S breaks type analysis. Crucially, the
legendre loop **never reads `S.gradients`** — its mere presence in the aggregate loaded from the
mutable Duplicated wrapper is enough. All the "exotic" suspects (FFTW-plan vectors with abstract
eltype, Type field, nested scratch) are innocent. This is the aggregate-size story again, error
#1's conclusion, now on the model side: `gradients` is by far the largest field of S
(~8×(24-byte data Vector + fat Spectrum) unrolled).

Next: `_mwe2_L5.jl` — (a) real-LTA subset counting (1/2/4/8 gradient fields) to find the size
threshold; (b) all-mimic gradients NamedTuple (8 MyLTA + 2 vecs) — if that fails, the
SpeedyWeather-free MWE is essentially done; conditional ablations either way.

## Error #2, L5 (2026-07-03) — ⭐⭐ SpeedyWeather-free reproducer; pure aggregate-size threshold

`_mwe2_L5.jl` (same harness; mimic S with `gradients` as the only knob):

| Step | Content | Result |
|------|---------|--------|
| RS1/RS2/RS4 | real gradients subset: 1 / 2 / 4 LTA fields | ✅ all |
| **RS8** | real gradients subset: 8 LTA fields (no eigen vecs) | ❌ **count threshold between 4 and 8** |
| G1 | mimic: 8 MyLTA (medium-fat MyLTASpectrum: 2 Int + 2 Vec) + 2 vecs | ✅ just under the limit |
| **G6** | mimic: 8 MyLTA with **field-exact spectrum mimic** (7 fields: 2 Int, arch, 2×Vec{UnitRange}, 2×Vec{Int}) | ❌ **SpeedyWeather-free reproducer** |
| **G7** | mimic: 16 MyLTA (medium-fat) + 2 vecs | ❌ second knob: count compensates fatness |
| G8 | mimic: 16 MyLTA field-exact | ❌ |
| G9 | fresh `gradient_arrays(Float32, spectrum)` (identity check) | ❌ → no aliasing/identity effect |

**Conclusion:** the failure is a cumulative **aggregate-size threshold** on the struct loaded
from the mutable Duplicated parent — LTA count × per-LTA nested-spectrum fatness trade off
against each other (8×exact ❌ vs 8×medium ✅ vs 16×medium ❌). No exotic ingredient (plans, Type
fields, nested scratch, aliasing) required. Error #1 and error #2 are now demonstrably the same
root cause in two guises: "aggregate too large for Enzyme's type analysis on Julia ≥ 1.11",
via tuple-of-views (error #1) or via getfield-of-fat-struct from a mutable parent (error #2).

Next: `_mwe2_L6.jl` — fully synthetic MWE (no SpeedyWeather import; sizes synthesized, T9-like:
lmax=mmax=10, ncoeffs=55, nlat_half=8; baseline = 16 exact-fat LTAs) + ablations the size story
may have obsoleted: no-scratch-struct variant, slim 7-field transform, immutable-wrapper
contrast, gradients-count bisect on the minimal config.

## Error #2, L6 (2026-07-03) — ⭐⭐⭐ FULLY SYNTHETIC MWE (Enzyme + Base only); scaffolding falls away

`_mwe2_L6.jl`: no SpeedyWeather import; all sizes/data synthesized (T9-like: lmax=mmax=10,
ncoeffs=55, nlat_half=8, nlayers=1).

| Step | Content | Result |
|------|---------|--------|
| F0 | synth baseline: full 21-field S, 16 exact-fat LTA, mut wrapper, shared-scratch parent | ❌ reproduces |
| F4 | NO scratch struct: separate `f_north`/`f_south` Duplicated args, `even`/`odd` allocated inside | ❌ scratch struct NOT needed |
| F5 | SLIM 7-field transform (only the 6 loop-read fields + unread `gradients`) | ❌ padding NOT needed |
| F45 | slim + no scratch combined | ❌ minimal shape holds |
| **F6** | minimal config + **IMMUTABLE wrapper** | ✅ **mut/imm contrast intact** |
| F8a | minimal, 8 exact-fat LTA | ❌ (8 suffices in minimal config too) |
| F9 | minimal, 16 LTA but SLIM (2-Int) nested spectrum | ✅ nested-vector fatness essential |
| F10 | minimal, no eigen vectors | ❌ eigen vecs NOT needed |

**Minimal failing shape after L6:** `mutable` wrapper → 7-field transform struct → plain CPU loop
(separate Duplicated args, locals inside); the only fat ingredient is the **unread** `gradients`
NamedTuple of 8 `MyLTA(Vector{Float32}, field-exact 7-field spectrum mimic)`. Immutable wrapper
with the identical content passes. Earlier "necessary ingredients" now ablated away: shared-parent
scratch getfield (L2's conclusion — superseded: it was just another size contribution),
KernelAbstractions (never needed for error #2), custom rules, eigen vectors, padding fields.

Next: `_mwe2_L7.jl` — LTA count bisect below 8, trivial 5-line loop, 2-/1-field transform
structs, plain-Matrix output, bare fat structs without the MyLTA wrapper, and (user question)
gradients as a struct / plain tuple with identical fields instead of a NamedTuple.

## Error #2, L7 (2026-07-03) — container kind irrelevant; loop complexity is the 3rd ingredient

`_mwe2_L7.jl` (base: L6's minimal shape — slim S, structured loop, 8 exact-fat LTAs, mut):

| Step | Content | Result |
|------|---------|--------|
| B0 | base re-check | ❌ (anchor holds) |
| **H8** | **gradients as a `struct` with the same 8 LTA fields** (user question) | ❌ **does NOT rescue** |
| **H8t** | gradients as a plain `Tuple` of the 8 LTAs | ❌ does NOT rescue |
| H1a | 4 LTAs | ✅ threshold 4<n≤8, matches L5's real-LTA counting |
| H2 | trivial elementwise loop instead of the structured one, same 8-LTA gradients | ✅ **loop complexity is load-bearing** |

**Answer to the struct-instead-of-NamedTuple question: NO.** NamedTuple, plain tuple, and struct
with identical fields all fail identically — Enzyme lowers them to the same aggregate; only the
content/size matters. Swapping the `gradients` NamedTuple for a struct in `SpectralTransform` is
therefore NOT a viable source fix.

The three final ingredients of error #2 (each ablation-verified):
1. getfield of the transform from a **mutable** Duplicated wrapper inside the differentiated code
   (immutable wrapper with identical content → fine);
2. a large nested **unread** field in the transform (≥8 array-wrapper structs with fat nested
   spectra; pure size threshold, container kind irrelevant);
3. a sufficiently complex differentiated loop (complex even/odd intermediates + strided view
   accumulation; a trivial elementwise loop → fine).

## ✅✅ FINAL error #2 MWE: [`MWE_enzyme_fat_aggregate.jl`](MWE_enzyme_fat_aggregate.jl)

> **NOTE (2026-07-03, later):** this section describes the FIRST frozen version (~150 lines,
> structured legendre loop). L8–L13 below reduced it much further; the file was REWRITTEN in
> place (~85 lines, plain structs, 7-line loop, exact-gradient check = 80.0). See the L13
> section for the final verification matrix. Historical record kept as-is below.

~150 lines, **Enzyme + Base only**, two-case contrast in one file. Fully verified (Enzyme 0.13.173):

| | mutable wrapper | immutable wrapper |
|--|--|--|
| Julia 1.12.6 | ❌ `EnzymeNoTypeError` | ✅ (sum \|df_north\| = 88.76) |
| Julia 1.10.11 | ✅ (133.91) | ✅ (96.75) |

Ready to post to https://github.com/EnzymeAD/Enzyme.jl/issues/3275 (or a follow-up issue) as the
dependency-free reproducer for the mutable-parent/aggregate-size shape (endpoint A.4/A.5 of the
plan below — DONE). Intermediate probe scripts `_mwe2_L4–L7.jl` can be deleted once the upstream
post is made (results all recorded above); `_mwe2_L1–L3.jl`, `_mwe2_e10.jl` and
`MWE_enzyme_mutable_transform.jl` are superseded by the final artifact.

## Error #2, L8–L9 (2026-07-03) — loop collapses to a bare view walk; bounds-provenance suspicion

`_mwe2_L8.jl`: how much of the structured loop is needed (user question: are `legendre!`,
`accumulate!`, `loss!` really that complex)?

| Step | Content | Result |
|------|---------|--------|
| B0 | baseline loop (anchor) | ❌ |
| **M1** | **view-walk only**: single ~15-line function; NO `accumulate!`, NO even/odd locals, NO `lon_offsets`/`conj`/`solid_angles`/`mmax_truncation`, NO `f_south`, NO `fill!` | ❌ **all of that machinery removable** |
| M1g | M1 body but 2-field struct (legendre + gradients) AND const loop bounds | ✅ |

`_mwe2_L9.jl` (micro-struct size compensation + view ablations): micro + 12/16/32 LTAs all ✅;
every ablation on the 7-field config ✅. **CONFOUND (my construction error):** all L9 walk
variants read the loop bounds from **consts**, whereas L8-M1 (the failing body) reads them from
S (`S.spectrum.lmax`, `S.grid.nlat_half`). Every L9 ✅ is therefore explained by a single sharper
hypothesis rather than by the nominal ablations: **the view extents must derive from loads out of
the Duplicated struct** (analysis-opaque runtime extents) — mirroring error #1, where literal
offsets collapsed the failing union. L9's nominal conclusions (both views needed, wrappers needed,
j-loop needed, fixed extent rescues) are all INVALID as stated pending `_mwe2_L10.jl`, which
re-runs them unconfounded (X-series: X1 = const-bounds control on the 7-field struct; X3/X3b =
which S-load; X2/X2b = struct shrink to 4/3 fields with S-bounds; X4/X5 = single-view and
bare-Matrix ablations with S-bounds; X6 = immutable contrast on the minimal failing config).

## Error #2, L10 (2026-07-03) — ⭐ ingredient #3 identified precisely: Duplicated-struct-derived loop bounds

`_mwe2_L10.jl` (all on 7-field struct + 8 LTA unless noted):

| Step | Content | Result |
|------|---------|--------|
| X0 | S-derived bounds (L8-M1 verbatim) — anchor | ❌ |
| **X1** | **identical body, CONST bounds** — the missing control | ✅ **hypothesis confirmed** |
| X3 | only `lmax` from `S.spectrum` (nlat const) | ❌ either load suffices |
| X3b | only `nlat_half` from `S.grid` (lmax const) | ❌ |
| X2 | S-bounds, 4-field struct (unread vector fields removed) | ✅ aggregate size still in play |
| X4s / X4l | S-bounds + only the specs view / only the legendre view | ❌ / ❌ **a single view suffices** |
| X5l | S-bounds + legendre as bare Matrix (Wrapped removed) | ✅ leg wrapper needed (or its fat spectrum's bytes) |
| X5s | S-bounds + output as bare Matrix | ❌ **output wrapper NOT needed** |
| X6 | immutable contrast on the failing config | ✅ contrast intact |

**Ingredient #3, final formulation:** the differentiated loop's trip counts / view extents must
derive from values **loaded out of the Duplicated (mutable-parent) struct chain** at runtime
(`S.spectrum.lmax` or `S.grid.nlat_half` — either alone suffices). With compile-time-constant
bounds the identical body passes. This is the precise analogue of error #1's runtime-length
ntuple (there, literal offsets collapsed the union and fixed it). L9's "both views needed /
j-loop needed / fixed extent rescues" claims were pure confound; single view + plain-Matrix
output fail fine.

## Error #2, L11 (2026-07-03) — all reductions hold independently

`_mwe2_L11.jl`:

| Step | Content | Result |
|------|---------|--------|
| Y0 | combo: plain-Matrix output + only legendre view + S-derived j bound (7-field, 8 LTA) | ❌ |
| **Y5** | Y0 minus the m/lm walk: opaque j loop over ONE **fixed-extent** view per ring | ❌ triangular walk NOT needed |
| **Y1** | no `f` input at all (2 Duplicated args: output + wrapper) | ❌ f NOT needed |
| Y2 / Y2b | 3-field struct (spectrum, legendre, gradients) + 16 / 24 LTA, lmax bound | ❌ / ❌ struct shrinks to 3 fields |
| Y6 | immutable contrast on Y0 config | ✅ |

Note Y5 nuance: with the bound itself opaque (S-load), even a **fixed-extent** view inside the
loop triggers — the runtime-*extent* variation of the views is not required, only the
opaque trip count around them. Next: `_mwe2_L12.jl` merges Y0+Y1+Y5+Y2 (Z0), probes a 2-field
struct with the bound read through the wrapper chain (`leg.spectrum.lmax`, Z8), the triviality
floor (no inner loop, Z2), plain-real output (Z5), plus immutable and `Const(w)` controls.

## Error #2, L12 (2026-07-03) — the floor

`_mwe2_L12.jl`:

| Step | Content | Result |
|------|---------|--------|
| Z0 | merged combo: 3-field struct, no f, fixed-extent view, opaque lmax bound, 16 LTA | ❌ |
| **Z8** | **2-field struct** (legendre + gradients), bound read through the wrapper chain (`leg.spectrum.lmax`) | ❌ smallest struct |
| Z1 | Z0 with 8 LTA | ✅ threshold rose to (8,16] at the leaner code — size story consistent |
| Z2 | no inner loop (single-element view read) | ✅ floor: the elementwise loop over the view is required |
| **Z5** | output = plain `Matrix{Float32}` | ❌ no complex needed |
| Z6 | IMMUTABLE contrast (smallest failing) | ✅ |
| Z7 | `Const(w)` control (smallest failing) | ✅ matches the real code: `Const(model)` always worked |

Minimal failing shape after L12: `mutable` wrapper → 2-field struct
(`legendre_polynomials::Wrapped`, `gradients::NamedTuple` of 16 `Wrapped`) → ~7-line loop whose
trip count loads `leg.spectrum.lmax` from the Duplicated chain, one fixed-extent `view` of
`leg.data` per iteration, elementwise `muladd` into a plain `Float32` output. `_mwe2_L13.jl`
checks the Z8+Z5 composition, whether `Wrapped` needs its `AbstractArray` interface at all, and
whether data lengths can shrink (55→8), then the final file gets frozen + verified on 1.10/1.12.

## Error #2, L13 (2026-07-03) — composition holds; MWE frozen at ~55 lines of trivial code

`_mwe2_L13.jl`: W0 anchor ❌; W1 (2-field + Float32 out) ❌; **W2 (`Wrapped` as a plain 2-field
struct, NO `AbstractArray` subtyping/interface) ❌**; W5 (all data lengths 55→8) ❌; W6 (bare +
short combined) ❌; W4 immutable contrast on exactly that shape ✅. Nothing about the array
wrapper beyond its 2-field layout matters; data lengths are irrelevant.

**FINAL MWE (rewritten in place):** [`MWE_enzyme_fat_aggregate.jl`](MWE_enzyme_fat_aggregate.jl)
— Enzyme + Base only, ~85 lines total (~55 without comments), plain structs, one ~7-line loop,
deterministic gradient check (sum |d legendre| = 80 exactly on success). Supersedes the earlier
150-line version (structured legendre loop, scratch machinery, complex numbers — all shown
unnecessary by L8–L13). Answer to "are `legendre!`/`accumulate!`/`loss!` needed in that
complexity": NO — they collapse to the 7-line `walk!`; the load-bearing loop property is solely
the Duplicated-chain-derived trip count around a view + elementwise accumulate (L10/X1).

Final verification matrix (Enzyme 0.13.173, exact gradient 80.0 on every SUCCESS):

| | mutable wrapper | immutable wrapper |
|--|--|--|
| Julia 1.12.6 | ❌ `EnzymeNoTypeError` | ✅ 80.0 |
| Julia 1.10.11 | ✅ 80.0 | ✅ 80.0 |

Probe scripts `_mwe2_L8.jl`–`_mwe2_L13.jl` are superseded by this file + the tables above and
can be deleted together with `_mwe2_L1–L7` once the upstream post is made.

## Error #2, L14 (2026-07-03) — field/operation necessity: almost everything falls; NOT raw bytes

`_mwe2_L14.jl`, one delta per case vs the frozen file (A0 anchor ❌ re-confirmed):

| Step | Content | Result |
|------|---------|--------|
| P1 | no `view` — direct `leg.data[l, j]` indexing | ❌ view NOT needed |
| P2 | plain `out[l,1] += lv[l]` (no muladd/constant) | ❌ |
| P3 | output as 1D `Vector{Float32}` | ❌ |
| P4 | plain `Reverse` (no `set_runtime_activity`) | ❌ |
| P5 | no `@inbounds` | ❌ |
| P6 | 5-field spectrum (drop `Arch` + `mmax`) | ✅ rescues at 16 (size-tradeoff again; keep 7 fields) |
| P7 | spectrum vectors all `Vector{Int}` (no `UnitRange`) | ❌ |
| **P8** | **gradients = 16 BARE spectra** (no Wrapped/data layer) | ❌ one nesting level gone |
| **P9/P9b** | gradients = 64 / **128 flat `Vector{Float32}`** | ✅/✅ **kills the raw-byte-count hypothesis** |

**Key insight (P8 vs P9b):** 128 flat Vector fields ≈ 3 KB of Vector aggregates PASS, while 16
nested struct-of-4-Vectors ≈ 1.8 KB FAIL → the trigger is the **number/depth of nested struct
fields each containing several Vector aggregates**, not total bytes. Also note P8: the gradient
spectra are **all-Int (inactive!)** — the fat unread field doesn't even need active data.
`_mwe2_L15.jl` merges the winners (K1), probes dropping the legendre wrapper too (bound from
`gradients.g1.lmax`, legendre as plain Matrix/Vector — K2/K3), + contrast/Const controls.

## Error #2, L15 (2026-07-03) — ABSOLUTE FLOOR: no custom array type, plain indexing

`_mwe2_L15.jl`:

| Step | Content | Result |
|------|---------|--------|
| K1 | all L14 winners combined (legendre still Wrapped for the bound) | ❌ |
| **K2** | + legendre = **plain `Matrix{Float32}`**, bound from `S.gradients.g1.lmax` | ❌ **no wrapper type left** |
| K3 | + legendre = plain `Vector{Float32}` (1D) | ✅ Matrix (2D indexing) required |
| K4 / K5 | immutable contrast / `Const(w)` control on K2 | ✅ / ✅ |

**FINAL MWE (rewritten in place again):**
[`MWE_enzyme_fat_aggregate.jl`](MWE_enzyme_fat_aggregate.jl) — now 2 struct types
(`Spectrum` = Ints + 4 `Vector{Int}`; `Transform` = `Matrix{Float32}` + NamedTuple of 16
`Spectrum`) + mut/imm wrappers + a **5-line plain nested loop** (`out[l] += S.legendre[l, j]`,
trip count from `S.gradients.g1.lmax`). No views, no muladd, no complex numbers, no custom
`AbstractArray`, no `@inbounds`, no `set_runtime_activity` — plain `Reverse`. The fat field is
entirely **inactive** (Int data) and unread except that one bound Int. Deterministic gradient
on success: sum |d legendre| = 40.0 exactly.

Final verification (Enzyme 0.13.173; exact gradient 40.0 on every SUCCESS):

| | mutable wrapper | immutable wrapper |
|--|--|--|
| Julia 1.12.6 (`--check-bounds=yes`) | ❌ `EnzymeNoTypeError` | ✅ 40.0 |
| Julia 1.12.6 (default flags) | ❌ (flag-independent) | ✅ 40.0 |
| Julia 1.10.11 | ✅ 40.0 | ✅ 40.0 |

Console outputs + the full error message are saved in
[`MWE_enzyme_fat_aggregate_output.txt`](MWE_enzyme_fat_aggregate_output.txt) (capture helpers:
`_mwe2_capture.jl`, `_mwe2_capture_verbose.jl`). Verbose capture note: the LLVM-level dump is
part of the exception's `showerror` (2.3 MB for this MWE, kept at
`_logs/MWE_enzyme_fat_aggregate_error_verbose_1.12.log`); its money line —
`Cannot deduce type of copy` on a **16-byte** `llvm.memcpy`, stacktrace `getproperty` →
`loss!` line 57 = the `w.transform` getfield — matches the real code's 24-byte copy at
`S = model.spectral_transform` exactly.

## Error #2, L16–L17 (2026-07-03) — the leading-Int fields are NOT ballast

`_mwe2_L16.jl`: the 5-field spectrum (lmax + 4 `Vector{Int}`, i.e. dropping `mmax` +
`architecture::Nothing`) **passes at 16, 24, 32 and even 48 spectra** (~4.6 KB of Vector
aggregates — triple the failing 7-field config's size). P6's rescue was therefore NOT size
compensation: at least one of the two "odd" fields is genuinely load-bearing.

`_mwe2_L17.jl` disambiguated (hypothesis confirmed):

| Step | Content | Result |
|------|---------|--------|
| **A16** | 6-field: `lmax` + `mmax` (no `Nothing`), 16 spectra | ❌ **two leading Ints suffice** |
| B16/B32 | 6-field: `lmax` + `Nothing` (no `mmax`), 16 / 32 spectra | ✅ / ✅ single Int never triggers |

**The `architecture::Nothing` ghost field is irrelevant; the TWO leading 8-byte Int fields are
load-bearing** — with two Ints the first Vector field sits at byte offset 16, exactly the size
of the un-deducible `llvm.memcpy`; with one Int (offset 8) the identical shape passes even at
double the count. This is a remarkably specific layout condition, valuable for the upstream
diagnosis. (The intermediate 6-field version was verified 1.12 mut ❌ / imm ✅ 40.0, 1.10 both ✅
before L18 shrank the struct further.)

## Error #2, L18 (2026-07-03) — vector-count dimension: ONE Vector per struct suffices

`_mwe2_L18.jl` (2 leading Ints + N `Vector{Int}` fields, 16/32 spectra):

| N vectors | 16 spectra | 32 spectra |
|--|--|--|
| 3 | ❌ | — |
| 2 | ✅ | ❌ |
| **1** | ✅ | ❌ |

Count × per-struct-vector-count trade off (4×16, 3×16, 2×32, 1×32 all fail; flat 128 vectors
without the 2-Int struct wrapper never fail). **True endpoint: the repeated unit is a 3-field
struct `(Int, Int, Vector{Int})`, 32 of them.** Every field of every struct in the MWE is now
provably necessary: `Transform.legendre` (differentiated data), `Transform.gradients` (the
trigger), `Spectrum.lmax`+`mmax` (the offset-16 layout condition, L17), `Spectrum.l_indices`
(the Vector aggregate unit), `MutModel.transform` (the mutable parent).

**FINAL-final MWE (32 × 3-field structs), fully verified** (Enzyme 0.13.173, exact gradient
40.0 on every SUCCESS):

| | mutable wrapper | immutable wrapper |
|--|--|--|
| Julia 1.12.6 (`--check-bounds=yes`) | ❌ `EnzymeNoTypeError` | ✅ 40.0 |
| Julia 1.12.6 (default flags) | ❌ (flag-independent) | ✅ 40.0 |
| Julia 1.10.11 | ✅ 40.0 | ✅ 40.0 |

Captures regenerated for this shape: `MWE_enzyme_fat_aggregate_output.txt` (console outputs,
full standard error, key verbose lines — the 16-byte `Cannot deduce type of copy` memcpy with
`getproperty` → `loss!` stacktrace) and `_logs/MWE_enzyme_fat_aggregate_error_verbose_1.12.log`
(836 KB full verbose dump). Probe scripts `_mwe2_L16–L18.jl` deletable with the rest after the
upstream post.

## Error #2, L19 (2026-07-03) — exact struct-count threshold bisected: 21 → 22

`_mwe2_L19.jl` (final 3-field shape, bisection between L18's 16 ✅ / 32 ❌): 24 ❌, 20 ✅, 22 ❌,
21 ✅ → **threshold is exactly 21 passes / 22 fails** (Julia 1.12.6, Enzyme 0.13.173). The MWE
was set to `NSPECTRA = 22` with the threshold denoted in the header and at the constant
(with a note to bump for margin if a future Enzyme/Julia version shifts the boundary).

22-spectra shape fully re-verified (exact gradient 40.0 on every SUCCESS):

| | mutable wrapper | immutable wrapper |
|--|--|--|
| Julia 1.12.6 (±`--check-bounds`) | ❌ `EnzymeNoTypeError` | ✅ 40.0 |
| Julia 1.10.11 | ✅ 40.0 | ✅ 40.0 |

Captures regenerated: `MWE_enzyme_fat_aggregate_output.txt` (50 lines) +
`_logs/MWE_enzyme_fat_aggregate_error_verbose_1.12.log` (792 KB; `Cannot deduce type of copy`,
16-byte memcpy, `getproperty` → `loss!`). `_mwe2_L19.jl` joins the deletable probe scripts.
**Campaign closed — the only remaining action is the upstream post (#3275 or follow-up).**

## Error #2, L20–L21 (2026-07-03) — API-knob test: maxtypedepth! useless, maxtypeoffset! FIXES it

Wrapper scripts `_mwe2_L20.jl` (MWE) / `_mwe2_L21.jl` (real `transform!` reproducer) set
`Enzyme.API.maxtypedepth!` / `maxtypeoffset!` before including the reproducer (defaults: 6 / 512):

| Setting | MWE (22 structs) | real transform! (m1) |
|--|--|--|
| depth 12 | ❌ | ❌ |
| depth 60 | ❌ | (not run) |
| depth 60 + offset 4096 | ✅ exact gradient 40.0 | ✅ |
| **offset 4096 only** | ✅ **exact gradient 40.0** | ✅ **all three cases** |

`maxtypedepth!` alone does nothing; **`maxtypeoffset!(4096)` is the effective knob** — consistent
with the failing copy being an aggregate whose relevant layout exceeds the 512-byte default
analysis window (22 × ~48-byte structs ≈ 1 KB unrolled). Unlike `looseTypeAnalysis!(true)`
(which makes Enzyme *guess* un-deducible types), `maxtypeoffset!` only enlarges the exact
analysis window — semantically safer as a SpeedyWeather-ext workaround, at some type-analysis
compile-time cost.

## Error #2, L22 (2026-07-03) — ⭐ maxtypeoffset!(4096) validated on the FULL barotropic testset

Driver `_driver_baro_offset.jl` (full `barotropic.jl`, `looseTypeAnalysis!(false)` +
`maxtypeoffset!(4096)`, Julia 1.12.6, `--check-bounds=yes`, 47 min):

| | Pass | Fail | Error | Broken |
|--|--|--|--|--|
| loose OFF, offset default (prior baseline) | 8 | — | **3** | 1 |
| **loose OFF, offset 4096 (this run)** | 16 | 4 | **0** | 1 |
| loose ON, offset default (prior baseline) | 16 | 4 | 0 | — |

**`maxtypeoffset!(4096)` eliminates all 3 error-#2 failures at full model scale** (the knob
generalizes from the isolated `transform!` to `Duplicated(model)`), reproducing the loose-ON
headline (16/4, 0 errors) with **exact** type analysis. The formerly-erroring
`Barotropic model parameters` testset (parameter AD, `Duplicated(model)`) now **compiles and
runs**, failing only on FD tolerance.

The 4 remaining failures are pre-existing FD-side discrepancies, NOT offset-induced wrong
gradients — verified by cross-checking Julia 1.10 (native AD, no workaround, trusted):
- `time_step!` barotropic.jl:219/220 — Enzyme gives ~4.0 where FD gives ~2.0 (factor-2). The
  **identical** pattern appears in the 1.10 log (Enzyme 4.0 / FD 2.0) — Enzyme is self-consistent
  across 1.10-native and 1.12-offset, so the discrepancy is on the FD/test side.
- `transform!(::Variables)` barotropic.jl:163 — the known marginal FD flake (atol=rtol=1e-3).
- `Barotropic model parameters` barotropic.jl:251 — also an FD-fail on trusted 1.10.

**Recommended ext change:** replace the gated `looseTypeAnalysis!(true)` with
`Enzyme.API.maxtypeoffset!(4096)` for Julia ≥ 1.11 — same test outcome, exact (non-guessing)
analysis, no correctness risk. Caveat: 4096 is tuned to the barotropic model; PrimitiveWet's
larger aggregates may need a bigger value — validate there before finalizing, and consider
loosening the barotropic.jl:163 tolerance to deflake. One open item (nice-to-have): a direct
1.10-vs-1.12-offset parameter-gradient value comparison to certify the parameter-AD gradient
numerically (state-AD is already self-consistent via the factor-2 cross-check).

## Error #2, L23 (2026-07-04) — PrimitiveWet dyn-core does NOT establish a clean 1.10 baseline

Ran `_driver_pwet.jl` (3 dyn-core testsets, all state-AD `Duplicated(vars)/Const(model)`) on
Julia 1.10.11 FIRST, to get the trusted native baseline before the 1.12+offset comparison.
Result: **not clean — three distinct pre-existing problems, none related to maxtypeoffset**
(log: `_logs/v1.10_primitivewet_dyncore_offset-era.log`):

1. **Non-finite AD gradients** — `dynamics_tendencies!` (primitivewet.jl:19) and
   `implicit_correction!` (:48) both FAIL `@test all(isfinite.(dvec))`. The Enzyme reverse pass
   returns NaN/Inf on 1.10 (native, no workaround). Root cause TBD: could be the
   `initialize_with_spinup!` state (5-day wet spinup → possible NaN/instability) or a genuine AD
   issue; needs a finiteness check of the spinup state itself.
2. **FiniteDifferences harness bug (NOT Enzyme)** — both testsets then throw
   `MethodError: cannot convert Array{ComplexF32,3} → SubArray{...}` inside
   `FiniteDifferences.NamedTuple_from_vec` (`to_vec.jl:267`) reconstructing the `Variables`.
   `to_vec` recorded `vorticity`/etc. as **SubArray-backed** `LowerTriangularArray`s (the
   `get_step` views) but rebuilds them from plain `Array`s; `LowerTriangularArrays`
   convert (lower_triangular_array.jl:14/636) has no Array→SubArray method. This is a
   test-infrastructure incompatibility with the view-backed `Variables`, independent of AD.
3. **Process SIGABRT (exit 134)** — after the first two testsets, Enzyme/LLVM crashes compiling
   the reverse pass of the 3rd testset (`transform!(::Variables)`) on 1.10 (LLVM.Attribute
   "ptr queue" dump, ~1.5 M-line bug-report dump; discarded, head kept in the archived log).

**Implication:** the maxtypeoffset validation on PrimitiveWet is BLOCKED until the pwet dyn-core
tests run cleanly on 1.10 (harness fix for #2 at minimum; investigate #1 and #3). The barotropic
offset result (L22) stands on its own. Decision needed from user on whether to fix the pwet
test harness / investigate the 1.10 crash before proceeding, or defer pwet and ship the
barotropic-validated offset change.

## Error #3, L24 (2026-07-04) — ⭐ crash root-caused: degenerate grid-field shadow in the _fourier! reverse rule

Isolated the crash with `_mwe3_transform_crash.jl` (just setup + the single
`autodiff(transform!(::Variables))`, no @testset/FD, flushed markers; spinup length as ARG).
Log: `_logs/pwet_transform_reverse_crash.log`. Findings, in order:

1. **Confirmed crasher = `transform!(::Variables)` reverse pass** (dies right after the
   "entering autodiff" marker; never returns).
2. **NOT a type-analysis failure** (unrelated to errors #1/#2, offset/loose knobs irrelevant):
   1.10 → `GC error (probable corruption)` + SIGABRT (exit 134); 1.12 → clean **`BoundsError`**.
   Same root cause, different manifestation (1.10 bounds-check bypassed in Enzyme codegen → OOB
   write corrupts the GC heap; 1.12 `--check-bounds=yes` catches it).
3. **NOT caused by the NaN state** (issue #1): with a 2-day spinup the primal state is fully
   finite (`all-finite=true`) yet it still crashes. (Issue #1 — non-finite gradients — is
   separately root-caused: the model spinup itself emits `NaN or Inf detected at time step 37`
   ≈ day 3; the 5-day-spinup config is numerically unstable. Independent of AD.)
4. **100% AD-specific** — a plain non-AD `transform!(deepcopy(vars), model)` on the same state
   runs fine (`PRIMAL transform!(::Variables) OK`). Only the reverse pass fails.
5. **Localized to the custom `_fourier!` EnzymeRule** (`SpeedyTransformsEnzymeExt.jl:118`, the
   `reverse` for `_fourier!(grid, f_north, f_south, S)`): it calls
   `_fourier!(dfnorthval, dfsouthval, grids.dval, S.val)` (forward FFT for the adjoint) →
   `_fourier_serial!` → `_apply_serial_fft!` (fourier.jl:91) → OOB
   `BoundsError: 0×1 view(::Array{Float32,3}, :, 0:0, 0) at index [1:20,1]`.
6. **Confirmed mechanism (instrumented the rule, shapes printed, reverted):**

   | array | size |
   |--|--|
   | `dfnorthval` / `f_north.val` | `(25, 4, 8)` ✓ (nfreq_max × nlayers × nlat_half) |
   | **`grids.dval`** (grid-field shadow) | **`(0, 1)`** ✗ degenerate/empty (should be N_points × 4) |
   | `size(grids.dval, 2)` → K=1 | forces the **serial** FFT path (`K>1` false) |

   **Root cause:** Enzyme's shadow of the intermediate grid field, `grids.dval`, is allocated
   **empty `(0,1)`** instead of matching the primal `grids.val` (N_points × 4). Two consequences:
   (a) `K = size(grids.dval, 2) = 1` routes to `_fourier_serial!` (barotropic nlayers=1 always
   uses serial and WORKS — so the serial path itself is fine; it's the size mismatch that kills
   it); (b) the `0`-sized `.data` makes `_apply_serial_fft!` index out of bounds.

**This is fixable in SpeedyWeather (no Enzyme upstream fix needed)** — unlike errors #1/#2.
The bug is that `grids.dval` reaches the `_fourier!` reverse rule with degenerate dims. Open
question for the fix: WHY is that shadow `(0,1)` — is the intermediate grid field constructed
inside `transform!` in a way Enzyme can't shadow (e.g. `similar`/reshaped/zero-sized alloc), or
is the rule being handed the wrong argument? Next step (fix development): trace how the grid
field reaching `_fourier!` is created in the `transform!(::Variables)` chain
(transform.jl:162 → spectral_transform.jl:464/407/475) and how its `Duplicated` shadow is set up.
Note: this is orthogonal to the maxtypeoffset ext change (which targets the type-analysis
errors); barotropic AD is unaffected (nlayers=1, and it passes).

### L24 fix diagnosis (2026-07-05) — narrowing the degenerate-shadow trigger

Trigger chain (from `transform!(::Variables)` PrimitiveEquation, transform.jl:149-151):
`transform!(grid_variables, prognostic_vars, scratch, S)` where
`grid_variables = get_prognostic_step(parent(vars.fused.grid), …)` — a **fused** field
(vor+div+temp+pres+humid packed, K≈17, unplanned) whose `.data` is already a view/reshape.
K unplanned → `_transform_chunked!` → `wrapped_view(field, :, chunk)` = `Field(view(data,:,chunk))`
→ a **view-of-a-view** into `_fourier!` (custom rule) → Enzyme grid-field shadow degenerates `(0,1)`.

Repro ladder (all 1.12, `--check-bounds=yes`):
- `_mwe4_chunk.jl`: plain contiguous field + forced chunking (`transform_batch=[1]`) → **✅ works**
  (finite, nonzero grad). ⇒ chunking alone is NOT the trigger; the fused/view data layout is needed.
- `_mwe5_pwet_fast.jl`: real pwet `transform!(::Variables)`, fresh `initialize!` (no 5-day spinup)
  → expected to reproduce (crash is structural, value-independent) — running.

Working hypothesis: Enzyme can't build a correct shadow for `Field{<:SubArray}` that wraps a
view of an already-viewed/reshaped fused-field, across the custom `_fourier!` rule boundary.

Confirmed constraints (2026-07-05):
- `_mwe5_pwet_fast.jl` (fresh `initialize!`, no spinup) reproduces the exact BoundsError → fast
  (~1–2 min) fix harness.
- **Fused prognostic field is K=17** (vor4+div4+temp4+pres1+humid4); scratch is sized
  `max(planned_K)=4`, so chunking is REQUIRED to fit scratch (not just FFTW alignment). Disabling
  `_needs_chunking` → primal `DimensionMismatch` in `_legendre!` (K=17 vs scratch 4). Can't plan
  K=17 either (`@assert max(planned_K) ≤ nlayers=4`). ⇒ the chunked path MUST be made AD-safe.

**Fix attempt #1 — ✅ WORKS.** `_transform_chunked!` (spectral→grid) materializes **contiguous**
chunk buffers (`Array(view.data)` + copy-back) instead of `wrapped_view`. `_mwe5` now returns a
finite, nonzero gradient (was BoundsError/GC-corruption). Confirms root cause = Enzyme can't
shadow the view-of-fused-parent chunk field across the custom `_fourier!` rule; contiguous copies
fix it. Applied symmetric fix to grid→spectral `_transform_chunked!` too.

Validation:
- ✅ **AD (1.12):** `_mwe5` returns finite nonzero gradient (was BoundsError). transform!(::Variables)
  reverse FULLY WORKS on 1.12 — the primary target (where the offset ext change applies).
- ⚠️ **AD (1.10):** the fourier-chunking crash is GONE (0 `_apply_serial_fft`/`_transform_chunked`
  frames), but a **separate pre-existing** `GC error (probable corruption)` is now exposed in the
  `UV_from_vordiv!` **KernelAbstractions-Enzyme** reverse (`_UV_from_vordiv_kernel!` / KA
  `EnzymeExt.cpu_aug_fwd`/`cpu_rev`, transform.jl:145). Runs earlier in fwd → later in reverse, so
  was masked by the fourier crash. NOT introduced by the fix (fix only touched `_transform_chunked!`);
  a distinct KA+Enzyme-1.10 bug. 1.12 does not hit it.
- ✅ **Primal equivalence** (`_mwe4b_primal.jl`): chunked-path (fix) vs batched reference,
  spectral→grid `max|Δ|=0.0` exact, grid→spectral `8.4e-8` (fp noise). Both methods fixed
  (spectral→grid + grid→spectral `_transform_chunked!`).
- ✅ **SpeedyTransforms unit tests: PASS** (0 failures/errors; 1858 AD-rule tests + MatrixSpectral
  agreement 112 + all primal transform round-trips) — no regression from the core-code change.
- Barotropic diff testset: regression-safe by construction (nlayers=1 → `_needs_chunking` false
  → changed code not reached).
- **perf caveat:** adds per-chunk `Array` copies to the fused transform (hot path). → REJECTED by
  user (allocations + slowdown); replaced by the L25 plan below.
- FD gradient-correctness check still blocked by the separate `to_vec` harness bug (issue #2).

## L25 PLAN (2026-07-05) — replace the materializing chunk fix with an allocation-free one

Constraint: primal `_transform_chunked!` must stay allocation-free (`wrapped_view`); only the
AD path may change. The bug to dodge: Enzyme mis-constructs the shadow of the chunk view
(view-of-already-viewed fused parent) when differentiating *through* the chunking code.

1. **Benchmark the rejected fix** (confirm cost, document): `_mwe4c_bench.jl`, primal chunked
   transform both directions (trunc 31, nlayers 8, batch [1,4] → 2 chunks), time+allocs,
   current tree (materializing) vs `git stash` baseline (views). → then revert the fix.
2. **Alternatives, ranked:**
   - **E (preferred): alias EnzymeRule for the view constructor** (`field_view`/`lta_view`):
     `augmented_primal` returns primal = `Field(view(x.val…))` AND shadow =
     `Field(view(x.dval…))` — WE build the shadow, Enzyme never has to; reverse is a no-op
     (shadow view aliases the parent shadow, accumulation flows automatically). Smallest fix,
     zero runtime cost, repairs every AD use of `wrapped_view`, not just chunking.
   - **A (user suggestion): custom rule one level up** — for `_transform_chunked!` (or the 4-arg
     `transform!`): `augmented_primal` = primal loop; `reverse` = per-chunk inner
     `autodiff(Reverse, transform!, …)` with **hand-built** `Duplicated(view(val), view(dval))`
     chunk args. Also sidesteps the shadow construction; reuses the existing `_fourier!` rules;
     no adjoint math re-derivation. More code than E, but confined to the known-broken site.
   - **D (fallback): preallocated chunk buffers** in `ScratchMemory` — kills allocations, keeps
     2 copies/chunk of memory traffic.
   - Rejected: full adjoint-math rule at `transform!` level (must re-derive quadrature/coslat
     adjoints — maintenance burden); restructuring `FusedParents` for one-level views (invasive).
3. **Diagnostic to include:** dump `typeof(field_chunk.data)` in the failing path — is there a
   `ReshapedArray` in the view chain? Pins the poison type for E's coverage and for the
   upstream Enzyme report.
4. **Order:** bench (1) → revert materializing fix → E → validate (`_mwe4` no-regression,
   `_mwe5` gradient, primal equivalence, SpeedyTransforms unit tests, zero-alloc check) →
   if E fights EnzymeRules aliasing semantics, fall back to A → same gates.

### L25 progress

1. ✅ **Bench confirms rejection** (`_mwe4c_bench.jl`, trunc 31 / nlayers 8 / 2 chunks, 300 reps):
   spec→grid 154.3 µs & 302 KB/call (materializing) vs 94.1 µs & 135 KB (views) — **+64% time,
   +124% allocs**; grid→spec allocs likewise ~2.2× (time noisy). Materializing fix REVERTED
   (`git restore`; primal back to allocation-light views).
2. ✅ **E implemented** (`SpeedyTransformsEnzymeExt.jl`): alias `augmented_primal`/`reverse` for
   `wrapped_view` — shadow built explicitly as `wrapped_view(x.dval, …)` (aliases parent shadow;
   reverse = no-op), + Const passthrough methods. Zero primal impact (rules only affect Enzyme's
   compilation).
3. Gates & gradient cross-checks (2026-07-05):
   - ✅ `_mwe5` 1.12: alias rule → finite, nonzero gradient (R1).
   - ⭐ **R1 (alias) vs R2 (materializing reference) DIFFER**: max|Δ|=153.46, rel=1.0.
   - ⭐ **Arbiter** (`_mwe4`, plain fields + forced chunking, where Enzyme differentiates the
     chunk loop natively = ground truth): native vs alias-rule gradients **BIT-IDENTICAL**
     (max|Δ|=0.0). → **alias rule CORRECT**; and max|native|=153.4596 = exactly the R1/R2 Δ →
     **the materializing version was silently WRONG in the fused context** (dropped the coeffs
     cotangent: `Array(view-of-fused-parent)` hits the same degenerate-shadow bug, silently as
     a no-op instead of crashing). The rejected fix was not just slow — it was incorrect; the
     boolean finite+nonzero gate had been insufficient.
   - ✅ 1.10 `_mwe5`: fourier-chunking crash GONE under the alias rule (0 fourier frames);
     remaining GC corruption = the separate pre-existing `UV_from_vordiv!` KA+Enzyme-1.10 bug
     (unchanged status: 1.12 fully works, 1.10 blocked by that distinct bug).
   - ✅ SpeedyTransforms unit tests PASS (1858 AD-rule tests + MatrixSpectral 112, 0 failures)
     with the alias rule loaded.
   - Primal untouched by E (rules affect only Enzyme compilation; `_transform_chunked!` back at
     HEAD views) → no re-bench needed.

### L26 (2026-07-05/06) — unit tests EXPOSE a deeper problem: the chunked reverse itself

Adding the requested EnzymeTestUtils unit tests (`wrapped_view Enzyme rules` testset) surfaced
failures that reframe L25:
- Pkg.test (julia **1.11.9** — bare default, not 1.10 as assumed): 9 failures — ETU consumer
  tests (shadow of argument 1) + the chunked-vs-batched integration test.
- ETU consumer probes turned out to be a red herring: the artificial broadcast-between-
  wrapped-views pattern hits an unrelated `EnzymeNoDerivativeError` on 1.12 (`_mwe6_rulecheck`).
- **The real finding** (`_mwe6b_consistency.jl`, rule active, plain LTA/Field, trunc5 NL4):
  chunked-vs-batched `dcoeffs` **INCONSISTENT on BOTH 1.10 and 1.12** (identical numbers:
  max|chunked|=30.2, max|Δ|=29.7; dfield both 0 ✓). Combined with L25's arbiter
  (chunked-native ≡ chunked-alias, bit-identical): **the CHUNKED reverse itself — natively,
  rule-irrelevant — disagrees with the batched reverse.** Linear op + verified-identical primals
  → one side is plainly wrong. L25's arbiter only proved rule==native, never correctness.
- ✅ `_mwe6c_fd.jl` FD ground truth: **batched CORRECT** (max|AD−FD| = 2.2e-5), **chunked WRONG**
  (29.7 of 30.2).
- ⭐ `_mwe6d_layers.jl` per-layer signature: **only the LAST chunk's gradient survives** —
  chunks-of-1: layers 1–3 ratio 0.0, layer 4 ratio 1.0 EXACT; chunks-of-2: layers 1–2 zero,
  3–4 exact. Mechanism: Enzyme reuses the rule-returned view shadows of the LAST loop iteration
  for ALL chunk reverses; the last-reversed chunk computes correctly and zeroes its dfield
  slice; every earlier chunk's reverse then reads that already-zeroed shadow → contributes 0.
  Per-iteration shadows returned by custom rules are not reliably cached across loop iterations
  (worth reporting upstream with the L25/L26 distillation).
- **Consequence: alias-rule approach (E) abandoned; plan A implemented instead** —
  kwarg-free cores `_chunked_spec2grid!`/`_chunked_grid2spec!` split out in
  spectral_transform.jl; custom rules FOR THE WHOLE CHUNK LOOP in the ext: forward = primal
  loop unchanged (allocation-free views); reverse = replay chunks with EXPLICITLY constructed
  `Duplicated(wrapped_view(x.val…), wrapped_view(x.dval…))` per chunk + inner
  `autodiff(Reverse, transform!…)` per chunk (batched path inside, no chunk-loop
  differentiation → dodges BOTH the degenerate-shadow crash (a) and the last-iteration shadow
  reuse (b)). wrapped_view alias rules deleted. Gates rerunning.

### L25 CLOSED (2026-07-05) — shipped: alias EnzymeRules for `wrapped_view`

Final state: error #3 fixed by ~40 lines of EnzymeRules in `SpeedyTransformsEnzymeExt.jl` only
(alias `augmented_primal`/`reverse` for `wrapped_view`, Duplicated + Const passthrough); source
of `_transform_chunked!` untouched (allocation-free views). Validated: bit-identical to native
Enzyme on the ground-truth config; real pwet `transform!(::Variables)` reverse works on 1.12;
unit tests green; 1.10 fourier crash gone (separate `UV_from_vordiv!` KA bug remains, own issue).
Bonus finding: the materializing alternative was silently WRONG in the fused context (dropped
coeffs cotangent) — plan A (rule for `_transform_chunked!`) not needed. CHANGELOG updated (#1143).

Remaining open threads (unchanged): pwet issue #1 (spinup NaN instability), #2 (`to_vec`
Array→SubArray harness bug, blocks FD checks), 1.10 `UV_from_vordiv!` KA+Enzyme GC bug,
maxtypeoffset ext decision + pwet offset validation, upstream posts (#3275 MWE + degenerate-
shadow rule bug: a custom-rule arg whose wrapper wraps a view-of-view gets a (0,1) shadow —
worth its own Enzyme.jl issue with `_mwe4`/`_mwe5` distilled).

**SpeedyWeather source-fix implications (error #2):**
- Container swap (struct/tuple for the gradients NamedTuple): ruled out (H8/H8t), and
  **confirmed on the REAL code (2026-07-03):** `gradient_arrays` was temporarily patched to
  return a `GradientArrays` struct with identical fields (drop-in; consumers only destructure);
  after full re-precompile, `MWE_enzyme_mutable_transform.jl` behaved identically
  (m1 mutable+`transform!` ❌ same `EnzymeNoTypeError`, m1b/i1 ✅). Patch reverted —
  a NamedTuple lowers to the same aggregate as a struct; the container kind cannot mitigate.
- Slimming `gradients` below the threshold: fragile, not a real fix (threshold is global to S).
- E3-d3 (2026-07-02) showed passing `Duplicated(S)` **directly** works — the failure needs the
  getfield-from-mutable-parent. Restructuring differentiated entry points to take S as an
  explicit argument instead of loading it from `model` inside would dodge it, but every
  `model.spectral_transform` load inside the differentiated region is a landmine.
- E8/E9: an **immutable** model wrapper with all 18 real fields passes — making the model structs
  immutable would fix parameter AD wholesale, at the cost of a breaking API change.
- Until then: the gated `looseTypeAnalysis!(true)` workaround remains the pragmatic option for
  parameter AD on Julia ≥ 1.11 (verified byte-identical results to 1.10 earlier).

# PLAN — how to continue (written at session end, 2026-07-02; updated 2026-07-03 after L2/L3)

## A. Error #2 SpeedyWeather-free MWE (the active campaign)

All probes: `julia +1.12 --project=SpeedyWeather --check-bounds=yes <file>`, script must call
`Enzyme.API.looseTypeAnalysis!(false)` after `using` (harmless if the ext no longer sets it).
Reproducer baseline: [`MWE_enzyme_mutable_transform.jl`](MWE_enzyme_mutable_transform.jl) (m1 ❌).
The `_fourier!`/`_legendre!` split probe is kept as [`_mwe2_e10.jl`](_mwe2_e10.jl) (both halves ❌).

1. ✅ **DONE (`_mwe2_L2.jl`)** — pinned the L1-vs-E10 delta: annotations and precompilation ruled
   out; the shared-parent getfield extraction (`scratch.north`/`.south`/`.column` from one
   Duplicated `ScratchMemory`, plus `S = w.spectral_transform` from the mutable wrapper) is the
   confirmed trigger. See "L-series round 2" above.
2. ✅ **DONE, negative result (`_mwe2_L3.jl`)** — swapped SpeedyTransforms types for local mimics
   (`MyTransform`/`MyScratch`/`MyLTA`, matching shape incl. `Ptr`-carrying fake FFTW plans and the
   fat LTA 2nd field) while keeping L2's confirmed getfield pattern: **did not reproduce**, either
   wrapper. Real type identity is still necessary, not just shape — the naive "just match shape"
   synthesis path from error #1 does not carry over.
3. ✅ **DONE (`_mwe2_L4.jl`)** — real `SpectralTransform` is the (only) trigger type; within it,
   the `gradients` NamedTuple is the (only) trigger field (R6). All exotic suspects innocent.
4. ✅ **DONE (`_mwe2_L5.jl`)** — pure aggregate-size threshold (real LTAs: 4 ✓/8 ❌; mimic:
   8×exact-fat ❌ = SpeedyWeather-free reproducer; 16×medium ❌; no aliasing effect).
5. ✅ **ENDPOINT REACHED (2026-07-03)** — `_mwe2_L6.jl`/`_mwe2_L7.jl` stripped the scaffold
   (no scratch struct, slim 7-field transform, no KA, no eigen vecs; container kind irrelevant;
   trivial loop does NOT trigger → structured loop is ingredient #3). Final artifact
   [`MWE_enzyme_fat_aggregate.jl`](MWE_enzyme_fat_aggregate.jl): Enzyme+Base only, verified
   1.12 mut ❌ / imm ✅ and 1.10 both ✅. **Remaining: post it upstream (#3275 or follow-up).**

## B. Upstream reporting

- Append to https://github.com/EnzymeAD/Enzyme.jl/issues/3275 (or open a second issue): the
  mutability finding — `getfield` of a large immutable struct (SpectralTransform, ~26 fields)
  from a **mutable** Duplicated parent, consumed by a transform, `EnzymeNoTypeError`
  ("cannot deduce type of copy", 24-byte Vector aggregate); include
  `MWE_enzyme_mutable_transform.jl` (SpeedyWeather-based) until A.4 delivers a dependency-free one.
- Also report the side finding: `_fourier!`'s EnzymeRule lacks an `augmented_primal` method for
  `Const` scratch memory (raw MethodError; SpeedyTransformsEnzymeExt.jl) — decide whether to add
  the missing rule methods or document Duplicated-scratch as a requirement.

## C. Shipped-configuration decisions (SpeedyWeather side)

- The Enzyme ext currently ships WITHOUT `looseTypeAnalysis!` — with the no-tuple source fix,
  state AD (`Duplicated(vars)`, `Const(model)`) is clean on 1.12, but **parameter AD
  (`Duplicated(model)`) still throws error #2** (testsets 5/6/7). Decide: re-add the gated flag
  (verified safe: 1.12 results byte-identical to 1.10) OR accept no-parameter-AD on ≥1.11 until
  error #2 is fixed upstream. CHANGELOG (#1148) currently mentions the flag — align text with the decision.
- `transform!(::Variables)` FD test (barotropic.jl:163) is marginal at `atol=rtol=1e-3`
  (FD noise, AD deterministic — see `_chk_transform` note above): loosen to ~2e-3 or test a
  mean-abs statistic to deflake.

## D. Remaining test work (from the original task list)

- primitivewet.jl: run the dyn-core group on 1.12+1.10 (driver `_driver_pwet.jl` exists;
  runs were killed for CPU); add the physics group afterwards.
- Decide `runtests.jl` wiring for `enzyme_rules.jl` (currently not included).
- Keep both Julia versions green: after any of the above, re-run `_driver_st_baro.jl` on 1.10
  and 1.12 (≈20–50 min each; expected: 0 Enzyme errors, 3 known FD-fails ± the flaky 4th).

---

# L27 (2026-07-06) — validate Plan A (chunk-loop custom rules) — SESSION IN PROGRESS

Picking up from L26: the working tree has **Plan A** implemented (custom Enzyme rules for the
whole chunk loop `_chunked_spec2grid!`/`_chunked_grid2spec!`; forward = primal loop, reverse =
per-chunk replay with explicit `Duplicated(wrapped_view(val), wrapped_view(dval))` + inner
`autodiff`). The `wrapped_view` alias rules (Plan E) were deleted. L26 note ended "Gates
rerunning" but no results were recorded — this session runs the decisive correctness gates that
plan A must pass, since L26 proved the previous approaches were silently WRONG (only the last
chunk's gradient survived).

Decisive gates (the ones that killed E):
1. `_mwe6c_fd.jl` — chunked vs batched vs FD ground truth (batched known-correct; chunked was
   WRONG = 29.7 of 30.2 under E/native). Plan A must make chunked match FD.
2. `_mwe6d_layers.jl` — per-layer ratio (under E: only last chunk's layer ≈ 1.0, rest 0.0).
   Plan A must give all layers ≈ 1.0.
3. `_mwe5_pwet_fast.jl` — real pwet `transform!(::Variables)` reverse finite+nonzero on 1.12.
4. SpeedyTransforms unit tests — no regression.

Note (env): `julia +1.12` juliaup shim fails under the tool sandbox ("failed to load a
configuration file"); use the direct binary
`~/.julia/juliaup/julia-1.12.6+0.aarch64.apple.darwin14/bin/julia` instead (1.10.11 analogous).
First-run compile of these gates > 10 min (Enzyme+FiniteDifferences+SpeedyWeather) — run
in background with generous timeout.

Results (filled in as they land):
- `_mwe6c_fd.jl` (1.12): ⏳ running
- `_mwe6d_layers.jl` (1.12): ⏳ pending
- `_mwe5_pwet_fast.jl` (1.12): ⏳ pending
- SpeedyTransforms unit tests: ⏳ pending

⚠️ **ROOT CAUSE FOUND — STACK OVERFLOW (2026-07-06):** both gates ran ~48 CPU-min with no result;
the `tee` logs revealed the reason:
`Warning: detected a stack overflow; program state may be corrupted`. The nested
`Enzyme.autodiff(Reverse, _transform_uc!, ...)` inside the `_chunked_spec2grid!` **reverse rule**
makes Enzyme **re-enter its own compilation** (autodiff-of-autodiff reentrancy) → stack overflow.
The 48 min was stack-overflow *recovery churning*, not compilation progress. **Plan A as written
(plain `Enzyme.autodiff` in the rule) is broken**, not merely slow.

**FIX (applied):** replace `Enzyme.autodiff` with **`Enzyme.autodiff_deferred`** in the two
per-chunk pullback helpers (function annotated as `Const(_transform_uc!)` / `Const(_transform_g2s!)`,
mode `set_runtime_activity(Reverse)` to match the outer call). `autodiff_deferred` is the API
designed to be invoked from *within* an already-differentiated context (custom rules, GPU kernels)
— it compiles the inner differentiation as a separate deferred unit instead of re-entering the
outer Enzyme pass. Standard remedy for nested-autodiff reentrancy.

Verified the API signature in isolation (`_check_deferred.jl`): a trivial `g!(y,x)=y.=2x` inside a
`@noinline` helper via `autodiff_deferred(set_runtime_activity(Reverse), Const(g!), Const, …)`
returns the correct gradient `[2,2,2]` and runs in seconds (no reentrancy). ✅

**User change (2026-07-06):** `@propagate_inbounds` added to all view constructors —
`wrapped_view` ×2 (array_utils.jl), `field_view` ×3 (RingGrids/field.jl), `lta_view` ×3
(LowerTriangularArrays). Lets a caller's `@inbounds` elide bounds checks on the underlying
`view()` (error #3 originally surfaced as a `BoundsError` in `_apply_serial_fft!`). Compatible
with Plan A (neutral/helpful; the reverse rule constructs these views on val+dval directly). Kept
in the tree; the deferred gate below recompiled with it.

**Env note:** the user's source edits invalidated the LTA/RingGrids/SpeedyWeather precompile
caches; the first re-run of the gate FAILED not in Enzyme but in **precompilation** with
`operation not permitted (EPERM)` writing pidfiles under `~/.julia/compiled/` — a tool-sandbox
restriction. Re-launched with the sandbox disabled so precompilation can write its cache; warm
cache lets later runs use the sandbox again.

Results (deferred fix + user's @propagate_inbounds):
- `_mwe6c_fd.jl` (1.12): ❌ ~40 CPU-min compile, no result (superseded — approach abandoned)

## L28 (2026-07-06) — ⭐ NEW APPROACH (user-requested): manual-adjoint rule for transform!

`autodiff_deferred` (L27) removed the stack overflow but STILL compiled ~40 CPU-min with no result
→ nested autodiff is non-viable regardless. User directed: "derive a completely new rule directly
for transform!". This sidesteps compile blowup entirely — the rule body is PRIMAL code (forward/
inverse transforms + scalings), so it compiles like normal primal.

**Adjoint derivation (from legendre.jl + docs/spectral_transform.md, FD-validated):**
The transform is linear; synthesis (spec→grid) = InvLegendre then InvFFT; analysis (grid→spec) =
FwdFFT then FwdLegendre WITH the `solid_angles` (ΔΩ = sinθΔθΔϕ) quadrature weight. Reusing the
existing `_fourier!` rule's `adjoint_scale` for the FFT adjoint:
- **spec→grid pullback** (field_bar → coeffs_bar):
  `dg = adjoint_scale .* FFT_fwd(field_bar)`; (if unscale_coslat) `dg .*= coslat⁻¹` per lat;
  `dg ./= solid_angles` per lat (cancels FwdLegendre's ΔΩ); `coeffs_bar = FwdLegendre(dg)`.
- **grid→spec pullback** (coeffs_bar → field_bar):
  `df = InvLegendre(coeffs_bar)`; `df .*= solid_angles` per lat; `field_bar = InvFFT(df ./ adjoint_scale)`.
solid_angles/coslat⁻¹ are length nlat (full); the Legendre loops index `[j]` for j in 1:nlat_half,
so slice the first nlat_half.

**FD validation (no Enzyme, `_adjoint_check.jl`, trunc5 NL4, batched S):**
spec→grid rel **7.2e-7** ✅; grid→spec rel **5.0e-6** ✅. Both hand adjoints correct.

**Chunking wrinkle:** for the fused K=17 (S.nlayers=4) case, `_fourier_serial!` asserts
`nlayers ≤ S.nlayers`, so the pullback must chunk over K (≤ largest planned batch) exactly like the
primal. Layers are independent so the chunked pullback is numerically identical to the whole-K one;
still validating a chunked S (transform_batch=[1]) against FD before wiring into the rule.

**Plan:** rule stays on `_chunked_spec2grid!`/`_chunked_grid2spec!` (Plan A's split — fires only
on the chunked path; batched path differentiates correctly+cheaply natively via the `_fourier!`
rule). Replace the nested-autodiff reverse bodies with the validated manual-adjoint pullback
(chunked). augmented_primal runs the primal chunk loop unchanged (allocation-free views).

### ✅✅✅ L28 SHIPPED (2026-07-06) — analytic-adjoint rule fixes error #3, compiles cheaply

Implemented in `SpeedyTransformsEnzymeExt.jl`: `spec2grid_pullback!` / `grid2spec_pullback!`
(chunked, allocating reverse-pass temporaries) + rewritten `augmented_primal`/`reverse` on
`_chunked_spec2grid!` / `_chunked_grid2spec!`. NO nested autodiff — the rule body is primal
transforms + scalings, so it compiles like normal primal code (the fatal flaw of Plan A/deferred
gone). `make_zero!` on the overwritten output shadow after propagating, matching the `_fourier!`
rules.

**All gates PASS (2026-07-06):**
| gate | result |
|--|--|
| `_adjoint_check.jl` (FD, no Enzyme) | spec→grid rel 7e-7, grid→spec rel 8e-6 — batched & chunked |
| `_mwe6c_fd.jl` (1.12) | chunked `max|AD−FD|=2.2e-5` **match=true**; batched match=true |
| `_mwe6d_layers.jl` (1.12) | chunk-1/chunk-2/batched: **all layers ratio 1.0, max|AD−FD|=0.0** (L26 last-chunk-only bug GONE) |
| `_mwe5_pwet_fast.jl` (1.12, --check-bounds=yes) | real fused K=17 `transform!(::Variables)` reverse **AUTODIFF RETURNED, finite, nonzero, n=131189** (error #3 crash FIXED) |
| `_chunked_rules_test.jl` (unit testset) | **8/8 on Julia 1.10 AND 1.12** |

Unit test rewritten: `spectral_transform_ad_rules.jl` "wrapped_view Enzyme rules" testset →
"chunked transform Enzyme rules" (removed the red-herring `field_chunk_double!`/`lta_chunk_double!`
ETU probes that tested the deleted alias rules; kept+extended the chunked-vs-batched integration
check to both spec→grid and grid→spec). CHANGELOG updated (analytic adjoint wording).

Scratch scripts kept for reference: `_adjoint_check.jl`, `_chunked_rules_test.jl`, `_check_deferred.jl`.

### L28b (2026-07-06) — UNIFIED rule + zero-alloc + inactive_type(SpectralTransform)

Per user direction, three further changes on top of the L28 chunked fix:

1. **Unified rule (chunked + batched).** The `transform!` 4-arg methods now forward to positional
   cores `_transform_grid!` / `_transform_spec!` (spectral_transform.jl); the Enzyme ext rules the
   CORES, so EVERY transform! (chunked or batched) is a single analytic-adjoint AD boundary.
   `_mwe6c_fd` batched now bit-identical to chunked (both 2.16e-5 vs FD; batched previously took the
   native `_fourier!`-rule path → 2.20e-5). Non-chunked no longer relies on native `_legendre!`
   differentiation.

2. **Zero-alloc pullback (user: "allocating in the rule is unacceptable").** The pullbacks reuse the
   passed `ScratchMemory` (`.north`/`.south`/`.column`) for the freq intermediates and accumulate
   straight into the input cotangent — no per-chunk `zeros`. Two new backward-compatible primal
   hooks: `reset::Bool=true` on the forward `_legendre!` (spec→grid pullback accumulates into
   `coeffs.dval`), and `add::Bool=false` threaded through the inverse `_fourier!` chain (grid→spec
   pullback accumulates into `field.dval`). Only remaining allocations are the inherent per-ring FFT
   plan outputs (same as the primal) + `adjoint_scale` once per call (same as the `_fourier!` rules).
   Re-validated: `_mwe6c_fd` chunked & batched match FD (2.16e-5); primal unchanged.

3. **`inactive_type(::Type{<:SpectralTransform}) = true`** in the ext (user: test dropping
   maxtypeoffset). Error #2's residual trigger was Enzyme building a SHADOW of the fat S aggregate
   when loading it from the `Duplicated` mutable model (the "cannot deduce type of copy" 24-byte
   memcpy at the `S = model.spectral_transform` getfield — confirmed via `_err2_full.jl`, "Failure
   within method ft!"). Marking S inactive (it is fixed geometry, never differentiated) stops the
   shadow construction. **`_maxoffset_drop_test.jl` at DEFAULT offset (512): m1 mutable-wrapper +
   transform! → SUCCESS** (was FAILED). So transform! rule + inactive_type together fix error #2 at
   the default offset in the MWE. Full-scale barotropic validation running (`_driver_baro_default_offset.jl`,
   maxtypeoffset=512) — note it also surfaces the PRE-EXISTING FD `to_vec` Array→SubArray harness bug
   (issue #2 in L23; orthogonal to Enzyme), so the metric is "0 EnzymeNoTypeErrors", not "0 fails".

**Unit tests (user: "test the new EnzymeRule" + "use test_reverse"):** `spectral_transform_ad_rules.jl`
"chunked transform Enzyme rules" testset now has THREE layers of coverage on the chunked rule:
1. **`test_reverse` (EnzymeTestUtils vs FD)** — both directions, chunked S. ⚠️ transform! returns the
   mutated output array, which makes ETU's FD Jacobian go out of bounds (`BoundsError` on a square
   transpose); WRAP it in a `nothing`-returning closure (like `_fourier!`) → works. Verified
   **2106/2106 on 1.12** (`_test_reverse_check.jl`).
2. **chunked-vs-batched consistency** — both directions match (rtol 1e-4).
3. **accumulation** — reverse accumulates into a pre-seeded cotangent (guards reset=false / add=true).
`runtests.jl` **un-gated on 1.12** (was `if VERSION < v"1.12"`).

**Full SpeedyTransforms suite on 1.12 (`Pkg.test`, first run before the test_reverse wrapper fix):**
1868 pass, 0 fail, only the 2 un-wrapped `test_reverse` calls errored (ETU BoundsError) — fixed by
the wrapper; re-run in progress. (`Pkg.test` on 1.10 is blocked by a pre-existing stale-Manifest /
StyledStrings env error, unrelated; isolated checks used there instead.)

### maxtypeoffset DROPPED (validated)
`SpeedyWeatherEnzymeExt.__init__`'s `maxtypeoffset!(4096)` REMOVED. `_driver_baro_default_offset.jl`
(full barotropic testset, maxtypeoffset=512 default): **0 EnzymeNoTypeErrors** (the 24 remaining
errors are all the pre-existing FD `to_vec` Array→SubArray harness bug — orthogonal). The fix is
exact type analysis via `inactive_type` + the analytic-adjoint rule, not a widened analysis window.
Caveat: PrimitiveWet parameter AD not separately validated (blocked by the to_vec harness bug +
spinup NaN; both pre-existing and out of scope). CHANGELOG entry moved to its own PR line (#1151).

### Files touched (PR #1151)
- `SpeedyTransforms/ext/SpeedyTransformsEnzymeExt.jl`: inactive_type(SpectralTransform); analytic
  adjoint pullbacks `spec2grid_pullback!`/`grid2spec_pullback!` (scratch-reuse, no per-chunk alloc);
  rules on `_transform_grid!`/`_transform_spec!`.
- `SpeedyTransforms/src/spectral_transform.jl`: positional cores `_transform_grid!`/`_transform_spec!`
  (kwarg forwarders); `_transform_chunked!` de-indirected.
- `SpeedyTransforms/src/legendre.jl`: `reset::Bool=true` kwarg on forward `_legendre!`.
- `SpeedyTransforms/src/fourier.jl`: `add::Bool=false` threaded through the inverse FFT chain.
- `SpeedyTransforms/src/array_utils.jl`, `RingGrids/src/field.jl`, `LowerTriangularArrays/src/lower_triangular_array.jl`:
  `Base.@propagate_inbounds` on the view constructors (user; qualified from bare `@propagate_inbounds`).
- `SpeedyWeather/ext/SpeedyWeatherEnzymeExt.jl`: removed the maxtypeoffset workaround.
- `SpeedyTransforms/test/spectral_transform_ad_rules.jl` + `runtests.jl`: testset + un-gate.
- CHANGELOG #1151.

Scratch/driver scripts kept: `_adjoint_check.jl`, `_chunked_rules_test.jl`, `_test_reverse_check.jl`,
`_maxoffset_drop_test.jl`, `_driver_baro_default_offset.jl`, `_err2_full.jl`, `_check_deferred.jl`.

### ⭐ UV_from_vordiv! 1.10 KA+Enzyme GC bug — RESOLVED (bonus, 2026-07-06)
The separate pre-existing GC-corruption in the `_UV_from_vordiv_kernel!` KA+Enzyme reverse on Julia
1.10 (L25) is GONE. `_mwe5_pwet_fast.jl` on **1.10**: full `transform!(::Variables)` reverse
`AUTODIFF RETURNED, finite, nonzero, n=131189`, no GC error. Almost certainly fixed by
`inactive_type(SpectralTransform)`: the kernel reads `S.vordiv_to_uv_x`/`vordiv_to_uv1/2`/`l_indices`;
marking S inactive removes the shadow handling of those S-derived kernel args that corrupted the KA
reverse on 1.10. So the fused `transform!(::Variables)` reverse now works on **both** 1.10 and 1.12.

### FINAL VERIFICATION MATRIX (2026-07-06)
| check | 1.10 | 1.12 |
|--|--|--|
| `_adjoint_check.jl` (hand adjoint vs FD, no Enzyme) | rel ~1e-6 | rel ~1e-6 |
| `_mwe6c_fd` chunked==batched==FD | (1.12 shown) | ✅ 2.2e-5 |
| `_mwe6d_layers` per-layer | — | ✅ ratio 1.0, err 0.0 |
| `_mwe5_pwet_fast` fused K=17 reverse | ✅ finite/nonzero | ✅ finite/nonzero |
| `_test_reverse_check` (ETU vs FD, wrapper) | ✅ 2106/2106 | ✅ 2106/2106 |
| full `SpeedyTransforms` suite (un-gated AD-rules) | (Pkg.test blocked by env manifest bug) | ✅ passed |
| barotropic full diff test @ default maxtypeoffset | — | ✅ 0 EnzymeNoTypeErrors (24 FD `to_vec` errors = pre-existing, orthogonal) |

**Remaining pre-existing / out-of-scope (NOT from this work):** FD `to_vec` Array→SubArray harness
bug (blocks FD-side of the diff tests for view-backed fused Variables); PrimitiveWet spinup NaN
instability; PrimitiveWet parameter-AD not separately validated.

### OPEN QUESTION (user, 2026-07-06): extend the analytic adjoint to the NON-chunked transform? [ANSWERED: yes, done — see L28b]

Non-chunked (batched) transform! is CORRECT (validated) but differentiates via **native
`_legendre!`** + the `_fourier!` rules. That native `_legendre!` reverse — with the fat
`S = model.spectral_transform` getfield from the MUTABLE model + the S-derived loop bounds — is
exactly error #2's ingredient #3 (the differentiated loop), which is why `maxtypeoffset!(4096)`
(SpeedyWeatherEnzymeExt.jl:12) is needed on Julia ≥1.11. **Hypothesis:** giving the WHOLE
`transform!` (chunked AND non-chunked) the analytic-adjoint rule makes it an AD boundary, removing
the differentiated legendre loop → could eliminate error #2 for transform! and the maxtypeoffset
dependence, plus cheaper compile. Cost: the rule treats S as `Const`, dropping parameter
sensitivity through the transform GEOMETRY (S.legendre_polynomials/solid_angles) — but physical
model parameters don't flow through S, so parameter AD for real params is unaffected; state AD
unaffected. Needs: restructure the rule onto a positional transform! core that both paths call,
then re-run the barotropic parameter-AD test with maxtypeoffset OFF to see if it's still needed.

⚠️⚠️ **SECOND RED FLAG — `autodiff_deferred` COMPILES but ~30+ CPU-min with no result
(2026-07-06):** the deferred fix removed the stack overflow (no more corruption warnings, memory
stable ~900 MB), but the chunked autodiff still sits at 30+ CPU-min in compilation with no output.
Nested autodiff — even deferred — is evidently too expensive to compile *at all* for the real
transform! (custom `_fourier!` rules + `_legendre!` + FFTW plans inside the inner differentiation,
compiled as a separate unit while the outer autodiff also compiles). **Conclusion forming:
nested-autodiff-in-a-rule (Plan A) is not viable on compile-time grounds, regardless of gradient
correctness.** Letting the run finish in the background only to learn whether the *gradient* is
correct (informs whether the chunk-replay LOGIC is right, which any alternative would reuse).

**Alternatives to weigh (each experiment ≈ 30 min compile, so choose deliberately):**
- **Native chunk-loop differentiation** — ruled out by L26 (last-chunk-only wrong, rule-irrelevant);
  `@propagate_inbounds` addresses bounds/crash, not the last-chunk correctness. Unlikely to help.
- **Preallocated contiguous chunk buffers + `copyto!`** (L25 option D-ish, but done right): copy
  chunk in → batched transform on contiguous buffer → copy out. Enzyme differentiates natively
  (copyto! adjoint = scatter-add; batched transform = correct+fast). No nested autodiff → normal
  compile cost. Risk: the `view(fused_parent,:,chunk)` read's adjoint scatter may still hit the
  degenerate-shadow bug for the view-of-view fused case (distinct from L24/L25's `Array(view)`
  which *allocated* — copyto! into a preallocated buffer may dodge it). Needs an experiment.
- **Manual adjoint in the reverse rule** (call the primal opposite-direction transform with
  adjoint scaling per chunk; reuse `_fourier!` rule's `adjoint_scale`): cheap compile + runtime,
  but re-derives the Legendre/quadrature adjoint → maintenance burden + correctness risk.
- Decision likely needs the user (domain expert on this AD code): compile-cost vs. allocation vs.
  hand-derived adjoint. Surface options.
