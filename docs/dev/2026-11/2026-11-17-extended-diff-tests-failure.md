# Extended differentiability tests failing on main (Byrne–O'Gorman longwave / Enzyme 0.13.220)

> Status: **completed** (uncommitted, awaiting review). Root cause: LLVM GVN bug llvm/llvm-project#219885, triggered by
> Enzyme's `!invariant.group` tags; MWE down to 65 lines of IR. Workaround (one non-inlined CPU column loop per
> parameterization) and warnings fix (inactive tracer `iterate` rule) applied; full extended suite passes.

Date of initial draft: 2026-11-17

Base revision: b06af911abf015a6cb2bcba991959908dbb87a03 (origin/main)

## Originating prompt

> Investigate the failures of the extended differentiabilty tests. First focus on the actual failure in
> the end. Try both to fix the issue in our code (possibly with a workaround) and setting up a MWE of the
> problem. Then also look into the warnings. Start a markdown log

## Revision log

- 2026-11-17: log started.
- 2026-11-17: user asked whether to cap the Enzyme version until the LLVM fix is backported to Julia; tested
  0.13.199/0.13.208, both abort like 0.13.219/0.13.220, so no cap (see "Restricting the Enzyme version? No").
- 2026-11-17: `@noinline` on the struct-taking `column_parameterization_cpu!` reverted (EnzymeNoTypeError in
  `time_step!` reverse); replaced by the by-name variant `_column_parameterization_cpu_by_name!` in the repo.
- 2026-11-17: warnings fix applied (inactive tracer `iterate` rule, `@test_broken` → `@test`, `trunc` → `truncation`).
- 2026-11-17: the session was resumed while the earlier instance kept running; both worked in parallel for ~10 min
  (duplicate `old_*`/`val_*` jobs). The earlier instance stopped; its log edits are merged here.

## Problem description

The workflow "Extended Tests" (`.github/workflows/CI_extended_tests.yml`, runs
`SpeedyWeather/test/differentiability/runtests.jl` on Julia 1.10, ubuntu-latest) fails on main.

Recent runs (public API `actions/workflows/CI_extended_tests.yml/runs`, job annotations):

| Commit | Date (UTC) | Result | Annotation |
|---|---|---|---|
| b06af911a Byrne–O'Gorman longwave as default, default CO₂ (#1279) | 2026-11-16 | **failure** | exit code **134** (SIGABRT), log line 324, after 54 min of the test step |
| caed82d64 Jin ocean albedo (#1371) | 2026-11-13 | success | |
| 3245c69fb Fix NaN reverse-mode gradients from two-band shortwave (#1401) | 2026-11-13 | success | |
| 069fb09f3 Two-band shortwave (#1262) | 2026-11-12 | failure | exit code 1 (test failure, fixed by #1401) |
| bc837bd60 Precipitable water output (#1370) | 2026-11-11 | failure | exit code 1 |

Job logs need GitHub authentication; the user pasted the relevant part of the b06af911a log.

### The failure (CI log, b06af911a)

The abort happens while *compiling* the reverse pass of the testset
"Differentiability: PrimitiveWet parameterization_tendencies! (parameter AD)"
(`primitivewet.jl:130–148`: `trunc = 5, nlayers = 1`, 6 h spin-up,
`autodiff(set_runtime_activity(Reverse), parameterization_tendencies!, Const, Duplicated(vars, dvars), Duplicated(m, dm))`).
Enzyme emits an invalid LLVM module and the LLVM verifier, run at the start of Enzyme.jl's `post_optimize!`
(`compiler/optimize.jl:550/573`), aborts the process:

```
Both operands to ICmp instruction are not of the same type!
  %inbounds5781 = icmp ugt i64 %arraylen5780, float 0xFFFFFFFFE0000000, !dbg !10662
Stored value type does not match pointer operand type!
  store float %.pre52369, i64* %errorbox5783, align 8, !dbg !10662, !noalias !3449
 floatLLVM ERROR: Broken module found, compilation aborted!
[4983] signal 6 (-6): Aborted
```

This is an array bounds check (`%arraylen`, `%inbounds`, `%errorbox` are the names Julia 1.10's codegen gives a
bounds check: index < length, else store the index into the error box) whose `i64` *index* has been replaced
by a `float` value — the comparison constant `float 0xFFFFFFFFE0000000` is the float with bit pattern
`0xFFFFFFFF`, i.e. an integer `-1` reinterpreted as a float. Some code path did a replace-all-uses of an
`i64` value with a `float` value. Release builds of LLVM do not assert on type-changing RAUW, so the module
is only caught by the verifier.

Previous testsets in the same file ("dynamics_tendencies!", "implicit_correction!") passed.

### Warnings (CI log)

Typical warnings the user sees (barotropic tests, also in `primitivewet.jl`):

```
┌ Warning: <unknown>:0:0: in function preprocess_julia_update_prognostic__74835 void ({ … } addrspace(11)*, {} addrspace(10)*, {} addrspace(10)*):
  Enzyme: Cannot deduce type of load   %unbox19.unpack27 = load i64, i64 addrspace(10)* %unbox19.elt26, align 8, !tbaa !87, … assumed {[-1]:Pointer}
└ @ LLVM ~/.julia/packages/LLVM/JGisQ/src/core/context.jl:175
```

and the same for `preprocess_julia_vertical_advection_…`. Further noise: the deprecation warning
"`trunc` is deprecated, use `truncation`" from `primitivewet.jl:136`, and many
"Time step changed … to match output frequency" infos.

## Background

Two things changed between the last green run (caed82d64, 2026-11-13) and the abort (b06af911a, 2026-11-16):

1. **Our code** (#1279): `ByrneOGormanLongwaveTransmissivity` becomes the default longwave transmissivity
   of `OneBandLongwave`, and `PrimitiveWetModel`/`PrimitiveDryModel` get
   `greenhouse_gases = (; co2 = CO2(spectral_grid))` by default. The new transmissivity reads the CO₂
   prognostic scalar (`vars.prognostic.greenhouse_gases.co2[]`, a `Ref`), takes a `log` of it and reads the
   humidity of the grid variables.
2. **Enzyme**: v0.13.220 (Enzyme_jll 0.0.310) was released 2026-11-15T22:45Z, so the Nov 16 CI run resolved
   it; the green Nov 13 runs had ≤ 0.13.219 (Enzyme_jll 0.0.307).

Julia itself did not change: 1.10.13 was released 2026-11-02, before the green runs.

## Setup

Work directory (outside the repo, scripts, envs, logs): `/p/tmp/maxgelbr/code/claude_enzyme_extended_2026-11-17/`

- Julia 1.10.13 official binary in `/p/tmp/maxgelbr/julia-1.10.13` (same as CI's `setup-julia '1.10'`).
- `env220/`: copy of the differentiability `Project.toml`, the five monorepo packages `Pkg.develop`ed by path
  (Julia 1.10 ignores `[sources]`), Enzyme 0.13.220 — mirrors CI.
- `env219/`: same, Enzyme pinned to 0.13.219 (Enzyme_jll 0.0.307).
- `run_subset.jl <files…>`: like `runtests.jl` but includes only the given test files.
- `job.sh`: SLURM wrapper (qos priority, 4 CPUs) running Julia 1.10.13 without extra flags, as CI does.

## Investigation log

### 2026-11-17

- Submitted `full220` (whole `runtests.jl`, Enzyme 0.13.220) and `pwet220` (`primitivewet.jl` only).
- `pwet220` reproduces the CI abort exactly — same instruction *names* (`%inbounds5781`, `%arraylen5780`,
  `%errorbox5783`), so the broken module is deterministic. Note: in a SLURM log, stdout (test summaries)
  written before the abort is lost (buffered), stderr (`@info`) survives.
- Standalone reproducer `repro_param_ad.jl` (only the parameter-AD testset: T6 L1 PrimitiveWet, 6 h spin-up,
  reverse `parameterization_tendencies!` with `Duplicated(model)`), aborts ~2.5 min after start:

  | Job | Enzyme | longwave transmissivity | greenhouse gases | Result |
  |---|---|---|---|---|
  | r220_default | 0.13.220 | ByrneOGorman (new default) | `(; co2)` (new default) | **abort** (same IR) |
  | r219_default | 0.13.219 | ByrneOGorman | `(; co2)` | **abort** |
  | r220_noghg | 0.13.220 | ByrneOGorman | `(;)` | **abort** |
  | r220_frierson | 0.13.220 | Frierson (old default) | `(; co2)` | pass (69.7 s) |
  | r220_frierson_noghg | 0.13.220 | Frierson | `(;)` | pass (68.7 s) |

  ⇒ **Not an Enzyme release regression.** Since the 0.13.219 run aborts too, the trigger is our code from #1279, specifically
  `ByrneOGormanLongwaveTransmissivity` in the column loop. The default CO₂ is irrelevant.
- Ablation of the transmissivity (copies of the scheme in `bo_variants.jl`, `repro_param_ad.jl bo=<variant>`):
  **every variant aborts**: exact copy, no humidity read (`noq`), humidity view hoisted, `ifelse` instead of
  `max`, no `max`, no CO₂ `log`/`max` (`nodry`), no `pressure_thickness` (`nodp`). Only the *type* of the bad
  replacement value changes: CI and `nodp` get `float 0xFFFFFFFFE0000000` (float bit pattern `0xFFFFFFFF`);
  the other variants get an array data pointer named `%arrayptr…_unwrap.pre`:
  ```
  %4000 = add float addrspace(13)* %arrayptr547511718_unwrap.pre, i64 -1
  store float addrspace(13)* %arrayptr547511718_unwrap.pre, i64* %errorbox5818
  ```
  A replacement value of arbitrary type points to a dangling `Value*` in some cache, not to a semantic
  problem in our code.
- IR dump (`repro_param_ad.jl dump` sets `Enzyme.Compiler.DumpPostEnzyme[]`, `DumpPrePostOpt[]`; stderr in
  `logs/r220_dump-*.ll`, 243 MB). The pre-`post_optimize!` module is saved as `ir/pre_post_optimize.ll` (93 MB, 86 functions).
  - In the module **after Enzyme's AD** and **before `post_optimize!`** the bounds check is well-typed:
    `%inbounds5781 = icmp ult i64 %3817, %arraylen5780`, `store i64 %.pre10513.unpack11601, i64* %errorbox5783`.
    ⇒ The AD transform is fine. The module is broken by the LLVM pass pipeline that Enzyme.jl runs in
    `post_optimize!` (`addOptimizationPasses!` + Julia legalization + machine passes, `compiler/optimize.jl`).
  - With the dump enabled (different heap layout) the same compile does not abort but **segfaults in
    `llvm::ScalarEvolution::getOperandsToCreate` ← `computeExitLimitFromICmp` ← `LoopDeletionPass`** inside the
    same `post_optimize!` run. Together with the arbitrary replacement types this suggests a stale `Value*` in an analysis
    cache (SCEV), i.e. some pass in the pipeline modifies IR without invalidating SCEV correctly.
  - Location (via `ir/dbgchain.py` on the `!dbg` of the bad instruction): function
    `diffejulia_column_parameterizations_cpu__11158` (46 160 lines, reverse pass of the fully inlined column loop),
    `getindex` ← `pressure_thickness` (`Δσ[k]`, vertical_coordinates.jl:124) ← `flux_to_tendency` (tendencies.jl:159)
    ← `surface_flux_to_tendency` (tendencies.jl:98) ← **`longwave_radiative_transfer!` (longwave_radiation.jl:263)**:
    `dTdt[ij, nlayers] += surface_flux_to_tendency(U / cₚ, pₛ, model)`. The index that gets replaced is
    `model.geometry.nlayers`. So the broken instruction is in the longwave radiative transfer, not in the new
    transmissivity. Swapping the transmissivity only perturbs the module enough to trigger the pipeline bug.
- Enzyme.jl's `addOptimizationPasses!` gained between 0.13.199 and 0.13.208: `SimpleLoopUnswitchPass(nontrivial =
  true)` + `JuliaLICMPass` in the first loop pipeline, `IRCEPass` (inductive range-check elimination, rewrites
  `icmp ult %idx, %len` bounds checks), `LoopInstSimplifyPass` + `LoopFullUnrollPass`, and a vectorizer block.
  0.13.208–0.13.220 have identical pipelines. Enzyme.jl's own Julia-side passes (`restore_alloca_type!`,
  `reinsert_gcmarker_pass!`, `safe_atomic_to_regular_store!`) report changes correctly, and LLVM.jl then invalidates all
  analyses, so they are not the stale-cache source.
- `patch_pipeline.jl` (env `PIPE=…`) re-evaluated `addOptimizationPasses!`/`post_optimize!` with textual edits.
  **This does not work**: Enzyme compiles from inside its `@generated` `thunk_generator`, which runs in the world age of
  its definition, so redefinitions from a script are invisible (backtraces still show the original
  `optimize.jl:550/573`). The six `pipe_*` results are void. (The differentiated *user* code is looked up in the
  caller's world, so overriding SpeedyWeather methods from a script does work.)
- Instead: `Enzyme_dev/` = copy of Enzyme 0.13.220, `Pkg.develop`ed into `env220dev/`, whose
  `addOptimizationPasses!`/`post_optimize!` read `ENV["ENZYME_PIPE"]` at runtime (`skip:<Pass>` drops a pass via a
  `@_maybe_add` macro; `verify` sets `verify_each`; `debuglog` sets `debug_logging`).

  | ENZYME_PIPE | Result |
  |---|---|
  | (none) | abort, same IR as CI (patched copy behaves like the release) |
  | `verify,debuglog` | **"Broken function found" right after `GVNPass` on `diffejulia_column_parameterizations_cpu__11181`** (38 354 instructions); verified fine after every earlier pass. GVN computed `MemoryDependenceAnalysis` and `PhiValuesAnalysis` fresh just before. Broken IR: `%3909 = add float %.pre52370, i64 -1`, `store float %.pre52370, i64* %errorbox5783` |
  | `skip:GVNPass` (only the GVN in `addOptimizationPasses!`) | Julia error (exit 1) instead of abort: `verify(mod)` in `JIT.move_to_threadsafe` finds `add float addrspace(13)* addrspace(11)* %_unwrap6553, i64 -1`. The later GVN runs in the legalization/machine pipelines break it instead |
  | `skip:IRCEPass`, `skip:SimpleLoopUnswitchPass`, `skip:JuliaLICMPass`, `skip:LoopDeletionPass`, `skip:IndVarSimplifyPass`, `skip:LoopFullUnrollPass` | abort (unchanged) |

  ⇒ **The culprit pass is LLVM's GVN** (Julia 1.10 ships LLVM 15). It replaces an `i64` value with an unrelated value
  of a different type. `%.pre…` names come from GVN's load PRE, which creates the new load with the type of the load it
  replaces, so a float `.pre` load can only reach `i64` users through corrupted bookkeeping: a stale entry keyed by
  `Value*` in GVN's value table or MemDep cache that now points at a re-allocated instruction. Together with the SCEV
  segfault under a different heap layout (dump run) this is a use-after-free-style LLVM bug, not a miscompile of our code.
  Same family as #1290 (LLVM segfault in Enzyme's reverse pass of the fused radiation loop, worked around in #1291
  by splitting the CPU loop) and as the 1.11.7 GVN segfault in SpeedyCalibration's `parameterization_tendencies!`
  gradient noted in October.

### Workaround in our code

Candidates as method overrides in the reproducer (`workarounds.jl`, knob `wk=`), all on the default model:

| Workaround | Result |
|---|---|
| `@noinline` on the generic `column_parameterization_cpu!` (every scheme's grid-point loop its own function) | pass (70.7 s) |
| `@inbounds` around the CPU grid-point loop (as the GPU kernel, `inbounds = true`) | pass (59.3 s) |
| `@noinline` only on the longwave loop | pass (79.8 s) |

Any change of the module's shape may dodge an LLVM bug of this kind, so the choice is about robustness. `@inbounds`
is nullified by `--check-bounds=yes`, which CLAUDE.md requires for unit tests. `@noinline` is structural: Enzyme
then emits one reverse function per scheme instead of one 38k-instruction function for all fused column loops.
Its runtime cost is one call per scheme per time step.

**First attempt, reverted**: `@propagate_inbounds` → `@noinline` on the generic `column_parameterization_cpu!` in
the repo (patch saved as `noinline_attempt1.patch` in the work directory). Only the plain-CPU path uses this method;
GPU and Reactant go through `column_parameterizations_kernel!`.

| Job (repo change, env220, Enzyme 0.13.220) | Result |
|---|---|
| `fix_param`: reproducer, default model | pass (73.1 s) |
| `fix_param_cb`: same with `--check-bounds=yes` | pass (73.6 s) |
| `fix_pwet`: `primitivewet.jl` | all 4 testsets pass (parameter AD 1m27s) |
| `fix_full`: full `runtests.jl` | **fails**: "forward vs reverse time_step! (PrimitiveWetModel)" errors in reverse mode with `EnzymeNoTypeError` in `column_parameterizations_cpu!` → `get_parameterizations(model)` |

As a real call, each parameterization struct taken from the `get_parameterizations(model)` NamedTuple has to be
materialized as an argument (untyped copy) under `Duplicated(model)`, which Enzyme's type analysis cannot follow.
The regression-guard comment in `primitivewet.jl:131-135` describes the same hazard. Reverted.

Second matrix, now also against `time_step!` (`repro_timestep.jl` = the PrimitiveWet case of `forward_mode.jl`:
forward and reverse `time_step!` with `Duplicated(model)`, T6 L2 Float64, adjoint identity check):

| Workaround | parameter-AD reproducer (CI abort) | `time_step!` fwd vs rev |
|---|---|---|
| none (main) | abort | pass (rel. diff 2e-15) |
| `@inbounds` CPU loop | pass | pass |
| `@noinline` generic loop | pass | `EnzymeNoTypeError` |
| `@noinline` longwave loop only | pass | pass |
| `valnoinline`: model path calls one `@noinline` function per parameterization *by name* (`Val(name)`, `getfield(model, name)` inside, so no struct is passed) | pass (69.5 s) | pass |

`valnoinline` keeps the structural benefit (one reverse function per scheme) without materializing structs. It
leaves `_column_parameterizations_cpu!(vars, ::NamedTuple, model)` untouched, which the convection tests call with
parameterizations that are not in the model. **Applied to the repo** (see Summary of changes). Every one of these
is a way of *not triggering* an LLVM bug (see root cause below), not a fix.

### MWE

1. **Standalone pipeline on saved IR** (`mwe_post_optimize_ir.jl`): parse `ir/pre_post_optimize.stripped.ll`
   (the module that entered `post_optimize!` in the dump run; `!enzymejl_allocart` metadata stripped because
   `RestoreAllocaType` dereferences the Julia object pointers stored there), verify the input (valid), run
   `Enzyme.Compiler.post_optimize!(mod, JIT.get_tm(); callconv = false)` in a fresh process without SpeedyWeather.
   **Reproduces the identical breakage** (`%inbounds5781 = icmp ugt i64 %arraylen5780, float 0xFFFFFFFFE0000000`,
   `store float %.pre34106, i64* %errorbox5783`, "Broken module found"). So the bug is deterministic given the IR,
   independent of the original process's heap.
   Pitfall: Enzyme's block names (`inverttop_amerge_amerge_…`) exceed LLVM's 1024-character
   `-non-global-value-max-name-size`; without raising it the parser truncates them, merges distinct labels and aborts
   with "Terminator found in the middle of a basic block!".
2. `mwe_split.jl`: same pipeline split right before the culprit GVN. `prefix` → `ir/pre_gvn.ll` (valid, 86 MB);
   `gvn` = **a single `GVNPass` in a fresh process breaks it** (`%3909 = add float %.pre52, i64 -1`,
   `store float %.pre52, i64* %errorbox5783`); `rest` = the remaining pipeline aborts as in CI. So no analysis
   state from earlier passes is needed: GVN alone on this IR is enough.
3. **Pure LLVM reproducer, no Julia process**: `opt` from `LLVM_full_jll` 15.0.7+14 (`LLVM version 15.0.7jl`, the
   LLVM of Julia 1.10):
   ```
   opt -passes=gvn -aa-pipeline=scoped-noalias-aa,tbaa,basic-aa -non-global-value-max-name-size=1000000 pre_gvn.ll -S -o /dev/null
   ```
   → "Both operands to a binary operator are not of the same type! … Broken module found". Still reproduces after
   `llvm-extract --func=diffejulia_column_parameterizations_cpu__11158` (`ir/colfun.ll`, 13 MB, 56k lines) and
   after stripping debug info and names (`ir/colfun_s.ll`, 12.5 MB).
4. Newer LLVM (`opt` from `LLVM_full_jll`, all Julia-patched builds) on `ir/colfun_s.ll`:

   | LLVM | Julia | GVN result |
   |---|---|---|
   | 15.0.7jl | 1.10 | broken: `i64` users get a `float` value (type mismatch) |
   | 16.0.6jl (`-opaque-pointers=0`) | 1.11 | broken: "Instruction does not dominate all uses! `<badref>` = load i64" (use of an erased load) |
   | 18.1.7jl | 1.12 | same as 16 |
   | 23.1.4jl (newest) | — | **segfault** in `GVNPass::AnalyzeLoadAvailability(LoadInst*, ReachingMemVal const&, Value*)` → `Instruction::getFunction()` on a dangling instruction |

   ⇒ A long-standing use-after-free in GVN's load analysis, i.e. a dependency that points at an erased load. It is
   present in every LLVM from 15 to 23, so not a regression. On opaque-pointer LLVMs it shows up as a dangling use or
   crash instead of a type mismatch. Not yet checked: a vanilla (non-Julia-patched) LLVM.
5. `llvm-reduce -j 16` with `ir/interesting.sh` (GVN output fails the verifier with a type mismatch;
   `llvm-reduce` only tests valid candidates) on `ir/colfun_s.ll`: 12.5 MB → 381 lines (`ir/reduced1.ll`,
   ~15 min), second round → **65 lines** (`ir/reduced2.ll`, copied to `mwe/gvn_invariant_group_llvm15.ll`;
   opaque-pointer variant `mwe/gvn_invariant_group_opaque.ll`). Two `load i64 … !invariant.group`, one plain load,
   one float store in a small loop nest; GVN rewrites `%4 = add nsw i64 %2, 0` with a non-`i64` operand:

   | `opt -passes=gvn -aa-pipeline=scoped-noalias-aa,tbaa,basic-aa` on the 65-line file | Result |
   |---|---|
   | 15.0.7jl (Julia 1.10) | broken module (type mismatch), crash while printing it |
   | 16.0.6jl (Julia 1.11) | same |
   | 18.1.7jl (Julia 1.12) | "Instruction does not dominate all uses!" |
   | 23.1.4jl | valid (but 23 still crashes on the full function and on upstream's reproducer, below) |
6. **Stripping `!invariant.group` from `ir/colfun_s.ll` (sed) makes GVN's output valid.** These tags come from
   Enzyme's C++ AD (`CacheUtility.cpp`, `GradientUtils`, `AdjointGenerator.h`: cache stores and re-loads of primal
   values): 5 925 in the post-AD module, 3 530 in the offending function. Neither Julia's codegen nor Enzyme.jl adds them.
7. **End-to-end confirmation**: `Enzyme_dev` option `ENZYME_PIPE=stripig` drops `!invariant.group` from non-pointer
   loads right before the `post_optimize!` pipeline (1 478 loads). On *unmodified* main (worktree `wt_main`, env
   `envdevmain`) the reproducer then **passes**; the same build without the option aborts as in CI.

### Root cause: LLVM bug llvm/llvm-project#219885 (fixed upstream by #219929, merged 2026-09-08)

Issue #219885 "`!invariant.group` on a non-pointer load crashes in GVN" has the same stack as our LLVM 23 crash
(`Instruction::getFunction()` ← `GVNPass::analyzeLoadAvailability(LoadInst*, ReachingMemVal const&, Value*)`).
Fix PR #219929 "[MemDep] Remove cached non-local defs referencing deleted instructions":
`MemoryDependenceResults::removeInstruction()` removed `NonLocalDefsCache` entries keyed by a removed non-pointer
instruction, but not entries whose cached *dependency result* is that instruction. The stale result stays reachable
after GVN erases the load, and a later non-local load query consumes it. LLVM 15.0.7's `removeInstruction` lacks
exactly the block the PR adds (checked in the llvmorg-15.0.7 source). The PR's regression test is named
`pr46054-md-nonlocaldefcache-cleanup.ll`, so the defect is at least as old as 2020. The Mono project hit it by
tagging non-pointer loads with `!invariant.group` and gated those tags off as a workaround.

Enzyme tags many non-pointer loads (`load i64`/`load float … !invariant.group`), so Enzyme-generated reverse passes
are an ideal trigger. Whether it fires depends on GVN's processing order and on which loads get erased. That is why
unrelated code changes switch it on or off: the new transmissivity in #1279 here, the fused radiation loop in #1290.

- The 15-line reproducer from #219885 (`mwe/upstream_219885.ll`, typed variant `mwe/upstream_219885_typed.ll`)
  crashes `opt` 23.1.4jl but not 15/16/18, whose GVN lacks the `ReachingMemVal` path it needs. Our 65-line file
  covers 15/16/18.
- The fix is on LLVM main only, not in any LLVM that Julia ships: 15 (1.10), 16 (1.11), 18 (1.12), 23.1.4jl.

### Restricting the Enzyme version? No

User question: since the bug is known and fixed upstream, restrict the Enzyme version until the fix is backported to
Julia? Reproducer on unmodified main (worktree `wt_main`), Julia 1.10.13, each Enzyme pinned in its own env
(`envmain013199`, `envmain013208`, `env219`, `env220`):

| Enzyme (release date) | Enzyme_jll | Result |
|---|---|---|
| 0.13.199 (2026-08-09, older post-AD pipeline without IRCE/unswitch/vectorizer) | 0.0.290 | abort, identical IR |
| 0.13.208 (2026-09-30) | 0.0.299 | abort, identical IR |
| 0.13.219 (2026-11-07, last green CI, with the old code) | 0.0.307 | abort |
| 0.13.220 (2026-11-15) | 0.0.310 | abort |

The defect is in Julia's bundled libLLVM, and every one of these Enzyme versions emits the triggering
`!invariant.group` tags and runs GVN after AD, so an Enzyme compat bound cannot avoid it. **No compat change made.**
Real fixes, all upstream: (a) backport llvm/llvm-project#219929 into Julia's LLVM fork plus Julia patch releases
(1.10.x for CI); or (b) Enzyme(.jl) stops emitting `!invariant.group` on non-pointer loads, or strips it before its
post-AD pipeline. `ENZYME_PIPE=stripig` above shows that (b) works end to end.

### Warnings: "Enzyme: Cannot deduce type of load … assumed {[-1]:Pointer}"

Emitted by Enzyme's `AdjointGenerator::visitLoadLike` when a load is not provably constant (runtime activity) and its
type is unknown. Enzyme then falls back to `defaultTypeTreeForLLVM(i64)` = Pointer. All warning sites (`update_prognostic!`,
`vertical_advection!`, the test helper `count_tracers`) iterate `model.tracers::Dict{Symbol, Tracer}`.
SpeedyWeatherEnzymeExt has *custom* forward/augmented/reverse rules for `iterate(::TRACER_DICT, …)`, which return the
boxed `Union{Nothing, Tuple{Pair{Symbol, Tracer}, Int}}`. The `i64` loaded from that box (the `Int` iteration state)
has no type information and is not provably inactive.

Experiment (worktree `wt_warn`, env `envwarn`, Enzyme 0.13.220): replace the three rules by
`EnzymeRules.inactive(::typeof(Base.iterate), ::SpeedyWeather.TRACER_DICT, args...) = nothing`. Enzyme's interpreter
blocks inlining for calls with inactive rules (`interpreter.jl:547/567`), just as for custom rules, so `iterate(::Dict)`
is never differentiated. Probe `warn_probe.jl` (parameters.jl, Barotropic `update_prognostic!` reverse, forward_mode.jl):

| Stage | main | inactive rule |
|---|---|---|
| parameters.jl (`count_tracers` testset) | 2 warnings | 0; the `@test_broken` "forward mode over a populated registry" (`EnzymeNoShadowError`) now **passes** |
| `update_prognostic!` Barotropic reverse | 1 | 0 |
| forward_mode.jl (Barotropic, ShallowWater, PrimitiveWet fwd vs rev) | 12, all 3 testsets pass | 0, all 3 testsets pass (9/9 each) |

**Applied to the repo** with the `@test_broken` turned into `@test` (see Summary of changes).

The other warning in the CI log, "`trunc` is deprecated", came from `primitivewet.jl:136` (`SpectralGrid(trunc = 5)`).
It is now `truncation = 6`, the same resolution. The "Time step changed … to match output frequency" lines are
`@info`s, not warnings, and were left alone.

## Summary of changes

Uncommitted on branch `mg/enzyme-failure` (base b06af911a). `SpeedyWeather` is already `0.25.0-DEV`, so no version
bump is needed. A CHANGELOG line is still to be added once there is a PR number.

1. `SpeedyWeather/src/parameterizations/tendencies.jl`, workaround for llvm/llvm-project#219885, **reversible by one
   constant**: `column_parameterizations_cpu!(vars, model)` branches on `const SPLIT_COLUMN_LOOPS = true` (constant-folded,
   no runtime cost). `true` unrolls over the *names* of the model's parameterizations
   (`_column_parameterizations_cpu_by_name!`), and each parameterization runs in its own `@noinline` function
   `_column_parameterization_cpu_by_name!(vars, ::Val{name}, model)` with `getfield(model, name)` inside: smaller Enzyme
   reverse functions, no struct copied into a non-inlined call. `false` is the original fused path
   `_column_parameterizations_cpu!(vars, get_parameterizations(model), model)`, which is unchanged. The constant and both
   new functions sit in one block delimited by `>>> WORKAROUND llvm/llvm-project#219885` / `<<< END WORKAROUND`, whose
   comment says how to check whether it is still needed and what to delete. Only the plain-CPU path is affected; GPU
   and Reactant use `column_parameterizations_kernel!`.
   - Revert check: `docs/dev/2026-11/2026-11-17-gvn-invariant-group.ll` (our 65-line reproducer, breaks LLVM 15/16/18)
     and `…-upstream.ll` (the #219885 reproducer, crashes LLVM 23.1.4jl). `opt -passes=gvn` on both with the `opt` of
     the matching `LLVM_full_jll` must give valid output. Decisive: `SPLIT_COLUMN_LOOPS = false` +
     `test/differentiability/primitivewet.jl`.
   - New unit test `SpeedyWeather/test/parameterizations/column_loops.jl`: both paths (and the default) give identical
     parameterization variables, tendencies and scratch for PrimitiveWet and PrimitiveDry, so the switch can be
     flipped without changing the model. To be deleted together with the workaround.
2. `SpeedyWeather/ext/SpeedyWeatherEnzymeExt.jl`: the custom forward/augmented-primal/reverse rules for
   `iterate(::TRACER_DICT, …)` are replaced by `EnzymeRules.inactive(::typeof(Base.iterate), ::TRACER_DICT, args...)`.
   This removes the "Cannot deduce type of load" warnings and makes forward mode over a populated registry work.
3. `SpeedyWeather/test/differentiability/parameters.jl`: comment updated; the `@test_broken` for forward mode over a
   populated tracer registry is now a `@test`.
4. `SpeedyWeather/test/differentiability/primitivewet.jl`: regression-guard comment updated;
   `SpectralGrid(trunc = 5)` → `SpectralGrid(truncation = 6)` (same resolution, no deprecation warning).

## Testing and verification

All on Julia 1.10.13, Enzyme 0.13.220 (as CI), unless noted.

| Test | Result |
|---|---|
| `repro_param_ad.jl` with the by-name workaround (override) | pass |
| `repro_timestep.jl` (forward_mode.jl PrimitiveWet case) with the by-name workaround (override) | pass, adjoint identity rel. diff 2e-15 |
| `warn_probe.jl` with the inactive iterate rule (worktree) | 0 warnings (main: 15), forward_mode.jl 3/3 testsets pass |
| `parameterizations/all_parametrizations.jl` + `convection.jl`, `--check-bounds=yes`, repo changes | 61/61 and 1148/1148 pass |
| Runic `--check` on the 4 changed files | clean |
| Full `differentiability/runtests.jl` (as CI), repo changes: job `final_full` | **pass, all 23 testsets**, 45.5 min; 0 "Cannot deduce type" and 0 deprecation warnings. Left: "Semi-implicit centering = 0.5 < 1 …" (a test's model configuration) and Enzyme's "fallback BLAS replacements for dgemm_64_" |

## Known limitations

- The code change does not fix the bug, it only avoids triggering it. Any later change to the column physics can
  trigger it again (as #1290 and #1279 did), until Julia ships an LLVM with llvm/llvm-project#219929.
- Not run locally: GPU tests, Reactant tests (`test/reactant/` loads Enzyme, but Reactant traces through its own
  `autodiff` path), the main test suite apart from the two parameterization files. The shared
  `SpeedyWeather/test/Manifest.toml` (Julia 1.12.2) is stale locally ("failed to find source of parent package:
  GPUArrays"), so the unit tests ran in the work-directory env on Julia 1.10.13.
- The 65-line MWE reproduces on Julia's LLVM 15/16/18, not on 23.1.4jl. 23 crashes on the full function and on
  upstream's 15-line reproducer instead. Not checked against a vanilla (non-Julia-patched) LLVM.

## Future work

- Ask for a backport of llvm/llvm-project#219929 to Julia's LLVM (1.10 LTS, 1.11, 1.12): MWEs in `mwe/`. Not filed;
  per the standing decision, upstream issues need the user's go-ahead.
- Alternatively suggest to Enzyme(.jl) to not tag non-pointer loads with `!invariant.group`, or to strip the tag
  before `post_optimize!`, as the Mono project did (validated here with `ENZYME_PIPE=stripig`).
- Once a fixed LLVM is available, the by-name `@noinline` split can stay (it costs one call per parameterization per
  time step) or be reverted.

## Follow-up 2026-11-17 (afternoon): performance of the workaround, Enzyme.jl MWE

### Performance

User question: does the `@noinline` per-parameterization split slow down the CPU model? Measured with the
`benchmark-regression` harness: `regression.jl benchmark` (`manual_benchmarking.jl --debug --length 2`, Julia 1.12.2),
on two frozen trees, base b06af911a and base + the uncommitted change to `tendencies.jl` and
`SpeedyWeatherEnzymeExt.jl`. They ran interleaved A-B-A-B in one SLURM job (40 of 128 cores reserved; an exclusive node
is impossible under the `priority` QoS cap of 64 CPUs/user). Paste-ready summary in
`SpeedyWeather/benchmark/regression/results/comment.md` (gitignored).

| T | L | base SYPD | candidate SYPD | change |
|---|---|---|---|---|
| 32 | 8 | 754 | 764 | +1% |
| 43 | 8 | 323 | 326 | +1% |
| 64 | 8 | 95 | 96 | +1% |
| 86 | 8 | 38 | 39 | +2% |
| 128 | 8 | 10 | 10 | +1% |
| 32 | 16 | 434 | 434 | +0% |
| 128 | 16 | 5.4 | 5.4 | +0% |
| geomean | | | | **+1%** (best of 2 each) |

The column physics alone (`column_parameterizations_cpu!`, BenchmarkTools median, two rounds each): T32 L8
4.54/4.55 → 4.47/4.48 ms (−1.7%), T64 L8 15.68/15.72 → 15.45/15.47 ms (−1.5%), T128 L16 81.7/82.8 → 81.2/82.8 ms
(−0.3%); 0 bytes allocated in both, so `Val(keys(get_parameterizations(model)))` is resolved statically.
**No regression**: the change is neutral to marginally faster.

Practicalities: compute nodes have no internet, so the benchmark envs created by `regression.jl` must be
resolved/instantiated on the login node first (`Pkg.resolve` fails offline for packages not yet in the depot). Two
trees with the same `HEAD` need separate `--workdir`s, as the env directory is keyed by the commit.

### Enzyme.jl MWE

- **`autodiff` reproducer with SpeedyWeather** (`mwe/mwe_speedyweather.jl`): the default PrimitiveWet T6 L1 model, no
  spin-up, `autodiff(Reverse, parameterization_tendencies!, Duplicated(vars), Duplicated(model))` → aborts
  ("Broken module found") on unmodified b06af911a.
- **Minimal failing SpeedyWeather configuration** (ablation via `repro_param_ad.jl` knobs, unmodified main, Enzyme 0.13.220):
  column parameterizations `boundary_layer + albedo + radiation` only, no shortwave; boundary layer with
  `SurfaceCondition`, no neutral wind, `ConstantSurfaceRoughness`, `BulkRichardsonDrag`; `OceanLandAlbedo(ocean =
  GlobalConstantAlbedo, land = LearnedLandAlbedo)`; longwave as default. Needed: the boundary layer (dropping it
  passes), `BulkRichardsonDrag`, `LearnedLandAlbedo`, longwave, **`Duplicated(model)`** (`Const(model)` passes) and
  the entry point `parameterization_tendencies!` (differentiating only `_column_parameterizations_cpu!` passes; the
  non-inlined `reset_variables!` that writes into `vars` before the column loop changes Enzyme's caching decisions).
- `llvm-reduce` keeping debug info (`ir/reduced_dbg.ll`, 314 lines): the loads involved are `model.geometry.nlayers`
  feeding a `SubArray` `reindex` in `surface_condition!` (rematerialized in the reverse pass with the same
  `!invariant.group`) and `model.geometry.npoints`, the bound of the longwave grid-point loop. SpeedyWeather's 3D grid
  variables and tendencies are `Field`s of `SubArray`s, layer views (`UnitRange`) into shared parent arrays.
- **No dependency-free `autodiff` reproducer found.** `mwe/mwe_standalone.jl` copies the minimal configuration's code
  (same loops, scheme code, `Duplicated` model and variables). It has switches for SpeedyWeather-like bundles + a
  `@generated` loop over a `NamedTuple` (`MWE_STRUCT`), `Field` wrappers (`MWE_FIELD`), packed `UnitRange` layer views
  (`MWE_PACKED`), the entry point with a non-inlined reset (`MWE_ENTRY`) and loop order/repetition (`MWE_ORDER`).
  ~40 combinations all compile fine, although the post-AD module has ~100 `i64` and ~600 `float` loads tagged
  `!invariant.group` and a 13–16k-line reverse function (real case: 448/1 677 and 46k lines).
- The reduced IR files do not break when run through the whole `Enzyme.Compiler.post_optimize!`: the earlier passes
  simplify them. They were reduced against GVN alone. The 93 MB `ir/pre_post_optimize.stripped.ll` does break through
  `post_optimize!` (`mwe_post_optimize_ir.jl`).
- Draft for an Enzyme.jl issue (not filed): `mwe/ENZYME_ISSUE_DRAFT.md`, with the SpeedyWeather reproducer, the
  65-line LLVM reproducer, the diagnosis and the stripping workaround.

### Reversible workaround (`SPLIT_COLUMN_LOOPS`)

User request: make the workaround easy to reverse once the LLVM fix reaches Julia. Implemented as a constant-folded
switch with both code paths kept (see Summary of changes, item 1). Validation (Julia 1.10.13, Enzyme 0.13.220):

| Test | Result |
|---|---|
| `test/parameterizations/column_loops.jl` + `all_parametrizations.jl`, `--check-bounds=yes` | 14/14 (paths bit-identical, Wet and Dry), 61/61 |
| reproducer `repro_param_ad.jl`, `SPLIT_COLUMN_LOOPS = true` (repo) | pass (65.9 s) |
| same, `SPLIT_COLUMN_LOOPS = false` (temporary worktree, otherwise identical) | **abort, "Broken module found"**: the switch restores the original behaviour, so the revert check works |
| full `differentiability/runtests.jl` with the final code (job `sw_full`) | running |
