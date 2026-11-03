# Reactant GPU: XLA fails to compile the PrimitiveWet time loop (concatenate fusion)

> Status: **in progress**. Cause: an XLA regression in Reactant 0.2.303 (openxla/xla#50225), reported as
> EnzymeAD/Reactant.jl#3681; Reactant capped at 0.2.302 (the XLA-flag workaround costs up to half the speed).
> Open: CI decisions (4b), lifting the cap once #3681 is fixed.

Date of initial draft: 2026-11-03

Base revision: 86c1669d3a3a78b561cd98aae1cb123ea3bd35c4 (`mg/fix-reactant`)

## Originating prompt

> ok i need to understand the time loop compile error. make a plan for the investigation

> I also need to understand why this didn't show up in either the reactant correctness tests nor the GPU tests

## Revision log

> okay execute the plan

Executing. Results are collected in [Results](#results) below as they come in.

## Results

Scripts: `…/snow-init-2026-11-03/xla_concat/` (`compile_check.jl <cpu|gpu> <configs…>`, env var `XLA_DEBUG` for
`xla_debug_options`, `REACTANT_DEBUG=1` for source locations; logs `1b_*.log`, `1c_*.log`, `1d_*.log`, `2a_*.log`).
Several Reactant processes share the A40 with `XLA_REACTANT_GPU_MEM_FRACTION=0.3`.

- **1a (H100, Slurm 2873655): fails too**, same `RET_CHECK`, Reactant 0.2.304. So the paper benchmark's PrimitiveWet
  configurations fail at compile time. Same job: CPU vs Reactant `land_fraction` identical, orography max 1.5 m
  difference everywhere, soil temperature differs at 328 of 6336 values by up to 56.5 K (as on the A40). So the
  soil-temperature difference is not from the land-sea mask (future work).
- **1b (A40, Reactant 0.2.299 + jll 0.0.428)**: env `xla_concat/env_r0299` (copy of `env_fix`, Reactant and
  Reactant_jll pinned, everything else preserved, Reactant precompiled on login02; 0.2.299 is the lowest
  SpeedyWeather's compat allows). *Was still compiling `wet` when paused*, log `1b_gpu_0299.log`.
- **1c (A40, CPU backend, 0.2.304): `wet` compiles** (399 s). So the Reactant Correctness CI (CPU) can't see this
  error. `wet_ci` was still running (log `1c_cpu_0304.log`). The run reported `vorticity finite: false`, most likely
  an artefact of the harness: it called `time_steps!` right after `initialize!(model)` without
  `initialize!(simulation; period)` (which the benchmark and `run!` do, e.g. scaling vor/div by the radius). Harness
  fixed afterwards; 1b/1c/1d still ran the old version, so their "finite" flags mean nothing, only OK/FAIL counts.
- **1d (A40, GPU, 0.2.304)**: `wet` FAIL (same `RET_CHECK`, reproduces). **`dry` OK** (compile 210 s), as on the
  H100 in benchmark 2871010: the trigger is PrimitiveWet-specific (humidity / moist physics). `wet_ci`,
  `wet_dynamics_only`, `wet_physics_only` were still running (log `1d_gpu_0304.log`).
- **2a (H100, Slurm 2873733): dump done**, `wet` failed as expected. `xla_concat/dump_h100/` has
  `…before_optimizations.txt` (656 KB), `…sm_9.0a_gpu_after_optimizations.txt` (2.4 MB, HLO after fusion, the
  input to the emitter that fails), buffer assignment etc. **Not analyzed yet.**
- **2b (started)**: the saved module has no source locations (`loc(unknown)`); Reactant only attaches them in debug
  mode (`Reactant.Ops.with_debug`, now `REACTANT_DEBUG=1`). Its 222 concatenates are mostly i64 index tensors for
  gathers/scatters (e.g. 56× `560x1 ×3 → 560x3`) and per-layer stacking (e.g. 9× `8 × (1x3168) → 8x3168`).
- XLA options available for 3b: `xla_disable_hlo_passes`, `xla_enable_hlo_passes_only`,
  `xla_gpu_experimental_enable_same_shape_multi_output_fusion`; no dedicated concatenate-fusion switch.

### Results after resuming (2026-11-03 afternoon)

- **1b: Reactant 0.2.299 compiles** T32 `wet` (444 s, A40). So the regression is in 0.2.300 to 0.2.304, most likely
  0.2.304's XLA/Enzyme-JAX update, and pinning 0.2.299 is a workaround.
- **1c**: CPU backend compiles `wet` and `wet_ci`. **1d** (GPU, 0.2.304): `wet` ✗, `dry` ✓, `wet_ci` ✗,
  `wet_dynamics_only` ✓, `wet_physics_only` ✗.
- **Paper benchmark 2873641 (H100, 0.2.304) finished**; only **T32** PrimitiveWet configurations fail to compile
  (14 `RET_CHECK`s: default setups #4, resolution #1 = T32 L8 default and #22 = T32 L8 matrix, Float32 #1, all four
  layer counts at T32, reverse mode T32 × 4 step counts). T43 to T512, the T64 grid suite and the reverse T64/T128/T256
  runs all compile. The bug depends on the shapes. It wrote its results into the `mg/fix-reactant` checkout
  (`benchmark_results.json`/`README.md` were identical between the branches).
- **All 21 SPEEDY.f90-like configurations failed in `initialize!`** with scalar indexing at `snow.jl:56`, a gap in
  6553d885: a prescribed soil temperature (`SeasonalLandTemperature`, `LandXYZ`) has no step dimension, so
  `data[:, 1, :]` is 3 indices on a 2D array; Base reshapes (`ReshapedArray` of a `ConcretePJRTArray`), then indexes
  element-wise. Fix on branch **`mg/fix-reactant-prescribed-snow-init`** (worktree, from `origin/main` 396ac47f,
  uncommitted): `data[:, 1]` if there are no steps; the CPU test "Initial snow depth …" now runs for a dynamic and a
  prescribed soil temperature.
- **2a: the failing fusion** (H100 dump, `xla_concat/fused_concatenate.7.hlo`): `%loop_concatenate_fusion.2 =
  (c64[2,33,560,1], f32[2,33,560]) fusion(…), calls=%fused_concatenate.7`, roots
  `complex(bitcast(concatenate.580) → [2,33,560,1], bitcast(param) → [2,33,560,1])` and `concatenate.580` itself.
  The two roots' epilogues differ (a rank-adding bitcast + complex vs identity), which is what the `RET_CHECK` rejects.
  `concatenate.580` stacks `[2,8,560]` blocks along dim 1 to 33 = 4 × 8 layers + 1: vorticity, divergence,
  temperature, humidity and surface pressure (PrimitiveWet only, hence dry ✓), 560 = spectral coefficients at T32,
  2 = probably the two leapfrog steps. Our code doesn't stack the prognostic variables, so a rewrite (Enzyme-JAX) or
  XLA built it. 51 multi-output fusions in the module contain a concatenate; only this one mixes root ranks.
- Reactant.jl's own pass configuration (`src/compiler/OptimizationPasses.jl`) differs between 0.2.299 and 0.2.304
  only in sharding/communication passes; the regression is in the binary Reactant_jll (Enzyme-JAX patterns or XLA).

### MWE, bisect and cause (2026-11-03 evening)

> the initialize fixes should be made here in `mg/fix-reactant`. Then I need a MWE of the fused kernel error. Set up
> the MWE (without Speedy), then determine which Reactant version brakes it. Prepare a GitHub issue text as well
> (don't post it, let me do that)

- Prescribed-snow fix moved to `mg/fix-reactant` (main checkout, uncommitted: `snow.jl`, `test/parameterizations/land.jl`).
  CPU `land.jl` tests 24/24 with `--check-bounds=yes`; SPEEDY.f90-like config at T43 on Reactant GPU: `initialize!`
  OK, time loop compiles (358 s), vorticity finite.
- **Workaround (3b): `xla_disable_hlo_passes = ["multi_output_fusion"]`**: T32 PrimitiveWet compiles (370 s) and stays
  finite. Not benchmarked yet.
- **The saved StableHLO module reproduces standalone** (`xla_concat/mwe/compile_module.jl`, only Reactant's XLA
  client): FAIL on 0.2.304 and 0.2.303, compiles on 0.2.302 and 0.2.299. Same input → **XLA regression**.
- **Bisect (MWE and module):** 0.2.302 / jll 0.0.435 OK, **0.2.303 / jll 0.0.436 FAIL**, 0.2.304 / 0.0.438 FAIL.
  XLA `058a1f6a` → `d5229a92` (Reactant.jl `6ad850aa` → `98b0a4eb`, Enzyme-JAX `75e50d4f` → `db6c53ad`, JAX
  `81ac7563` → `3e6335a8`). Culprit candidate in that range: openxla/xla#50225 (`e6914da6`, 2026-10-26) "Support
  multi-output fusion in the concatenate emitter", which added the `RET_CHECK` and `UseConcatenateFusion`'s multi-root rule.
- `loop_concatenate_fusion.2` (2a) is *not* the failing fusion (its root ranks differ, so the new rule forbids it;
  0.2.299 builds and compiles it too). Applying the new rule to all fusions of the dump leaves one with differing
  epilogues: `input_concatenate_fusion.10` (`fused_concatenate.17`), roots `select(divide(concatenate.588, …))` and
  `concatenate.591`, both `f32[1,2,3168]`.
- **MWE** (`xla_concat/mwe/mwe_final.jl`, copy `docs/dev/2026-11/reactant-issue-concat-fusion-ret-check-mwe.jl`):
  `f(a, b, m) = (cat(a, b; dims = 2) ./ m, cat(2a, b .+ 1; dims = 2))` with `(100, 1, 1)` inputs. Compiles: single
  output, both outputs with the same epilogue, no size-1 dimension, scalar epilogue. Disabling `multi_output_fusion`
  fixes it too. Guess (not verified in XLA's source): the two epilogue indexing maps differ only in form for the
  size-1 dimension.
- **Issue draft**: `docs/dev/2026-11/reactant-issue-draft-concat-fusion-ret-check.md`, posted by the user as
  [EnzymeAD/Reactant.jl#3681](https://github.com/EnzymeAD/Reactant.jl/issues/3681).

### Speed of the workaround and decision

> measure the speed impact

> add the workaround for this issue to the `mg/fix-reactant` as well

> okay no, then we can't use this workaround. Pin the reactant version and reference the github issue
> https://github.com/EnzymeAD/Reactant.jl/issues/3681 in the comment

H100, `xla_concat/speed/speed_check.jl` (timed as the benchmark, FFT transform, `2 × n_timesteps` steps, median of 5
runs, all finite), Slurm 2884389 (0.2.304) and 2884390 (0.2.302), PrimitiveWet L8, SYPD:

| Truncation | 0.2.304 default | 0.2.304 no `multi_output_fusion` | 0.2.302 default | 0.2.302 no `multi_output_fusion` |
|---|---|---|---|---|
| T32 | RET_CHECK | 15344 | **29733** | 15540 |
| T43 | 17182 | 8758 (−49%) | 16980 | 8701 |
| T64 | 9543 | 5201 (−45%) | 9792 | 5225 |
| T128 | 2094 | 1435 (−31%) | | |
| T256 | 369 | 302 (−18%) | | |

(blank: not yet in when this was written; see the job logs `xla_concat/speed/speed_*.log`.) Disabling the pass costs
up to half the speed (it halves compile time too), on both versions; 0.2.302 and 0.2.304 are equally fast where both
compile. **Decision: no workaround, Reactant capped at 0.2.302.** The workaround (committed in 75eaf319) is removed
again (uncommitted on `mg/fix-reactant`); `Reactant = "0.2.299 - 0.2.302"` with a comment referencing #3681 in the
`[compat]` of LowerTriangularArrays, RingGrids, SpeedyTransforms, SpeedyWeatherInternals, SpeedyWeather and
`SpeedyWeather/test/reactant` (all already `+DEV`/`-DEV`). A scratch env with the capped SpeedyWeather resolves to
Reactant 0.2.302 / Reactant_jll 0.0.435. The benchmark env's (untracked) Manifest still has 0.2.304 and needs a
`Pkg.resolve()`.

## Handoff 2026-11-03

**Where things are.** `X=/p/tmp/maxgelbr/code/Speedy2/claude_reactant_gpu_2026-10-13/snow-init-2026-11-03/xla_concat`.
Three runs were left running detached (`setsid`) on login02's A40 and survive the session; each writes
`$X/<name>.log` and, when done, `$X/<name>.exit`:

| Run | Step | What it decides |
|---|---|---|
| `1b_gpu_0299` | Reactant 0.2.299 on GPU, `wet` | regression in 0.2.300 to 0.2.304? pinning 0.2.299 rescues the paper benchmark? |
| `1c_cpu_0304` | CPU backend, `wet_ci` (`wet` done: compiles) | – |
| `1d_gpu_0304` | GPU 0.2.304: `wet_ci`, `wet_dynamics_only`, `wet_physics_only` (done: `wet` fails, `dry` compiles) | which part of the model triggers it |

Read the results with
`grep -E "^Reactant|^[a-z_]+ +(OK|FAIL)" $X/1b_gpu_0299.log $X/1c_cpu_0304.log $X/1d_gpu_0304.log; cat $X/*.exit`.
Only OK/FAIL counts in these three (harness fixed after they started).

The paper benchmark (Slurm 2873641, `mg/benchmark-reverse` code, Reactant 0.2.304) was still running; its
PrimitiveWet configurations will fail with this error. The main checkout was switched to `mg/fix-reactant` while it
runs: switch back to `mg/benchmark-reverse` before it finishes, as it writes `benchmark_results.json`/`README.md`
into the checked-out branch at the end.

**Next steps, in order.**

1. Collect 1b/1c/1d results (above).
   - If 0.2.299 compiles: pin Reactant 0.2.299 + Reactant_jll 0.0.428 in the benchmark env (copy of the pins in
     `$X/env_r0299`) and rerun the paper benchmark (step 3a). Optionally bisect 0.2.300 to 0.2.303 (needs downloads).
   - If 0.2.299 fails too: the regression is older than 0.2.299 or caused by our code since 2026-10-18 → step 2 is the
     main path.
2. Analyze the dump (step 2a/2b): in `$X/dump_h100/*sm_9.0a_gpu_after_optimizations.txt` find fusions whose
   computation contains `concatenate` and whose root is a `tuple` (multi-output), e.g.
   `grep -n "fusion(.*kind=k" …after_optimizations.txt | grep -i concat` and the `fused_computation` bodies they
   call. The failing one has roots with different shapes/indexing after the concatenate. For its Julia source, rerun
   with `REACTANT_DEBUG=1` (stack traces as locations, tracing is slow) and the same `XLA_DEBUG` dump options, on the
   H100 via `$X/submit_dump.sh` (add `export REACTANT_DEBUG=1`).
3. Try workarounds (step 3b) with the fixed harness, which now also checks finiteness properly:
   `XLA_DEBUG='(; xla_gpu_experimental_enable_same_shape_multi_output_fusion = false)'`, and
   `XLA_DEBUG='(; xla_disable_hlo_passes = ["<pass>"])'` with fusion pass names taken from the dump (e.g. the
   multi-output fusion pass). Run as
   `XLA_REACTANT_GPU_MEM_FRACTION=0.3 julia --project=$X/../env_fix $X/compile_check.jl gpu wet`.
4. Then steps 2c (reproducer), 4a (issue) and 4b (CI decisions) as planned.

**Other open items.** CPU vs Reactant initial soil temperature (up to 56.5 K at 328 values, not the land-sea mask,
which is identical): pre-#1284, look at `LandBucketTemperature`/`SeasonalLandTemperature` initialization on Reactant.

## Problem description

Found while verifying the initialization fixes of
[2026-11-03-reactant-gpu-snow-init-scalar-indexing.md](2026-11-03-reactant-gpu-snow-init-scalar-indexing.md)
(6553d885). With those fixes, `initialize!` of a T32 L8 `PrimitiveWetModel` completes on Reactant GPU (A40, login02),
but compiling the time loop as the benchmark does,

```julia
@compile compile_options = SpeedyWeather.reactant_compile_options() SpeedyWeather.time_steps!(simulation, 10)
```

fails inside XLA's GPU backend (`xla::gpu::ConcatenateFusion` → `ConcatenateFusionKernelEmitter::EmitEntryFunction`,
reached from `ThunkEmitter::EmitWhile`, i.e. inside the time loop's while body):

```
INTERNAL: RET_CHECK failure (external/xla/xla/codegen/emitters/concatenate_kernel_emitter.cc:200)
ComputeEpilogueInputToOutputIndexing(fusion_spec_.fusion_hero(i), fusion_spec_.fusion_root(i), &mlir_context_)
== epilogue_indexing UseConcatenateFusion should have ensured all roots share the same epilogue indexing
```

This is an internal consistency check in XLA: the fusion pass built a (multi-output) concatenate fusion whose roots
don't share one epilogue indexing, which the emitter requires. So it is an XLA bug, triggered by something in our
program. If it also fails on the H100, every PrimitiveWet configuration of the paper benchmark fails again, now at
compile time.

## Background

What is known (2026-11-03):

- **Versions.** Last good reactant-gpu benchmark (H100, committed README/JSON) used **Reactant 0.2.295**
  (Reactant_jll 0.0.424), including the whole PrimitiveWet suite. PrimitiveWet also compiled and ran on the A40 on
  2026-10-13 (0.2.29x). The failure is on **Reactant 0.2.304** (benchmark env, Reactant_jll 0.0.438).
- **0.2.304** was released 2026-11-01 18:45 with "update XLA/PJRT deps" (Reactant.jl#3633) and an Enzyme-JAX
  update. Prime suspect, since the failure is in XLA.
- **Scope so far.** In the failed benchmark 2871010 (H100, 0.2.304) Barotropic, ShallowWater (all resolutions) and
  PrimitiveDry T32L8 (182 s) compiled; PrimitiveWet never got to compiling (failed in `initialize!`). The init fixes
  can't be the trigger: the compiled time loop doesn't trace `initialize!` code.
- **Artefacts.** The module Reactant handed to XLA (after Reactant/Enzyme-JAX passes) is saved at
  `/p/tmp/maxgelbr/code/Speedy2/claude_reactant_gpu_2026-10-13/snow-init-2026-11-03/xla_concat/a40_timeloop_post_xla_compile.mlir`
  (927 KB, 222 `stablehlo.concatenate`). Reproducer: `diagnose_snow_diff.jl` in the same folder (env `env_fix`,
  A40, about 10 min).
- **Locally available versions.** Reactant 0.2.292 to 0.2.299 and 0.2.304, Reactant_jll 0.0.418 to 0.0.428 and
  0.0.438. 0.2.300 to 0.2.303 would have to be downloaded.
- **Relevant XLA debug options** in the jll (`strings libReactantExtra.so`): `xla_gpu_experimental_enable_same_shape_multi_output_fusion`,
  `xla_gpu_unsupported_enable_triton_multi_output_fusion`. They can be passed via `Reactant.CompileOptions(;
  xla_debug_options = (; …))`, the mechanism we already use for `xla_gpu_enable_triton_gemm`.
- **Techniques that worked in October** (memory/plan docs): `Reactant.MLIR.IR.DUMP_MLIR_ALWAYS`, rerunning a dumped
  pipeline on the saved module in seconds, module reduction with placeholder `custom_call`s
  (`claude_reactant_gpu_2026-10-13/dft_impl_2026-10-18/reduce_module.jl`).

### Why CI didn't catch this, or the three init bugs

| | Reactant Correctness (GitHub Actions) | GPU tests (Buildkite, CUDA) |
|---|---|---|
| Reactant backend | **CPU** | GPU |
| PrimitiveWet on Reactant | yes: init + time steps, but with `convection = nothing` and constant transmissivity (`test/reactant/setup.jl`) | **commented out** in `test/GPU/reactant.jl` (`#= … =#`, "the compiled_call just takes too long in the CI with consumer GPUs"); only `BarotropicModel` runs |
| Init bugs (view broadcast, traced branches) | invisible: eager broadcasts on Reactant CPU take a host path without tracing | would have failed in `initialize!`, but the test is disabled |
| XLA compile error | probably invisible: XLA's CPU backend has its own fusion pipeline and emitters (to verify, step 1c) | would show if it is GPU-generic and Reactant ≥ 0.2.304, but the test is disabled |
| Reactant version | compat `"0.2.299"` = [0.2.299, 0.3), so the newest 0.2.x at run time | same (CUDA test env has no Reactant bound) |

So the GPU-only bugs had no test that could see them. The time at which a release lands also decides which version
CI uses: 0.2.304 was released on the evening #1284 was merged, so #1284's PR runs most likely used ≤ 0.2.303
(to verify from the CI logs). The benchmark environment isn't pinned either and picked up 0.2.304 untested.

Why the PrimitiveWet GPU test was disabled: the Buildkite CUDA job has a 40 min limit, the suite took 39 min on the
CI's A40 in October with `reactant.jl` about 18 min (XLA compile of the loop program, host-bound).

## Investigation plan

Each step says what it decides. Steps 1a to 1d are cheap and independent, and decide whether the rest is needed.

### Step 1: scope (≤ 30 min each)

- **1a. H100 vs A40.** Slurm job 2873655 (`h100_check.jl`: T32L8 PrimitiveWet, compile + 10 steps, plus CPU vs
  Reactant land-sea mask/orography/soil temperature) and the running paper benchmark 2873641 (its 4th config of the
  first suite is PrimitiveWet). *Decides urgency:* if the H100 fails too, the paper benchmark is blocked → go to the
  workaround in 3a right away.
- **1b. Reactant version.** Same script on the A40 with 0.2.299 (jll 0.0.428) in a copy of `env_fix`. If it compiles,
  bisect 0.2.300 to 0.2.303. Optionally hold Reactant.jl fixed and swap only the jll (if compat allows) to separate
  the Julia side (pass list, lowering) from the jll (XLA, Enzyme-JAX). *Decides* regression vs. old bug, and gives the
  version range for the issue. Precompile each env on login02 first (shared depot).
- **1c. CPU backend.** Same compile with `Reactant.set_default_backend("cpu")`. *Decides* whether the Reactant
  Correctness CI could ever see it.
- **1d. Model scope.** On the failing version: PrimitiveDry T32L8 on the A40 (compiled on the H100); PrimitiveWet with
  the correctness-test config; `dynamics_only = true`; physics only; then set parameterizations to `nothing` one at
  a time (convection, large-scale condensation, land, ocean, sea ice, radiation, surface fluxes). *Decides* which
  component's code gets fused into the failing concatenate.

### Step 2: localize in XLA (if 1 shows a deterministic failure)

- **2a. Dump.** Compile with `XLA_FLAGS="--xla_dump_to=<dir> --xla_dump_hlo_as_text"` (or via `xla_debug_options`) and
  find the failing fusion: the last fusion emitted before the error, with a `concatenate` hero and several roots.
  Keep its HLO text.
- **2b. Map to source.** Follow the operands of that fusion back to `stablehlo` ops in the saved module and their
  `loc(...)` (Julia file:line). *Decides* whether the concatenate comes from our code (e.g. a `vcat`/`cat`, or a
  slice-assign pattern in a kernel or transform) or from an Enzyme-JAX rewrite (`concat_*`, `*_pad` patterns).
- **2c. Minimal reproducer.** Extract the fusion's computation as a standalone StableHLO module and compile it with
  Reactant on GPU, then reduce (October's `reduce_module.jl`). Goal: a module of a few ops that fails on 0.2.304 and
  compiles on the last good version.

### Step 3: workarounds (evaluate the cheapest first, measure SYPD vs the committed results)

- **3a. Pin Reactant** in the benchmark env to the last good version (0.2.295 or whatever 1b finds). Immediate
  unblock for the paper benchmark; needs the init fixes to work on that version too (should, they're plain Julia).
- **3b. XLA flag** in `reactant_compile_options()` that keeps XLA from building this fusion (candidates: the
  multi-output fusion options above, others from the jll's `DebugOptions`). Check the perf impact.
- **3c. Exclude an Enzyme-JAX pattern** via `excluded_passes` if 2b shows a rewrite creates the concatenate.
- **3d. Change our code** if one obvious spot produces the pattern and the change is cheap.

### Step 4: upstream and CI

- **4a. Issue** at Reactant.jl (and openxla/xla if it reproduces as plain HLO) with the reproducer from 2c and the
  version range from 1b.
- **4b. CI gaps** (needs the user's decisions):
  1. A cheap Reactant **GPU** PrimitiveWet test: `initialize!` only costs no time-loop compile and would have caught
     all three init bugs; a short compile at T32 for the time loop, if the CI budget allows.
  2. Reactant version: the compat floor alone lets every new 0.2.x into CI and benchmarks unannounced. Options: an
     upper bound bumped deliberately (as with the 0.2.296 cap in October), or a scheduled CI run against the newest
     release.
  3. Pin the benchmark env's Reactant (Manifest is gitignored, so the version is whatever was resolved last).

## Testing and verification

- The failure reproduces deterministically with `diagnose_snow_diff.jl` / `h100_check.jl` before any change.
- A workaround or fix counts as verified when: the T32L8 PrimitiveWet time loop compiles and 10 steps stay finite on
  the A40 and the H100; prognostics stay close to CPU (≈1e-4, as in October); a short reactant-gpu benchmark subset
  (benchmark100, PrimitiveWet row) gives SYPD within ±20% of the committed 27700.
- For 3b/3c additionally: the CPU Reactant correctness suite (`SpeedyWeather/test/reactant`, no `--check-bounds`).

## Known limitations

- Bisecting inside XLA (between two jll builds) would need building Reactant_jll; the plan stops at release
  granularity plus a minimal reproducer.
- login02's A40 is shared: compile times there are not representative, use it for pass/fail only. Only one Reactant
  process fits on it at a time (XLA preallocates 75% of the memory). H100 via Slurm `gpushort`: 2 jobs per user, the
  paper benchmark occupies one.
- Long runs must be detached with `setsid` (session teardowns killed `nohup` runs twice on 2026-11-03).

## Future work

- Unrelated finding of the same verification: the initial soil temperature differs between CPU and Reactant at
  164 of 3168 points by up to 56.5 K (CPU 285 K = ocean value, Reactant land values). Code older than #1284; step 1a
  compares the land-sea mask and orography to see where it comes from.
