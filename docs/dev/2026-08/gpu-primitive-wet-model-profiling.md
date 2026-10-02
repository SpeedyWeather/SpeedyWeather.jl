# GPU profiling of the full `PrimitiveWetModel` time step

> Status: **in progress**. Phases A–C measured on both an H100 and an A40: the GPU is idle
> 70–78 % of every time step, `dynamics` issues 97 % of the launches, and an unplanned `K=2`
> FFT batch size accounts for ~60 % of them. Phase D is blocked (`ncu` counters are
> admin-restricted). Paused before confirming the `K=2` fix — see *Next steps*.

Date of initial draft: 2026-08-07

Base revision: `f417c407` (`mg/optimize-legendre-gpu`)

## Originating prompt

> Make a plan how to profile the GPU performance of the complete `PrimitiveWetModel` to
> identify performance bottlenecks on GPU.

## Revision log

- **2026-08-07, initial draft.** Five-phase plan: measurement harness with correctness guards,
  NVTX phase annotation under `nsys`, trace analysis (GPU-busy fraction, kernel table, sync
  hunting), `ncu` deep dives on the top kernels, synthesis into a ranked bottleneck list.
- **2026-08-07, "Detail Phase B - D more for the executing agent. For Phase A we will restrict
  us to 8 layers for now."** Phase A restricted to `nlayers = 8` (the `L = 16` comparison is
  deferred; layer scaling is the subject of the layer-count investigation). Phases B–D expanded
  into step-by-step instructions: driver skeleton with exact phase ranges, `nsys`/`ncu` command
  lines and SLURM scripts, analysis recipes with decision rules, and known failure modes.

- **2026-08-07, execution started.** Scripts for Phases A–D written under
  `SpeedyWeather/test/GPU/modelbench/` and smoke-tested on CPU before spending GPU job time:
  - `MODELBENCH_CPU=1` added to the Phase B driver so the B3 equivalence check and a driver
    smoke test can run on a login node; profiling itself is still GPU-only.
  - **B3 passes bitwise on CPU** (T31, 10 steps, all of vorticity/divergence/temperature/
    humidity identical) — the driver does mirror `time_step!`.
  - Bug caught by that check: the state snapshot listed surface pressure as `:pres`, but the
    prognostic variable is `:pressure`, so it was being silently skipped. Both scripts now
    assert every expected name exists instead of `hasproperty`-skipping, so an upstream rename
    fails loudly rather than quietly shrinking coverage.
  - Phase C report names verified against nsys 2025.3.2: `nvtx_kern_sum` gives the per-phase
    kernel breakdown directly. `nvtx_pushpop_sum` turned out to be **empty** for these traces;
    `nvtx_gpu_proj_sum` carries both the projected GPU time *and* the host range span, so it is
    the only NVTX report actually needed.
  - **Trap in `nvtx_gpu_proj_sum`, worth remembering:** `Total Proj Time` is the range's *span
    on the GPU timeline*, which includes the idle gaps between kernels — it is **not** GPU busy
    time, and `Total GPU Ops` is a count, not a duration. Deriving "GPU busy" from those
    columns overstates it wildly (99 % instead of the true ~20 %). True per-phase busy time
    comes from summing kernel durations in `nvtx_kern_sum`. `phase_table.jl` now reports both
    (`busy ms` and `span ms`) so the distinction cannot be lost again.
- **2026-08-07, Phases A–C measured; paused.** Findings recorded below. The `ncu` preflight
  failed with `ERR_NVGPUCTRPERM`, blocking Phase D. Work paused at the user's request before
  running the `K=2` confirmation experiment; `verify_batch_k2.jl` is written and ready.
  - *One conclusion drawn and then retracted within the session:* the ~200× slowdown under
    nsys was first attributed to `--trace=osrt`, and Phase B was re-run on the A40 without it.
    The re-run showed the same overhead, disproving that. Corrected in the methodological note
    below; the mistaken claim had also propagated into the plan's next steps and was fixed
    there too. The A40 traces remain useful as a cross-GPU check, just not for host spans.
  - *Two GPUs, discovered mid-run:* SLURM nodes have H100s, the login node has an A40. All
    principal findings were reproduced on both.

## Problem description

The GPU Legendre transform work (`docs/dev/2026-08/legendre-gpu-optimization.md`) made the
Legendre kernels 5–37× faster, so they are no longer the bottleneck of `transform!`. What is
*not* known is where the full `PrimitiveWetModel` time step now spends its time on GPU, and
whether that changes with resolution. Candidate suspects, none yet quantified:

- **Kernel-launch overhead / GPU starvation**: a time step launches many small kernels
  (parameterizations are column-wise, dynamics has ~36 kernel files); at low truncation the GPU
  may be idle most of the time waiting for the host to launch work.
- **Parameterizations vs dynamics vs transforms**: no per-phase breakdown exists for GPU at all.
- **Hidden per-step synchronization**: `nan_detection!` (`SpeedyWeather/src/output/feedback.jl:118`)
  does a scalar `@allowscalar vor[2, end]` read every step when `feedback.debug == true`
  (the default) — a forced device-to-host copy + sync that also breaks CPU/GPU pipelining.
  There may be others (implicit `Mem.copy!`s, scalar reductions, output hooks).
- **Individual slow kernels**: kernels far from memory-bandwidth roofline, like the Legendre
  kernels were before optimization.

The goal of this work is a measured, per-phase and per-kernel account of a steady-state GPU
time step at several truncations, ending in a ranked list of bottlenecks with an estimate of
what fixing each would buy — i.e. the input for the *next* optimization plan, not the
optimizations themselves.

## Background

### What a time step looks like

`time_step!(vars, time_stepping, model::PrimitiveEquation)`
(`SpeedyWeather/src/time_stepping/time_integration.jl:68`) is a short, stable sequence — the
natural phase boundaries for profiling:

1. `reset_tendencies!`
2. `parameterization_tendencies!` (+ `greenhouse_gases_time_step!`, `ocean_timestep!`,
   `sea_ice_timestep!`, `land_timestep!`)
3. `dynamics_tendencies!` (includes the grid→spectral transforms of tendencies)
4. `horizontal_diffusion!` + `implicit_correction!` (via `diffusion_and_implicit!`)
5. `update_prognostic!` (leapfrog + tracer/ocean/land stepping, filters)
6. `transform!(vars, model)` (fused/batched spectral→grid: Fourier via CUDA-graph cache +
   Legendre)
7. per-step host-side wrapper: `progress!` (incl. `nan_detection!`), `output!`, `callback!`

### Existing measurement infrastructure and prior findings

- `SpeedyWeather/benchmark/manual_benchmarking.jl` — whole-model SYPD table in
  `benchmark/README.md`, GPU columns included. **Caveat**: the layer-count investigation
  (`docs/dev/2026-08/layer-count-performance-anomaly.md`) showed this harness silently counts
  post-NaN no-op steps as work, times over as little as ~1 s, and includes setup — its numbers
  identify configurations but cannot be trusted for regression math without the fixes proposed
  there. Every profiling run here must check `model.feedback.nans_detected` afterwards.
- `SpeedyWeather/test/GPU/legendrebench/` — the sweep/`ncu` pattern from the Legendre work
  (`bench_legendre.jl`, `ncu_target.jl`); this plan reuses that structure at model scope.
- `SpeedyWeather/test/GPU/graphbench/` — CUDA-graph cache analysis; the Fourier graph cache is
  bounded and keyed by device pointer (fix in `SpeedyTransformsCUDAExt.jl`).
- The benchmark env (`SpeedyWeather/benchmark/Project.toml`) has `CUDA`, `BenchmarkTools` and a
  path-dev'd `SpeedyWeather`.

### Tooling on this cluster

- GPU nodes via SLURM: `--partition=gpu --qos=gpushort --account=flai --gres=gpu:1`
  (see `submit_bench.sh`).
- **The GPU model is not fixed across the partition.** The first profiling job landed on an
  **NVIDIA H100 80GB HBM3** (driver 580.126.09), whereas the Legendre optimization work
  (`legendre-gpu-optimization.md`) was carried out and reported against an **A40**. The two are
  very different roofline targets:

  | | A40 | H100 80GB HBM3 |
  |---|---:|---:|
  | peak HBM bandwidth | ~0.70 TB/s | ~3.35 TB/s |
  | L2 cache | 6 MB | 50 MB |

  **Both are reachable, and cheaply:** SLURM GPU nodes carry H100s, while the *login node*
  (`login02`) has an A40 attached that can be used directly with no queue wait. Since the GPU
  partition is often saturated (all 4 GPUs allocated on all 12 nodes when this was written),
  the practical workflow is: iterate on the A40 locally, and use SLURM for the H100 numbers.
  Profiling both is worth the small extra cost — the bottleneck *ranking* is expected to
  differ, and the A40 is the hardware the Legendre work was tuned on.

  Consequences: (a) "% of peak bandwidth" figures may not be compared across the two — the
  Legendre kernels' reported 37–41 % of A40 peak would be a far smaller fraction of H100 peak;
  (b) the L2-residency argument in the Legendre plan's Step 0 (an 8.4 MB coefficient array
  overflowing a 6 MB L2) does **not** hold on H100's 50 MB L2, so that regression may simply
  not reproduce here; (c) launch overhead matters relatively *more* on H100, since the same
  fixed per-launch cost buys much faster kernels. **Every result must record which GPU it ran
  on** (`nvidia-smi` output is captured at the top of each job log).
- `nsys` 2025.3.2 at `/usr/local/bin/nsys`; `ncu` at `/usr/local/cuda-13.0/bin/ncu`.
  Both are host tools wrapping the `julia` process — this does **not** conflict with the rule
  "don't load system CUDA modules for Julia" (no module load needed; CUDA.jl keeps using its
  own artifacts).
- `CUDA.@profile external=true` brackets the capture window via `cudaProfilerStart/Stop`, so
  `nsys profile --capture-range=cudaProfilerApi` records only the steady-state steps, excluding
  compilation, warmup and graph capture.
- NVTX ranges (`NVTX.jl`) give named phase bars in the timeline and enable per-phase GPU-time
  attribution in `nsys stats` (`nvtx_gpu_proj_sum` report).

## Summary of changes

All new scripts live in `SpeedyWeather/test/GPU/modelbench/` with their own `Project.toml`
(contents specified in Phase B1), mirroring `legendrebench/`. No model source is modified in
phases A–D (see the NVTX decision below).

### Phase A — measurement harness with correctness guards (`bench_model.jl`)

A trustworthy macro-timing baseline, replacing the flawed harness for this purpose:

- Configurations: `trunc ∈ {31, 63, 127, 255}`, `nlayers = 8` (fixed for now), `Float32`,
  default `PrimitiveWetModel`, `output=false`, no callbacks. T31 probes the
  launch-overhead-bound regime, T255 the compute-bound regime. The per-launch vs per-kernel
  cost separation comes from the truncation sweep plus the Phase C launch counts (the layer
  sweep originally planned for this is deferred).
- Protocol per configuration: `initialize!`; spin up ≥ 32 steps via `run!` (JIT, FFT plans,
  CUDA-graph capture all settle); then time N ≥ 100 steps wallclock with
  `CUDA.synchronize()` only at start/end of the timed window (not per step, to keep pipelining
  realistic); report ms/step and SYPD.
- **Guards** (lessons from the layer-count investigation): after the timed window assert
  `model.feedback.nans_detected == false` — otherwise print `UNSTABLE` instead of a number;
  never include `initialize!`/`run!` setup in the window; record step count actually executed.
- Also record ms/step with `feedback = nothing` vs default feedback, and `dynamics_only = true`
  vs full physics — two cheap A/B switches that already bound the cost of the per-step NaN sync
  and of the entire parameterization suite before any profiler is attached.

### Phase B — phase-annotated timeline (`nsys_target.jl` + `submit_nsys.sh`)

**Goal:** one `.nsys-rep` per truncation in which each captured steady-state time step is
decomposed into named NVTX phases.

**NVTX decision**: annotation lives in the driver script, *not* in model source, so phases
A–D need no source change, no version bump, no new package dependency. If profiling becomes
routine, a follow-up PR can add `NVTX.@annotate` permanently to the ~7 phase functions
(Oceananigans-style; no-op overhead when no profiler is attached) — out of scope here.

**B1 — environment.** Create `SpeedyWeather/test/GPU/modelbench/Project.toml` with `[deps]`
`CUDA`, `NVTX`, `JSON3`, `Printf`, `Statistics`, `SQLite` (for Phase C2), `SpeedyWeather`, and
`[sources] SpeedyWeather = {path = "../../.."}` (mirroring `legendrebench`). Instantiate once:
`julia --project=SpeedyWeather/test/GPU/modelbench -e 'using Pkg; Pkg.instantiate()'`
(first CUDA precompile takes minutes; fine on the login node).

**B2 — driver** `nsys_target.jl <trunc> <nlayers=8> <n_warmup=32> <n_capture=20>`:

- Build the model exactly as in Phase A; `simulation = initialize!(model)`.
- Then `SpeedyWeather.initialize!(simulation; steps = n_warmup + n_capture, output = false)`.
  This is the key trick: `initialize!(simulation; ...)`
  (`SpeedyWeather/src/models/simulation.jl:53`) does the clock setup, the radius scaling
  (`scale_prognostic!`) and the initial `transform!` — after it, plain
  `SpeedyWeather.time_step!(simulation)` calls are exactly `run!`'s loop body, and
  `SpeedyWeather.finalize!(simulation)` at the end unscales. The driver never re-implements
  scaling.
- Warmup: `for _ in 1:n_warmup; SpeedyWeather.time_step!(simulation); end; CUDA.synchronize()`.
  32 steps cover JIT compilation, FFT planning and all Fourier CUDA-graph captures (the cache
  is bounded, so no captures happen later).
- Captured window:

  ```julia
  CUDA.@profile external = true begin
      for _ in 1:n_capture
          instrumented_time_step!(simulation)
      end
      CUDA.synchronize()
  end
  ```

  `cudaProfilerStart/Stop` from `external = true` is what `--capture-range=cudaProfilerApi`
  keys on, so warmup never enters the report.
- `instrumented_time_step!` replicates the two nested functions
  `time_step!(simulation, time_stepping)` (`time_integration.jl:14`) and
  `time_step!(vars, time_stepping, model::PrimitiveEquation)` (`time_integration.jl:68`),
  wrapping each phase in `NVTX.@range`. Skeleton (the executing agent must copy the calls from
  the *current* source at implementation time, not from this sketch, and
  `@assert model.dynamics && !model.dynamics_only` up front so the replicated branch is the one
  the real code takes):

  ```julia
  function instrumented_time_step!(simulation)
      (; variables, model) = simulation
      (; time_stepping) = model
      vars = variables
      NVTX.@range "step" begin
          NVTX.@range "reinitialize"      SpeedyWeather.reinitialize!(model, variables)
          NVTX.@range "reset_tendencies"  SpeedyWeather.reset_tendencies!(vars, time_stepping)
          NVTX.@range "parameterizations" begin
              SpeedyWeather.greenhouse_gases_time_step!(vars, model)
              SpeedyWeather.parameterization_tendencies!(vars, model)
          end
          NVTX.@range "ocean_land" begin
              SpeedyWeather.ocean_timestep!(vars, model)
              SpeedyWeather.sea_ice_timestep!(vars, model)
              SpeedyWeather.land_timestep!(vars, model)
          end
          NVTX.@range "dynamics"           SpeedyWeather.dynamics_tendencies!(vars, model)
          NVTX.@range "diffusion_implicit" SpeedyWeather.diffusion_and_implicit!(vars, model)
          NVTX.@range "leapfrog"           SpeedyWeather.update_prognostic!(vars, model)
          NVTX.@range "transform"          SpeedyWeather.transform!(vars, model)
          NVTX.@range "particle_advection" SpeedyWeather.particle_advection!(vars, model)
          NVTX.@range "clock"              SpeedyWeather.time_step!(vars.prognostic.clock, time_stepping)
          NVTX.@range "feedback_output" begin
              SpeedyWeather.progress!(model.feedback, vars, model)
              SpeedyWeather.output!(model.output, simulation)
              SpeedyWeather.callback!(model.callbacks, vars, model)
          end
      end
      return nothing
  end
  ```

- After the window: `SpeedyWeather.finalize!(simulation)`, then
  `@assert !model.feedback.nans_detected`, and print host-side ms/step of the captured window
  (for the Phase A cross-check).

**B3 — equivalence check** (`--check` flag on the driver, run once, and re-run after any
rebase): build two simulations from identically-configured models, calling `Random.seed!` with
the same seed before each construction (initial conditions may draw random numbers). Run one
with `run!(sim1, steps = k)` (k ≈ 10); run the other through
`initialize!(sim2; steps = k, output = false)`, k × `instrumented_time_step!(sim2)`,
`finalize!(sim2)`. Compare `Array(...)` of the prognostic spectral state (vorticity,
divergence, temperature, humidity, pressure) with `==`. Bitwise equality is expected — the
kernels are deterministic (the forward-Legendre atomics were removed in the Legendre work). If
it fails, the driver has diverged from `time_step!`; fix the driver, never loosen the
comparison.

**B4 — SLURM + nsys invocation.** `submit_nsys.sh` modeled on `submit_bench.sh`
(same partition/qos/account, `--time=02:00:00` suffices):

```bash
module load julia
NSYS=/usr/local/bin/nsys
mkdir -p SpeedyWeather/test/GPU/modelbench/reports
for trunc in 31 63 127 255; do
    $NSYS profile \
        --trace=cuda,nvtx,osrt \
        --capture-range=cudaProfilerApi --capture-range-end=stop \
        --force-overwrite=true \
        -o SpeedyWeather/test/GPU/modelbench/reports/model_T${trunc}_L8 \
        julia --project=SpeedyWeather/test/GPU/modelbench \
        SpeedyWeather/test/GPU/modelbench/nsys_target.jl $trunc 8
done
```

Notes: ~20 captured steps keep reports small; if a report is still unwieldy, drop `osrt` from
`--trace` (loses OS-thread context, keeps everything needed for C1–C4). Julia must NOT be run
with system CUDA modules loaded; `nsys`/`ncu` as host wrappers are fine.

### Phase C — trace analysis (per-phase and per-kernel accounting)

Runs on the login node (no GPU needed once the `.nsys-rep` files exist).

**C1 — canned reports.** For each report:

Implemented as `analyze_nsys.sh`, which loops over the Phase B traces and runs (report names
verified against nsys 2025.3.2 via `nsys stats --help-reports`):

```bash
/usr/local/bin/nsys stats --force-export=true \
    --report nvtx_gpu_proj_sum --report nvtx_pushpop_sum --report nvtx_kern_sum \
    --report cuda_gpu_kern_sum --report cuda_api_sum \
    --report cuda_gpu_mem_time_sum --report cuda_gpu_mem_size_sum \
    --format csv --output <base> <base>.nsys-rep
```

Interpretation:

- `nvtx_gpu_proj_sum` — GPU kernel/memcpy time *projected* onto each NVTX range (via launch
  correlation), plus instance counts → per-phase GPU-busy time and launches. The trustworthy
  per-phase number.
- `nvtx_pushpop_sum` — host-side wall span of each range (`NVTX.@range` is a push/pop range,
  so this is the right report, not `nvtx_sum`).
- `nvtx_kern_sum` — kernels grouped by enclosing NVTX range; gives the per-phase kernel
  breakdown directly, without hand-mapping the global kernel table back to phases.
  (`cuda_gpu_kern_sum:nvtx-name` is an equivalent alternative.)
- Approximate per-phase GPU-idle fraction = 1 − projected-GPU-time / host-span. This is only
  approximate because kernels launched in phase *n* can still be running while the host is
  already in phase *n+1* (an NVTX range is a host-side interval and launches are async); where
  that matters, use C2.

**C2 — sqlite post-processing** (`analyze_nsys.jl`, only if C1 leaves ambiguity).
`nsys stats --force-export=true` already writes a `.sqlite` next to the report. Relevant
tables: `CUPTI_ACTIVITY_KIND_KERNEL` / `_MEMCPY` / `_MEMSET` (device intervals),
`NVTX_EVENTS` (range text, start/end), `CUPTI_ACTIVITY_KIND_RUNTIME` (API calls with
`correlationId` linking a launch/memcpy to its device activity). Compute exactly: merged-union
GPU-busy time per phase and overall, launches per phase per step, and for every
device-to-host copy the API call and NVTX phase it originated from.

**C3 — headline tables + decision rules.** Per configuration:

1. **Phase table**: phase | host ms/step | GPU ms/step | launches/step | idle %. The headline
   result — where the step goes, and whether each phase is kernel-bound (GPU busy) or
   launch/host-bound (GPU idle while the CPU sits in the range).
2. **Kernel table** (`cuda_gpu_kern_sum`): top ~15 kernels by cumulative time with count and
   mean duration. KernelAbstractions kernels appear as `gpu_<kernel_function_name>`; map back
   to source with `grep -rn "function <name>" SpeedyWeather/src SpeedyTransforms/src
   RingGrids/src LowerTriangularArrays/src`.
3. **Sync & transfer hunt** (`cuda_api_sum` + C2): per-step counts of `cuMemcpyDtoH*`,
   `cuStreamSynchronize`/`cuCtxSynchronize`/`cuEventSynchronize`, attributed to phases.
   Expected: exactly one scalar DtoH per step from `nan_detection!` inside `feedback_output`;
   anything else is a finding.
4. **Graph coverage**: `cuGraphLaunch` count per step from `cuda_api_sum`; fraction of
   per-step GPU time spent in graph-launched kernels vs individually-launched ones.

Rules of thumb for classification: a phase is *launch-bound* if idle ≳ 30 % and its mean
kernel duration ≲ 20 µs; *bandwidth/compute-bound* if a few long kernels dominate with low
idle; *sync-stalled* if idle clusters immediately after a DtoH or explicit sync. Compare T31
against T255 to locate the launch-bound → compute-bound transition, and report the overall
GPU-busy fraction per configuration.

**C4 — optional second-level ranges.** If a dominant phase is internally opaque (most likely
`parameterization_tendencies!` or `dynamics_tendencies!`), rerun Phase B once with sub-ranges.
Preferred: call the phase's internal per-component functions from the driver with their own
`NVTX.@range`s, re-validated by B3. If the internals don't decompose cleanly at call level,
temporary `NVTX.@range` lines in model source in a **local scratch commit that is never
pushed** are acceptable — note it in the revision log and drop the commit afterwards.

### Phase D — `ncu` deep dives on the top kernels

**D1 — kernel selection.** From the C3 kernel tables take kernels with cumulative share ≥ 5 %
at T127 or T255 **and** mean duration ≥ 20 µs; cap at 5. Launch-bound phases contribute no
candidates by construction — their remedy is fusion/graphs, not kernel tuning, and `ncu`
(which serializes and replays kernels) cannot measure launch overhead anyway.

**D2 — permission preflight.** Hardware counters may be admin-restricted
(`NVreg_RestrictProfilingToAdminUsers`). Before investing in long runs, submit a tiny job:

```bash
/usr/local/cuda-13.0/bin/ncu --set basic -k "regex:gpu_" --launch-count 1 \
    julia --project=SpeedyWeather/test/GPU/modelbench \
    SpeedyWeather/test/GPU/modelbench/nsys_target.jl 31 8 2 2
```

If it errors with `ERR_NVGPUCTRPERM`, Phase D is blocked pending cluster support (report to
the user); Phases A–C are unaffected.

**D3 — profiled runs.** Reuse `nsys_target.jl` with a short window (`ncu` replay is slow;
keep total profiled launches ≤ ~10, e.g. `n_warmup=4 n_capture=6`). One job per kernel:

```bash
/usr/local/cuda-13.0/bin/ncu --set full \
    -k "regex:<kernel_name>" \
    --launch-skip 6 --launch-count 3 \
    -f -o SpeedyWeather/test/GPU/modelbench/reports/ncu_<kernel>_T<trunc> \
    julia --project=SpeedyWeather/test/GPU/modelbench \
    SpeedyWeather/test/GPU/modelbench/nsys_target.jl <trunc> 8 4 6
```

`--launch-skip` skips the first in-window instances so the profiled ones are steady-state.
The kernel-name regex alone is usually selective enough; if one kernel name is launched from
several phases, add `--nvtx --nvtx-include "<phase>/"` to restrict to the Phase B range
(consult `ncu --help` for the include syntax if no kernels match). Note `CUDA.@profile
external=true` gates only `nsys`; `ncu` profiles regardless — harmless.

**D4 — what to record.** Export text for archiving
(`ncu --import <file>.ncu-rep --page details > <file>.txt`) and record per kernel: duration,
grid/block dims, registers/thread, achieved occupancy, SM "Speed of Light" %, Memory SOL %,
DRAM throughput in GB/s **and the GPU it was measured on** (as % of that GPU's peak — ~3.35
TB/s on H100, ~0.70 TB/s on A40; see *Tooling* above), L2 hit rate, and the top stall reason. Classify: *bandwidth-bound* (Mem SOL ≥ 60 %), *compute-bound* (SM SOL ≥ 60 %),
*latency/occupancy-bound* (both ≤ 40 %). Append one short section per kernel to this document:
classification + one-line diagnosis + the obvious remedy direction.

### Phase E — synthesis

Append to this document: the per-phase tables, the kernel table, and a **ranked bottleneck
list**, each entry with (a) measured cost share per configuration, (b) diagnosis
(launch-bound / bandwidth-bound / sync-stall / algorithmic), (c) estimated ceiling from fixing
it (Amdahl), (d) proposed remedy. That list seeds the follow-up optimization plan(s); no
optimization is implemented under this plan.

## Findings

### Phase A — macro timing (2026-08-07)

`nlayers = 8`, `Float32`, default `PrimitiveWetModel`, no output, no callbacks. Every window is
≥ 4 s of steady-state stepping after 32 spin-up steps; every row passed the stability guard
(`nans_detected == false` **and** all spectral prognostic variables finite), so no row is a
blown-up run being rewarded for failing early.

**H100 80GB HBM3** (SLURM job 1729291) and **A40** (login node), ms per time step:

| trunc | H100 default | H100 nofeedback | H100 dynamics_only | A40 default | A40 nofeedback | A40 dynamics_only |
|------:|-------------:|----------------:|-------------------:|------------:|---------------:|------------------:|
| 31    | 2.760  | 2.317 (−16.1 %) | 2.573 (−6.8 %) | 3.590  | 3.214 (−10.5 %) | 3.294 (−8.2 %) |
| 63    | 5.120  | 3.784 (−26.1 %) | 4.938 (−3.5 %) | 6.256  | 5.052 (−19.2 %) | 6.183 (−1.2 %) |
| 127   | 10.485 | 7.879 (−24.9 %) | 9.882 (−5.8 %) | 13.298 | 9.863 (−25.8 %) | 12.837 (−3.5 %) |
| 255   | 23.833 | 20.490 (−14.0 %) | 24.687 (+3.6 %) | 30.600 | 25.541 (−16.5 %) | 30.214 (−1.3 %) |

(The single positive entry, H100 T255 `dynamics_only`, is within run-to-run noise; read it as
"indistinguishable from zero", not as physics being free-plus.)

#### A1. The feedback path costs more than all of the physics

Disabling `Feedback` saves **10–26 %** of every time step, at every resolution, on both GPUs.
Disabling the *entire* parameterization suite — convection, radiation, large-scale
condensation, surface fluxes, boundary layer, vertical diffusion, land, ocean, sea ice — saves
**1–8 %**. A diagnostic convenience is more expensive than the model physics.

The mechanism is host round-trips per time step, in `progress!`
(`SpeedyWeather/src/output/feedback.jl:90`):

- `max_speed(vars)` → `extrema(vars.grid.u)` — a device reduction whose result is read on the
  host; cost grows with grid size.
- `temperature_range(vars)` → `extrema(vars.grid.temperature)` — a second one.
- `nan_detection!` → `GPUArrays.@allowscalar vor[2, end]` — a scalar device-to-host read, every
  step, whenever `feedback.debug` (default `true`).

Each forces a synchronization that also drains the host's ability to run ahead queueing
kernels, so the true cost exceeds the transfers themselves. Note `verbose` defaults to
`isinteractive()` (false in a script), so **this cost is paid even when nothing is printed** —
`max_speed`/`temperature_range` are gated on `show_umax`/`show_temperature_range` (both default
`true`) and the progress-meter counter, not on `verbose`.

That the cost peaks at T63–T127 and eases at T255 fits: the fixed sync latency is amortized
over more real work as resolution grows, while the `extrema` reductions themselves grow.

#### A2. The model is latency-bound on GPU, not bandwidth- or compute-bound

The H100 is only **1.22–1.30× faster than the A40**, and that ratio is flat across all four
truncations:

| trunc | 31 | 63 | 127 | 255 |
|---|---:|---:|---:|---:|
| A40 / H100 speedup | 1.30 | 1.22 | 1.27 | 1.28 |

For reference the hardware ratios are ~**4.8×** in peak HBM bandwidth and ~**1.8×** in FP32
throughput. A speedup that matches neither, and does not improve as the problem grows, points
at a limiter both GPUs share — per-launch and per-sync latency — rather than at bandwidth or
arithmetic. This is independent corroboration of the launch-bound hypothesis, obtained before
any profiler was attached, and it argues that fusion / graph capture / sync removal will pay
better than tuning individual kernels.

#### A3. Consequence for the benchmark suite

`SpeedyWeather/benchmark/README.md` numbers are produced with the default `Feedback`, so the
published GPU SYPD figures carry this 10–26 % overhead. That is arguably the honest
"as-shipped" number, but it means the benchmark is partly measuring its own instrumentation.
Worth stating explicitly in the benchmark README either way.

### Phase B/C — per-phase and per-kernel accounting (2026-08-07)

Traces: `reports/model_T{31,63,127,255}_L8.nsys-rep` (H100, SLURM job 1729291) and
`reports/a40_T*_L8.nsys-rep` (A40, login node, `--trace=cuda,nvtx` only). 20 captured
steady-state steps each. The B3 equivalence check passed on GPU as well as CPU, so the driver
profiles the same computation `run!` does.

Tables below are from the **H100** traces. Kernel *durations* and *counts* from CUPTI are
reliable; host-side spans in those traces are inflated by `osrt` (see the note below), so
"GPU busy %" is computed against the **native** ms/step from Phase A, not against the traced
host span.

#### C1. The GPU is idle 70–78 % of every time step

| trunc | native ms/step (Phase A) | GPU kernel busy ms/step | GPU busy % | kernels/step |
|------:|-------------------------:|------------------------:|-----------:|-------------:|
| 31  | 2.760  | 0.619 | 22 % | 301 |
| 63  | 5.120  | 1.167 | 23 % | 545 |
| 127 | 10.485 | 2.571 | 25 % | 1073 |
| 255 | 23.833 | 7.020 | 29 % | 2237 |

(Caveat: kernel durations are measured under the profiler and may differ slightly from native,
so the busy fractions are approximate — but not by enough to change the conclusion.)

This is the quantitative form of the Phase A2 inference: the model does not keep the GPU fed.
Note also that **kernel count per step grows with truncation** — 301 → 2237 — which it should
not, since the code path is identical at every resolution. That is the thread that leads to the
root cause below.

#### C2. Per-phase breakdown (GPU busy time, per step)

T255, H100 (`busy` = summed kernel durations inside the range; `%ofbusy` = share of the step):

| phase | busy ms | kernels/step | % of busy |
|---|---:|---:|---:|
| dynamics | 5.895 | 2180 | 84.0 % |
| transform | 0.736 | 11 | 10.5 % |
| parameterizations | 0.213 | 7 | 3.0 % |
| feedback_output | 0.051 | 4 | 0.7 % |
| diffusion_implicit | 0.048 | 8 | 0.7 % |
| ocean_land | 0.024 | 9 | 0.3 % |
| leapfrog | 0.026 | 8 | 0.4 % |
| reset_tendencies | 0.026 | 10 | 0.4 % |
| **TOTAL (step)** | **7.020** | **2237** | 100 % |

`dynamics` dominates at every truncation (70 % at T31 → 84–87 % at T127/T255) and issues
**97 % of all kernel launches**. Every other phase launches ≤ 11 kernels per step.

`feedback_output` shows only 0.05 ms of GPU time yet costs 10–26 % of the wallclock step
(Phase A1) — exactly the signature of a host-side stall rather than device work.

#### C3. Root cause: an unplanned FFT batch size `K=2` degenerates into per-ring launches

The kernel table is dominated by launches that do almost no work:

| kernel | T31 | T63 | T127 | T255 | avg duration |
|---|---:|---:|---:|---:|---:|
| `gpu_broadcast_kernel_linear` | 96 | 192 | 384 | 768 | **1.23 µs at every resolution** |
| `preprocess_kernel` (cuFFT) | 8 | 44 | 132 | 404 | ~1.3 µs |
| `regular_fft_factor` (cuFFT) | — | 16 | 56 | ~192 | ~1.6 µs |

A kernel whose duration is constant at 1.23 µs from T31 to T255 is not doing resolution-scaled
work; it exists to be launched. The counts are exactly **`4 × nlat_half`** (nlat_half = 24, 48,
96, 192 for these truncations), which is precisely what `_fourier_serial!` emits: it loops

```julia
for (k, k_grid) in zip(1:nlayers, eachlayer(field))     # SpeedyTransforms/src/fourier.jl:188
    for j_north in 1:nlat_half                          # ... one FFT per ring ...
```

over `nlayers × nlat_half` and applies an FFT to the northern *and* southern ring — i.e.
`2 × K × nlat_half` launches, which for `K = 2` is `4 × nlat_half`. ✔ matches the measurement
at all four truncations.

That path is taken when no batched plan exists for the requested batch size
(`SpeedyTransforms/src/fourier.jl:5`):

```julia
K = size(field, 2)
if K > 1 && haskey(S.rfft_plans_batched, K)
    _fourier_batched!(...)     # one FFT per ring, all layers at once
else
    _fourier_serial!(...)      # host loop over layers x rings x hemispheres
```

and the planned set (`SpeedyWeather/src/dynamics/spectral_grid.jl:339`) is

```julia
default_transform_batch(::Type{<:AbstractArchitecture}, nlayers) =
    Int[1, nlayers, 2nlayers, 4nlayers + 1, 6nlayers + 1, 9nlayers + 1]   # = [1,8,16,33,49,73]
```

**`K = 2` is not in it.** One transform per time step requests `K = 2`, misses the plan, and
degenerates into ~768 one-microsecond launches plus ~600 tiny cuFFT kernels at T255 — roughly
**60 % of all kernel launches in the entire time step**, doing a negligible share of the work.

This independently corroborates an observation in `layer-count-performance-anomaly.md`, which
recorded "exactly the same single serial fallback (`K=2`, two calls per two steps)" for every
layer count but treated it as incidental. It is not incidental — it is the bulk of the launch
traffic and therefore, in a launch-bound regime, the bulk of the time step.

**Status of this finding:** the mechanism is identified from kernel counts matching
`4 × nlat_half` at four truncations plus the code path, and is consistent with the earlier
independent measurement. What is **not yet done** is the direct confirmation and the fix
measurement — see *Next steps*.

#### C4. Sync and transfer traffic (T31, H100)

| CUDA API call | calls/step |
|---|---:|
| `cuLaunchKernelEx` | 197 |
| `cuLaunchKernel` | 104 |
| `cuMemcpyDtoDAsync_v2` | 96 |
| `cuGraphLaunch` | 3 |
| `cuMemcpyDtoHAsync_v2` | 3 |
| `cuStreamSynchronize` | 6 |

Three device-to-host copies and six stream synchronizations per step — the `feedback_output`
cost of Phase A1. Only **3 graph launches per step**: the CUDA-graph capture covers the Fourier
step of `transform!` and essentially nothing else, so the ~2200 launches in `dynamics` are all
individual. The 96 device-to-device async copies per step at T31 are also worth a look.

#### C5. Phase D is blocked on this cluster

`ncu` fails with `ERR_NVGPUCTRPERM` on the SLURM H100 nodes: hardware performance counters are
restricted to admin users. Enabling them needs cluster support
(`NVreg_RestrictProfilingToAdminUsers=0`).

This matters less than it would have: `ncu` answers "why is this kernel slow", and the finding
is that the kernels are *not* slow — the step is launch- and sync-bound. `ncu` cannot measure
launch overhead. Phases A–C stand on their own.

### Methodological note: nsys imposes a large fixed per-step overhead (and it is *not* `osrt`)

Under nsys a time step takes ~0.5–0.7 s instead of 2.8–30 ms. **An initial reading blamed
`--trace=osrt` and Phase B was re-run on the A40 with `--trace=cuda,nvtx` only. That was
wrong** — the numbers are essentially identical with and without `osrt`:

| trunc | H100 traced, **with** osrt | A40 traced, **without** osrt | native (A40) |
|------:|---------------------------:|-----------------------------:|-------------:|
| 31  | 566.9 ms/step | 578.4 ms/step | 3.59 ms |
| 63  | 484.3 | 577.0 | 6.26 |
| 127 | 515.3 | 603.1 | 13.30 |
| 255 | 565.0 | 670.9 | 30.60 |

Dropping `osrt` changed nothing. Note also that the traced cost is nearly **flat in
truncation** (578 → 671 ms/step while the native step grows 8.5× and the launch count grows
7.4×), so it is neither per-unit-work nor per-launch, but a roughly constant per-step cost of
the profiling machinery itself.

Consequences for reading these traces:

- **Usable:** kernel durations, kernel counts, launch counts, API call counts, per-phase kernel
  time (`nvtx_kern_sum`) — everything Phases C1–C4 rest on.
- **Not usable:** traced wallclock, host range spans, and any GPU-idle fraction derived from
  them, in *either* trace set. This is why C1 computes GPU busy % against the **native**
  Phase A ms/step rather than against anything in the trace.

The cause of the fixed overhead has not been chased (it is not needed for the findings). If a
clean host-span attribution is ever wanted, that is the thing to investigate first.

### Phase E — ranked bottleneck list (preliminary)

Ordered by expected payoff. Items 1 and 2 are measured; item 3's mechanism is identified but
its fix is not yet measured.

1. **Unplanned `K=2` FFT batch → per-ring serial launches** (Phase C3). ~60 % of all kernel
   launches per step, doing a negligible fraction of the work, in the phase that holds 84 % of
   GPU time. *Diagnosis:* launch-bound. *Proposed remedy:* add `2` to
   `default_transform_batch` — a one-line change; `transform_batch` is already a `SpectralGrid`
   kwarg, so it is testable without touching source. *Expected ceiling:* if launch cost is
   ~5–7 µs and ~1370 of 2237 launches disappear at T255, that is ~7–10 ms off a 23.8 ms step,
   i.e. plausibly **1.4–1.7×**. Must be measured, not assumed — see *Next steps*.
2. **Per-step host round-trips in `progress!`** (Phase A1). Costs **10–26 %** of every step at
   every resolution on both GPUs — more than the entire parameterization suite. *Diagnosis:*
   sync-stall. *Remedy options:* run `max_speed`/`temperature_range`/`nan_detection!` every N
   steps instead of every step; keep the reductions on device and only copy when actually
   displayed; or make them opt-in. Note `verbose = isinteractive()` means a script pays this
   cost while printing nothing. *Expected ceiling:* up to the full 10–26 %.
3. **`dynamics` launches ~2200 kernels/step** (Phase C2), of which the `K=2` fallback explains
   the majority but not all. Once item 1 lands, re-profile: the residual per-launch traffic is
   the next target, and the candidates are kernel fusion in the tendency chain and extending
   CUDA-graph capture beyond the Fourier step (only **3 graph launches per step** today, C4).
4. **96 device-to-device async copies per step at T31** (C4). Not yet attributed to a source
   location; cheap to chase once the above are done.

Explicitly **not** a bottleneck, contrary to the usual expectation: individual kernel
efficiency. The two Legendre kernels are now the only large ones (`gpu_forward_legendre_kernel`
1.23 ms at T255) and the model is nowhere near bandwidth-bound — the H100, with ~4.8× the
bandwidth of the A40, delivers only ~1.28× the throughput.

## Next steps (not yet done)

1. **Run `verify_batch_k2.jl`** (written, not yet executed — the A40 was still tracing). It
   (a) instruments `_fourier!` to record every requested `K` and whether a batched plan existed,
   directly confirming the single `K=2` fallback and the predicted `4 × nlat_half` launch count,
   and (b) A/B-times the model with `transform_batch = [1, 2, L, 2L, 4L+1, 6L+1, 9L+1]` against
   the default. This converts finding C3 from "identified" to "confirmed and quantified".
2. **Analyse the A40 traces** (`reports/a40_T*_L8.nsys-rep`, all four written to disk; run
   `analyze_nsys.sh` with its glob adjusted to `a40_T*`, then `phase_table.jl`). Purpose is a
   **cross-GPU check of the kernel counts and per-phase shares**, in particular that the
   `4 × nlat_half` broadcast-launch pattern reproduces on the A40. It will *not* give cleaner
   host spans — see the methodological note; that hope was based on a mistaken `osrt`
   diagnosis.
3. **Decide on the `progress!` fix** (item 2 above) — a behavioural change, so it wants a
   decision rather than a unilateral edit.
4. **Ask cluster support about `NVreg_RestrictProfilingToAdminUsers`** if kernel-level counter
   analysis is ever wanted; not needed for the findings above.

## Testing and verification

- Phases A–D change no model source, so no unit tests are needed; the deliverables are scripts
  + measurements.
- The Phase B driver's bitwise equivalence check against `run!` guards against profiling a
  different computation than the real model runs.
- Every reported number carries the `nans_detected` guard; unstable configurations are reported
  as such, never as timings (T255 with 16 layers is a known suspect, see the layer-count
  investigation).
- Phase A is re-run once at the end to confirm the instrumented driver's per-phase times sum to
  ≈ the uninstrumented ms/step (NVTX overhead sanity check).

## Documentation changes

Findings and tables are appended to this document (status → **in progress** → **completed**).
No user-facing docs change. If the follow-up NVTX-in-source PR happens, it documents itself.

## Known limitations

- Single-GPU, CUDA/A40 only; conclusions about launch overhead may shift on other backends
  (Metal/AMDGPU) and on GPUs with different launch latency.
- The driver-script approach duplicates ~20 lines of `time_step!`; it can silently drift from
  the real time step if that function changes — mitigated by the equivalence check, but the
  check must be re-run after rebases.
- `nsys` osrt tracing adds some host overhead; per-phase *GPU* times are robust, per-phase CPU
  spans are upper bounds.
- Default-configuration physics only; bottleneck ranking may differ with other
  parameterization choices (e.g. different convection schemes).

## Future work

- The actual optimization plan(s) that Phase E seeds (kernel fusion, extending CUDA-graph
  capture beyond the Fourier step, removing per-step syncs, tuning individual kernels).
- Optional permanent `NVTX.@annotate` instrumentation PR.
- Folding the Phase A harness guards (NaN check, setup exclusion, longer windows) back into
  `SpeedyWeather/benchmark/manual_benchmarking.jl`, as already recommended by the layer-count
  investigation.
