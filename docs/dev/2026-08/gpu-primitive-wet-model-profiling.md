# GPU profiling of the full `PrimitiveWetModel` time step

> Status: **in progress**. Scripts for Phases A–C written and smoke-tested; GPU runs submitted.

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
  (see `submit_bench.sh`); GPUs are Nvidia A40 (peak ~696 GB/s HBM, 6 MB L2).
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
DRAM throughput in GB/s and as % of the A40 peak (~696 GB/s), L2 hit rate, and the top stall
reason. Classify: *bandwidth-bound* (Mem SOL ≥ 60 %), *compute-bound* (SM SOL ≥ 60 %),
*latency/occupancy-bound* (both ≤ 40 %). Append one short section per kernel to this document:
classification + one-line diagnosis + the obvious remedy direction.

### Phase E — synthesis

Append to this document: the per-phase tables, the kernel table, and a **ranked bottleneck
list**, each entry with (a) measured cost share per configuration, (b) diagnosis
(launch-bound / bandwidth-bound / sync-stall / algorithmic), (c) estimated ceiling from fixing
it (Amdahl), (d) proposed remedy. That list seeds the follow-up optimization plan(s); no
optimization is implemented under this plan.

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
