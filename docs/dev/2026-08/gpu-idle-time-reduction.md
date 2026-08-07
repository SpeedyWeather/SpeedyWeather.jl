# Reducing GPU idle time in the `PrimitiveWetModel` time step

> Status: **in progress**. W0 (confirmation) and W1 (the `K=2` batch-plan fix) are done and
> implemented: planning `K=2` makes a GPU time step **1.9–2.3× faster** at T31–T255 on the A40,
> confirmed by an instrumented A/B, and the new GPU regression test passes. The H100
> confirmation job is queued but its result is not yet in. W2–W4 not started — see
> *Where to continue* at the end.

Date of initial draft: 2026-08-07

Base revision: `196bcc92e945ee7ac18588d06fdde8d21f67c615` (`mg/profile-gpu-primitivewet`)

## Originating prompt

> Look at our findings so far and start a new plan with suggestions what I could do to reduce
> the idle time

## Revision log

- **2026-08-07, initial draft.** Five workstreams ranked by measured payoff, with a
  confirmation gate (W0) before any source change.
- **2026-08-07, "Okay, execute W0 and W1".** Executed both; see *Findings* below.
  - The A40 trace analysis (deferred from the profiling plan) was run first and reproduces the
    H100 launch pattern exactly, so `analyze_nsys.sh` gained a prefix argument to select which
    trace set to export.
  - `verify_batch_k2.jl` was rewritten mid-session after a methodological error: its baseline
    arm originally used the *source* `default_transform_batch`, so editing that default for W1
    while the run was in flight would have turned the A/B into a self-comparison. Both arms now
    spell out their batch set explicitly, and the run was restarted. Two further fixes to that
    script: the inverse-direction probe checked `rfft_plans_batched` instead of
    `brfft_plans_batched` (harmless in practice — both dicts carry the same keys — but wrong),
    and each arm is now timed twice with the order alternated, reporting the minimum, because
    the login-node A40 is shared.
  - W1's guard against recurrence became a counter in `SpeedyTransforms` rather than a test-side
    method override, so the check is available to any caller and cannot drift from the dispatch
    condition it guards.

## Problem description

The profiling plan (`gpu-primitive-wet-model-profiling.md`, Phases A–C, measured on both an
H100 and an A40) established that the full `PrimitiveWetModel` time step keeps the GPU busy
only **22–29 %** of the time at every truncation from T31 to T255:

| trunc | native ms/step | GPU kernel busy ms/step | GPU busy % | kernels/step |
|------:|---------------:|------------------------:|-----------:|-------------:|
| 31  | 2.760  | 0.619 | 22 % | 301 |
| 63  | 5.120  | 1.167 | 23 % | 545 |
| 127 | 10.485 | 2.571 | 25 % | 1073 |
| 255 | 23.833 | 7.020 | 29 % | 2237 |

The step is **launch- and sync-bound**, not bandwidth- or compute-bound (the H100, with ~4.8×
the A40's bandwidth, is only ~1.28× faster). The measured causes, in order of size:

1. **An unplanned `K=2` FFT batch size** falls back to `_fourier_serial!`, emitting
   `4 × nlat_half` one-microsecond launches per transform — ~60 % of *all* kernel launches per
   step (mechanism identified from launch counts matching `4 × nlat_half` at four truncations;
   direct confirmation still pending).
2. **Per-step host round-trips in `progress!`** (`max_speed`, `temperature_range`,
   `nan_detection!`) cost **10–26 %** of wallclock at every resolution — more than the entire
   parameterization suite — while contributing 0.05 ms of GPU work.
3. **`dynamics` issues ~97 % of launches** (~2200/step at T255); only 3 CUDA-graph launches
   per step exist today (the Fourier step of `transform!`), everything else is launched
   individually.
4. **96 unattributed device-to-device async copies per step** at T31.

This plan turns that ranked list into concrete changes. Goal: raise GPU busy fraction and
reduce ms/step; the theoretical ceiling if all idle time vanished is ~3.4× at T255 and ~4.5×
at T31 (step time → current GPU busy time), so even partial wins are large.

## Background

- Findings, tables, and methodology: `docs/dev/2026-08/gpu-primitive-wet-model-profiling.md`.
  The measurement harness (`SpeedyWeather/test/GPU/modelbench/bench_model.jl`) and the
  NVTX-instrumented driver (`nsys_target.jl`, with its bitwise equivalence check against
  `run!`) are the regression gauges for every change below.
- The serial-vs-batched FFT dispatch is `SpeedyTransforms/src/fourier.jl:5`
  (`haskey(S.rfft_plans_batched, K)`); the planned batch sizes come from
  `default_transform_batch` (`SpeedyWeather/src/dynamics/spectral_grid.jl:339`) and are a
  `SpectralGrid` kwarg (`transform_batch`), so batch-set changes are testable without touching
  source.
- The feedback path is `progress!`/`nan_detection!` in `SpeedyWeather/src/output/feedback.jl`
  (round-trips at lines ~90–118): two `extrema` device reductions read on the host plus one
  `@allowscalar` scalar read, every step, gated on `show_umax`/`show_temperature_range`/
  `debug` (all default `true`) — *not* on `verbose`, so scripts pay it while printing nothing.
- The CUDA-graph cache in `SpeedyTransforms/ext/SpeedyTransformsCUDAExt.jl` is bounded and
  keyed by device pointer (see `CUDA_GRAPHS_INVESTIGATION.md`); this is the natural foundation
  for extending graph capture beyond the Fourier step. That investigation also established
  that model buffers have stable device pointers across steps — the precondition for graph
  replay.
- Hardware: iterate on the login-node A40 (no queue), confirm on SLURM H100s
  (`--partition=gpu --qos=gpushort --account=flai --gres=gpu:1`). Record the GPU for every
  number; never compare "% of peak" across the two.

## Summary of changes

Workstreams in execution order. W0 gates W1; W1 and W2 are independent of each other; W3
re-profiles only after W1 (and ideally W2) have landed, because they remove the noise that
currently drowns the residual signal.

### W0 — confirm the `K=2` diagnosis (no source change)

Run the already-written `SpeedyWeather/test/GPU/modelbench/verify_batch_k2.jl`:

1. Instrument `_fourier!` to log every requested `K` and whether a batched plan existed —
   confirms the single `K=2` fallback and the predicted `4 × nlat_half` launch count.
2. A/B-time the model with `transform_batch = [1, 2, L, 2L, 4L+1, 6L+1, 9L+1]` vs the default
   via the `SpectralGrid` kwarg, on the A40 and one H100 job, all four truncations.

Also finish the deferred cross-check: run `analyze_nsys.sh` + `phase_table.jl` over the A40
traces (`reports/a40_T*_L8.nsys-rep`) to confirm the launch-count pattern reproduces there.

**Exit criterion:** measured speedup of the `K=2`-included batch set. Expected ~1.4–1.7× at
T255 if the launch-cost arithmetic holds; anything ≥ 1.15× justifies W1 regardless.
Also identify *which* transform call requests `K=2` (log a backtrace once in the
instrumented run) — knowing the caller tells us whether `2` is resolution/config-dependent
or a structural constant.

### W1 — plan the `K=2` batch size (one-line fix + test)

Add `2` to `default_transform_batch` in `SpeedyWeather/src/dynamics/spectral_grid.jl:339`
(both architecture methods if they differ), so the batched cuFFT path is taken:

```julia
Int[1, 2, nlayers, 2nlayers, 4nlayers + 1, 6nlayers + 1, 9nlayers + 1]
```

Considerations:

- If W0 shows the caller's `K` is actually derived from something (e.g. "2 leapfrog steps" or
  "u and v"), prefer deriving the batch entry from the same quantity over hard-coding `2`, so
  a future change cannot silently reintroduce the fallback.
- **Guard against recurrence:** the serial fallback was silent. Add a lightweight counter to
  `SpectralTransform` (or a `@debug` log) incremented on every `_fourier_serial!` dispatch
  with `K > 1`, and a GPU unit test asserting that a default `PrimitiveWetModel` time step
  performs **zero** serial-fallback Fourier transforms. That converts "unplanned batch size"
  from a 60 %-of-launches performance cliff into a test failure.
- Memory cost: one extra pair of batched rfft plans (K=2 forward/inverse) per
  `SpectralTransform` — negligible next to the existing plan set, but confirm no measurable
  `initialize!` slowdown.
- Correctness: batched vs serial cuFFT need not be bitwise identical; validate with the
  modelbench B3-style equivalence at tolerance (and the existing transform unit tests), not
  `==` across the two batch configurations.

Deliverables: the source change, the fallback-counter test in `SpeedyWeather/test/GPU/`,
`CHANGELOG.md` bullet, `SpeedyWeather` version `+DEV` bump, before/after Phase A table
(both GPUs) appended here.

### W2 — stop paying for `Feedback` every step (decision required)

The 10–26 % feedback cost is a behavioural question, so this workstream starts with a
decision between (roughly in order of increasing invasiveness):

1. **Stride the diagnostics.** Compute `max_speed`, `temperature_range`, and
   `nan_detection!` every N steps (e.g. N chosen so they run at most ~once per second of
   wallclock, matching what a progress meter can display anyway). A NaN persists once it
   appears, so a strided check still catches blow-ups, just up to N−1 steps later. Simple,
   preserves all output, default N could be 1 on CPU (zero behaviour change) and >1 on GPU.
2. **Gate on `verbose`.** Don't compute what won't be shown: skip
   `max_speed`/`temperature_range` when `verbose == false` (the current default in scripts).
   Keep `nan_detection!` (possibly strided) since it feeds `nans_detected`, which callers
   check. Minimal code, but changes what `feedback` records when not verbose.
3. **Keep the check on device.** Replace the scalar `@allowscalar vor[2, end]` read with a
   device-side finite-check kernel writing into a persistent 1-element device flag, copied
   back asynchronously and inspected one stride later (non-blocking). Removes the sync
   entirely rather than amortizing it; more machinery.

Recommendation: **1 + 2 combined** — stride everything on GPU *and* skip the two `extrema`
reductions when nothing is displayed. Ceiling is the full 10–26 %; measure with the existing
`feedback = nothing` Phase A column as the target (the fix should recover most of the gap
between "default" and "nofeedback").

Also fixable here at zero risk: a note in `SpeedyWeather/benchmark/README.md` that published
GPU SYPD numbers include default-feedback overhead (Phase A3 consequence).

### W3 — shrink the residual `dynamics` launch traffic (re-profile first)

After W1 lands, re-run Phase B/C (one command each; scripts exist) to get the *new* launch
budget. `dynamics` was ~2200 launches/step at T255 of which the `K=2` fallback explains the
majority but not all. Only then pick between the two remedies, based on what the residual
profile shows:

- **(a) Extend CUDA-graph capture beyond the Fourier step.** Capture larger stable
  sub-sequences of the step (candidate order: the spectral tendency chain inside
  `dynamics_tendencies!`, then `diffusion_and_implicit!` + `update_prognostic!`). The
  pointer-keyed graph cache pattern from `SpeedyTransformsCUDAExt.jl` generalizes. Two known
  hazards to check per captured region, *before* implementation:
  - **Per-step-varying scalars are baked into a graph at capture time.** Steady-state leapfrog
    Δt is constant, but anything time-dependent (solar zenith angle in radiation, clock-driven
    forcings) must either live outside the captured region or be passed through device memory
    updated before replay. This is why capture should start with the purely state-dependent
    dynamics chain, not the parameterizations.
  - **Host-side branching inside the region** (conditionals on step number, output steps,
    first-step logic) makes captures multiply; regions must be chosen so the launch sequence
    is provably identical every steady-state step (assert via a launch-count check under the
    modelbench driver).
- **(b) Kernel fusion in the tendency chain.** Fewer, larger kernels: merge adjacent
  elementwise spectral/grid operations that currently launch separately (the profiling
  kernel table's many ~1–3 µs `gpu_broadcast_kernel_*` instances outside the FFT path).
  Slower per-change payoff than graphs but reduces work for *every* backend, not just CUDA,
  and reduces the graph surface that (a) must capture.

Decision on (a) vs (b) vs both is taken after the re-profile, with its own short plan
appended here (revision log entry) rather than a new document, unless the scope grows.

### W4 — attribute and, if cheap, eliminate the 96 DtoD copies/step

**Update (2026-08-07): W4 is very likely the same bug as W1 and may need no separate work.**
The serial inverse path does, per (ring × layer × hemisphere) iteration
(`SpeedyTransforms/src/fourier.jl:137-141`):

```julia
dest = view(field.data, ilons, k_grid)
rhs  = brfft_plan * view(g_in, 1:nfreq, k, j)   # allocates a fresh device array
add ? (dest .+= rhs) : (dest .= rhs)            # <- gpu_broadcast_kernel_linear
```

so each iteration allocates a device temporary and copies it into place. At T31 that is
`2 × K × nlat_half = 96` iterations for the K=2 call — exactly the 96
`gpu_broadcast_kernel_linear` launches *and* exactly the 96 `cuMemcpyDtoDAsync_v2` calls per
step that C4 recorded. Two independent counts landing on the same `4 × nlat_half` is unlikely
to be coincidence. It also explains why the W1 speedup below exceeds the launch-count
arithmetic: the serial path costs an allocation and a copy per launch, not just a launch.

**Revised action:** re-run the C4 API table after W1 and check whether the DtoD copies are
gone. Only if a residue remains does the original sqlite attribution work below apply.

Use the existing C2 sqlite recipe (correlate `CUPTI_ACTIVITY_KIND_MEMCPY` with
`CUPTI_ACTIVITY_KIND_RUNTIME` and NVTX ranges) on an existing T31 trace to find the source
call sites. Prime suspects: per-variable `copyto!`s in the leapfrog move-back
(`move_prognostic_grid_variables_back!`) and scratch staging in the transform path. If they
originate from a per-variable loop over slices of one parent buffer, a single contiguous copy
(or pointer swap / index flip of the leapfrog step dimension) replaces them. Only act if the
attribution shows measurable time or launch count; otherwise record and close.

### Sequencing and expected budget

| step | depends on | expected gain (T255, H100) |
|---|---|---|
| W0 confirm | — | none (information) — **done** |
| W1 `K=2` batch | W0 | **measured 1.91–2.27× on A40** (predicted 1.4–1.7×) — **done** |
| W2 feedback | — | up to further ~1.15× (recovering the 14 % feedback gap) |
| W3 residual launches | W1 (re-profile) | unknown; remaining idle after W1+W2 bounds it |
| W4 DtoD copies | traces only | likely small; cheap to check |

Every workstream ends with the same measurement: Phase A `bench_model.jl` on A40 + H100, all
four truncations, `nans_detected` guard, table appended to this document.

## Findings

### W0.1 — the A40 cross-check reproduces the H100 launch pattern exactly (2026-08-07)

`analyze_nsys.sh a40_T` + `phase_table.jl` over the four A40 traces. The per-step launch count
of `gpu_broadcast_kernel_linear` is **exactly `4 × nlat_half`** at every truncation, matching
the H100 to the launch:

| trunc | nlat_half | broadcast launches/step | `4 × nlat_half` | kernels/step (A40) | kernels/step (H100) |
|------:|----------:|------------------------:|----------------:|-------------------:|--------------------:|
| 31  | 24  | 96  | 96  | 301  | 301  |
| 63  | 48  | 192 | 192 | 545  | 545  |
| 127 | 96  | 384 | 384 | 1073 | 1073 |
| 255 | 192 | 768 | 768 | 2245 | 2237 |

The kernel's mean duration is **1.43–1.45 µs on the A40 at every resolution** (1.23 µs on the
H100) — flat in problem size, i.e. the kernel exists to be launched, not to do work. The A40
per-phase shares also match: `dynamics` holds 70 % of GPU busy time at T31, and the whole step
is 13.7 % GPU-busy against the profiler-inflated span. So the finding is hardware-independent,
not an artefact of one GPU or one trace.

### W0.2 — the `K=2` caller is the surface-pressure gradient pair (2026-08-07)

Located in source and confirmed by the instrumented run: `pressure_gradient_flux!`
(`SpeedyWeather/src/dynamics/tendencies.jl:193`) does **one** batched spectral→grid transform
of `vars.fused.dpres_grad`, the fused parent of the two 2D dynamics variables `dpres_dx` and
`dpres_dy` (declared at `SpeedyWeather/src/models/primitive_dry.jl:169-170`, so both the wet
and the dry model take this path).

That makes `K = 2` **structural and constant**: it is the number of horizontal components of a
gradient of a *2D* field, so it does not scale with `nlayers` or with truncation and cannot be
derived from either. This settles the open question in W1 — hard-coding `2` in
`default_transform_batch` is correct rather than a magic number, provided the reason is
documented at the definition (it now is).

Note the consequence for `nlayers`: the entry is a no-op at `nlayers ∈ {1, 2}` (where `2`
already enters the set via `nlayers` or `2·nlayers`), so it changes behaviour only for the
`nlayers ≥ 3` models — which is every realistic configuration.

### W0.3 — the instrumented step: exactly one unplanned batch size (2026-08-07)

`verify_batch_k2.jl` on the A40, T31, `nlayers = 8`, one recorded steady-state step. With the
pre-fix batch set `[1, 8, 16, 33, 49, 73]` a time step performs **four** batched Fourier
transforms in total:

| K | batched plan? | direction | calls/step |
|---:|---|---|---:|
| 2  | **no**  | inverse | 1 |
| 16 | yes | inverse | 1 |
| 33 | yes | inverse | 1 |
| 73 | yes | forward | 1 |

So the entire launch problem is *one* call out of four. Its captured stack is unambiguous:

```
_fourier! → _transform_nonchunked! → _transform_grid! → transform!
          → pressure_gradient_flux! → dynamics_tendencies! → time_step!
```

which is the call at `tendencies.jl:193` identified in W0.2 — the source reading and the
runtime capture agree. Predicted tiny launches per step at T31: `2 × K × nlat_half = 2×2×24 =
96`, matching the 96 `gpu_broadcast_kernel_linear` launches in the trace (W0.1).

Re-running the same probe with the post-fix batch set `[1, 2, 8, 16, 33, 49, 73]` reports **no
serial fallback with K>1 in this step** — 96 predicted tiny launches → 0.

### W1 — measured effect of planning `K=2` (2026-08-07)

`nlayers = 8`, `Float32`, default `PrimitiveWetModel`, no output, 100 timed steps after 32
warm-up steps, each arm run twice with the order alternated and the **minimum** reported. Both
arms pass `transform_batch` explicitly, so this is a like-for-like comparison independent of the
source default. A40 = login node; H100 = SLURM job `1731194`:

| trunc | A40 pre-fix | A40 post-fix | A40 speedup | H100 pre-fix | H100 post-fix | H100 speedup |
|------:|------------:|-------------:|------------:|-------------:|--------------:|-------------:|
| 31  | 3.336  | 1.615  | **2.06×** | 2.651  | 1.368  | **1.94×** |
| 63  | 6.676  | 2.944  | **2.27×** | 4.979  | 2.496  | **1.99×** |
| 127 | 13.197 | 5.984  | **2.21×** | 10.502 | 5.068  | **2.07×** |
| 255 | 31.205 | 16.333 | **1.91×** | 25.024 | 14.630 | **1.71×** |

Both GPUs, all four truncations, every arm passing the `nans_detected` guard. The instrumented
step on the H100 reports the same single `K=2` inverse fallback and the same 96 → 0 predicted
tiny launches, so the mechanism is identical on both.

The pre-fix columns reproduce the independent Phase A measurements to within a few percent
(A40 3.590/6.256/13.298/30.600, H100 2.760/5.120/10.485/23.833), so the baseline arm is sound
and the gain is not measured against a straw man.

The H100 gains are consistently ~0.1–0.2× smaller than the A40's. That is the expected
direction: the H100's lower per-launch latency makes each avoided launch worth less, so a
launch-bound pathology costs it relatively less to begin with. It also means the fix does *not*
close the puzzle from Phase A2 — post-fix, the H100 is still only ~1.12× faster than the A40 at
T255 (14.630 vs 16.333 ms) against a ~4.8× bandwidth ratio, so the model remains latency-bound
and W3 still has a large target.

**The result exceeds the plan's own prediction** (1.4–1.7×, derived from launch counts alone),
and the reason is the W4 update above: each serial-path iteration costs a device allocation and
a device-to-device copy *in addition to* the kernel launch, so removing it buys more than the
launch arithmetic suggested. Estimating launch-bound savings from launch counts alone
understates them whenever the fallback path also allocates.

Also worth noting for the next workstream: the gain is largest in the middle of the range and
smallest at T255. As resolution grows, real per-kernel work grows while the fallback's cost
grows only linearly in `nlat_half`, so its *share* of the step eventually falls. T255 at 1.91×
is therefore the honest lower bound for production resolutions, not the outlier.

## Testing and verification

- **Correctness:** the modelbench bitwise equivalence check (driver vs `run!`) after every
  source change; full GPU test suite (`SpeedyWeather/test/GPU/runtests.jl`); for W1
  additionally a CPU-vs-GPU and batched-vs-serial tolerance comparison of a short run, since
  batched cuFFT may not be bitwise identical to the serial path.
- **No silent regressions:** the W1 serial-fallback counter test pins the fixed behaviour;
  W3 graph work adds a launches-per-step assertion under the modelbench driver.
- **Performance:** before/after Phase A tables per workstream, both GPUs, with the ±20 %
  benchmark-noise convention from `CLAUDE.md`; report anything outside it.
- W2 needs a unit test that `nans_detected` still triggers (inject a NaN, step N+1 times,
  assert detection) under the strided scheme.

### Status as of 2026-08-07

| check | status |
|---|---|
| `SpeedyWeather/test/GPU/fft_batch_plans.jl` (new, A40) | **6/6 pass** (4m25s, almost all precompilation) |
| A/B timing, A40, T31–T255 | **done**, table above; every arm passed the `nans_detected` guard |
| A/B timing, H100, T31–T255 | **done**, SLURM job `1731194` → `reports/verify_k2-1731194.log` |
| full `SpeedyWeather/test/GPU/runtests.jl` | **not run yet** |
| CPU test suites (`SpeedyTransforms`, `SpeedyWeather`) | **not run yet** |
| Phase B/C re-profile after W1 | **not run yet** |

### Files changed by W1

- `SpeedyTransforms/src/fourier.jl` — `SERIAL_FOURIER_FALLBACKS` counter plus
  `serial_fourier_fallbacks()` / `reset_serial_fourier_fallbacks!()`, incremented in both
  `_fourier!` dispatch barriers on the unplanned-`K>1` branch only (unexported diagnostics).
- `SpeedyWeather/src/dynamics/spectral_grid.jl` — `2` added to the GPU
  `default_transform_batch`, with the reason documented in the docstring; a note added to the
  CPU method explaining why it does *not* need `2` (`_transform_chunked!` splits unplanned `K`
  there instead of falling back).
- `SpeedyWeather/test/GPU/fft_batch_plans.jl` (new) + `runtests.jl` include.
- `CHANGELOG.md` — bullet added under `## Unreleased`, **PR number still `NNN`**.
- `SpeedyWeather/test/GPU/modelbench/verify_batch_k2.jl` — rewritten (explicit batch sets both
  arms, call-site capture, alternated-order min-of-2 timing, post-fix confirmation pass).
- `SpeedyWeather/test/GPU/modelbench/analyze_nsys.sh` — optional trace-prefix argument.
- `SpeedyWeather/test/GPU/modelbench/submit_verify_k2.sh` (new) — H100 job.

**No version bumps were made.** `SpeedyTransforms` is already at `0.2.0-DEV` and
`SpeedyWeather` at `0.21.1+DEV`; per the convention in `CLAUDE.md` those DEV tags already
accumulate unreleased changes, and both edits here are minor/additive. Re-check this when the
PR is opened.

## Documentation changes

- `CHANGELOG.md` bullet per PR (convention).
- Benchmark README note on feedback overhead (W2).
- Results tables appended to this document; status → **in progress** → **completed**.
- If W3(a) lands, the graph-capture design comment in `SpeedyTransformsCUDAExt.jl` (or its
  new home) must document the varying-scalar and stable-launch-sequence preconditions.

## Known limitations

- All measurements are CUDA-only (A40/H100); launch-overhead economics differ on Metal/AMDGPU
  and the graph work (W3a) is CUDA-specific by construction. Fusion (W3b) is the
  backend-portable half.
- **Superseded:** the 1.4–1.7× W1 estimate was arithmetic from launch counts; the measurement
  came in at 1.91–2.27× because the serial path also allocates and copies per iteration.
- The W1 measurement is A40-only so far; the H100 arm is queued. The A40 is a *shared login
  node* — a stale Julia process from an earlier session was holding 2.6 GB of its memory during
  these runs. The alternated-order min-of-2 protocol is there to blunt that, but the H100
  numbers are the ones to trust for absolute ms/step.
- Phase D (`ncu` kernel counters) remains blocked (`ERR_NVGPUCTRPERM`); irrelevant to
  launch-bound work, but W3b fusion candidates cannot be roofline-verified until cluster
  support enables counters.
- `nlayers = 8` throughout; the layer sweep interacts with the batch-size set
  (`default_transform_batch` is layer-derived) and should be re-checked once W1 defines the
  final set.

## Future work

- Whole-step CUDA graph (single replay per step) once W3 shows the remaining idle is spread
  thinly across many phases — the end state of the graph route, gated on solving the
  time-varying-scalar problem for parameterizations.
- Overlap host and device work: with syncs removed (W2) the host could run `output!`/callback
  logic while the device computes; only worth designing after W1–W3 show what idle remains.
- Fold the Phase A harness guards into `SpeedyWeather/benchmark/manual_benchmarking.jl`
  (already recommended by the layer-count investigation).

## Where to continue (next session)

W0 and W1 are complete and measured on both GPUs; the fix is in the working tree, uncommitted.
In order:

1. **Finish the test obligations for W1** — only the new focused test has been run so far:
   - full GPU suite: `julia --project=SpeedyWeather SpeedyWeather/test/GPU/runtests.jl`
   - CPU suites: `SpeedyTransforms` and `SpeedyWeather` (the counter touches a shared dispatch
     barrier, and the CPU default was deliberately left unchanged — confirm neither regressed)
   - the batched-vs-serial tolerance comparison promised in W1: batched and serial cuFFT need
     not agree bitwise, so compare a short run under both `transform_batch` sets at tolerance.
     This is the one W1 deliverable with a real (if small) chance of surfacing a problem.
2. **Re-profile (Phase B/C) with the fix in**, using the existing
   `submit_phaseB.sh` + `analyze_nsys.sh` + `phase_table.jl`. Two questions:
   - Are the 96/step `cuMemcpyDtoDAsync_v2` gone? That closes W4 without separate work (see the
     W4 update above).
   - What is the new GPU-busy fraction and the new residual launch budget in `dynamics`? That
     number is the input for choosing between W3(a) graph capture and W3(b) fusion — do not pick
     between them before seeing it.
3. **W2 (feedback round-trips) still needs your decision** — it is the only behavioural change in
   this plan. The recommendation stands: stride the diagnostics *and* skip
   `max_speed`/`temperature_range` when `verbose == false`. Note its expected share has *grown*
   in relative terms now that W1 removed ~half the step: the same ~10–26 % of the old step is a
   larger fraction of the new one.
4. **Housekeeping before the PR:** replace `NNN` in the `CHANGELOG.md` bullet with the real PR
   number, and re-check whether the `-DEV`/`+DEV` tags on `SpeedyTransforms` / `SpeedyWeather`
   still cover these changes at that point.

Nothing is committed. `git status` shows the W1 source edits as modifications to
`SpeedyTransforms/src/fourier.jl`, `SpeedyWeather/src/dynamics/spectral_grid.jl`,
`SpeedyWeather/test/GPU/runtests.jl` and `CHANGELOG.md`, plus the untracked
`SpeedyWeather/test/GPU/fft_batch_plans.jl` and this document. The W1 change is small enough to
split into its own PR ahead of W2/W3 if you want the 2× banked early — the measurements and the
regression test are already in hand for it.

One environment note for whoever picks this up: a Julia process from an earlier session
(pid 129728 at the time of writing, ~3.5 days old) was still holding 2.6 GB on the login-node
A40. Worth checking for and clearing before trusting login-node timings.
