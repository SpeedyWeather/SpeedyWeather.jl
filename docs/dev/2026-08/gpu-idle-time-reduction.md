# Reducing GPU idle time in the `PrimitiveWetModel` time step

> Status: **in progress**. W0, W1, W2 and W4 are **done**; W3 is **re-scoped and not started**.
> W1 (planning the `K=2` FFT batch) and W2 (striding the feedback diagnostics) together make a
> GPU time step **1.67× faster at T255 and 2.2–2.45× at T31–T127** on the H100, measured, with
> the full test suite behind them. W4 closed for free: the 96–768 device-to-device copies per
> step were W1's bug. The post-fix re-profile then **invalidated W3's premise** — with CUDA
> graph nodes expanded the step is **66–78 % GPU-busy**, not the 22–29 % this plan was written
> against, so the model is no longer launch-bound and neither of W3's two proposed remedies
> targets what is left. See finding W3.1 and *Where to continue* at the end.

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
    condition it guards. *(Superseded 2026-08-08: the counter was removed, see below.)*
- **2026-08-08, "continue the GPU idle-time reduction".** Short session; work started on the
  W1 test obligations and on the W2 diagnosis. Outcomes:
  - the stale login-node A40 process was cleared (see *Known limitations*), so future login-node
    timings are no longer taken against a GPU with 2.6 GB held by a dead benchmark;
  - `verify_batch_equivalence.jl` was written (the batched-vs-serial tolerance comparison that W1
    still owed) but **has not been run**;
  - **the `SERIAL_FOURIER_FALLBACKS` counter was removed at the user's request**, reverting
    `SpeedyTransforms/src/fourier.jl` to `HEAD`. W1 is now a *one-line* change to
    `default_transform_batch` plus its test. The recurrence guard is instead a direct assertion
    that the batched `K = 2` plans exist, which is the same condition `_fourier!` dispatches on;
    `SpeedyTransforms` needs no version bump any more, only `SpeedyWeather`;
  - W2 gained a measured mechanism (finding W2.0 below) which changes its character: the largest
    part is a *bug*, not a behavioural trade-off, and so needs no decision.
- **2026-08-08 (second session), "for W2b I want to stride the `nan_detection!` and progress meter
  calculations, make this a parameter called `interval`, set the default to 50 steps. Consider
  reusing the `check_iterations` property of `ProgressCore` as well."** W2 implemented as one
  striding mechanism covering both halves:
  - the user's decision on the open W2b question is **option 1 (stride)**, with the stride exposed
    as `Feedback.interval`, default `50`, and the same value used on CPU and GPU (the original
    plan sketched `1` on CPU — dropped, a single default is simpler and 50 steps of latency on a
    NaN abort is harmless on either device);
  - **W2a is subsumed rather than implemented as sketched.** The plan had proposed *skipping*
    `max_speed`/`temperature_range` whenever the progress meter is disabled. Instead they are
    strided by `max(interval, check_iterations)`, which reuses ProgressMeter's own adaptive
    counter for the enabled case (unchanged behaviour) and falls back to `interval` when the
    counter is pinned at `1` because the meter is disabled. Cost when nothing is displayed drops
    from 1 sync/step to 1/50 rather than to 0; the residue is negligible and, unlike skipping, it
    keeps `FEEDBACK_UMAX`/`FEEDBACK_TMIN`/`FEEDBACK_TMAX` populated for any consumer;
  - one addition not in the plan: `nan_detection!` also runs unconditionally on the **last** time
    step of a run. Without it a run shorter than `interval` would be checked only at step 0, which
    would silently hollow out the ~30 `@test model.feedback.nans_detected == false` assertions in
    the test suite, none of which runs 50 steps.

- **2026-08-08 (fourth session), "continue the gpu-idle-time-reduction".** No new source behaviour;
  the session closed the measurement debts that W1 and W2 had left open. The previous session had
  *produced* a post-W1+W2 Phase A run and a post-W1 Phase B trace set but analysed neither, so both
  were sitting on disk unread. Outcomes:
  - the combined W1+W2 Phase A result is in (finding W1W2.1): **1.67–2.45× faster** than the
    pre-plan baseline on the H100, **2.6–2.8×** on the A40, and the `nofeedback` gap that motivated
    W2 has collapsed from 14–26 % to −0.3…3.0 %, i.e. into the noise;
  - the Phase B/C re-profile is in (finding W3.0) and it **closes W4 outright**: the 96–768
    `cuMemcpyDtoDAsync_v2` per step are gone, along with the per-step device allocations. Every
    host-side per-step count is now *resolution-independent*;
  - one methodological defect found in the re-profile and fixed in the harness: `nsys` defaults to
    `--cuda-graph-trace=graph`, under which the kernels inside a CUDA graph are **not collected at
    all**. Post-W1 the K=2 FFT moved from the serial path (counted) into a graph (uncounted), so the
    post-fix "GPU busy ms" and "kernels/step" understate the real GPU work and must not be compared
    against the pre-fix ones. `submit_phaseB.sh` gained a granularity argument and the re-profile was
    re-run at `node` granularity (job `1732864`);
  - **that re-profile re-scopes W3 entirely** (finding W3.1): with graph nodes expanded the step is
    **66–78 % GPU-busy**, not 22–29 %. The model is no longer launch-bound, and 44 % of the remaining
    GPU time at T255 is inside cuFFT — which neither W3(a) nor W3(b) as written is aimed at;
  - mid-session the user reported the **progress bar never being drawn**. Diagnosed and fixed
    (finding W2.3); it is a pre-existing `main` bug, not a W2 regression. The first draft of the fix
    also computed a bar width that left room for SpeedyWeather's long progress string; **the user
    asked for that helper to be dropped**, so the fix is now only the field type and its default;
  - the CPU `SpeedyWeather` suite has **one failure, `dynamics/dispatch.jl`, that this plan did not
    cause**: reverting only `feedback.jl` to pre-W2 in the same environment still reproduces it, and
    the JET report names only KernelAbstractions and ModelParameters. See the status table;
  - **a caution for future sessions**: an A/B experiment that reverts a tracked file with
    `git checkout <rev> -- <path>` leaves that revert *staged*, and it was committed from the index
    as `ef75b832 "feedback reverse"` before it could be undone. Use a scratch copy for such
    experiments, or `git restore --staged` immediately after.

- **2026-08-08 (third session), "continue the GPU idle-time reduction"** plus, mid-session, *"Does the
  `max(interval, iterations)` make sense? And also calling `progress!` / `next!` still at every step?"*
  Both questions were answered by measurement and the first one found a real defect in W2.1:
  - `max(interval, check_iterations)` was **replaced by `enabled ? check_iterations : interval`**.
    The `max` rested on the claim that `check_iterations` is huge while the meter is enabled
    (32681, measured in the previous session). That measurement came from a bare `next!` loop with
    no work per iteration and does not describe a model run: `calc_check_iterations` converges to
    `feedback_dt / step_time`, so in a real *enabled* `PrimitiveWetModel` run it is **3–10**, and
    `max(50, 3…10) = 50` would have made the displayed max speed and temperature range 5–15× staler
    than the progress bar that shows them. See finding W2.2;
  - `progress!` / `ProgressMeter.next!` **stays at every step**: measured 4.6 ns/call on a disabled
    meter, and `p.counter` is the clock every stride in `progress!` is taken modulo, so it cannot be
    strided itself.

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
- **Guard against recurrence:** the serial fallback was silent. Assert in a GPU unit test that
  a default `SpectralTransform` and a default `PrimitiveWetModel` actually carry batched
  `K = 2` plans — `haskey(S.rfft_plans_batched, K)` is the exact condition `_fourier!`
  dispatches on, so the test fails the moment the fallback returns. (An earlier draft added a
  `SERIAL_FOURIER_FALLBACKS` counter to `SpeedyTransforms/src/fourier.jl` for this; it was
  **removed on 2026-08-08** — a global mutable counter in a hot dispatch barrier is a large
  permanent surface for a test-only concern, and checking the plan key tests the same
  condition with no source change at all.)
- Memory cost: one extra pair of batched rfft plans (K=2 forward/inverse) per
  `SpectralTransform` — negligible next to the existing plan set, but confirm no measurable
  `initialize!` slowdown.
- Correctness: batched vs serial cuFFT need not be bitwise identical; validate with the
  modelbench B3-style equivalence at tolerance (and the existing transform unit tests), not
  `==` across the two batch configurations.

Deliverables: the source change, the batch-plan test in `SpeedyWeather/test/GPU/`,
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
| W1 `K=2` batch | W0 | **measured 1.91–2.27× on A40, 1.71–2.07× on H100** — **done** |
| W2 feedback | — | **measured: feedback cost 14–26 % → −0.3…3.0 %** (finding W1W2.1) — **done** |
| W3 residual launches | W1 (re-profile) | **re-scoped, see W3.1**: the step is now 66–78 % GPU-busy, so the whole remaining idle is worth ≤1.36× at T255 and the launch-reduction premise is gone |
| W4 DtoD copies | traces only | **closed at zero cost** — the 96–768 DtoD copies/step were W1's bug and are gone (finding W3.0) |

Combined W1 + W2 on the H100: **2.45× at T31, 2.39× at T63, 2.20× at T127, 1.67× at T255.**

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

### W1.1 — batched and serial cuFFT agree *bitwise* (2026-08-08)

`verify_batch_equivalence.jl` on the A40, `Float32`. The plan expected a small round-off difference
(different batch shape → possibly different cuFFT algorithm) and defined a relative pass criterion
against the CPU-vs-GPU difference the suite already tolerates. The measured difference is **exactly
zero** on every comparison:

| comparison | max abs difference |
|---|---:|
| one `K=2` transform, forward, T31/63/127/255 | `0.000e+00` at every truncation |
| one `K=2` transform, inverse, T31/63/127/255 | `0.000e+00` at every truncation |
| 100-step `PrimitiveWetModel` T31 L8: vorticity, divergence, temperature, pressure, humidity | `0.000e+00` for all five |

with the CPU-vs-GPU reference difference of the same run at `8.7e-03` relative (divergence) — i.e.
the change is not merely inside the tolerated band, it does not move the state at all. That is the
expected outcome in hindsight: for the same transform length cuFFT selects by length, and batching
changes only how many independent transforms are issued per call, not the algorithm applied to each.

Scope note: the model-run arm compares the five array-valued prognostics at the top level.
Namespaced state (`ocean`, `land`, `tracers`) is a `NamedTuple` and is skipped by the comparison
loop; it is only reachable through those five, so a difference there would have shown up here, but
it is not directly asserted.

### W2.0 — the feedback stride degenerates to *every step* exactly when nothing is printed (2026-08-08)

`progress!` (`SpeedyWeather/src/output/feedback.jl:90-99`) *already* strides the two device
reductions:

```julia
every_nsteps = feedback.progress_meter.core.check_iterations
feedback.show_umax && mod(counter, every_nsteps) == 0 && max_speed(vars)
feedback.show_temperature_range && mod(counter, every_nsteps) == 0 && temperature_range(vars)
```

`check_iterations` is ProgressMeter's adaptive "how often is it worth checking the clock"
counter. It is only ever updated inside `updateProgress!`, whose **first line is
`!p.enabled && return`** (ProgressMeter 1.11.0, `src/ProgressMeter.jl:247`). `enabled` is
`feedback.verbose`. So with `verbose = false` the counter never leaves its initial value `1`,
`mod(counter, 1) == 0` is always true, and both `extrema` reductions — each a device→host
sync — run on **every** time step. Measured directly, independent of SpeedyWeather:

| `Progress(…; enabled)` | `check_iterations` after 20 000 `next!` |
|---|---:|
| `false` (`verbose = false`, the default in scripts and benchmarks) | **1** — i.e. every step |
| `true`  (interactive REPL) | 32681 — i.e. effectively never |

> **Read the second row with care** (see W2.2): 32681 is what a *bare* `next!` loop converges to,
> because `check_iterations` tracks `feedback_dt / step_time`. A loop that also integrates the
> atmosphere converges to 3–10 instead. The first row is the finding; the second is an artefact of
> how it was measured, and building on it as if it described a model run was an error.

This inverts the intent: the diagnostics are cheap-by-striding when they are displayed, and
un-strided when they are *not*. It explains the whole Phase A3 gap (10–26 % of wallclock for
0.05 ms of GPU work) without appealing to any deliberate design choice, and it means the first
part of W2 is a **bug fix, not a behavioural change** — no decision needed for it:

- **W2a (no decision needed):** do not compute `max_speed`/`temperature_range` when the progress
  meter is disabled. Their only consumer is `progress_string` via `ProgressMeter.speedstring`,
  which is only reached when printing, so nothing observable is lost. Gate on
  `feedback.progress_meter.core.enabled` (equivalently `verbose`) rather than on
  `check_iterations`, and keep the existing stride for the verbose case.
- **W2b (decision still open):** `nan_detection!` is *not* strided at all — it runs every step
  and its `@allowscalar vor[2, end]` read is a full device sync that stops the host running
  ahead, which is the specific thing a launch-bound step cannot afford. Striding it delays the
  early-abort in `time_stepping!` (`time_integration.jl:39,74`) by up to N−1 steps; a NaN
  persists, so it is still detected, just later. Options 1/3 of the original W2 list apply here
  only.

### W2.1 — the implemented stride (2026-08-08)

`Feedback` gains one option:

```julia
"[OPTION] Interval in time steps between NaN checks and progress meter diagnostics ..."
interval::Int = 50
```

and `progress!` becomes (the display stride as **revised** by W2.2):

```julia
(; counter, n, check_iterations, enabled) = feedback.progress_meter
interval = max(1, feedback.interval)

every_nsteps = enabled ? check_iterations : interval
if mod(counter, every_nsteps) == 0
    feedback.show_umax && max_speed(vars)
    feedback.show_temperature_range && temperature_range(vars)
end

progress!(feedback)

last_step = counter == n - 1
feedback.debug && (mod(counter, interval) == 0 || last_step) &&
    nan_detection!(feedback, vars, model)
```

Three deliberate asymmetries:

- **`enabled ? check_iterations : interval` for the display quantities, `interval` alone for the NaN
  check.** `check_iterations` answers "how often can the meter possibly *show* a new value", which
  is the right bound for `max_speed`/`temperature_range` and meaningless for NaN detection. Branching
  on `enabled` keeps the enabled case exactly as it was — and bounds its cost in *wallclock*
  (one synchronization per displayed frame, i.e. per `feedback_dt`, at any resolution) — while the
  disabled case, where `check_iterations` is frozen at W2.0's `1`, falls back to the step-count
  bound `interval`.
- **`max(1, interval)`** so `interval = 0` means "every step" instead of a `DivideError` in `mod`.
- **`|| last_step`** so the final state of every run is NaN-checked regardless of length. Without
  it, `interval = 50` would reduce a 10-step test run to a single check at step 0 — i.e. a check of
  the initial conditions only — and the ~30 `nans_detected == false` assertions across the suite
  would become nearly vacuous. Cost is one extra sync per `run!`, not per step.

`counter` runs `0 … n-1` inside `progress!` (it is incremented by `ProgressMeter.next!`, which
does so whether or not the meter is enabled), so step 0 and the last step are both checked and
`n` is the total step count from `clock.n_steps`.

### W2.2 — `check_iterations` in a *model* run is 3–10, not 32681 (2026-08-08)

`calc_check_iterations` (ProgressMeter 1.11.0, `src/ProgressMeter.jl:246`) is

```julia
iterations_per_dt = (p.check_iterations / (t - p.tlast)) * p.dt
round(Int, clamp(iterations_per_dt, 1, p.check_iterations * 10))
```

whose fixed point is `dt / step_time` — the number of iterations that fit in one displayed frame.
It is therefore **a function of how slow the loop body is**, and a bare `next!` loop (3 µs/iteration
→ 32681) says nothing about a loop that integrates the atmosphere. Measured on a real enabled
`PrimitiveWetModel`, 300 steps, CPU:

| run | `check_iterations` after 300 steps |
|---|---:|
| T31, `nlayers = 8`, `verbose = true` | **10** |
| T63, `nlayers = 8`, `verbose = true` | **3** |
| any resolution, `verbose = false` | 1 (never updated, W2.0) |

Consequence for W2.1 as first written: `max(interval, check_iterations)` = `max(50, 3…10)` = `50`,
so the enabled meter would have kept printing at `feedback_dt` = 0.1 s while the max speed and
temperature range inside it changed only every 50 steps — 5–15× staler than the bar carrying them,
which is precisely the behaviour change the `max` was introduced to avoid. Branching on `enabled`
instead gives each rule the regime it was derived for, and makes the enabled-case cost
*wallclock*-bounded (one sync per displayed frame) rather than step-bounded.

The companion question — whether `progress!` (`ProgressMeter.next!`) should also be strided — is
answered no: **4.6 ns/call** with the meter disabled (`updateProgress!` returns on its first line,
so no `time()`, no device access, and `lock_if_threading` does not lock single-threaded), and
`p.counter` is the clock that every `mod(counter, …)` above is taken against, so striding `next!`
would stride the strides.

### W2.3 — the progress bar was never drawn: `barlen = 0`, not `nothing` (2026-08-08)

Reported by the user mid-session while running a model on CPU: the progress meter printed its text
(`100% Time: 0:00:10 (2000-03-01, 703.93 years/day, 88 m/s, [-81, 26] ˚C)`) but no bar.

Not caused by W2 — `main` has the identical code. `Feedback.progress_bar_length` was declared
`::Int = 0` while its own docstring said "nothing = full window width". ProgressMeter's `barlen` is
`Union{Int, Nothing}`, where `nothing` means "fit to the terminal" and an integer is taken literally;
`barstring` draws nothing for `barlen == 0`. So the field could never express the behaviour it
documented, and every run got a zero-length bar. Fixed by retyping it `Union{Int, Nothing}` with
default `nothing`, passed straight through to ProgressMeter.

One consequence to be aware of: ProgressMeter's automatic width reserves 14 characters for the speed
field, whereas SpeedyWeather's `progress_string` is 54 characters in a typical run and up to 60 in
the worst case (measured: date 10, speed 22, maximum wind 10, temperature range 18). On a narrow
terminal the line will therefore wrap. An earlier draft of this fix sized the bar itself from those
four measured widths; it was **removed at the user's request** in favour of the plain pass-through.
Setting `progress_bar_length` explicitly remains the escape hatch, and `0` still means "no bar".

Unrelated but noticed while reading the same struct: `Feedback.show_time` is declared and documented
but never read anywhere in `src/`. Left alone, recorded here.

### W1W2.1 — combined effect of W1 + W2 on the step time (2026-08-08)

Phase A (`bench_model.jl`, `nlayers = 8`, `Float32`, no output, every arm passing the
`nans_detected` guard). H100 = SLURM job `1732490` → `reports/phaseA_gpu_w1w2.json`; A40 =
login node, T31/T63 only → `reports/phaseA_gpu_w2b_a40.json`. Baselines are the pre-plan Phase A
runs (`phaseA_gpu.json`, `phaseA_gpu_a40.json`) measured with the same harness:

| trunc | H100 before | H100 after | speedup | A40 before | A40 after | speedup |
|------:|------------:|-----------:|--------:|-----------:|----------:|--------:|
| 31  | 2.760  | **1.126**  | **2.45×** | 3.590 | **1.277** | **2.81×** |
| 63  | 5.120  | **2.138**  | **2.39×** | 6.256 | **2.386** | **2.62×** |
| 127 | 10.485 | **4.777**  | **2.20×** | — | — | — |
| 255 | 23.833 | **14.280** | **1.67×** | — | — | — |

W2 is verified *within* the run rather than by dividing two harnesses: the `nofeedback` variant is
the same model with `feedback = nothing`, so the gap between it and `default` **is** the feedback
cost, and it is what W2 was aimed at:

| trunc | feedback cost before (H100) | after (H100) | before (A40) | after (A40) |
|------:|----------------------------:|-------------:|-------------:|------------:|
| 31  | 16.1 % | **3.0 %**  | 10.5 % | **1.0 %** |
| 63  | 26.1 % | **−0.2 %** | 19.2 % | **0.9 %** |
| 127 | 24.9 % | **−0.3 %** | — | — |
| 255 | 14.0 % | **1.8 %**  | — | — |

Negative entries are `nofeedback` measuring *slower* than `default` — i.e. the remaining cost is
below the run-to-run noise of the harness, which is the intended end state. The stride therefore
recovered essentially all of the Phase A3 gap, not merely part of it.

The `dynamics_only` variant tells the complementary story: the whole parameterization suite costs
1.8–10.2 % of the post-fix step (was 1.2–8.2 % of a step twice as long), so in absolute terms it is
unchanged — W1 and W2 removed overhead, not physics.

Note the shape of the speedup: **largest at T31, smallest at T255**, the same gradient as W1 alone.
Both fixes remove per-step *overhead* whose size grows more slowly than the real work, so their
share shrinks as resolution grows. T255's 1.67× is the honest production figure.

### W3.0 — the post-W1+W2 re-profile: every per-step host count is now resolution-independent (2026-08-08)

Phase B re-run on the H100 with the fix in (SLURM job `1732491`, traces `reports/w1_T*_L8`, the
pre-W1 `model_T*` set preserved), exported with `analyze_nsys.sh w1_T` and tabulated by
`phase_table.jl`. 20 captured steady-state steps per truncation. Host-side CUDA API calls **per
step**:

| API call | T31 → | T63 → | T127 → | T255 → | after (all truncations) |
|---|---:|---:|---:|---:|---:|
| kernel launches (`cuLaunchKernel` + `…Ex`) | 301 | 545 | 1073 | 2237 | **97.2** |
| `cuMemcpyDtoDAsync_v2` | 96 | 192 | 384 | 768 | **0** |
| `cuMemAllocFromPoolAsync` | 101 | 197 | 389 | 773 | **1.2** |
| `cuStreamGetCaptureInfo` | 2333 | 3485 | 5789 | 10397 | **1162** |
| `cuStreamSynchronize` | 6 | 6 | 6 | 6 | **0.5** |
| `cuMemcpyDtoHAsync_v2` | 3 | 3 | 3 | 3 | **0.2** |
| `cuGraphLaunch` | 3 | 3 | 3 | 3 | **4** |

Every count that used to scale with `nlat_half` is now flat, and identical at all four truncations —
which is what a launch sequence determined by the *program* rather than by the *grid* looks like.
Three consequences:

- **W4 is closed with no separate work.** The 96/step `cuMemcpyDtoDAsync_v2` at T31 (768 at T255) are
  at exactly zero. The W4 update's prediction — that they were the serial inverse path's
  per-iteration `dest .= rhs` temporaries and not an independent phenomenon — is confirmed. The
  per-step device allocations (`cuMemAllocFromPoolAsync`, 773/step at T255) vanished with them,
  which is the allocation half of the same line and the reason W1 beat its own launch-count estimate.
- **The two per-step host round-trips of W2 are visible as such**: `cuMemcpyDtoHAsync_v2` 3 → 0.2 and
  `cuStreamSynchronize` 6 → 0.5 per step, i.e. one every other step rather than six every step.
- **`cuStreamGetCaptureInfo` is now the most frequent CUDA call in the step** at 1162/step — twelve
  per kernel launch, CUDA.jl checking whether it is inside a capture. At 73 ns each it is only
  0.085 ms/step, so it is not worth attacking on its own, but it scales with launch count and is a
  reason to prefer *fewer, larger* launches over cheaper ones in W3.

**Caveat that must be carried into W3 — these traces cannot see inside CUDA graphs.**
`nsys` defaults to `--cuda-graph-trace=graph`, which records a graph as a single GPU operation and
**does not collect the kernels inside it**. Pre-fix, the K=2 FFT ran on the serial path and every one
of its kernels was counted; post-fix it runs batched inside a graph and none of them are. Hence:

- the post-fix kernel table contains **no FFT kernels at all** (pre-fix: 964/step at T255, 2.695 ms),
  and the "69 distinct kernels" is 69 *non-graph* kernels;
- the post-fix "GPU busy" of 0.267 ms/step (T31) and 2.797 ms/step (T255) **excludes all graph work**
  and is therefore a lower bound. The pre-fix 0.619/7.020 ms included the serial FFT, so
  **the pre/post GPU-busy and kernels-per-step numbers are not comparable** and the apparent fall in
  busy fraction (29 % → 19 % at T255) is largely an artefact of this. The host-side launch counts
  above are unaffected and are the trustworthy half of this trace set.

`submit_phaseB.sh` now takes the granularity as a second argument, and the re-profile at
`node` granularity ran as job `1732864` (prefix `w1g_T`, so `w1_T` and `model_T` both survive).
That is the trace set W3 must be sized from — see finding W3.1, which supersedes everything
GPU-side in this one.

The API totals here invite a wrong conclusion, and W3.1 corrects it:
at T255 kernel launches cost 3.60 ms/step and graph launches 2.02 ms/step, which leaves ~8 ms/step
of the 14.5 ms step outside the CUDA API entirely. That looks like host-side Julia time. It is not —
`cuGraphLaunch` blocks while the graph's GPU work runs, and this trace cannot see that work. Do not
size W3 from this paragraph; use W3.1.

### W3.1 — with graph nodes expanded, the step is 66–78 % GPU-busy: the launch-bound era is over (2026-08-08)

Re-profile at `--cuda-graph-trace=node` (SLURM job `1732864`, traces `reports/w1g_T*_L8`), which
expands each CUDA graph into its constituent kernels and so counts the FFT work the previous trace
set dropped. Against the *native* Phase A step time (`phaseA_gpu_w1w2.json`, same H100, no profiler):

| trunc | GPU busy ms/step | native ms/step | GPU busy % | kernels/step (incl. graph-internal) | host launches/step |
|------:|-----------------:|---------------:|-----------:|------------------------------------:|-------------------:|
| 31  | 0.745  | 1.126  | **66 %** | 319  | 97.2 |
| 63  | 1.569  | 2.138  | **73 %** | 615  | 97.2 |
| 127 | 3.735  | 4.777  | **78 %** | 1283 | 97.2 |
| 255 | 10.465 | 14.280 | **73 %** | 2843 | 97.2 |

against the 22–29 % that opened this plan. **The premise of W3 as written no longer holds.** The step
is not launch-bound any more: 97 host launches per step drive 319–2843 kernels, because 96 % of the
kernels now run inside four replayed graphs, and the GPU is busy for roughly three quarters of the
step. The remaining idle is 22–34 %, not the 71–78 % the plan was scoped against, and the ceiling
from removing *all* of it has fallen from ~3.4× to ~1.36× at T255.

Where the GPU time actually goes now (T255, of 10.465 ms busy):

| phase | busy ms/step | share | kernels/step |
|---|---:|---:|---:|
| `dynamics` (includes the transforms it calls) | 5.574 | 53 % | 1414 |
| `transform` | 4.552 | 44 % | 1387 |
| `parameterizations` | 0.213 | 2 % | 7 |
| everything else | 0.126 | 1 % | 35 |

so the spectral transform *is* the model at this resolution. Inside it the cost is not one big kernel
but many small cuFFT ones — `preprocess_kernel` 606/step (0.836 ms), `regular_fft_factor<2>` 224/step
(0.380 ms), `postprocess_kernel` 196/step (0.284 ms) — plus three Bluestein variants (lengths 729,
625, 1024; 0.455 ms together), which is cuFFT's fallback for lengths it cannot factor cheaply. The
two Legendre kernels are 1.231 ms (forward) and 0.739 ms (inverse). This is *GPU work*, not launch
overhead: it is already inside graphs, so making the launches cheaper cannot touch it.

**Consequence for W3: neither (a) nor (b) as scoped.** Extending graph capture (a) attacks a launch
cost that is now 0.30 ms/step of `cuLaunchKernelEx` at T255 — 2 % of the step. Elementwise fusion (b)
targets the `gpu_broadcast_kernel_*` family, which no longer appears anywhere near the top of the
kernel table. The two candidates that the measurement actually points at are:

- **the FFT itself** — 44 % of GPU busy at T255, spent in hundreds of ~1.4 µs cuFFT kernels per step.
  Worth checking whether the ring lengths driving the Bluestein path can be avoided (grid choice, or
  padding the transform length to a factorable one), and whether the per-ring plans can be merged
  further so cuFFT issues fewer, larger kernels;
- **the residual 22–34 % idle**, which is now small enough that it needs attribution before anything
  is built: at T31 (66 % busy, the worst case) it is 0.38 ms/step, and the profiler's own inflation
  is of the same order, so the first job is to establish that it is real.

Caveat on this trace set: `node` granularity has "significant runtime overhead" by nsys's own
documentation, so the busy *totals* may be slightly inflated. They are consistent with the
`graph`-granularity set on everything both can see (identical host launch counts, identical non-graph
kernel counts, `dynamics` busy 1.720 → 5.574 ms explained entirely by the graph-internal FFTs), so
the direction is not in doubt even if the last few percent are.

## Testing and verification

- **Correctness:** the modelbench bitwise equivalence check (driver vs `run!`) after every
  source change; full GPU test suite (`SpeedyWeather/test/GPU/runtests.jl`); for W1
  additionally a CPU-vs-GPU and batched-vs-serial tolerance comparison of a short run, since
  batched cuFFT may not be bitwise identical to the serial path. **Done, and the tolerance was
  not needed — the two paths agree bitwise; see finding W1.1.**
- **No silent regressions:** the W1 batch-plan test pins the fixed behaviour;
  W3 graph work adds a launches-per-step assertion under the modelbench driver.
- **Performance:** before/after Phase A tables per workstream, both GPUs, with the ±20 %
  benchmark-noise convention from `CLAUDE.md`; report anything outside it.
- W2 needs a unit test that `nans_detected` still triggers (inject a NaN, step N+1 times,
  assert detection) under the strided scheme. **Done** — `"NaN detection with strided feedback"`
  in `SpeedyWeather/test/output/feedback.jl`, 5/5 passing on CPU. It covers both branches of the
  new condition (the strided check and the last-step check) plus a healthy run, because only the
  last-step branch protects the many short-run `nans_detected == false` assertions elsewhere in
  the suite.

### Status as of 2026-08-08

| check | status |
|---|---|
| `SpeedyWeather/test/GPU/fft_batch_plans.jl` (new, A40) | **8/8 passed** as part of the full GPU suite below |
| A/B timing, A40, T31–T255 | **done**, table above; every arm passed the `nans_detected` guard |
| A/B timing, H100, T31–T255 | **done**, SLURM job `1731194` → `reports/verify_k2-1731194.log` |
| full `SpeedyWeather/test/GPU/runtests.jl` (A40) | **passed**, exit 0, no failures across all 26 testsets |
| CPU suite `SpeedyTransforms` | **passed**, exit 0 |
| CPU suite `SpeedyWeather` | **one pre-existing failure, unrelated**: `dynamics/dispatch.jl` (the JET no-runtime-dispatch guard) fails, all other 50 testsets pass. Bisected by reverting *only* `feedback.jl` to the pre-W2 revision in the same environment — it still fails, and the JET report names only KernelAbstractions' CPU `__run`/`__thread_run` and `ModelParameters.withunits`, nothing in SpeedyWeather. It passes in a freshly-created worktree, so it is dependency-version-sensitive; `SpeedyWeather/test/Manifest.toml` is gitignored, so different checkouts of the same source disagree. Not caused by this plan; worth its own issue |
| batched-vs-serial comparison | **done and stronger than required — bitwise identical**, finding W1.1 |
| Phase A after W1+W2 (H100, T31–T255; A40, T31/T63) | **done**, finding W1W2.1 |
| Phase B/C re-profile after W1 | **done twice**: `graph` granularity (job `1732491`, finding W3.0) and `node` granularity (job `1732864`, finding W3.1). Only the second can see inside CUDA graphs |
| W2 source change | **implemented**, findings W2.1/W2.2, CPU feedback tests 11/11 |
| progress bar fix | **implemented**, finding W2.3, covered by the `"Progress bar length"` testset |

### Files changed by W1

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
- `SpeedyWeather/test/GPU/modelbench/verify_batch_equivalence.jl` (new, 2026-08-08) — the
  batched-vs-serial tolerance comparison. **Untested: written but never run.**

### Files changed by W2

- `SpeedyWeather/src/output/feedback.jl` — new `Feedback.interval` option (default `50`) and the
  strided `progress!` shown in finding W2.1.
- `SpeedyWeather/test/output/feedback.jl` — new `"NaN detection with strided feedback"` testset
  with a `NaNInjector` callback: a NaN written mid-run is still detected under `interval = 4`; a
  NaN in a run shorter than `interval` is caught by the last-step check; a healthy run is not
  flagged. Passes (5/5) on CPU.
- `SpeedyWeather/benchmark/README.md` — note that published numbers include feedback overhead, and
  that pre-stride GPU numbers paid it every step.
- `CHANGELOG.md` — bullet added under `## Unreleased`, PR [#1185].

### Files changed by the progress-bar fix (finding W2.3)

- `SpeedyWeather/src/output/feedback.jl` — `progress_bar_length` retyped
  `Union{Int, Nothing}` and defaulted to `nothing`.
- `SpeedyWeather/test/output/feedback.jl` — `"Progress bar length"` testset.
- `CHANGELOG.md` — bullet under `## Unreleased`, PR [#1185].

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

W0, W1, W2 and W4 are done and measured; the source changes are committed on
`mg/profile-gpu-primitivewet` and carry PR number [#1185]. What is left:

1. **W3 needs a fresh decision, because finding W3.1 removed its premise.** The step is 66–78 %
   GPU-busy, so the whole remaining idle is worth at most ~1.36× at T255 and ~1.5× at T31, and both
   remedies the plan proposed (more CUDA graphs / elementwise fusion) aim at launch traffic that is
   already down to 97 host launches per step. Where the time now is, at T255: cuFFT 44 % of GPU busy
   spread over ~1400 small kernels per step, the two Legendre kernels 19 %, parameterizations 2 %.
   Before building anything:
   - **check the busy fraction without a profiler.** `node` granularity carries "significant runtime
     overhead" by nsys's own documentation, and the whole re-scoping rests on that one number. CUDA
     events around the step in `nsys_target.jl`, or simply comparing the summed kernel time against
     the native Phase A step, would confirm it independently and cheaply;
   - then pick a target *inside* the transform rather than around it. The concrete lead is the
     Bluestein fallback (lengths 729, 625, 1024 at T255, 0.455 ms/step together): cuFFT uses it for
     lengths it cannot factor cheaply, and the lengths are the reduced grid's ring lengths, so a grid
     or padding choice may avoid it. The `preprocess`/`postprocess` kernel pairs (802 launches,
     1.12 ms/step) are the other visible overhead of issuing many short transforms.
2. **Finish the A40 arm.** Post-fix Phase A on the A40 covers only T31 and T63
   (`reports/phaseA_gpu_w2b_a40.json`). T127 and T255 are one login-node command and make the
   two-GPU comparison complete.
3. **Re-run the GPU suite.** `SpeedyWeather/test/GPU/runtests.jl` last passed on the A40 *before* W2
   and the progress-bar fix landed. Both touch `Feedback`, which every GPU run exercises.
4. **Housekeeping before the PR:** re-check the `-DEV`/`+DEV` tags on `SpeedyWeather`
   (`0.21.1+DEV`) — `Feedback.progress_bar_length` changed type from `Int` to
   `Union{Int, Nothing}`, which is a (small) public API change and may argue for a minor bump
   rather than `+DEV`. `SpeedyTransforms` is untouched by the final form of W1.
5. **Not ours, but file it:** `SpeedyWeather/test/dynamics/dispatch.jl` fails in this checkout and
   passes in a fresh worktree of the same source. `SpeedyWeather/test/Manifest.toml` is gitignored,
   so the JET guard silently depends on whichever KernelAbstractions/ModelParameters versions a
   given checkout happens to have resolved. Either pin the manifest for that test or relax the
   guard to SpeedyWeather-owned modules.

Environment notes for whoever picks this up:

- the login-node A40 is shared, and a stale 2.6 GB Julia process from an earlier session was
  cleared on 2026-08-08 — the W1 A40 numbers were taken under that contention, later ones are not,
  so do not read a small A40 shift between sessions as a code effect. Check
  `nvidia-smi --query-compute-apps=pid,used_memory --format=csv` for strays before trusting any
  login-node timing; other users' processes are not yours to kill.
- `submit_phaseB.sh` takes `[prefix] [graph_granularity]`. Existing trace sets:
  `model_T*` (pre-W1, H100), `a40_T*` (pre-W1, A40), `w1_T*` (post-W1+W2, H100, `graph`),
  `w1g_T*` (post-W1+W2, H100, `node`). Use a new prefix; the script passes
  `--force-overwrite=true`.
- **use `--cuda-graph-trace=node` for any question about GPU work.** At the default `graph`
  granularity the kernels inside a CUDA graph are not collected at all, which cost this
  investigation one wrong conclusion (see the correction at the end of finding W3.0).
