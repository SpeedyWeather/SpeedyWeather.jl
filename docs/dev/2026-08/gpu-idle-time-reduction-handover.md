# GPU idle time reduction — continuation handover

> Status: **in progress**. M1 is complete and did what it was written to do: it reversed the priority
> order. Both readings of the W3.1 phase table were wrong — cuFFT is **69 % of GPU busy at T255**
> (not 44 %, not 19 %), and non-transform dynamics is **7.9 %**. See finding **M1.1**. L1 and L4
> together address at most ~25 % of busy; the per-ring FFT loop is the target.
>
> **Session of 2026-08-14 ended here. Resume at *State and next actions* at the bottom.** M2, M3 and
> all implementation leads are not started; the GPU environment is now working (finding E1.1) and the
> "NaN at step 1" scare was a harness-usage artifact, fully explained and closed (finding E1.2).

Date of initial draft: 2026-08-14

Base revision: `c343adf6b07270b74b1d8505e1e87302a1c86b2f` (`main`, v0.22.0)

## Originating prompt

> I want to continue the GPU idle time reduction investigation. Where do you see further potential?

> based on these findings write a handover file, I can execute on a GPU node. First the "What i
> would do first, in order". Then from the leads include 1) and 4), and 5) as a maybe

## Revision log

- **2026-08-14, initial draft.** Written against `main` after re-deriving the numbers in finding
  W3.1 of [`gpu-idle-time-reduction.md`](gpu-idle-time-reduction.md). No measurements were taken for
  this draft; every figure below is arithmetic on that plan's published tables, and is labelled as
  such. Two corrections to the parent plan's *interpretation* motivated it (see *Problem
  description*): the phase table's rows are disjoint, not nested, and "GPU busy" is not utilization.
- **2026-08-14, session 2 (execution).** Ran M1 on the existing `w1g_T*` traces. It confirmed the
  disjoint-rows claim but **refuted both of the *Problem description*'s FFT shares** (44 % and
  18.7 %): the true figure is 69 % at T255, because the W3.1 kernel table is truncated to 15 rows
  and hides a long tail of per-ring-length cuFFT variants. Findings M1.1, E1.1, E1.2 added; parent
  plan corrected; *State and next actions* appended. M2–M4 and all leads remain unstarted. The
  *Problem description* above is left as originally written, with its errors, so the correction in
  M1.1 is legible — **do not read its two FFT percentages as current.**

## Problem description

The parent plan closed W1, W2 and W4 and re-scoped W3 on the basis of finding W3.1: with CUDA graph
nodes expanded, the `PrimitiveWetModel` step is 66–78 % GPU-busy, so removing *all* remaining idle is
worth at most ~1.36× at T255. That conclusion is correct as stated but leads to the wrong next step,
for two reasons.

**(1) The phase attribution is ambiguous, and the two readings disagree by a factor of five.** The
W3.1 phase table is headed "`dynamics` (includes the transforms it calls)", but its four rows sum to
exactly the busy total:

```
5.574 (dynamics) + 4.552 (transform) + 0.213 (parameterizations) + 0.126 (else) = 10.465 ms = GPU busy
```

Disjoint rows summing exactly to the total is not what nesting produces. If the rows really are
disjoint, then **non-transform dynamics is 53 % of GPU busy at T255 across 1414 kernels/step** and
is completely unexamined, while `transform` is 44 %. The parent plan's *Where to continue* section
instead states "cuFFT 44 % of GPU busy", which the same trace's kernel table contradicts: cuFFT is
0.836 + 0.380 + 0.284 + 0.455 = **1.955 ms = 18.7 %** of busy, and the two Legendre kernels a
further 1.970 ms = 18.8 %. Which reading is right decides whether the target is inside `transform!`
or outside it. M1 settles it from traces that already exist.

**(2) "GPU busy" is not "GPU utilized", so the ~1.36× ceiling is not the real ceiling.** Per step the
model issues **8 · nlat_half cuFFT invocations** — 4 transforms × 2 hemispheres × nlat_half rings,
the loop at [fourier.jl:228-243](../../../SpeedyTransforms/src/fourier.jl#L228-L243) — i.e. 1536 at
T255 and 192 at T31. Each moves one ring × K layers:

| quantity (T255, H100) | value | source |
|---|---:|---|
| mean cuFFT kernel duration | ~1.4 µs | W3.1 (0.836 ms / 606 `preprocess_kernel`) |
| bytes per invocation (nlon ≈ 400, K = 33, Float32) | ~50–100 KB | derived |
| time for that traffic at HBM peak (3.35 TB/s) | ~0.05 µs | derived |
| total FFT traffic per step, all four transforms, read+write | ~0.3 GB | derived |
| ⇒ FFT time at peak bandwidth | ~0.09 ms | derived |
| measured cuFFT time per step | 1.955 ms | W3.1 |

The FFT therefore runs at roughly **5 % of achievable bandwidth**; ≳95 % of each kernel's duration is
fixed dispatch cost. This is the same pathology W1 fixed, one level up: W1 removed the per-**layer**
loop, the per-**ring** loop is still there. A step that is 73 % *busy* at 5 % *efficiency* has more
headroom than the idle-gap arithmetic suggests, and the headroom is reached by making the busy time
productive rather than by closing the idle gap.

A third point, not a correction but a scoping gap: the re-scoping arithmetic was done at T255 only.
Host-side cost is resolution-*independent* by finding W3.0 (97.2 launches/step at every truncation,
0.30 ms of `cuLaunchKernelEx` + 0.085 ms of `cuStreamGetCaptureInfo`). At T31 that fixed 0.385 ms is
**34 % of the 1.126 ms step**, against 2 % at T255. Conclusions about launch overhead do not transfer
downward in resolution, and low resolution is where ensembles run.

## Background

- Parent plan and all measurements: [`gpu-idle-time-reduction.md`](gpu-idle-time-reduction.md).
  Methodology and harness: [`gpu-primitive-wet-model-profiling.md`](gpu-primitive-wet-model-profiling.md).
- Legendre kernel history, including the forward-transform L2 mechanism L4 builds on:
  [`legendre-gpu-optimization.md`](legendre-gpu-optimization.md).
- ~~**Harness availability — check this before anything else.**~~ **Resolved 2026-08-14: the harness
  is present** at `SpeedyWeather/test/GPU/modelbench/` (untracked in the working tree, hence the
  original doubt), complete with `reports/` and all four `w1g_T*` traces. The paragraph below is
  **stale — ignore it.** ~~The modelbench harness
  (`bench_model.jl`, `nsys_target.jl`, `analyze_nsys.sh`, `phase_table.jl`, `submit_phaseB.sh`,
  `reports/`) is not present in `main` under `SpeedyWeather/test/GPU/modelbench/`; it lived on
  `mg/profile-gpu-primitivewet`, which is not in this clone either. Locate it on the cluster or
  restore the branch before starting.~~ Every task below assumes it. Existing trace sets referenced:
  `model_T*` (pre-W1, H100), `a40_T*` (pre-W1, A40), `w1_T*` (post-W1+W2, H100, `graph`
  granularity), `w1g_T*` (post-W1+W2, H100, `node` granularity).
- **Always use `--cuda-graph-trace=node`** for any question about GPU work; at the default `graph`
  granularity the kernels inside a CUDA graph are not collected at all. This cost the parent
  investigation one wrong conclusion.
- Hardware: iterate on the login-node A40 (shared — check
  `nvidia-smi --query-compute-apps=pid,used_memory --format=csv` for strays first), confirm on the
  SLURM H100s (`--partition=gpu --qos=gpushort --account=flai --gres=gpu:1`). Record the GPU for
  every number.

## Summary of changes

Two blocks. **M1–M4 are measurements and must be run in order** — each can change what follows.
**L1, L4, L5 are implementation leads**, to be started only after M1–M2 have confirmed the target.

### M1 — resolve the phase-attribution ambiguity (no new job, ~30 min)

Re-tabulate the existing `w1g_T255_L8` trace as *kernel name × NVTX range*, rather than the
per-phase totals `phase_table.jl` currently prints. The question to answer is a single yes/no:

> Are the cuFFT and Legendre kernels counted inside the `dynamics` row, or only inside the
> `transform` row?

Method: the sqlite recipe from the profiling plan's C2 section — join `CUPTI_ACTIVITY_KIND_KERNEL`
against the NVTX range table and group by `(range, kernel name)`. Sanity checks that make the answer
unambiguous without trusting the join:

- the `transform` row has 1387 kernels/step at T255 and `dynamics` has 1414. If the four transforms
  contribute ~1536 cuFFT invocations plus 4 Legendre launches, they cannot all fit in the 1387 row,
  which would already imply that the two transforms called from inside `dynamics_tendencies!`
  ([tendencies.jl:157](../../../SpeedyWeather/src/dynamics/tendencies.jl#L157) and
  [tendencies.jl:193](../../../SpeedyWeather/src/dynamics/tendencies.jl#L193)) are attributed to
  `dynamics`, and the two called from
  [transform.jl:151](../../../SpeedyWeather/src/time_stepping/transform.jl#L151) and
  [transform.jl:162](../../../SpeedyWeather/src/time_stepping/transform.jl#L162) to `transform`;
- count how many kernels in the `dynamics` row are *not* cuFFT/Legendre. That number is the true
  size of the dynamical core outside the spectral transform.

**Exit criterion / decision gate:**

- If non-transform dynamics is ≳2000 kernels and ≳3 ms/step at T255 → **stop and re-plan**: that is
  a larger target than L1 and L4 combined and nothing in the parent plan has looked at it.
- If non-transform dynamics is a few hundred kernels and ≲1 ms/step → the transform is confirmed as
  the target and L1/L4 are correctly prioritised.

Record the resulting kernel × range table in this document as finding **M1.1**, and correct the
parent plan's W3.1 table header and its "cuFFT 44 %" sentence in the same edit.

### M2 — confirm the busy fraction without nsys, and measure host slack (1 login-node session)

The whole re-scoping rests on one `node`-granularity number that nsys itself documents as carrying
"significant runtime overhead". Two independent, cheap checks:

**M2a — CUDA.jl's own CUPTI profiler.** Run ~20 steady-state steps under
`CUDA.@profile trace=true` and sum the device-side kernel durations, then divide by the profiler-free
step time from `bench_model.jl`. First verify the FFT and Legendre kernels actually appear in the
table — if they do not, CUDA.jl is not expanding graph nodes either and this check is void, in which
case fall back to M2b alone.

**M2b — host-slack probe (no profiler at all).** Insert a host-side busy-wait of `x` µs per time step
(a spin loop, not `sleep`, to avoid scheduler granularity) and sweep `x ∈ {0, 50, 100, 200, 400}` µs
at T31 and T255. Then:

- if step time rises by ~`x` immediately, the host is on the critical path and the *idle* is host-bound;
- if step time is flat until some `x₀` and rises linearly after, `x₀` **is** the host slack per step —
  a direct, profiler-free measurement of how far the host runs ahead of the device.

This is the number that decides whether the residual 22–34 % idle is worth attacking at all, and at
which resolution.

**Exit criterion:** a busy-fraction figure from at least one method that is independent of the
`w1g_T*` traces, plus a per-resolution host-slack value. Record as finding **M2.1**. If M2a and W3.1
disagree by more than ~10 percentage points, the W3.1 re-scoping must be revisited before L1 starts.

### M3 — `FullGaussianGrid` sweep (zero code, one login-node command + one H100 job)

The per-ring FFT loop exists because `OctahedralGaussianGrid` (the default,
[spectral_grid.jl:9](../../../SpeedyWeather/src/dynamics/spectral_grid.jl#L9)) has a different `nlon`
per ring, and cuFFT bakes the transform length into the plan. On a full grid all rings share `nlon`,
so a single `plan_many` over `(nlon, nlat·K)` could replace 2·nlat_half invocations — a ~380×
reduction in invocations at T255 — at the cost of 1.91× more grid points (154 368 → 294 912 at T255).

That specialised path is *not* what M3 measures. M3 measures the **penalty side only**, today,
with no source change, so that the payoff side is bounded before anyone writes the code:

```julia
spectral_grid = SpectralGrid(truncation=..., nlayers=8, architecture=SpeedyWeather.GPU(),
                             Grid=FullGaussianGrid)     # note: capital-G `Grid` kwarg
```

Run the standard Phase A `bench_model.jl` sweep (T31/T63/T127/T255, `nlayers = 8`, Float32, no
output, `nans_detected` guard) for `FullGaussianGrid` against the `OctahedralGaussianGrid` baseline
in `reports/phaseA_gpu_w1w2.json`.

Interpretation: today both grids pay one cuFFT call per ring, so the full grid should come out
*slower* by roughly the extra grid-point work in the grid-space parts of the step. That measured
slowdown is the price a single-call FFT path would have to beat. Since the FFT is ~19 % of GPU busy
and parameterizations are 2 %, a slowdown of less than ~15 % at T255 makes the specialised full-grid
FFT path clearly worth writing; more than ~40 % probably kills it.

**Exit criterion:** the two-grid Phase A table, recorded as finding **M3.1**, plus a go/no-go on a
follow-up plan for a uniform-`nlon` batched FFT path. Note explicitly that `dealiasing` recomputes
`nlat_half` per grid type, so the two arms are matched on truncation, not on point count.

### M4 — start L1 (below), then re-measure

M4 is not a separate measurement: it is the instruction that L1 is the first *implementation* task,
and that it ends with the same Phase A sweep plus a `node`-granularity re-profile under a new trace
prefix (`submit_phaseB.sh <prefix> node`). Do not start it before M1's decision gate.

---

### L1 — reduce transforms per step from four to two

**What exists today.** One `PrimitiveWetModel` step performs exactly four batched transforms
(finding W0.3 measured this directly):

| # | direction | K (nlayers = 8) | call site | destination |
|---|---|---:|---|---|
| 1 | inverse | 4L+1 = 33 | [transform.jl:151](../../../SpeedyWeather/src/time_stepping/transform.jl#L151) | `vars.fused.grid` ← `vars.fused.prognostic` |
| 2 | inverse | 2L = 16 | [transform.jl:162](../../../SpeedyWeather/src/time_stepping/transform.jl#L162) | `vars.fused.uv_grid` ← `vars.fused.spectral_scratch`, `unscale_coslat = true` |
| 3 | inverse | 2 | [tendencies.jl:193](../../../SpeedyWeather/src/dynamics/tendencies.jl#L193) | `vars.fused.dpres_grad` ← surface-pressure gradient pair |
| 4 | forward | 9L+1 = 73 | [tendencies.jl:157](../../../SpeedyWeather/src/dynamics/tendencies.jl#L157) | spectral tendencies |

Because the per-ring FFT cost is overhead-dominated and nearly independent of `K` (see *Problem
description*), **the FFT cost of a step scales with the number of transforms, not with the number of
layers transformed**. Merging 2 into 1 removes a quarter of the ring invocations; merging 3 as well
removes half.

There is already an explicit TODO for the first merge at
[transform.jl:158-164](../../../SpeedyWeather/src/time_stepping/transform.jl#L158-L164), which names
two blockers. Both are real, and the second is the one that matters.

**L1.1 — `unscale_coslat` on a layer subrange (easy, required by every route).** The unscaling is not
part of the FFT; it is a separate KA kernel applied to the Fourier-space scratch after the inverse
Legendre transform, at [legendre.jl:197-219](../../../SpeedyTransforms/src/legendre.jl#L197-L219) and
called from [legendre_ka.jl:117](../../../SpeedyTransforms/src/legendre_ka.jl#L117). It already takes
an `nlayers` bound ("scale only the layers this transform wrote"). Generalise that to a layer
*offset + length*, or pass a length-`K` device mask and make the kernel
`g[m,k,j] *= mask[k] ? coslat⁻¹[j] : 1`. Thread the new argument through the `transform!` chain:
[spectral_transform.jl:475](../../../SpeedyTransforms/src/spectral_transform.jl#L475),
[:518](../../../SpeedyTransforms/src/spectral_transform.jl#L518),
[:526](../../../SpeedyTransforms/src/spectral_transform.jl#L526),
[:534](../../../SpeedyTransforms/src/spectral_transform.jl#L534), plus the CPU
[legendre.jl:44/94](../../../SpeedyTransforms/src/legendre.jl#L44) and the matrix transform
[matrix_transform.jl:214/238](../../../SpeedyTransforms/src/matrix_transform.jl#L214). Keep the
boolean as the default so no caller changes.

**L1.2 — contiguity, the actual blocker.** A single `transform!` call needs its source coefficients
contiguous in one `LowerTriangularArray` and its destination contiguous in one `Field`. The
prognostic fuse parent is 4-D (it carries the leapfrog step axis: `SpectralXYZT`,
[primitive_dry.jl:145-147](../../../SpeedyWeather/src/models/primitive_dry.jl#L145-L147)) whereas the
scratch parent holding U and V is 3-D (`SpectralXYZ`,
[barotropic.jl:87-88](../../../SpeedyWeather/src/models/barotropic.jl#L87-L88)). That is the
"dimensions don't quite align" in the TODO. Two routes:

- **(iii) staging copies — recommended, low risk, do this first.** Keep the variable system
  untouched. Allocate one staging pair (a `K = 6L+1` spectral buffer and its grid counterpart), copy
  the prognostic step slice and the U/V scratch into it with one kernel, run **one** transform, then
  scatter the grid side back into `vars.fused.grid` and `vars.fused.uv_grid`. Extra traffic at T255:
  ~4 MB spectral + ~20 MB grid, read+write ≈ 50 MB ⇒ ~0.02 ms at HBM peak, against a predicted
  saving of ~0.49 ms. The copies are two extra kernel launches per step, both inside the existing
  graph region. If the measured win matches the prediction, route (i) can be done later purely to
  remove the copies.
- **(i) fuse-family surgery — zero-copy, invasive, only if (iii) proves the win.** The grid side is
  already aligned (`u`, `v` are `GridXYZT(pg)` at
  [barotropic.jl:84-85](../../../SpeedyWeather/src/models/barotropic.jl#L84-L85), the same dims as
  the `:grid` family), so moving them from `fuse = :uv_grid` to `fuse = :grid` is mechanical. The
  spectral side requires giving `scratch.a`/`scratch.b` a leapfrog axis so they can join
  `:prognostic`, which changes every user of `vars.scratch.a/.b` (notably `UV_from_vordiv!` at
  [transform.jl:145](../../../SpeedyWeather/src/time_stepping/transform.jl#L145)) and costs one extra
  copy of those two fields in memory (~4 MB at T255 L8, negligible). Note
  `_assert_fuse_alignment(fused.prognostic, fused.grid)`
  ([variables.jl:260-262](../../../SpeedyWeather/src/variables/variables.jl#L260-L262)) requires the
  member order of the two parents to correspond, so `a↔u`, `b↔v` must be declared in matching order.

`6 * nlayers + 1` is **already in the default GPU batch set**
([spectral_grid.jl:355-356](../../../SpeedyWeather/src/dynamics/spectral_grid.jl#L355-L356)), so the
merged transform gets a batched cuFFT plan with no change to `default_transform_batch`. Assert that
in the test, exactly as `fft_batch_plans.jl` does for `K = 2`.

**L1.3 — fold the `K = 2` surface-pressure gradient into the merged call (optional, do last).** Its
source `vars.fused.dpres_grad_spec` and destination `vars.fused.dpres_grad` are 2-D
([primitive_dry.jl:169-172](../../../SpeedyWeather/src/models/primitive_dry.jl#L169-L172)); a 2-D
variable occupies one slot on the fused layer axis, which the prognostic parent already does for
surface pressure, so the shapes are not the obstacle. The obstacle is *ordering*: the gradient is
currently computed inside `pressure_gradient_flux!` during `dynamics_tendencies!`, i.e. after the big
inverse transform, whereas merging requires computing ∇(ln pₛ) from the current prognostic pressure
*before* it. That input is available at that point, so the reordering is legal, but it is a change to
the dycore's call order and carries more risk than L1.2. Take it only if L1.2's measurement confirms
the model (per-transform, not per-layer, cost).

**Predicted payoff** (arithmetic on W3.1's cuFFT total of 1.955 ms/step at T255; the T31 column
assumes ~200 of that resolution's 319 kernels/step are cuFFT — 8·nlat_half = 192 invocations — at the
measured ~1.4 µs, i.e. ~0.28 ms of the 0.745 ms busy. **Verify that split in M1** before relying on
the T31 column):

| | T255 | T31 |
|---|---:|---:|
| L1.2 alone (4 → 3 transforms) | ~0.49 ms ≈ 3.4 % of step | ~0.07 ms ≈ 6 % of step |
| L1.2 + L1.3 (4 → 2 transforms) | ~0.98 ms ≈ 6.9 % of step | ~0.14 ms ≈ 12 % of step |

Plus two fewer Legendre kernel launches per step and correspondingly fewer graph nodes. The relative
gain is larger at low resolution, the opposite gradient to W1 and W2 — worth stating in the results
table, since it makes L1 an ensemble/low-resolution win first and a T255 win second.

**Deliverables:** source change; a GPU test asserting the merged `K = 6L+1` batched plan exists and
that a step performs the expected number of transforms; a bitwise/tolerance equivalence run against
`main` in the style of `verify_batch_equivalence.jl` (the coslat subrange is the part most likely to
be silently wrong — check `u`, `v` on the grid explicitly, not just the prognostics); Phase A table
both GPUs; `node`-granularity re-profile; `CHANGELOG.md` bullet; version bump per `CLAUDE.md`.

### L4 — the forward Legendre transform is 1.7× the inverse

At T255 the two Legendre kernels cost 1.231 ms (forward) and 0.739 ms (inverse), together 18.8 % of
GPU busy — the same order as the entire FFT. The inverse was optimised to 37–41 % of peak bandwidth
by [`legendre-gpu-optimization.md`](legendre-gpu-optimization.md); the forward was not, and that plan
records the likely mechanism: moving to a 2-D `(row, layer)` launch put all layers in flight
simultaneously, so the accumulation array no longer fits L2 (8.4 MB against the A40's 6 MB; the H100
has 50 MB, so the effect should be *weaker* there — check this, it is a falsifiable prediction).

The forward mega-transform runs at `K = 9L+1 = 73`, the largest batch in the step, which is exactly
the regime where an L2-residency effect would bite hardest.

Experiments, in order, all on the A40 first (where the effect should be strongest) then the H100:

1. **Establish the scaling.** Time the forward transform alone (existing `SpeedyTransforms`
   benchmarks) at fixed truncation for `K ∈ {8, 16, 33, 49, 73}` and plot ms **per layer**. Flat ⇒
   no residency effect and L4 is dead. Rising with `K` ⇒ the hypothesis holds and the knee locates
   the working-set limit.
2. **Test the fix cheaply before writing a kernel.** `_transform_chunked!` already exists
   ([spectral_transform.jl:545](../../../SpeedyTransforms/src/spectral_transform.jl#L545)) and splits
   a call into layer chunks. Force the forward mega-transform through it at the chunk size the knee
   suggests and re-measure. This costs more FFT invocations (chunking multiplies the ring loop by the
   number of chunks), so it is a *diagnostic*, not a candidate fix — read the Legendre kernel time
   from the trace, not the total.
3. **If confirmed**, the fix is inside `forward_legendre_kernel!`
   ([legendre_ka.jl](../../../SpeedyTransforms/src/legendre_ka.jl), the kernel whose comment block
   begins "Parallelised over the *output* coefficients"): tile the layer axis so a workgroup's
   working set stays L2-resident, or order the thread table so that layers are the slowest-varying
   index rather than the fastest. Both are launch-configuration changes, not algorithmic ones.

Constraint carried over from the Legendre plan: **do not add type parameters to `SpectralTransform`**.
Its type name is within ~10 characters of the 1024-character cap enforced by
`SpeedyTransforms/test/type_name_length.jl`, and exceeding it silently disables nested Enzyme AD.

**Predicted payoff:** bringing the forward kernel to the inverse's efficiency would recover ~0.4–0.5
ms/step at T255 (≈3.5 % of step). Unlike L1 this is real work being done inefficiently, so the gain
should grow with resolution — the opposite gradient to L1, which is why the two are complementary.

### L5 — *maybe*: let the independent ring FFTs run concurrently

Filed as an experiment, not a plan. The strongest form of the idea is not generic multi-streaming of
the model; it is this specific observation:

**The 2 · nlat_half cuFFT invocations inside one transform are mutually independent, but the CUDA
graph currently serialises them.** `run_graph!`
([SpeedyTransformsCUDAExt.jl:299-336](../../../SpeedyTransforms/ext/SpeedyTransformsCUDAExt.jl#L299))
captures the ring loop from a single stream, so the captured graph is a linear dependency chain of
~2·nlat_half nodes (the extension's own comment confirms the node count). A graph captured with a
fork/join across several streams would express the independence, and the driver could then run many
of those ~1.4 µs kernels concurrently — with **no numerical change and no algorithmic change**, since
it is the same kernels on the same data.

If the FFT is dispatch-bound at ~5 % of peak bandwidth (see *Problem description*), concurrency is
precisely the lever that converts that into throughput, and it is complementary to L1: L1 reduces the
number of invocations, L5 overlaps whatever remains.

Minimal experiment before committing to anything:

1. Outside SpeedyWeather, replay a hand-built graph of N independent small cuFFT calls captured
   (a) from one stream and (b) from 4–8 forked streams, and compare total GPU time. If (b) is not
   materially faster, stop — the driver is not going to overlap them inside the model either.
2. Only if it is: change the capture in `run_graph!` to fork the ring loop across a small pool of
   streams and join before the graph ends. Hazards to check — the gather/scatter kernels that bracket
   the ring loop must remain ordered against all FFT nodes; the graph cache is keyed by device
   pointer and `add` mode and must stay so; capture failure must still fall back to `loop!()`.

Explicitly out of scope for L5: multi-streaming the dynamical core or the parameterizations. Those
have genuine data dependencies and the payoff is unquantified until M1 says how much time they use.

## Findings

### M1.1 — phase attribution resolved; cuFFT is 69 % of GPU busy (H100, 2026-08-14)

Re-read of the existing `w1g_T*_L8` traces (H100, `--cuda-graph-trace=node`, 20 captured steps).
The harness *is* present at `SpeedyWeather/test/GPU/modelbench/` — the *Background* warning is stale.
No `sqlite3` binary on the login node; Python's stdlib `sqlite3` was used instead.

**(a) The NVTX ranges are disjoint siblings, not nested.** Dumping one T255 step's host ranges:

```
step                 dur 14854.58us
  reset_tendencies   start+     0.81us  dur   120.33us
  parameterizations  start+   121.51us  dur   246.12us
  ocean_land         start+   368.20us  dur   119.74us
  dynamics           start+   488.33us  dur  7506.18us
  diffusion_implicit start+  7995.79us  dur   115.33us
  leapfrog           start+  8111.53us  dur    98.58us
  transform          start+  8210.45us  dur  6631.56us
```

`dynamics` ends at +7994 µs and `transform` starts at +8210 µs: no overlap. **The parent plan's W3.1
table header "`dynamics` (includes the transforms it calls)" is wrong** — the rows are disjoint, as
the exact sum to the busy total implied. The `dynamics` row contains the *two* transforms called from
`dynamics_tendencies!`, and the `transform` row the *two* from `transform.jl`.

**(b) But the disjoint reading's conclusion is wrong too.** Attributing kernels by name rather than
by range — which is boundary-independent, and necessary because the host runs far enough ahead of the
device that projecting GPU timestamps onto host range boundaries smears work across ranges (it parks
1.66 ms/step and 454 kernels under `feedback_output` at T255) — gives:

| category | T31 | T63 | T127 | T255 |
|---|---:|---:|---:|---:|
| cuFFT | 0.439 (59.0 %) | 1.154 (73.6 %) | 2.881 (77.1 %) | **7.180 (68.6 %)** |
| Legendre | 0.037 (4.9 %) | 0.086 (5.5 %) | 0.337 (9.0 %) | 1.970 (18.8 %) |
| transform gather/scatter | 0.038 (5.1 %) | 0.061 (3.9 %) | 0.146 (3.9 %) | 0.485 (4.6 %) |
| **other (dycore + physics)** | 0.231 (31.0 %) | 0.268 (17.1 %) | 0.371 (9.9 %) | **0.831 (7.9 %)** |
| GPU busy ms/step | 0.745 | 1.569 | 3.735 | 10.465 |
| kernels/step | 319 | 615 | 1283 | 2843 |

**Both readings in the *Problem description* were wrong.** The parent plan said cuFFT was 44 % of
busy; the correction in this document's *Problem description* said 18.7 %, counting only the
`preprocess`/`postprocess`/two named kernels visible in W3.1's truncated top-15 kernel table. The
actual total over *all* cuFFT kernels is **69 %** at T255 — the truncated table was hiding a long
tail of per-ring-length FFT variants. Sub-breakdown at T255:

| cuFFT kernel class | ms/step | kern/step | mean µs |
|---|---:|---:|---:|
| `regular_fft*` | 3.321 | 1150 | 2.89 |
| `vector_fft*` | 1.185 | 338 | 3.51 |
| other `fft` variants | 1.100 | 328 | 3.35 |
| `preprocess_kernel` | 0.836 | 606 | 1.38 |
| `regular_bluestein_fft` | 0.455 | 112 | 4.06 |
| `postprocess_kernel` | 0.284 | 196 | 1.45 |
| **total** | **7.180** | **2730** | 2.63 |

Kernel duration distribution: min 1.15 µs, median 2.53 µs, p75 3.36 µs, max 8.86 µs.

**Decision gate result.** The gate asked whether non-transform dynamics is ≳2000 kernels and ≳3
ms/step (→ re-plan) or a few hundred kernels and ≲1 ms/step (→ transform confirmed). It is **93.2
kernels/step and 0.831 ms/step** — the low branch, decisively, and *resolution-independent in count*
(93.2 kernels/step at every truncation, the dycore is fully batched). The transform is confirmed as
the target, but **not the part of it L1 and L4 address**:

- **L1** (4 → 2 transforms) removes at most half the ring invocations. Real, but bounded.
- **L4** (forward Legendre) targets 1.970 ms = 18.8 % at T255 — and note it *collapses* at low
  resolution (4.9 % at T31), so it is a T255-and-up lead only.
- The **per-ring FFT loop is 69 % of GPU busy** and is the item currently filed under *Future work*.
  At a 2.63 µs mean for kernels moving ~50–100 KB, this remains dispatch-bound, so the
  *Problem description*'s core argument survives its own arithmetic being wrong — it is in fact
  understated.

Note the FFT share is *non-monotonic* (59 → 74 → 77 → 69 %): T255 is where Legendre finally becomes
significant, so the FFT's share dips even as its absolute cost quadruples.

### E1 — environment notes and a false alarm (2026-08-14)

Two things cost time on 2026-08-14 and are recorded so the next session does not repeat them.

**E1.1 — CUDA.jl needed its runtime pinned.** On this login node `using CUDA` failed with
`CUDA runtime not found` / "JLLs were precompiled without an NVIDIA driver present", even though the
A40 and driver 580.126.09 (CUDA 13.0) are visible to `nvidia-smi`. `Base.compilecache` on
`CUDA_Runtime_jll` alone did **not** fix it. What worked:

```julia
julia --project=SpeedyWeather/test/GPU/modelbench -e 'using CUDA; CUDA.set_runtime_version!(v"12.9")'
# then restart Julia
```

This writes `LocalPreferences.toml` in the modelbench project (and adds `CUDA_Runtime_jll` to
`[extras]` in its `Project.toml`). Both are **new, uncommitted** files as of this writing — decide
whether to commit them. After the restart `CUDA.functional()` is `true`, runtime 12.9.0.

**E1.2 — the "NaN at time step 1" is a harness artifact, not a model bug.** Driving the model with a
bare `time_step!` loop:

```julia
m = PrimitiveWetModel(spectral_grid); s = initialize!(m)
for i in 1:5; SpeedyWeather.time_step!(s); end     # <-- NaN at step 1
```

produces NaN in the prognostics at step 1 for `PrimitiveWetModel` **and** `PrimitiveDryModel`, on
**CPU and GPU alike** and in **both Float32 and Float64**, while `BarotropicModel` and
`ShallowWaterModel` stay finite. That combination — CPU as well as GPU, F64 as well as F32 — rules
out a GPU or precision bug, and the checkout is *pristine*: `HEAD` is exactly `origin/main`
(`c343adf6`, v0.22.0) with **no** tracked file under any `src/` modified. The GPU CI passing is
therefore correct and consistent.

The cause is that `run!` performs a prologue that a raw `time_step!` loop skips —
[simulation.jl:78](../../../SpeedyWeather/src/models/simulation.jl#L78) calls
`transform!(variables, model, initialize = true)`, and `run!` also radius-scales vorticity and
divergence for the dynamical core (undone by `unscale!` at
[simulation.jl:95](../../../SpeedyWeather/src/models/simulation.jl#L95)). Stepping without it starts
the dycore from an unscaled, untransformed state. Confirmations:

- `run!(s, period=Day(1))` on `PrimitiveWetModel` T32 L8 CPU → all finite, `maxabs(vor) = 2.57e-5`;
- inserting `SpeedyWeather.initialize!(s; period=Day(10), output=false)` before the raw loop → 40 GPU
  steps at T32 L8, all finite. (`maxabs(vor) = 163.8` there because that state is still
  radius-*scaled*; `run!` unscales before returning. Do not mistake this for an instability.)

**The harness was checked and is already correct** — `bench_model.jl:78` and `nsys_target.jl:139`
both call `SpeedyWeather.initialize!(simulation; steps, output = false)` before their
`time_step!` loops, and `nsys_target.jl:104-128` additionally asserts bitwise equivalence against
`run!`. Only the ad-hoc one-liner used during the environment smoke test omitted it. **No existing
measurement is invalidated**, and no harness fix is needed. The rule to remember: any *new* driver
calling `time_step!` directly must call that prologue first.

## Testing and verification

- **Correctness gate for L1 and L4:** the modelbench bitwise-equivalence check (driver vs `run!`),
  then `SpeedyWeather/test/GPU/runtests.jl` in full, then the CPU `SpeedyWeather` and
  `SpeedyTransforms` suites. Note the parent plan's caveat: `dynamics/dispatch.jl` fails in some
  checkouts for dependency-version reasons unrelated to this work (`test/Manifest.toml` is
  gitignored) — confirm any failure reproduces on a clean `main` before attributing it.
- **L1 specifically:** batched-vs-merged tolerance comparison of a ≥100-step run, asserting `u` and
  `v` on the grid as well as the five prognostics. W1.1 found the K-batching change to be bitwise
  identical; the `unscale_coslat` subrange in L1.1 is the part that can plausibly break, and it will
  break *silently* (a coslat factor applied to the wrong slots looks like a plausible wind field).
- **L4:** the existing `SpeedyTransforms` transform unit tests plus a CPU-vs-GPU comparison at the
  new launch configuration.
- **Performance:** Phase A on both GPUs, all four truncations, `nans_detected` guard, with the ±20 %
  noise convention from `CLAUDE.md`. Every table appended to this document with the GPU named.
- **Regression guard:** after L1, assert the transform count per step (the modelbench driver can
  count `transform!` entries) so a future change cannot silently split the merged call back apart —
  the same failure mode as the `K = 2` fallback W1 fixed.

## Documentation changes

- `CHANGELOG.md` bullet per PR under `## Unreleased`.
- Findings M1.1, M2.1, M3.1 and the L1/L4 results tables appended to this document; status updated.
- **Correct the parent plan in the same PR as M1**: the W3.1 table header ("includes the transforms
  it calls") and the "cuFFT 44 % of GPU busy" sentence in its *Where to continue* section.
- If L5 is pursued, the concurrency preconditions must be documented at the capture site in
  `SpeedyTransformsCUDAExt.jl`, next to the existing pointer-stability and node-count comments.

## Known limitations

- ~~Every number in this document is arithmetic on the parent plan's published tables~~ — **partly
  superseded 2026-08-14**: the M1.1 table is measured (H100, from the `w1g_T*` traces). Everything
  *outside* M1.1/E1 is still arithmetic, including the ~5 %-of-peak-bandwidth figure and the
  per-transform payoff predictions in L1. M1.1 did not measure bandwidth, only durations and counts;
  the dispatch-bound conclusion still rests on the 2.63 µs mean against ~50–100 KB per invocation,
  which is an inference, not a counter reading. Phase D would settle it.
- The 66–78 % busy fraction that motivates the whole re-scoping still rests on one
  `node`-granularity nsys trace (M2 addresses this). **M1.1 does not help here** — it re-reads the
  same trace, so it inherits any profiler-overhead bias in the *busy fraction*. It is however
  unaffected in the *ratios between kernels*, which is what M1.1 actually claims.
- ~~The modelbench harness is not in `main`~~ — **refuted 2026-08-14, the harness is present**; see
  *Background*. Effort estimates hold.
- `nlayers = 8` throughout. L1's payoff is per-transform and therefore roughly independent of
  `nlayers`, whereas the FFT's *work* is not — a layer sweep would separate the two and has been
  outstanding since the parent plan.
- CUDA-only. L1 is backend-neutral (fewer transforms is fewer transforms everywhere); L5 is
  CUDA-graph-specific and would need an equivalent on Metal/AMDGPU.
- Phase D (`ncu` counters) remains blocked by `ERR_NVGPUCTRPERM`. It is the direct way to settle
  whether the FFT and Legendre kernels are bandwidth- or occupancy-limited, which is the question
  under both L4 and L5; unblocking it would make both cheaper to decide.

## Future work

- **Uniform-`nlon` batched FFT path**, gated on M3: one `plan_many` per transform instead of
  2·nlat_half calls, on a full grid. Largest single ceiling identified so far (~1.8 ms/step at T255)
  and the only route that removes the per-ring loop rather than shrinking it.
- **Custom KernelAbstractions FFT** over all rings in one launch, which is the same win without
  giving up the reduced grid: one workgroup per (ring, layer), mixed radix in shared memory —
  octahedral ring lengths are `16 + 4j ≤ 784 = 2⁴·7²`, all shared-memory-resident in Float32. It
  would also remove cuFFT's `preprocess`/`postprocess` kernel pairs (802 launches, 1.12 ms/step at
  T255) and the Bluestein fallback. High effort, backend-neutral, and it subsumes L5.
- **Whole-step CUDA graph, re-scoped to low resolution.** The parent plan demoted this on T255
  arithmetic (0.30 ms of a 14.28 ms step). The same fixed 0.385 ms is 34 % of the step at T31, where
  ensembles run. Worth reopening as a low-resolution optimisation once M2b has measured the host
  slack.
- Roofline the transform once `ncu` is unblocked; fold the Phase A guards into
  `SpeedyWeather/benchmark/manual_benchmarking.jl`.

## State and next actions

> Written at the end of the 2026-08-14 session. Read this first when resuming.

### What is done

| item | status |
|---|---|
| **M1** — phase-attribution ambiguity | **done**, finding M1.1. Answer: *neither* published reading was right; cuFFT is **69 %** of GPU busy at T255. Decision gate took the "transform confirmed" branch, but pointed past L1/L4 to the per-ring FFT loop. |
| Correcting the parent plan | **done**, in `gpu-idle-time-reduction.md`: the W3.1 table header, the two "FFT 44 %" claims, and the *Where to continue* paragraph all now carry the M1.1 correction. |
| GPU environment | **working**, see E1.1 — required `CUDA.set_runtime_version!(v"12.9")`. |
| "NaN at step 1" | **closed, not a bug**, see E1.2. Harness was already correct; no measurement invalidated. |
| **M2, M3, M4, L1, L4, L5** | **not started.** |

Hardware available: A40 on the login node (shared — check for strays first; one ~490 MiB stray was
present all session and was left alone), H100 via SLURM (`--partition=gpu --qos=gpushort
--account=flai --gres=gpu:1`, harness scripts `submit_phaseA.sh` / `submit_phaseB.sh` already carry
these). All M1 numbers above are **H100**, re-read from the existing `w1g_T*` traces.

### Next actions, in order

1. **Re-read M1.1, then decide whether M3 is still the right next measurement.** M1 changed the
   picture enough that the M2/M3/M4 ordering written before it is no longer obviously right:
   - **M2** (independent busy-fraction + host slack) is now *less* urgent for choosing a target —
     M1.1 identifies the target from kernel durations alone, which no profiler-overhead argument
     touches. It is still needed before anyone attacks *idle*, and M2b's host-slack number is still
     the gate on the low-resolution whole-step-graph item under *Future work*.
   - **M3** (`FullGaussianGrid` sweep) is now the *highest-value* measurement, because it prices the
     one route that removes the per-ring loop rather than shrinking it, and that loop is 69 % of
     busy. It is also zero-code: one login-node command plus one H100 job. **Recommend starting
     here.** Note the go/no-go thresholds in the M3 section were written against a 19 % FFT share and
     should be *loosened* now that the share is 69 % — a full-grid slowdown well above the old
     "~40 % kills it" line could still be worth paying.
2. **Run M3** exactly as specified (`Grid=FullGaussianGrid`, capital-G kwarg; four truncations;
   matched on truncation, not point count, because `dealiasing` recomputes `nlat_half`). Record as
   M3.1.
3. **Run M2b** (host-slack spin-loop sweep) — cheap, profiler-free, and it is the only thing that
   tells us whether the low-resolution whole-step graph is worth reopening. M2a is optional; if
   CUDA.jl does not expand graph nodes it is void anyway.
4. **Re-prioritise L1/L4/L5 against M1.1 before writing any code.** As measured at T255:
   L1 targets the *number* of transforms (bounded by ~half the ring invocations), L4 targets 18.8 %,
   L5 targets concurrency across the ring loop. Given cuFFT is 69 %, **L5 and the two *Future work*
   FFT items now look stronger relative to L1/L4 than this document originally assumed** — L5's
   step 1 (the standalone one-stream vs forked-stream cuFFT replay experiment) is cheap and is
   arguably the better next implementation probe than L1. L4 additionally *shrinks* at low
   resolution (18.8 % at T255 → 4.9 % at T31), so it is a high-resolution-only lead.

### Gotchas for the next session

- `module load julia` does **not** persist between separate shell invocations here; load it in every
  command. Julia is 1.12.2 by default.
- There is **no `sqlite3` binary** on the login node. Use Python's stdlib `sqlite3` (verified
  working) for any trace re-analysis.
- Attribute trace work by **kernel name, not by NVTX range**. The host runs far enough ahead of the
  device that projecting GPU timestamps onto host range boundaries misfiles ~1.7 ms/step at T255
  into `feedback_output`. This is the trap that produced the original "44 %" error.
- The harness *is* present at `SpeedyWeather/test/GPU/modelbench/` — the *Background* section's
  warning that it is missing is **stale**; ignore it.
- `trunc=` still works but is deprecated in favour of `truncation` (`truncation = trunc + 1`), so the
  harness emits a depwarn. Existing baselines were taken with `trunc=`, so keep using it for
  comparability rather than silently shifting resolution by one.
- Uncommitted and undecided: `LocalPreferences.toml` and the `[extras]` addition to `Project.toml` in
  the modelbench project (both from E1.1).
