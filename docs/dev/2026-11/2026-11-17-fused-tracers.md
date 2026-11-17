# Fuse the tracers: one batched transform and one kernel per operation for all tracers

> Status: **in progress**. Phase 0 ran on an NVIDIA H100 and says go; phase 1 is implemented (all tracers processed
> together by default on GPU) and tested on CPU, GPU and with Enzyme, but the Reactant correctness test fails and
> the benchmarks of the implementation are still running (see "State at the end of the session").

Date of initial draft: 2026-11-17

Base revision: 174afe29 (`mg/unroll-tracer-loops`) plus the uncommitted part 2 of
`2026-11-16-unroll-tracer-loops.md` (tracers as a NamedTuple fixed at model construction), which this plan builds on.

## Originating prompt

> okay next, could we actually look into making a plan whether it's possible to also fuse multiple tracers in case
> they are defined with our fuse system for variables and advect them together? This should only be done if it's
> really faster on GPU

> decisions: 1) all tracers 2) yes fall back to one-by-one in this case 3) go with your proposal 4) for now just
> split in chunks when it's more than the largest allocated transform (set by `transform_batch`) 5) make the edits
> now, then we hand over to a gpu machine with nvidia gpu

## Revision log

- 2026-11-17: initial draft.
- 2026-11-17: open questions decided: fuse all tracers automatically; fall back to the per-tracer path when any tracer
  is deactivated; the proposed go/no-go thresholds; split transforms into chunks of whole tracers when they exceed the
  transform's largest batch (`S.nlayers`, from `transform_batch`), no tracer-dependent transform construction for now;
  the phase 0 edits are made here, the NVIDIA runs on a GPU machine.
- 2026-11-17: phase 0 implemented. 0a as an opt-in suite (`--tracers`, `:benchmark510`) of the benchmark, so that a
  regular benchmark run doesn't take about 40 runs longer; 0b as the script
  `SpeedyWeather/benchmark/fused_tracers_prototype.jl`. Both are meant to run in the same session on the GPU machine
  (one handover instead of two), the gate of 0a is then applied to the results of both. The prototype merges the
  `+A*div` and flux kernels into one kernel over all tracers.
- 2026-11-17: phase 0 run on an NVIDIA H100 80GB (the cluster's GPU nodes, as the stored `gpu-nvidia` results; Julia
  1.12.2, CUDA.jl 6.4.1), the full benchmark with the flags of the stored run (`--repeats 10 --length 2 --hi-res`) plus
  `--tracers`, and the prototype; results below. **Go**: 8 tracers add far more than 20% at T ≤ 128, the prototype
  saves 24–72% of the time step with 4 to 8 tracers at T ≤ 128 and 55–68% at T256. The CPU part of the gate failed
  on the Apple M3 (up to 7% slower), so as planned the fused path is the default on GPU only. Reactant: the tracer suite
  only ran for the ShallowWaterModel (also in a partial run in a separate worktree, to keep the stored results); every
  PrimitiveWetModel construction on Reactant fails with scalar indexing in `free_troposphere_layers`
  (`parameterizations/radiation/clouds.jl`), with or without tracers, a bug unrelated to this plan.
- 2026-11-17: phase 1 implemented. Deviations from the plan: (1) the switch is a model option `fuse_tracers`
  (default `default_fuse_tracers(architecture)`: `true` on GPU, `false` on CPU and Reactant), not only the
  architecture, so that both paths can be tested against each other on any architecture, and the fused path is taken
  if `fuse_tracers` and all tracers are active (`fused_tracers(vars, model)`); (2) the fluxes are `DynamicsVariable`s
  (as the dycore's `uq`, `vq`) named `<tracer>_u_flux`, `<tracer>_v_flux` in the namespaces `:tracers` (spectral) and
  `:grid_tracers` (grid), ordered `[tendencies; u fluxes; v fluxes]` by a `variables` method for the NamedTuple of
  tracers; the per-tracer path uses them too, so the tracers don't declare the scratch `a`, `b` any more; (3) the
  forward transform is one batch of 3N·L layers (in chunks), the prototype showed it best for few tracers and about
  equal to three batches for more; (4) the spectral flux divergence (`divergence!`) has the same layer limit as the
  transforms, it is chunked likewise (`stacked_chunks`); (5) vertical advection uses a stacked kernel (`k =
  mod1(kn, L)` and the stencil shifted to the tracer's layers) instead of reshaping, reusing the
  `reconstructed_at_face` methods; (6) moving the grid step back and (un)scaling the grid tendencies always act on
  the tracer parents (one broadcast each, Reactant-safe as for the dycore), the per-tracer loops there are removed.

## Problem description

Every tracer is processed on its own, in each of the per-tracer functions that `foreach_active_tracer` calls once
per active tracer. Per tracer and time step that is:

| Operation | Function | Transforms | Kernels |
|---|---|---|---|
| spectral → grid | `transform_tracer!` | 1 (inverse, L layers) | |
| horizontal advection | `horizontal_advection_tracer!` → `horizontal_advection!` (the "sequential" version in `tendencies_sequential.jl`) | 3 (forward, L layers each: `A_tend_grid`, `u*A`, `v*A`) | `+A*div`, flux `(u*A, v*A)`, spectral flux divergence |
| vertical advection (primitive equation models) | `vertical_advection_tracer!` | | 1 |
| horizontal diffusion | `horizontal_diffusion_tracer!` | | 1 |
| time step | `update_prognostic_tracer!` | | 1 (leapfrog + filters) |
| move grid step back, scale grid tendencies | `move_grid_tracers_back!`, `scaling.jl` | | 1 each |

So N tracers cost 4N transforms and about 7N kernel launches per time step, all for arrays of only L layers. The
dynamical core does the opposite: it batches its variables into fused parents (`vars.fused.prognostic`,
`vars.fused.grid`, `vars.fused.spectral_tendencies`, …) and transforms them in a few calls with 2L to 9L+1 layers
each ("mega-batched" transforms).

The idea: allocate the tracer variables in fused parents too, and replace the N calls per operation with one call
on the fused parent, i.e. one transform with N·L layers instead of N transforms with L layers, and one kernel over
N·L layers instead of N kernels. This only pays off if it is measurably faster on GPU, so phase 0 measures before
anything is implemented.

## Background

### Why it should be faster on GPU

- **Small transforms are limited by kernel launches and latency, not by work.** In the stored NVIDIA results
  (`SpeedyWeather/benchmark/README.md`, `gpu-nvidia`), the PrimitiveWetModel at T32 runs at the same speed with 4, 8,
  12 or 16 layers (8527, 8601, 8371, 8556 SYPD): at T32 a wider batch is practically free. One time step there takes
  about 0.8 ms (SYPD 8349 at Δt = 2400 s).
- **Wider batches stay sub-linear at high resolution.** Tripling the layers from 8 to 24 slows the model down by only
  1.3× at T128 (423 → 321 SYPD), 1.6× at T256 (78 → 50) and 1.9× at T512 (11 → 5.7).
- **GPU graphs don't remove the launches.** CUDA/HIP graphs (`gpu_graphs`, `SpeedyTransformsCUDAExt.jl`) only cover the
  Fourier part of a transform (gather, per-ring cuFFTs, scatter). The Legendre part and all elementwise kernels are
  launched one by one. The graphs are also cached per device pointer of the field, so every tracer gets its own
  cached graph today.
- On GPU the default transform up to T100 is the `MatrixSpectralTransform` (dense matrix multiplies), for which more
  columns per call is the most efficient case.

What isn't known is the actual cost of a tracer today, i.e. how much of a time step N tracers add. That decides
whether there is anything worth saving.

### Why it's less clear on CPU and Reactant

- CPU transforms are not latency-bound and split batches that weren't planned into chunks of the largest planned
  batch (`L` by default), so a fused transform does the same work in the same chunks. Fusing should be roughly
  neutral there; the plan only needs it not to be slower.
- Reactant compiles the whole time step into one XLA program, which may already overlap or merge independent
  operations. The benefit has to be measured separately.

### The fuse system

- A variable declared with `fuse = :symbol` is allocated as a view into one parent per `(namespace, fuse)` group.
  Members are stacked along the layer axis (axis 2) of the parent; the step axis stays the 3rd. `FusedParent` stores
  the `slot_map` (member name → layer range), and `vars.fused.<symbol>` holds the parent.
- Variables of different variable types can share a group, but the same name may appear only once per group. Fuse
  symbols must be unique across namespaces.
- `_assert_fuse_alignment` checks that two parents (e.g. `:prognostic` and `:grid`) list the same members in the same
  order, so that a batched transform maps slot k of one to slot k of the other.
- The per-variable views keep working, so all existing per-tracer code (and `set!`, output, `copy!`,
  `on_architecture`, Reactant's `copy_variables!`) works unchanged on fused tracer variables. The dycore already
  relies on this.

### Transform batch sizes and scratch memory

On GPU the transform's scratch memory is sized to the largest planned batch, `max_transform_batch(L) = 9L + 1` (the
PrimitiveWetModel's tendency batch), and FFT plans exist for the batch sizes in `default_transform_batch`. Other batch
sizes are planned on first use and cached (`ensure_batched_plans!`), but a batch wider than the scratch memory
doesn't fit. A fused tracer transform of N·L layers therefore fits into one call for up to 9 tracers; more tracers
have to be split into chunks of whole tracers. The CPU already splits unplanned batch sizes into chunks.

### Deactivated tracers

A deactivated tracer must stay frozen. Its spectral and grid variables must not change, bitwise, which the tests from
part 1 check. A fused operation processes all tracers in the parent, so it would also recompute the grid of a
frozen tracer and compute its tendencies. The time step could skip the update with a per-tracer mask, but the grid
would still be re-transformed once from the frozen spectral state, which differs from the grid that was frozen (that
one was computed before the last time step). See the open questions.

## Phase 0: measure (go/no-go)

No changes to `src/`. Two steps, the second only if the first shows enough potential.

### 0a. What do tracers cost today?

A new benchmark suite in `SpeedyWeather/benchmark/define_benchmarks.jl` (worth keeping anyway, also to catch tracer
regressions): SYPD with N = 0, 1, 2, 4, 8 and 16 tracers (all active) for

- ShallowWaterModel at T32, T64, T128 (L = 1),
- PrimitiveWetModel at T32, T64, T128, T256 (L = 8),

with the default transform per architecture. Run on `gpu` (NVIDIA, the machine of the stored results), `cpu-arm`,
and `reactant-nvidia-gpu`, plus `metal` locally as an early indication.

The difference between N tracers and none is the most that fusing can save. **Stop** if 8 tracers add less than
about 20% to the time step at T ≤ 128 on NVIDIA: then fusing can't give a speedup worth the code.

### 0b. Prototype the transforms and the horizontal advection

The 4 transforms per tracer are most of the expected saving. A benchmark script (not in `src/`) allocates N tracer
variables in fused parents for a real model and times, per time step, the per-tracer path (as today) against:

- one inverse transform of the fused prognostic tracers (N·L layers),
- the `+A*div` and flux kernels as one launch over N·L layers,
- the three forward transforms as three N·L batches, and alternatively as one 3N·L batch.

For N = 1, 2, 4, 8 at the resolutions of 0a, on NVIDIA GPU and CPU.

**Go** for phase 1 if, on NVIDIA, the fused path reduces the total time step by at least 20% with 4 to 8 tracers at
T ≤ 128, is not slower at T256, and is not more than 5% slower on CPU (otherwise phase 1 uses the fused path on GPU
only). The thresholds are a proposal (see open questions).

The results of phase 0 go into this document whatever the outcome.

### Implementation of phase 0

- **0a**, `SpeedyWeather/benchmark/`: `BenchmarkSuite` has a field `ntracers` (tracers `tracer1`, `tracer2`, … from
  `benchmark_tracers`, zero initial conditions). Its table gets the columns "Tracers" and "vs no tracers" (the wallclock
  time relative to the same configuration without tracers) when the number of tracers varies. The suite
  `:benchmark510` is opt-in with `--tracers`, as `--SPEEDY`, and its numbers are also stored in the `tracers` field of
  the architecture's record in `assets/benchmark_results.json`.
- **0b**, `SpeedyWeather/benchmark/fused_tracers_prototype.jl`: for each configuration and N = 1, 2, 4, 8 it builds the
  model with N tracers, runs 4 steps for a realistic state and times (median of 7 samples of at least 0.2 s each):
  the transform and horizontal advection of the tracers as the model does it now; the prototype with the forward
  transforms as three N·L batches and as one 3N·L batch (both split into chunks of whole tracers that fit `S.nlayers`);
  and one full time step. It also checks that the prototype gives the model's grid tracers and advection tendencies.
  `--quick` runs only T32 with 1 and 2 tracers.

### Results of phase 0

Tested on CPU (Apple M3, 4 threads, Julia 1.12.7):

- A mini version of the suite (ShallowWaterModel T32 with 0, 2, 4 tracers): 1.84× and 2.76× the time without
  tracers. Tracers are expensive for a model as cheap as the 1-layer ShallowWaterModel.
- The prototype gives results bitwise identical to the model's in all 28 cases (and no compared array is zero).
- On CPU the prototype is about neutral for the ShallowWaterModel (1.00× to 1.05×) and up to 15% slower for the
  PrimitiveWetModel (0.85× to 0.99×); the worst case, T64 with 8 tracers, is 7% of the time step slower. That's beyond
  the 5% of the gate, so phase 1 would use the fused path on GPU only. Presumably the per-tracer arrays of L layers
  stay in cache, the stacked ones of N·L layers don't; not investigated.
- Tracers are expensive on CPU too: 8 tracers' transforms and horizontal advection alone take 40% of a PrimitiveWetModel
  T64 time step (16.6 of 42.0 ms), and the ShallowWaterModel T32 time step with 8 tracers takes about 4× as long as
  without (0.73 vs. about 0.18 ms).

<details><summary>Prototype on CPU (Apple M3, 4 threads)</summary>

| Model | T | L | Tracers | per tracer [ms] | fused, 3 batches [ms] | fused, 1 batch [ms] | speedup | time step [ms] | saved of time step | max. rel. difference |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 32 | 1 | 1 | 0.070 | 0.071 | 0.097 | 0.98× | 0.245 | -0% | 0.0e+00 |
| ShallowWaterModel | 32 | 1 | 2 | 0.135 | 0.133 | 0.133 | 1.02× | 0.315 | 1% | 0.0e+00 |
| ShallowWaterModel | 32 | 1 | 4 | 0.264 | 0.262 | 0.260 | 1.02× | 0.454 | 1% | 0.0e+00 |
| ShallowWaterModel | 32 | 1 | 8 | 0.533 | 0.508 | 0.535 | 1.05× | 0.727 | 4% | 0.0e+00 |
| ShallowWaterModel | 64 | 1 | 1 | 0.298 | 0.292 | 0.322 | 1.02× | 1.122 | 0% | 0.0e+00 |
| ShallowWaterModel | 64 | 1 | 2 | 0.599 | 0.586 | 0.611 | 1.02× | 1.368 | 1% | 0.0e+00 |
| ShallowWaterModel | 64 | 1 | 4 | 1.159 | 1.135 | 1.136 | 1.02× | 1.959 | 1% | 0.0e+00 |
| ShallowWaterModel | 64 | 1 | 8 | 2.341 | 2.274 | 2.250 | 1.04× | 3.144 | 3% | 0.0e+00 |
| ShallowWaterModel | 128 | 1 | 1 | 1.596 | 1.604 | 1.604 | 1.00× | 5.763 | -0% | 0.0e+00 |
| ShallowWaterModel | 128 | 1 | 2 | 3.242 | 3.174 | 3.173 | 1.02× | 7.416 | 1% | 0.0e+00 |
| ShallowWaterModel | 128 | 1 | 4 | 6.393 | 6.339 | 6.319 | 1.01× | 10.720 | 1% | 0.0e+00 |
| ShallowWaterModel | 128 | 1 | 8 | 12.884 | 12.785 | 12.647 | 1.02× | 17.290 | 1% | 0.0e+00 |
| PrimitiveWetModel | 32 | 8 | 1 | 0.391 | 0.447 | 0.431 | 0.91× | 5.658 | -1% | 0.0e+00 |
| PrimitiveWetModel | 32 | 8 | 2 | 0.785 | 0.886 | 0.885 | 0.89× | 6.213 | -2% | 0.0e+00 |
| PrimitiveWetModel | 32 | 8 | 4 | 1.563 | 1.711 | 1.740 | 0.91× | 7.162 | -2% | 0.0e+00 |
| PrimitiveWetModel | 32 | 8 | 8 | 3.343 | 3.474 | 3.462 | 0.97× | 9.124 | -1% | 0.0e+00 |
| PrimitiveWetModel | 64 | 8 | 1 | 2.029 | 2.147 | 2.154 | 0.95× | 25.038 | -0% | 0.0e+00 |
| PrimitiveWetModel | 64 | 8 | 2 | 4.160 | 4.359 | 4.349 | 0.96× | 27.528 | -1% | 0.0e+00 |
| PrimitiveWetModel | 64 | 8 | 4 | 8.256 | 9.152 | 9.082 | 0.91× | 32.236 | -3% | 0.0e+00 |
| PrimitiveWetModel | 64 | 8 | 8 | 16.554 | 19.457 | 19.588 | 0.85× | 42.030 | -7% | 0.0e+00 |
| PrimitiveWetModel | 128 | 8 | 1 | 12.355 | 12.587 | 12.746 | 0.98× | 116.176 | -0% | 0.0e+00 |
| PrimitiveWetModel | 128 | 8 | 2 | 25.258 | 26.307 | 26.023 | 0.97× | 129.822 | -1% | 0.0e+00 |
| PrimitiveWetModel | 128 | 8 | 4 | 49.631 | 53.067 | 52.709 | 0.94× | 155.195 | -2% | 0.0e+00 |
| PrimitiveWetModel | 128 | 8 | 8 | 99.123 | 105.525 | 105.364 | 0.94× | 208.535 | -3% | 0.0e+00 |
| PrimitiveWetModel | 256 | 8 | 1 | 82.535 | 83.316 | 83.617 | 0.99× | 619.271 | -0% | 0.0e+00 |
| PrimitiveWetModel | 256 | 8 | 2 | 164.205 | 166.340 | 166.505 | 0.99× | 722.753 | -0% | 0.0e+00 |
| PrimitiveWetModel | 256 | 8 | 4 | 328.678 | 333.626 | 333.745 | 0.99× | 882.095 | -1% | 0.0e+00 |
| PrimitiveWetModel | 256 | 8 | 8 | 666.807 | 740.590 | 689.685 | 0.97× | 1326.712 | -2% | 0.0e+00 |

</details>

Run on an NVIDIA H100 80GB HBM3 (AMD EPYC 9554, Julia 1.12.2, CUDA.jl 6.4.1, default transform):

- **0a, the tracer suite** (table below, from the full `gpu-nvidia` benchmark run): still running at the end of the
  session, see the handover. <!-- 0A-SUMMARY -->
- **0b, the prototype** (table below): the fused transforms and horizontal advection are 1.4× to 7.6× faster than
  per tracer and save 24–72% of the whole time step with 4 to 8 tracers at T ≤ 128 and 55–68% at T256, matching the
  model to ≤ 2e-6 (relative). A single tracer is already 1.4–2.8× faster (one forward transform instead of three).
  At T128 and T256 the transforms are so latency-bound that a transform of 64 layers takes about as long as one of 8.
- **Decision: go** for phase 1, on GPU (the CPU part of the gate failed on the M3, see above).
- Reactant on the H100 (ShallowWaterModel only, see the revision log): 1, 2, 4, 8 and 16 tracers take 1.4×, 1.9×,
  2.8×, 4.8× and 7.7× the time without tracers at T32 (3.3× with 8 tracers at T128), i.e. XLA does not merge the
  per-tracer work by itself.

<!-- 0A-TABLE -->

<details><summary>Prototype on NVIDIA H100</summary>

| Model | T | L | Tracers | per tracer [ms] | fused, 3 batches [ms] | fused, 1 batch [ms] | speedup | time step [ms] | saved of time step | max. rel. difference |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ShallowWaterModel | 32 | 1 | 1 | 0.110 | 0.110 | 0.067 | 1.64× | 0.569 | 7% | 0.0e+00 |
| ShallowWaterModel | 32 | 1 | 2 | 0.221 | 0.124 | 0.069 | 3.21× | 0.700 | 22% | 7.5e-07 |
| ShallowWaterModel | 32 | 1 | 4 | 0.441 | 0.126 | 0.098 | 4.49× | 0.966 | 35% | 7.5e-07 |
| ShallowWaterModel | 32 | 1 | 8 | 0.905 | 0.135 | 0.130 | 6.94× | 1.505 | 51% | 7.6e-07 |
| ShallowWaterModel | 64 | 1 | 1 | 0.351 | 0.348 | 0.257 | 1.36× | 1.058 | 9% | 0.0e+00 |
| ShallowWaterModel | 64 | 1 | 2 | 0.701 | 0.534 | 0.285 | 2.46× | 1.421 | 29% | 1.7e-06 |
| ShallowWaterModel | 64 | 1 | 4 | 1.403 | 0.629 | 0.394 | 3.56× | 2.147 | 47% | 1.7e-06 |
| ShallowWaterModel | 64 | 1 | 8 | 2.803 | 0.596 | 0.569 | 4.93× | 3.594 | 62% | 1.7e-06 |
| ShallowWaterModel | 128 | 1 | 1 | 2.887 | 2.886 | 1.623 | 1.78× | 8.484 | 15% | 0.0e+00 |
| ShallowWaterModel | 128 | 1 | 2 | 5.779 | 2.909 | 1.655 | 3.49× | 11.394 | 36% | 0.0e+00 |
| ShallowWaterModel | 128 | 1 | 4 | 11.545 | 2.970 | 2.363 | 4.88× | 17.209 | 53% | 0.0e+00 |
| ShallowWaterModel | 128 | 1 | 8 | 23.093 | 3.078 | 3.151 | 7.50× | 28.843 | 69% | 0.0e+00 |
| PrimitiveWetModel | 32 | 8 | 1 | 0.137 | 0.135 | 0.077 | 1.79× | 1.093 | 6% | 0.0e+00 |
| PrimitiveWetModel | 32 | 8 | 2 | 0.274 | 0.130 | 0.074 | 3.68× | 1.298 | 15% | 2.2e-07 |
| PrimitiveWetModel | 32 | 8 | 4 | 0.549 | 0.131 | 0.127 | 4.31× | 1.738 | 24% | 6.6e-07 |
| PrimitiveWetModel | 32 | 8 | 8 | 1.095 | 0.158 | 0.185 | 6.91× | 2.242 | 42% | 8.4e-07 |
| PrimitiveWetModel | 64 | 8 | 1 | 0.598 | 0.596 | 0.216 | 2.77× | 1.627 | 23% | 0.0e+00 |
| PrimitiveWetModel | 64 | 8 | 2 | 1.196 | 0.510 | 0.285 | 4.20× | 2.193 | 42% | 1.4e-06 |
| PrimitiveWetModel | 64 | 8 | 4 | 2.394 | 0.452 | 0.523 | 5.30× | 3.428 | 57% | 1.9e-06 |
| PrimitiveWetModel | 64 | 8 | 8 | 4.795 | 0.631 | 0.852 | 7.60× | 5.907 | 70% | 1.3e-06 |
| PrimitiveWetModel | 128 | 8 | 1 | 3.101 | 3.110 | 1.855 | 1.67× | 7.005 | 18% | 0.0e+00 |
| PrimitiveWetModel | 128 | 8 | 2 | 6.245 | 3.274 | 1.977 | 3.16× | 10.136 | 42% | 0.0e+00 |
| PrimitiveWetModel | 128 | 8 | 4 | 12.441 | 3.518 | 2.902 | 4.29× | 16.393 | 58% | 0.0e+00 |
| PrimitiveWetModel | 128 | 8 | 8 | 24.794 | 4.013 | 4.020 | 6.18× | 28.913 | 72% | 0.0e+00 |
| PrimitiveWetModel | 256 | 8 | 1 | 8.084 | 8.075 | 4.777 | 1.69× | 18.421 | 18% | 0.0e+00 |
| PrimitiveWetModel | 256 | 8 | 2 | 16.124 | 8.946 | 5.482 | 2.94× | 26.540 | 40% | 0.0e+00 |
| PrimitiveWetModel | 256 | 8 | 4 | 32.147 | 10.368 | 8.607 | 3.73× | 42.688 | 55% | 0.0e+00 |
| PrimitiveWetModel | 256 | 8 | 8 | 64.490 | 12.920 | 12.957 | 4.99× | 75.365 | 68% | 0.0e+00 |

</details>

### Handover to the NVIDIA GPU machine

On the machine of the stored `gpu-nvidia` results, from the repository root, on the branch with these changes:

1. Instantiate the benchmark environment: `julia --project=SpeedyWeather/benchmark -e 'using Pkg; Pkg.instantiate()'`.
2. Check the prototype works (a few minutes):
   `julia --project=SpeedyWeather/benchmark SpeedyWeather/benchmark/fused_tracers_prototype.jl gpu --quick`.
   The last column ("max. rel. difference") must be small (≲ 1e-5 in Float32); "Inf" means a compared array was zero.
3. **0b**, the prototype:
   `julia --project=SpeedyWeather/benchmark SpeedyWeather/benchmark/fused_tracers_prototype.jl gpu --output fused_tracers_gpu-nvidia.md`.
4. **0a**, the benchmark with the tracer suite: `cd SpeedyWeather/benchmark && julia --project=. manual_benchmarking.jl gpu --tracers`.
   This is a full benchmark run: it replaces the `gpu-nvidia` section of `README.md` and its record in
   `assets/benchmark_results.json`, now with the tracer suite. Don't use `--suites benchmark510` for this, a partial run
   replaces the stored results of all other suites of that architecture with nothing. If there is time, also
   `reactant-nvidia-gpu --tracers`.
5. Copy both tables (the tracer suite from the `gpu-nvidia` section of the README, the prototype's from
   `fused_tracers_gpu-nvidia.md`) into "Results of phase 0" above and apply the gate: phase 1 only if the tracer suite
   shows that 8 tracers add at least 20% at T ≤ 128, and the prototype's "saved of time step" is at least 20% with 4 to
   8 tracers at T ≤ 128, with no slowdown at T256. Record the decision in the revision log, then commit and push.

## Phase 1: implementation (phase 0 said go)

As planned below, with the deviations in the revision log. The code: `dynamics/tracers.jl` (variables, `fused_tracers`,
`transform_stacked!`), the fused drivers next to their per-tracer counterparts (`transform_tracers!`,
`horizontal_advection_tracers!`, `vertical_advection_tracers!`, `horizontal_diffusion_tracers!`,
`update_prognostic_tracers!`), `move_grid_tracers_back!` and the scaling on the parents.

### State at the end of the session of 2026-11-17 (handover)

Tests of the implementation:

- CPU, `test/dynamics/fused_tracers.jl`: fused and per-tracer paths **bitwise identical** for all four models
  (incl. chunked transforms with 4 and 11 tracers), deactivation switches to per-tracer and back, stacked kernels
  equal the per-tracer kernels: 94/94 passed, in 12 min (mostly compilation, 4 model types with 3 to 11 tracers).
- GPU (A40), `test/GPU/tracers.jl`: fused (default on GPU) vs per-tracer to `rtol = 1e-4`, 58/58 passed (then the
  redundant 3-tracer ShallowWater case was dropped).
- Existing tracer tests (`tracer_advection.jl`) pass.
- JET (`dispatch.jl`): all `@test_opt` in this file, **also at the base revision** `ac3c5befb`, fail in the freshly
  instantiated test environment because of runtime dispatch inside KernelAbstractions (`__thread_run`) and
  ModelParameters. With `ignored_modules = (Base, KernelAbstractions, ModelParameters)` all pass, including the new
  fused-tracer testset (ShallowWater, PrimitiveWet). Environment issue, not this change.
- Enzyme forward vs reverse (`forward_mode.jl`, Julia 1.10.12, Enzyme 0.13.220): Barotropic with a tracer and the new
  PrimitiveWet with 2 fused tracers pass (11/11, 12/12).
- **Reactant correctness: FAILS** in the existing single-tracer case (per-tracer path, `test_model(PrimitiveDryModel;
  tracers = (:tracer,))`): the grid tracer tendency matches (1.6e-5) but the spectral tracer tendency differs (max rel.
  0.39), and with it the time stepping. Suspected cause: the per-tracer path now writes the tendency and both fluxes of
  a tracer as views of the same fused buffer, one in-place update after another (`transform!` into each, then
  `divergence!` adding into the tendency) — the pattern that Reactant mis-handles within one trace (see the comment in
  `scale_tendencies!`); before, the fluxes were the separate scratch arrays `a`, `b`. To do: confirm (e.g. the base
  revision passes with the same Reactant 0.2.308), then fix, e.g. separate flux buffers on Reactant or one batched
  update; the new fused case (`tracers = (:tracer1, :tracer2), fuse_tracers = true`) hasn't run yet as the script
  stopped at the first failure. Run with the benchmark environment (Reactant 0.2.308) and a script that includes only
  the tracer cases of `test_correctness.jl`.

Benchmarks of the implementation, SLURM jobs still running (logs in `/p/tmp/maxgelbr/fused-tracers-phase1/`), each in
its own throwaway worktree at `df8f45ecd` (partial `--suites benchmark510 --tracers` runs, results in that worktree's
`SpeedyWeather/benchmark/README.md`, not to be committed):

| Job | Worktree `/p/tmp/maxgelbr/code/…` | What |
| --- | --- | --- |
| `p1-gpu-fused` | `SpeedyWeather-phase1` | H100, fused (default), `--repeats 10` |
| `p1-rea-pertracer` | `SpeedyWeather-phase1-rea` | Reactant H100, per tracer (default) |
| `p1-rea-fused` | `SpeedyWeather-phase1-fuseall` | Reactant H100, fused (`default_fuse_tracers` patched to `true`) |
| `p1-cpu-pertracer` | `SpeedyWeather-phase1-cpu` | x86 CPU node, per tracer (default) |
| `p1-cpu-fused` | `SpeedyWeather-phase1-fuseall-cpu` | x86 CPU node, fused (patched) |

Phase 0 jobs still running (logs in `/p/tmp/maxgelbr/fused-tracers-phase0/`): `ft-bench-gpu`, the full H100 benchmark
with the tracer suite, in the worktree `SpeedyWeather-phase0-gpu` at `ac3c5befb`: when done, copy its
`SpeedyWeather/benchmark/README.md` and `assets/benchmark_results.json` (only the `gpu-nvidia` section/record changes)
into the branch and the tracer table into "Results of phase 0" (`<!-- 0A-… -->` above); `ft-proto-cpu`, the prototype
on an x86 CPU node (`fused_tracers_cpu-x86.md`), a second CPU data point for the CPU default.

To decide then: the CPU default (`p1-cpu-*` and the x86 prototype against the 5% of the gate) and the Reactant
default (`p1-rea-*`, ShallowWater only, after the Reactant failure above is fixed). The worktrees `SpeedyWeather-phase0-*`
and `SpeedyWeather-phase1*` can be removed afterwards (`git worktree remove`).

Other findings: (1) every PrimitiveWetModel construction on Reactant fails (scalar indexing in
`free_troposphere_layers`, `clouds.jl:238`), independent of tracers; (2) with Leapfrog, the first `run!` after
deactivating a tracer still copies its 1st leapfrog step to the 2nd (as without fused tracers), frozen from then on;
(3) compute nodes have no internet: a new worktree needs `RingGrids/Artifacts.toml` (gitignored) copied for the assets;
(4) the benchmark manifest was updated (`Pkg.update`, CUDA stayed at 6.4.1), the old one is in
`/p/tmp/maxgelbr/benchmark-Manifest.toml.bak-20261117`.

### 1. Fused tracer variables

`variables(::Tracer, name, model)` declares fuse symbols, so that all tracers are always allocated in fused parents.
This costs no extra memory, and per-tracer code keeps working on the views.

| Variables | Fuse group |
|---|---|
| prognostic spectral tracers | `:tracers_spectral` |
| grid tracers | `:tracers_grid` |
| spectral tracer tendencies, spectral fluxes `u*A`, `v*A` | `:tracer_tendencies_spectral` |
| grid tracer tendencies, grid fluxes `u*A`, `v*A` | `:tracer_tendencies_grid` |

The fluxes become per-tracer scratch members (e.g. named `Symbol(name, :_u_flux)`) instead of the shared
`vars.scratch.a`/`b`, so that one forward transform covers the tendency and both fluxes of all tracers. That's 2N
extra grid and 2N extra spectral fields of L layers each. Alignment checks as for the dycore's pairs.

### 2. Fused drivers

At each of the sites that now call `foreach_active_tracer`, call one fused function instead when the fused path
applies (all tracers active, see open questions, and architecture preference from phase 0):

- `transform!`: one batched inverse transform `vars.fused.tracers_spectral` → `vars.fused.tracers_grid`.
- horizontal advection: one `+A*div` kernel and one flux kernel over N·L layers, the forward transform(s), one spectral
  flux divergence over N·L layers.
- vertical advection: one launch over N·L layers.
- horizontal diffusion, time step (`update_prognostic!` is elementwise and works on the fused parent as it is),
  moving the grid step back and scaling the grid tendencies: one call each on the parents.

`foreach_active_tracer` and the per-tracer functions stay for the other case (some tracers deactivated).

### 3. Kernels over stacked tracers

The tracers are stacked along the layer axis, but `u`, `v`, the divergence, `w` and the diffusion coefficients have L
layers. Two ways, picked per kernel:

- index with `k = mod1(k_stacked, L)` for the shared fields (horizontal advection, diffusion);
- for vertical advection, whose stencil must not cross from one tracer into the next: reshape the parent data
  `(npoints, N·L, nsteps)` → `(npoints, L, N·nsteps)` (free, column-major) and fold the tracer index into the step
  index the existing kernels already take (`ξ[ij, k, s]`). The existing advection schemes are then reused unchanged.

### 4. Batch sizes

A helper `transform_tracers!` splits the fused transform into chunks of whole tracers that fit the transform's
largest batch `S.nlayers` (set by `transform_batch`), i.e. one call for up to 9 tracers by default (more for the
3N·L batch), as `transform_slabs!` in the prototype. The transform's construction doesn't depend on the tracers for
now. (Possible later: the `tracers` field comes before `spectral_transform` in the model structs, so the default
transform could plan the tracer batch sizes.)

### 5. Versioning, changelog

`SpeedyWeather` only, already at `0.25.0-DEV`. Changelog line under `## Unreleased`. If phase 0 adds the benchmark
suite, that goes into the same PR or a separate one.

## Testing and verification

- **Same results:** the fused path against the per-tracer path, for BarotropicModel, ShallowWaterModel,
  PrimitiveDryModel and PrimitiveWetModel with 3 tracers over 20 steps. Bitwise on CPU if the transforms do the same
  chunks (to be confirmed), otherwise to round-off; on GPU to round-off.
- Deactivating one of several tracers switches to the per-tracer path, and the deactivated tracer stays bitwise frozen
  (existing test).
- More tracers than fit into one transform batch (e.g. 10 with L = 8) run in chunks with the same results.
- JET runtime-dispatch test with 3 tracers; Enzyme forward vs reverse with 2 tracers; Reactant correctness case with 2
  tracers; the GPU tests (`test/GPU/`) with 3 tracers on CUDA, AMDGPU and Metal.
- The phase 0 benchmark suite again, against the per-tracer results, on NVIDIA GPU, CPU and Reactant GPU.

## Documentation changes

- `docs/src/tracers.md`: a short note that tracers are advected together (batched) when all are active, which is
  faster on GPU, and that deactivating a tracer switches to per-tracer processing.
- Benchmark README: the new tracer suite.

## Known limitations

- With a deactivated tracer, all tracers are processed one by one (if the recommendation for open question 2 is taken).
- The tracers are batched among themselves, not with the dycore's variables.
- More tracer memory: 2N grid and 2N spectral flux fields instead of two shared ones.

## Future work

- Merge the tracers into the dycore's batches (`:prognostic`/`:grid`, `:spectral_tendencies`/`:grid_tendencies`), so
  that tracers add no transform calls at all. Needs fuse groups across namespaces and has to respect the different
  step choice of 2D models.
- Retire the sequential `horizontal_advection!` in `tendencies_sequential.jl` once tracers no longer use it.

## Open questions

Decided 2026-11-17 (see revision log):

1. All tracers are fused automatically, no per-tracer `fuse` option.
2. With any tracer deactivated, all tracers are processed one by one (the existing per-tracer path).
3. The go/no-go thresholds as proposed: 0a stops below a 20% cost of 8 tracers at T ≤ 128 on NVIDIA; 0b goes at a
   ≥ 20% shorter time step with 4 to 8 tracers at T ≤ 128, not slower at T256, at most 5% slower on CPU.
4. Transforms wider than the transform's largest batch (`S.nlayers`, from `transform_batch`) are split into chunks of
   whole tracers; no tracer-dependent transform construction for now.
5. The phase 0 edits are made on this machine; the NVIDIA runs are done on a GPU machine (see "Handover").
