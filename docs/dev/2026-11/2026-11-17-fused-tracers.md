# Fuse the tracers: one batched transform and one kernel per operation for all tracers

> Status: **planned**. Measurement-gated: phase 0 decides with GPU benchmarks whether phase 1 (the implementation)
> is done at all.

Date of initial draft: 2026-11-17

Base revision: 174afe29 (`mg/unroll-tracer-loops`) plus the uncommitted part 2 of
`2026-11-16-unroll-tracer-loops.md` (tracers as a NamedTuple fixed at model construction), which this plan builds on.

## Originating prompt

> okay next, could we actually look into making a plan whether it's possible to also fuse multiple tracers in case
> they are defined with our fuse system for variables and advect them together? This should only be done if it's
> really faster on GPU

## Revision log

- 2026-11-17: initial draft.

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

## Phase 1: implementation (only if phase 0 says go)

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
largest batch, i.e. one call for up to 9 tracers by default (two or three for the 3N·L batch). Optionally the model's
default `spectral_transform` plans the tracer batch sizes: the `tracers` field comes before `spectral_transform` in
the model structs, so its default can depend on the number of tracers (see open questions).

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

Recommendations in **bold**:

1. **Fuse all tracers automatically**, or only those declared with a `fuse` option per tracer
   (`Tracer(fuse = :group)`)? Automatic needs no new user-facing option and costs no memory. Per-tracer groups would
   allow several groups but add an API without a known use.
2. **Deactivated tracers: switch to the per-tracer path when any tracer is deactivated**, or process all tracers fused
   and mask the time step? The switch keeps the current behaviour exactly, and the per-tracer code already exists.
   Masking keeps one path, but recomputes the frozen tracer's grid once (from the frozen spectral state) and wastes work
   on frozen tracers.
3. **Go/no-go thresholds:** 0a stops below a 20% cost of 8 tracers at T ≤ 128 on NVIDIA; 0b goes at a ≥ 20% shorter
   time step with 4 to 8 tracers at T ≤ 128, not slower at T256, at most 5% slower on CPU. Adjust as you see fit.
4. **More than 9 tracers: split the transforms into chunks first**, and only make the model's default transform plan
   the tracer batch sizes if phase 0b shows that one wide batch is clearly faster than chunks.
5. **Phase 0 runs:** I can run `cpu-arm` and `metal` on this machine. The NVIDIA and Reactant GPU runs need the GPU
   machine of the stored results. Who runs them?
