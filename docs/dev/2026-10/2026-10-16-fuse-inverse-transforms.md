# Fusing the U, V and ∇ln pₛ spectral → grid transforms (3 → 2 per time step)

> Status: **in progress**. Implemented on `mg/fuse-more-inverse` (not committed); CPU verification
> under way, GPU runs and benchmarks outstanding.

Date of initial draft: 2026-10-16

Base revision: `40f13c90` (`mg/fuse-more-inverse`, identical to `main`)

## Originating prompt

> Could we fuse the three inverse transforms per timestep? Analyse and make a plan with pro and contra.

## Revision log

> I only want A, partial

- 2026-10-16: initial draft with three options: A (partial, 3 → 2 calls), B (full, 3 → 1 by changing
  the variable layout), C (full, through staging buffers). The recommendation was B in two phases
  with A first.
- 2026-10-16: scoped to Option A only on request. B and C are kept below under *Alternatives
  considered*.

> okay implement it, go ahead, commit nothing yet

- 2026-10-16: implemented as planned. Deviations and details:
  - `pressure_gradient_flux!` is now `pressure_gradient_flux!(vars, time_stepping)`. The call in
    `benchmark/benchmark_suite.jl` is updated. `benchmark/CUDA/benchmark_suite.jl` is left alone: it
    already used an outdated signature (`lf` argument, `vars.grid.temp`) before this change. The
    `pressure_gradient_flux!` row of the benchmark tables now measures only the flux; the gradient
    and its transform are part of `transform!`.
  - The unit test is `test/dynamics/pressure_gradient_transform.jl`. It only covers
    `PrimitiveWetModel`, since dry and wet share the code path, which halves its compile time. It
    sets random vorticity, divergence and pressure first, because the initial surface pressure has no
    gradient. It compares with `≈` (`rtol = 1e-5`) for both transforms; bitwise equality on CPU is
    checked by the 1-day regression run instead.
  - `test/GPU/vertical_integration.jl` needed no change: it copies `dpres_dx/dy` with `.=` between
    two fields that now both carry steps.

> call the fuse groups "uv_dpres_grad_grid" and "uv_dpres_grad_spectral"

- 2026-10-16: renamed `:uv_grid` → `:uv_dpres_grad_grid` and `:spectral_scratch` →
  `:uv_dpres_grad_spectral` on request, also for the 2D models, which declare these groups.

> Drop 2

- 2026-10-16: removed `2` from the GPU `default_transform_batch`. After this change no model emits
  K = 2 when nlayers > 1, and for nlayers = 1 the `2 * nlayers` entry already covers it. The other
  entries are each still used by one model. PrimitiveDry's prognostic batch (3L+1) and ShallowWater's
  tendency batch (3L+2) remain unlisted and are planned on first use, as before.

## Problem description

A primitive-equation time step does **one** grid → spectral transform: all tendencies go in one batch
of 9L+1 columns (wet) or 6L+1 columns (dry). It does **three** spectral → grid transforms:

| # | call site | spectral input → grid output | K (wet / dry) | `unscale_coslat` |
|---|---|---|---|---|
| 1 | `transform!(vars, model)`, [transform.jl:149-151](../../../SpeedyWeather/src/time_stepping/transform.jl) | `:prognostic` (vor, div, temp, pres, [humid]) → `:grid` | 4L+1 / 3L+1 | no |
| 2 | same function, [transform.jl:162-166](../../../SpeedyWeather/src/time_stepping/transform.jl) | `:spectral_scratch` (U, V in `scratch.a/b`) → `:uv_grid` (u, v) | 2L | yes |
| 3 | `pressure_gradient_flux!` in `dynamics_tendencies!`, [tendencies.jl:199-203](../../../SpeedyWeather/src/dynamics/tendencies.jl) | `:dpres_grad_spec` → `:dpres_grad` (∇ln pₛ) | 2 | yes |

Calls 2 and 3 both unscale by 1/cos(lat), so they can be a single call of width 2L+2. Three things
prevent that today:

1. **Step dimension.** `:uv_grid` (`XYZT`) has the leapfrog step axis, but `dpres_dx/dy` are 2D. A
   fuse group that has a member with a step axis rejects 2D members
   ([dimensions.jl:294](../../../SpeedyWeather/src/variables/dimensions.jl)).
2. **Separate fuse groups.** U, V sit in `:spectral_scratch` and the pressure gradient in its own
   `:dpres_grad_spec` group.
3. **Call site.** Call 3 happens at the start of the next step's `dynamics_tendencies!`, not in
   `transform!`. It still reads the same state: spectral pressure at the current step, which nothing
   changes between the end-of-step `transform!` and the next `dynamics_tendencies!`.

Call 1 is not unscaled, has a step axis on both sides, and is not part of this plan (see
*Alternatives considered*).

## Background

### What removing the K = 2 call buys, per transform

- **`MatrixSpectralTransform`.** This is the default on GPU for `truncation <= 100` (`WhichTransform`)
  and on Reactant. One call is one real `gemm` with `backward_stacked` (npoints × 2 nharmonics), and
  the whole matrix is read on every call whatever K is. On the H100 (timings from
  [2026-10-02-matrix-transform-blas-forward.md](2026-10-02-matrix-transform-blas-forward.md)) a T64
  call costs 0.134 / 0.177 / 0.285 ms for K = 8 / 33 / 73. So about 0.12 ms per call does not depend
  on K. Merging the K = 2 call into the K = 2L call saves one full read of the matrix.
- **`SpectralTransform` on CPU (the CPU default).** Wide calls are split into the pre-planned FFT batch
  sizes `[1, L]` (`_needs_chunking`, `_transform_chunked!`,
  [spectral_transform.jl:452-517](../../../SpeedyTransforms/src/spectral_transform.jl)). Today the
  2L call becomes L, L and the K = 2 call becomes 1, 1. The fused 2L+2 call becomes L, L, 1, 1: the
  same chunks, the same plans, and per column the same arithmetic. **No speedup is expected, and
  results should be bitwise identical to `main`.**
- **`SpectralTransform` on GPU (default for `truncation > 100`).** It saves one Legendre launch, one
  unscale launch and one graph replay of about 2·nlat_half + 4 nodes (the K = 2 gather, per-ring
  cuFFTs and scatter). That is about 200 of the 1283 kernels per step at T127 counted in
  [gpu-idle-time-reduction.md W3.1](../2026-08/gpu-idle-time-reduction.md). These are small-batch
  cuFFTs. Expect a few percent at most; this has to be measured.

### Estimated effect for the matrix transform (H100)

The rows below are **estimates, not measurements**. Spectral → grid times at K = 2, 16 and 18 are
linearly interpolated from the measured K = 8 / 33 / 73. Step times come from the matrix-transform
SYPD on the same H100 with the default Δt (2400 / 1200 / 900 s), using ms/step = Δt·1000 / (365·SYPD).
The model is `PrimitiveWet`, L = 8.

| T | today: K = 16 + 2 | fused: K = 18 | saving | ms/step | ≈ share of the step |
|---|---|---|---|---|---|
| 32 | 0.110 ms | 0.058 ms | 0.05 ms | 1.13 | 5 % |
| 64 | 0.272 ms | 0.151 ms | 0.12 ms | 1.16 | 10 % |
| 86 | 0.656 ms | 0.389 ms | 0.27 ms | 2.25 | 12 % |

### Fusion constraints in the variable system

- A fuse group is keyed by (namespace, fuse symbol) and can span variable types (`ScratchVariable` and
  `DynamicsVariable` can share one). In a 3D parent (`:spectral_scratch`, `LMZ`) a 2D member takes
  one slot with no step axis. In a 4D parent (`:uv_grid`, `XYZT`) a `GridXYT` member takes one slot
  with the step axis, like `pressure` in `:grid`.
- **There is a precedent on the tendency side.** Intermediates such as `uT_anomaly` and
  `kinetic_energy` are declared with the tendency step dimension only so that they fit into
  `:grid_tendencies`.
- Nothing except `transform.jl` uses `parent(vars.fused.spectral_scratch)` or
  `parent(vars.fused.uv_grid)` as a whole. `move_prognostic_grid_variables_back!` copies the
  `:uv_grid` parent step 2 → 1 ([leapfrog.jl:134](../../../SpeedyWeather/src/time_stepping/steppers/leapfrog.jl)),
  so it will also copy the two `dpres` slots.
- All six functions that read `dpres_dx/dy` already get `time_stepping`: `pressure_gradient_flux!`,
  both `vertical_integration!` methods, `surface_pressure_grid_tendency!`, `vertical_velocity!` and
  `vordiv_grid_tendencies!`.

## Summary of changes

1. **Variables**, [primitive_dry.jl:166-169](../../../SpeedyWeather/src/models/primitive_dry.jl):
   - `dpres_dx`, `dpres_dy`: `Grid2D()`, `fuse = :dpres_grad` → `GridXYT(pg)`, joining u, v.
   - `dpres_dx_spec`, `dpres_dy_spec`: `fuse = :dpres_grad_spec` → joining `scratch.a/b` (still
     `Spectral2D`).
   - The two joined groups are renamed: `:uv_grid` → `:uv_dpres_grad_grid` (now [u, v, dpres_dx,
     dpres_dy], 2L+2 slots per step) and `:spectral_scratch` → `:uv_dpres_grad_spectral` (now [a, b,
     dpres_dx_spec, dpres_dy_spec], 2L+2 slots). The names are declared in `barotropic.jl`, so the 2D
     models use them too, holding only [u, v] and [a, b].
   - The `:dpres_grad` and `:dpres_grad_spec` groups disappear.
2. **`transform!(vars, ::PrimitiveEquation)`**, [transform.jl](../../../SpeedyWeather/src/time_stepping/transform.jl):
   - after `UV_from_vordiv!`, compute `∇!(dpres_dx_spec, dpres_dy_spec, pres, S)` with
     `pres = get_prognostic_step(vars.prognostic.pressure, time_stepping, S)`. For leapfrog that is
     the same step as `DynamicalCore()` picks today;
   - the existing `:uv_dpres_grad_spectral` → `:uv_dpres_grad_grid` call with `unscale_coslat = true`
     now covers 2L+2 columns;
   - update the comment there and the TODO at line 160.
3. **`dynamics_tendencies!`**, [tendencies.jl](../../../SpeedyWeather/src/dynamics/tendencies.jl):
   - `pressure_gradient_flux!` keeps only the flux `(u, v)·∇ln pₛ`. Its `∇!` and transform go, and
     its now unused `S` argument can go too;
   - the six readers of `dpres_dx/dy` (lines 194, 242, 292, 400, 453, 564) select the current step:
     `dpres_dx = get_prognostic_step(vars.dynamics.dpres_dx, time_stepping, DynamicalCore())`, and
     the same for `dpres_dy`. The kernels keep receiving 2D views, so their bodies do not change;
   - the docstring of `surface_pressure_grid_tendency!` (line 396) says ∇ln pₛ comes from
     `pressure_gradient_flux!`; change that to `transform!`.
4. **GPU FFT plans**, [spectral_grid.jl:373-375](../../../SpeedyWeather/src/dynamics/spectral_grid.jl):
   add 2L+2 to the GPU `default_transform_batch`, drop 2 (unused now, see the revision log), and
   update its docstring. The scratch is large enough because 2L+2 ≤ 9L+1. CPU is unchanged.
5. `CHANGELOG.md` entry. `SpeedyWeather` is already `0.23.0+DEV`, and no other submodule changes.

What follows from this:

- **2 spectral → grid transforms per primitive-equation step instead of 3**, not counting tracers and
  the random pattern.
- The memory change is negligible: `dpres_dx/dy` gain a second step, which is two extra 2D grid
  columns.
- ∇ln pₛ is now ready right after `transform!`, so parameterizations could read it if they ever need
  to.

## Testing and verification

- **New unit test** `SpeedyWeather/test/dynamics/pressure_gradient_transform.jl`. For `PrimitiveDryModel`
  and `PrimitiveWetModel` (T31, L = 4), with both `SpectralTransform` and `MatrixSpectralTransform`:
  - after `initialize!`, compare u, v with references from separate `transform!(…, U/V; unscale_coslat = true)`
    calls, and `dpres_dx/dy` with `transform!` of `∇(pres)` unscaled. Use `==` for the CPU FFT
    transform and `isapprox` for the matrix transform;
  - assert the member order of `vars.fused.uv_dpres_grad_grid` and `vars.fused.uv_dpres_grad_spectral`
    (slot map keys),
    and that `vars.fused` has no `:dpres_grad` or `:dpres_grad_spec`.
- **Bitwise regression on CPU**: a 1-day `PrimitiveWetModel` run with the default transform, comparing
  all prognostic variables `==` against `main`. Same chunks and plans, so this should hold.
- **Existing tests**: `variables/fuse.jl`, `dynamics/time_stepping.jl`, `dynamics/dispatch.jl` (JET);
  on an H100 with the CUDA env, `test/GPU/runtests.jl` for both transforms.
  `test/GPU/vertical_integration.jl` copies `dpres_dx/dy` with `.=` between two fields that now both
  have steps, so it should work unchanged; check it.
- **Differentiability (extended tests)**:
  - the `transform!(::Variables)` AD test now also produces ∇ln pₛ, so add `dpres_dx/dy` to its
    outputs (they live in `vars.dynamics`, not in `outgroup = :grid`);
  - the `dynamics_tendencies!` AD test (`ingroup = :grid`) should still pass, but it no longer covers
    the ∇ln pₛ transform, because ∇ln pₛ now arrives as an input;
  - run x86 through workflow_dispatch, since the Enzyme segfaults there do not reproduce locally.
- **Benchmarks**:
  - a micro-benchmark of `transform!(vars, model)` alone, before and after, on CPU and the H100 with
    both transforms; this is the direct check of the estimates above;
  - then `manual_benchmarking.jl` on CPU (expected within ± 20 % noise of the README) and on GPU (the
    MT rows T32–T86 should gain ≈ 5–12 %; the FFT rows T128+ are measured, not predicted).

## Documentation changes

- `default_transform_batch` docstring (adds 2L+2).
- Comments in `transform.jl` and `pressure_gradient_flux!`, and the docstring of
  `surface_pressure_grid_tendency!`.
- No user-facing docs page lists the fuse groups, so `docs/src` needs no change.

## Known limitations

- Still two spectral → grid calls. The prognostic batch (4L+1) stays separate because it has no
  coslat unscaling and its spectral side carries the step axis.
- No gain on CPU with the default `SpectralTransform`: the chunking into planned batches makes fused
  and separate calls do identical work.
- `vars.dynamics.dpres_dx/dy` now have a step axis. That is internal but visible through `vars`, and
  their previous-step copy is unused.
- ∇ln pₛ is computed in `transform!`, so it also runs when `model.dynamics = false` (negligible cost).
  A callback that changes spectral pressure without transforming again sees a stale ∇ln pₛ, which is
  the same contract as for every other grid variable.
- The gains in this plan are interpolated from stand-alone transform timings and have not been measured
  in the model yet.

## Alternatives considered (not pursued)

- **Full fusion through the variable layout (3 → 1 call).** U, V and ∇ln pₛ would join the
  `:prognostic` / `:grid` parents with the step axis, and the coslat unscaling would move to grid
  space after the transform. The estimated matrix-transform saving is about 9–23 % (H100, T32–T86),
  about twice A's. It changes the variable layout of all four models, removes `:uv_grid` and
  `:spectral_scratch`, and makes u, v and ∇ln pₛ differ from `main` at ULP level. Not pursued, on
  request. This plan is a strict subset of that one, so it is not wasted work if it is ever picked up.
- **Full fusion through staging buffers.** Copy every input into a buffer of width 6L+3, transform
  once, copy the outputs back. No layout change, but it needs an extra grid buffer and an extra pass
  every step, and it repeats slot bookkeeping that `FusedParent` already does.

## Future work

- Fold tracers (one call of width L each) into a batched call. This needs fusion across namespaces,
  since tracers live in `:tracers`.
- On ARM CPUs, where the FFTW alignment problem behind the chunking does not occur, check whether
  batched plans for the model's call widths are faster than the chunks.
