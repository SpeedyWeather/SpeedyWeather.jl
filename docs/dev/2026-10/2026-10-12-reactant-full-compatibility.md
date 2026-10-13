# Full Reactant compatibility incl. the primitive equation models

> Status: **in progress**. All models run on Reactant (CPU backend), CPU-vs-Reactant correctness tests pass; awaiting review, GPU (Reactant+CUDA) untested.

Date of initial draft: 2026-10-12

Base revision: c0acb2ca (main at 7c1beecb incl. #1302 merged into `mg/new-reactant`)

## Originating prompt

> We want to achieve full compatability with Reactant. Reactant was updated just a few hours ago to
> be able to deal with complex valued kernels. This was the one feature we were missing. We already
> started a lot of this work on `mg/reactant-primitive` / #985. However this branch is stale and out
> of date now. We want to start in a fresh branch but take a lot of the changes from
> `mg/reactant-primitive`. Start by taking over the `test/reactant` changes that add correctness
> tests. Then inspect the Reactant extensions. Then the rest of the code changes. Apply what's
> necessary. Then, try to run the Reactant correctness tests and fix upcoming issue.

## Revision log

> I just merged a fix for Reactant MatrixSpectralTransform into `main`. Merge that into this branch
> as well and take note #1302

Merged `main` (squash merge of #1302, identical to the already merged branch).

> in case upstream changes in Reactant are advisable, do them first here in the Reactant extension
> and mark them with #TODO: move upstream to Reactant

Fixes that belong into Reactant are in our Reactant extensions, marked `# TODO: move upstream to Reactant`.

> Next, we want to test if the differentiation now works. We need to revise the differentiation.jl
> test for that.

`test/reactant/differentation.jl` compares the gradient of one time step (raised, Enzyme reverse)
against Enzyme on the CPU, currently `@test_broken` because of two upstream issues (complex constants
can't be raised, no reverse rule for a masked reduction). Found and fixed (in the extension) that
`deepcopy` of concrete Reactant arrays breaks the aliasing of views into fused parents, which made
`make_zero(::Variables)` return shadows holding the primal values.

> Upstream in Reactant `main` the issue with complex constant should be fixed. Apply the kernel
> workaround for SpeedyTransforms

`_UV_from_vor_kernel!` and `_UV_from_vordiv_kernel!` are branch-free: the special cases (first
element, last two rows) zero the coefficients, the neighbours lm±1 wrap around with `rem` (the only
index clamp Reactant can raise) and are multiplied by zero there. Raised, their branches became a
masked reduction without a reverse rule in Enzyme-MLIR. Results are bitwise identical on the CPU,
the raised gradients agree with finite differences. The kernels themselves are 20-30% slower on
1 CPU thread (5-15% on Reactant), within noise for a full time step (Barotropic, ShallowWater T63,
PrimitiveDry T31).

> Now we have to make the benchmarks runnable with Reactant. Also inspect how we do that in
> Terrarium.jl as a reference

As in Terrarium.jl's benchmarks, architecture-specific hooks (`build_runner`, device
synchronization, memory footprint, compiled calls) with the Reactant methods in
`benchmark/reactant_runner.jl`, included only when Reactant is loaded. Unlike Terrarium (one
compiled loop with a fixed number of steps) the time step is compiled once per configuration,
reused for the warm-up and the timed run, and the compile time is reported in its own column.
Timing covers only the time loop (as the README states), synchronizing Reactant by reading the
clock (`KernelAbstractions.synchronize` is a no-op for Reactant). Reactant runs use the
`MatrixSpectralTransform`, configurations too large for it (T>150) and resulting duplicates are
dropped. `Base.summarysize` misses Reactant's device buffers, they are added. A failing
configuration is recorded as failed instead of aborting the suite. Fixed: the Reactant models were
built with a progress meter as repeated keyword arguments silently take the last value.

> Can you find from where the warning about type-promotion from ifelse comes from exactly and fix it?

Reactant's `ifelse` warning came from `max`/`min` with integer literals in broadcasts over traced
arrays (sea ice `filter!`, Reactant traces the literal as an `Int64` number). Literals of the arrays'
number format there and in the land-sea mask clamping, Held-Suarez forcing and 10m wind output.

> yes look into them

Two Reactant failures found on the way, not covered by tests before:
- Held-Suarez forcing: its initialization broadcast Float64 constants into Float32 Reactant arrays.
  Reactant then swaps in the Float64 result buffer without converting, any later compiled function
  using the array fails (RuntimeProgramInputMismatch, a segfault in a plain Reactant MWE). The
  constants are now in the forcing's number format (`oftype`). The Reactant bug itself (eager
  broadcast into a concrete array with a different element type) is not worked around, to report
  upstream.
- 10m wind output: `pressure(nlayers, pₛ, coordinate)` in a broadcast, Reactant traces `nlayers` so
  `pressure(k::Integer, ...)` doesn't apply. The bottom layer height is now computed in a kernel.

## Problem description

`mg/reactant-primitive` (#985) got the primitive equation models close to running with Reactant but
diverged from `main` (300 commits, ~80 files). Reactant v0.2.292 now supports complex numbers in
kernels (Reactant.jl#3438), the last missing feature. The goal is to run every model on Reactant
with results matching the CPU.

## Background

- Reactant compiles `time_step!(simulation)` once (`@compile`) and calls it every time step, see
  `SpeedyWeatherReactantExt.time_stepping!`. Everything inside is traced: no data-dependent control
  flow, no transfers to the CPU, no scalar indexing. Numbers in the `Clock` are traced (concrete
  Reactant numbers), other numbers in the model are compile-time constants.
- Outside of `@compile` (`initialize!`), Reactant arrays are concrete, operations on them are
  compiled on the fly (KernelAbstractions kernels, broadcasting) or wrapped in `@maybe_jit`.
- KernelAbstractions kernels are lowered natively on the CPU/CUDA backends (not raised), so
  control flow inside kernels is fine.

## Summary of changes

From `mg/reactant-primitive`, only the changes required on top of `main` were taken over. Purely
"more parametric" changes (#1041, `Bool`/`Int`/`Second` fields as type parameters) were dropped as
plain numbers are not traced by `@compile`.

Model and numerics (all architectures)
- `Geometry`, `SigmaCoordinates`, `SigmaPressureCoordinates`: separate type parameters for vectors
  of different length (Reactant's kernel arrays carry their size as a type parameter).
- `ImplicitPrimitiveEquation`: new `initialize_traceable!` (broadcasting + Gauss-Jordan kernel)
  used for `ReactantDevice`; `Δt` is a `scalar_ref` (traced on Reactant). On Reactant
  `reinitialize!` recomputes the operators every step from a temperature profile that is only
  updated if the time step changed (no branching), equivalent to `main`'s behaviour.
- Leapfrog: the step `lf` is selected inside the kernel instead of a view with a traced index.
- Time step scaling without `oftype` (traceable), `secondofday` for any `AbstractDateTime`,
  zenith functions accept `Dates.AbstractDateTime`.
- `AbstractCPU` dispatch instead of `CPU`/`GPU` (geopotential, vertical integration),
  `set_scalar!` (SpeedyWeatherInternals) for single-element writes.
- `inv(x)` → `1 / x` in kernels (Reactant.jl#2560, Reactant+CUDA only).
- Logical indexing replaced by broadcasting (ocean, land), `@maybe_jit interpolate!` for boundary
  conditions, land-sea mask `mask!` for array masks, `extrema` → `minimum`/`maximum` (Reactant.jl#2387).
- RingGrids interpolation: no `undef` allocations, no reshape of matrices (`flatten_layers`),
  generic `average_on_poles` returning `1 × nlayers` matrices (no `vec`).

Reactant extensions
- `SpeedyWeatherReactantExt`: output and feedback (progress meter) between compiled steps
  (`do_output!` returns false, `progress!` returns inside `@compile`), `Feedback` not traced
  (its IO buffers can't be), `reinitialize!` for the implicit solver, `scalar_ref`.
- `RingGridsReactantExt`: `FieldReactantStyle` broadcasting (out-of-place broadcasts return
  fields), bounds checks, branchless ring search, `update_locator!` without safety checks.
- `LowerTriangularArraysReactantExt`: fix `BroadcastStyle` for the 5-parameter type.
- Upstream candidates (`# TODO: move upstream to Reactant`): broadcasting into views of concrete
  arrays (silently wrong results otherwise), `Dates.year` for traced days, `deepcopy` of concrete
  arrays keeping views aliased with their (copied) parent.

## Testing and verification

- `SpeedyWeather/test/reactant`: CPU-vs-Reactant correctness (tendencies after one step,
  prognostic and grid variables after 100 steps following a 10-day spin-up) for
  `BarotropicModel`, `ShallowWaterModel`, `PrimitiveDryModel`, `PrimitiveWetModel`
  (constant longwave transmissivity, no convection to avoid chaotic amplification), NetCDF and
  Zarr output, `MatrixSpectralTransform` outside of compilation, `@maybe_jit`.
- Default models (incl. Betts-Miller convection, Frierson longwave, default feedback with progress
  meter) run without NaNs; one-step tendencies agree to ~1e-4 (norm-relative) with convection active.
- CPU test suites of SpeedyWeatherInternals, LowerTriangularArrays, RingGrids, SpeedyTransforms pass.
- New unit tests: `initialize_traceable!` vs `initialize!`, `secondofday`, `set_scalar!`,
  no-op conversions of `Field`/`LowerTriangularArray`.

## Known limitations

- Reactant on GPU (CUDA) is not tested in this PR, `SpeedyWeather/test/GPU/reactant.jl` is still
  not included in the GPU test suite.
- Particles are not supported on Reactant.
- Differentiation with Reactant (Enzyme reverse mode, kernels raised to StableHLO) is blocked
  upstream until a Reactant release includes the fix for raising kernels that store complex
  constants (`div[I] = 0` in `_divergence_kernel!`, fixed on Reactant's main). With a computed zero
  there instead (tried outside of this PR only) the gradient of a `BarotropicModel` time step agrees
  with Enzyme on the CPU to a relative error of 1e-6. ShallowWater and primitive equation models
  are untested.
- The prescribed seasonal climatologies (e.g. `SeasonalOceanClimatology`) are untested inside
  the compiled time step.

## Future work

- Move the `# TODO: move upstream to Reactant` fixes into Reactant, remove the `inv` and `extrema`
  workarounds once Reactant.jl#2560 and #2387 are fixed.
- Re-enable `@trace` for the time loop (checkpointing) in `time_stepping!`.
- Tighten the tolerances of the correctness tests, differences are orders of magnitude smaller
  than the current tolerances.
