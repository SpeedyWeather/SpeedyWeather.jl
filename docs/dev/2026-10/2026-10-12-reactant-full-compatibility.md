# Full Reactant compatibility incl. the primitive equation models

> Status: **in progress**. All models run on Reactant (CPU backend), CPU-vs-Reactant correctness tests pass. On Reactant+CUDA the PrimitiveWetModel is wrong (upstream layout bug, see "GPU (Reactant + CUDA)"); making the implicit kernels raisable as a workaround.

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

(2026-10-13, GPU session) Ran the PrimitiveWetModel with Reactant on CUDA (A40, login02). It goes
wrong because of an upstream bug (kernel operand layouts, EnzymeAD/Reactant.jl#3377), see
"GPU (Reactant + CUDA)". `raise = true` avoids that bug but failed on the implicit solver kernels.
That's why they are being rewritten.

> We need to continue the Reactant investigation that got aborted. You were just in the process of
> rewriting the implicit correction so that the kernel can be raised.

> make regular updates to the progress markdown file, so we don't lose any progress

> We might now replace the `while` loop again with a for loop as the AMD issue should be fixed

Implicit kernels raise with `@inbounds` in the closures and a zero-factor Gauss-Jordan elimination.
The `while` loop is replaced by `for` again. On the GPU, the raised correction kernel runs a matvec as
a TF32 Triton GEMM (see "GPU (Reactant + CUDA)").

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

## GPU (Reactant + CUDA)

Investigation on login02 (NVIDIA A40), Reactant 0.2.292 / Reactant_jll 0.0.418, environment
`SpeedyWeather/test/GPU/CUDA`. Scripts and HLO dumps are in
`/p/tmp/maxgelbr/code/Speedy2/claude_reactant_gpu_2026-10-13/`.

Findings (2026-10-13)
- The PrimitiveWetModel initializes and steps on CUDA but is NaN by step 2 (land points first).
  Root cause upstream: the `enzymexla_compile_gpu` custom call that KernelAbstractions kernels lower
  to has no `operand_layouts`, so XLA's layout assignment can hand a kernel an operand with a
  non-default layout (seen: `c64[2,33,560]{0,2,1} bitcast(...)` into `UV_from_vordiv!`), and the
  kernel reads scrambled data. It shows up in `transform!(vars, model)` when the leapfrog step copy,
  `UV_from_vordiv!` and the batched matrix transform are compiled together. Fixed by the open PR
  EnzymeAD/Reactant.jl#3377 (verified with an MWE). Issue draft:
  `reactant-issue-draft-kernel-operand-layout.md`, MWE `reactant-issue-kernel-operand-layout-mwe.jl`.
- Float32 matmuls use TF32 on the A40 (~3e-4 relative transform error).
  `Reactant.with_config(dot_general_precision = PrecisionConfig.HIGHEST)` fixes that.
- Workaround until #3377 is released: raise the kernels (`@compile raise = true`), XLA then owns
  the layouts. `@jit raise = true` of the transform is correct, but `@compile raise = true time_step!`
  failed raising the implicit solver kernels:
  - `implicit_invert_S_kernel!` raised to wrong results (the `if k != pivot` around the row
    elimination). Now branch-free: the factor is zero for the pivot row (CPU bitwise identical).
  - `implicit_primitive_leapfrog_kernel!` failed to raise ("'affine.load' op operand cannot be used as
    a symbol") in the `ntuple` closure of step 2 (geopotential integration with a `while` loop).
    `inbounds = true` of `@kernel` doesn't reach into closures, so the bounds checks (early exits)
    remained. Fixed with `@inbounds` inside the closures. The `while` loop (#1193, workaround for
    AMDGPU.jl#1015) is back to `for r in k:nlayers`: #1015 was fixed in LLVM (llvm-project#215829) and
    backported to the AMDGPU LLVM backend jll (Yggdrasil#14465). Not verified on an AMD GPU here, and
    the AMDGPU compat (`2.1.2` in SpeedyWeatherInternals) still allows versions without the fix.
    Variants tested on the CPU backend (`leapfrog_variants.jl`): `while` without
    `@inbounds` (HEAD) reproduces the error; `while`, `for r in k:nlayers` and a masked `for r in 1:nlayers`,
    all with `@inbounds`, raise and agree with the CPU (rel 8e-8 to 9e-8), and are bitwise identical to the
    source kernel on the CPU.
  - Kernel test (`implicit_kernels_raise.jl`: model-like views into fused parents, vs the CPU).
    Working tree: zero-factor Gauss-Jordan, `for r in k:nlayers` + `@inbounds` in the closures.
    Both kernels raise. On the CPU backend both match the CPU (rel 1e-7). On the GPU backend
    `implicit_invert_S_kernel!` matches (1.6e-6), but the raised `implicit_primitive_leapfrog_kernel!`
    differs by 9.7e-5, also with `dot_general_precision = HIGHEST`.
  - Why (`leapfrog_raise_diag.jl`): raising turns the vertical matvecs of steps 1 and 2 (L·ΔD, R·G_T)
    into 16 `stablehlo.dot_general`s without a `precision_config`. Reactant's `dot_general_precision`
    only applies to dots it traces itself, not to dots created by the raising/optimization passes.
    The error (1.5e-4 vs a Float64 reference, CPU 1.7e-7, spread over most elements) is consistent
    with TF32. `NVIDIA_TF32_OVERRIDE=0` has no effect. No Enzyme-JAX pass sets dot precision.
    In XLA's optimized HLO the step-1 dot (real and imaginary parts merged, `f32[1120,8]`) is a
    `__triton_nested_gemm_fusion`; the triangular vector dots became reductions.
  - Workaround: `xla_gpu_enable_triton_gemm = false`, either as a compile option
    (`@compile compile_options = Reactant.CompileOptions(; raise = true, xla_debug_options = (; xla_gpu_enable_triton_gemm = false)) f(...)`,
    the macros don't take `xla_debug_options` directly) or via `XLA_FLAGS=--xla_gpu_enable_triton_gemm=false`. The raised kernel then matches the CPU (rel 1.6e-7,
    same distance to Float64 as the CPU). Upstream fix candidate: dots created by raising should carry
    Reactant's `dot_general_precision`.

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
- On Reactant the implicit operators are recomputed every time step (`reinitialize!`, kept for now):
  about 0.5 ms per step (PrimitiveDry T31L8 7.3 vs 6.7 ms, T63L16 within noise) and ~8 s compile
  time, although the time step only changes in the first steps (Leapfrog: Δt/2, Δt, 2Δt).
  Alternatives: a `@trace if` that only captures arrays (failed when it traced the model), or
  reinitializing outside the compiled time step during the spin-up steps in `time_stepping!`.
- Tighten the tolerances of the correctness tests, differences are orders of magnitude smaller
  than the current tolerances.
