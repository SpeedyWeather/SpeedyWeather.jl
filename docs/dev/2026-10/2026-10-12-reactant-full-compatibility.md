# Full Reactant compatibility incl. the primitive equation models

> Status: **in progress**. See "Handoff 2026-10-16" first. All models run on Reactant 0.2.293 (`raise = true`, central compile options), on the CPU backend and on CUDA; the time loop is one compiled program per output interval and can be differentiated through several time steps. Everything since f5082066 is uncommitted; three upstream Reactant issues are drafted, not posted.

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

> Don't try to run the potential fix, we want to concentrate of getting Reactant to work with `raise = true`

> no, use reactant#main which includes the constant complex fix

> but still rewrite kernels to branchless form when this doesn't come with performance penalties on GPU

`@compile raise = true time_step!` on CUDA next failed on a possibly undefined variable in
`pseudo_adiabat!` (fixed) and on complex-constant stores (Reactant.jl#3450). For the latter, a
separate environment uses Reactant `main` with `Reactant_jll` 0.0.419. `_divergence_kernel!` is
branch-free (no GPU penalty). The branch-free `zero_last_degree_kernel!` is 3.6-4.3x slower at large
sizes.

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
- `ImplicitPrimitiveEquation`: one `initialize!` for all architectures, the operators are built via
  broadcasting on the device (traceable), S⁻¹ via `invert_S!`: LU decomposition on the CPU
  (`invert_S_lu!`, as before) or, for `ReactantDevice`, Gauss-Jordan elimination in a kernel
  (`invert_S_gauss_jordan!`); `Δt` is a `scalar_ref` (traced on Reactant). On Reactant
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
- New unit tests: S⁻¹ inverts S, `invert_S_gauss_jordan!` vs `invert_S_lu!`, `secondofday`, `set_scalar!`,
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
- Full PrimitiveWetModel, `@compile raise = true time_step!` (`reactant_wet_compare_raise.jl`, RAISE,
  PRECISION=highest, NOTRITON). With the implicit kernels raisable, two more failures:
  - `pseudo_adiabat!` (Betts-Miller convection): `temp_parcel_dry` was assigned only in `if !saturated`
    but read in the `else` of `if saturated`. Julia can't prove it's defined there, so an
    `UndefVarError` check (`llvm.intr.trap`) stays in the column kernel. Fixed: defined on all paths
    (bitwise identical). The `StepRange` of `(nlayers - 1):-1:1` also leaves an overflow-check function
    with a trap in the module, but it is never called.
  - Complex-constant stores (`div[I] = 0` in `_divergence_kernel!`, `data[I] = 0` in
    `zero_last_degree_kernel!`): `cannot raise op to stablehlo llvm.mlir.constant`. Reactant.jl#3450,
    fixed in EnzymeAD/Enzyme-JAX#3488 (merged 2026-10-12 19:29), in `Reactant_jll` 0.0.419 (built
    2026-10-13 07:02), but Reactant 0.2.292 and `main` pin `Reactant_jll` 0.0.418. Test environment:
    `claude_reactant_gpu_2026-10-13/env_reactant_main` (Reactant `main` d2285670 developed from
    `reactant_main/` with jll compat widened, `Reactant_jll` 0.0.419).
- Branch-free kernels on CUDA (no Reactant, A40, `branchfree_kernels_bench.jl`, results bitwise equal):
  - `_divergence_kernel!` (first element masked, neighbours via `previous_lm`/`next_lm`, last element
    via `ifelse`): equal or 1-2% faster once bandwidth-bound (>= ~20 μs, e.g. T255L32 50.4 vs 49.9 μs),
    up to ~3 μs slower in a few latency-bound cases (6-10 μs). Applied.
  - `zero_last_degree_kernel!` with `data[I] = ifelse(l == lmax, 0, data[I])` reads and writes all
    elements instead of writing only the last degree: equal up to ~5 μs, but 3.6-4.3x slower from
    T170L64, T255L32, T340L16 upwards (T255L64 16 vs 63 μs, T340L64 27 vs 115 μs). The model only
    calls it on single-layer tendencies (`η`, surface pressure in `horizontal_diffusion.jl`), where the
    kernel is launch-bound and equally fast. A penalty-free branch-free alternative for multi-layer
    arrays: launch over (m, k) only, with `lm = nonzeros(lmax, m)`.
- With Reactant `main` + `Reactant_jll` 0.0.419, `@compile raise = true time_step!` failed on both
  branch-free kernels: `ifelse` on a `Complex` compiles to a select of an LLVM aggregate
  (`llvm.select` on `!llvm.array<2 x f32>` with an aggregate `llvm.mlir.constant dense<0.0>`), which
  can't be raised (EJ#3488 covers scalar constants only). New `componentwise_ifelse`
  (SpeedyWeatherInternals.Utils, exported, tested) selects real and imaginary part separately and is
  used for the complex selects in `_divergence_kernel!` and `zero_last_degree_kernel!`.
  Possible upstream issue: raise aggregate selects/constants (or split them into scalars).
  Both kernels then raise (`branchfree_kernels_raise.jl`, Reactant CPU backend): `zero_last_degree!`
  bitwise equal to the CPU, `divergence!` max abs diff 4e-6 (values O(10), Float32 rounding).
- 2026-10-13 17:30: **`@compile raise = true time_step!` of the PrimitiveWetModel compiles on CUDA**
  (325 s; Reactant `main`, `Reactant_jll` 0.0.419, HIGHEST precision, Triton GEMM off). After
  `initialize!` all 112 arrays agree with the CPU (< 1%). Step 1 has NaNs in the spectral tendencies:
  `tendencies.divergence` 256 (= 32 x 8 layers), `tendencies.temperature` 256, `tendencies.pressure`
  29, none in `tendencies.grid.temperature`. Everything in grid space is NaN afterwards.
  Ruled out: the partial launch of `implicit_invert_S_kernel!` over `1:truncation`, raised on the
  GPU it keeps the untouched last row of `S⁻¹` (`invert_S_partial_raise.jl`).
  NaN pattern after step 1: divergence/temperature at 32 lm in all layers (m = 0 for l = 1-33 except
  l = 2 and 8, plus one at m = 32), pressure at m = 0, l = 3-32 except 8.
- Stage bisection (`reactant_wet_bisect_stages.jl`: the stages of `time_step!` each `@jit raise = true`
  on the GPU, compared with the CPU after each): all physics stages are NaN-free. The first NaNs
  appear in `dynamics_tendencies!`, in `dynamics.div_mean` (15 coefficients, m = 0,
  l = 3-5, 10-12, 17-20, 25-29), and from there in `tendencies.pressure`. `implicit_correction!`
  spreads them to all layers and to divergence and temperature at the same lm (62 lm then).
  `div_mean` comes from `_vertical_integration_spectral_kernel!` in `vertical_integration!`, which
  only reads the (NaN-free) prognostic divergence. Isolated (`vertical_integration_raise.jl`): the
  NaNs appear only when raised, with Triton on or off.
  Also seen: shortwave and cloud diagnostics differ by 1.7% (`cloud_top` up to 71%) after
  `parameterization_tendencies!`, to check separately.
- **Root cause (upstream raising bug):** a kernel that accumulates a `ComplexF32` in a loop and stores
  it gets the store packed into an i64 by LLVM (`zext` of the real-part bits, `shl` of the imaginary
  part, `or`). Raised, the zero-extension becomes `stablehlo.convert i32 -> i64` (sign-extending)
  followed by `add`. A negative real part decrements the bits of the imaginary part: ±0 (all m = 0
  modes) becomes NaN, otherwise the imaginary part is off by 1 ulp. MWE
  `mwe_complex_store_signext.jl` (Reactant CPU backend, so not GPU-specific): summing a complex
  matrix over its second dimension into a vector or a matrix column fails. Separate real
  accumulators (`complex(re, im)` stored) or accumulating in memory (`out[i] += ...`) work, as does a plain copy.
  - Fixed in `_vertical_integration_spectral_kernel!` and `implicit_primitive_leapfrog_kernel!`
    (real and imaginary part accumulated separately; all operators are real). The implicit kernel
    with m = 0 test data had NaN at all 264 m = 0 coefficients when raised. Now NaN-free, matches the CPU
    to 1e-7, and the new kernel is bitwise identical to the old one on the CPU (incl. signed zeros,
    `implicit_m0_raise.jl`). `test/dynamics/implicit.jl` passes.
  - To report upstream (Reactant.jl / Enzyme-JAX raise of `zext`). Issue draft:
    `reactant-issue-draft-complex-store-signext.md` with `reactant-issue-complex-store-signext-mwe.jl`
    (minimal kernel: weighted sum with a runtime loop bound, `-1 + NaN*im` instead of `-1 + 0im`).
- 2026-10-13 18:35: **PrimitiveWetModel with `@compile raise = true` runs on CUDA, 10 steps without
  NaNs** (`reactant_wet_compare_raise.jl`, Reactant `main` + `Reactant_jll` 0.0.419, HIGHEST precision,
  Triton GEMM off, compile 328 s). Max relative difference to the CPU after step 1 / step 10: vorticity
  3.4e-6 / 7.7e-5, divergence 9.3e-4 / 2.2e-4, temperature 2.4e-6 / 1.8e-6, humidity 2.6e-5 / 1.6e-5,
  pressure 9.6e-8 / 2.7e-8. Still different: shortwave and cloud diagnostics
  (`surface_shortwave_down` 1.7% from the first `parameterization_tendencies!`, `cloud_top` at some
  points). Not specific to raising or the GPU: the unraised Reactant run on the CPU backend
  (previous session, `compare_cpu.log`) shows identical numbers (1.66e-2, 0.714). To investigate
  separately (shortwave/cloud scheme on Reactant vs CPU).
- CUDA (no Reactant) after the split accumulators: GPU tests `vertical_integration.jl` (6/6) and
  `primitive_wet.jl` (2/2, 2/2) pass. Old vs new kernels are bitwise identical on the GPU
  (`split_accumulator_bench.jl`). Vertical integration is equally fast. The implicit kernel is equally fast up to
  L48, but 1.25-1.53x slower at L64 (T127L64 1.26 vs 1.92 ms, T255L64 4.7 vs 5.8 ms; the kernel
  spills to local memory there). The first split version used tuples of (re, im) pairs for
  `geopotential` and `G`. Variants (`implicit_variants.jl`, raised with m = 0 data in
  `implicit_variants_raise.jl`): keeping `geopotential`/`G` as complex `ntuple`s (variant a) still
  gives NaNs when raised (527, all m = 0), so these intermediates are packed too. Separate real
  `ntuple`s for real and imaginary parts (variant b) are NaN-free, bitwise identical. GPU timing
  (`implicit_variants_bench.jl`, μs, old / source (pairs) / variant b): T31L8 13.4 / 12.9 / 11.7,
  T63L8 13.6 / 13.1 / 11.9, T127L16 38.7 / 41.3 / 35.4, T255L32 940 / 926 / 891, T127L64 1277 / 1938 / 1851,
  T170L64 2184 / 2820 / 3319, T255L64 4638 / 5784 / 7072. Variant b is 5-13% faster than the old kernel up to L32.
  Every split version is 1.25-1.5x slower at L64, where the kernel spills its `ntuple`s to local
  memory. Open decision: accept, keep the old kernel for non-Reactant architectures, or restructure
  for large nlayers.

> Eventually we'll have to add a Reactant test to the GPU test suite. We can't afford a full correctness
> test, but we must add a test that checks that the PrimitiveWetModel steps without NaNs. Use `raise=true` there

- `test/GPU/reactant.jl` rewritten and included again in `test/GPU/runtests.jl` (CUDA only):
  PrimitiveWetModel (default resolution, `MatrixSpectralTransform`), `@compile raise = true
  time_step!`, 10 steps through `time_stepping!`, `@test` no NaNs in vorticity, divergence,
  temperature, humidity and pressure. Barotropic and unraised runs removed (unraised is wrong on CUDA
  until #3377). Needs a Reactant release with `Reactant_jll` >= 0.0.419; with 0.2.292 the raise fails
  on complex-constant stores.
- First run of that test: NaNs in all prognostic variables, although `reactant_wet_compare_raise.jl`
  runs the same 10 steps cleanly. The difference: `initialize!(simulation)` transforms the initial
  conditions with `@maybe_jit` (`Reactant.@jit` without raising), so on CUDA the time stepping starts
  from a grid state corrupted by the layout bug. The comparison script redoes that transform with
  `raise = true`.

> The GPU test isn't how we would typically run our models. What if we make `raise=true` the default
> for `@maybe_jit` can we then just run a standard model script?

> A new Reactant version will be released in the next days, don't worry about this

- `raise = true` by default: `_jit` (`@maybe_jit`, SpeedyWeatherInternalsReactantExt) and the
  `@compile` of `time_step!` in `time_stepping!` (SpeedyWeatherReactantExt, also its docstring example).
  `test/GPU/reactant.jl` is now a standard model script: `initialize!(model)`, `run!(simulation, steps = 10)`,
  NaN checks. Not changed: `test/reactant/test_correctness.jl` and `benchmark/reactant_runner.jl`
  compile `time_step!` themselves without `raise`.
- With `@maybe_jit` raising, `initialize!(model)` failed: `_interpolate_kernel!` (RingGrids, via
  `@maybe_jit interpolate!` in `initialize!(::EarthOrography)`) couldn't be raised (`llvm.ptrtoint`).
  Traced by `@jit`, the pole averages are 0-dimensional arrays (numbers on CPU/GPU), so the ternary
  `(A_northpole, A_northpole) : (A[...], A[...])` mixed types. Fixed with `A_northpole[]`/`A_southpole[]`
  (works for numbers too). After that, `initialize!(model)` (167 s) and `initialize!(simulation)` (100 s)
  run raised on the GPU without NaNs (`initialize_raise.jl`).

## Handoff 2026-10-13

State at the end of the session (branch `mg/new-reactant`, last commit `e7f6942f`).

Uncommitted in the working tree:
- `raise = true` default in `_jit` (SpeedyWeatherInternalsReactantExt) and `time_stepping!`
  (SpeedyWeatherReactantExt, incl. docstring).
- RingGrids `_interpolate_kernel!`: `A_northpole[]`/`A_southpole[]`.
- `test/GPU/reactant.jl` as a standard model script (`initialize!`, `run!(steps = 10)`, NaN checks), included
  again in `test/GPU/runtests.jl` for CUDA.
- The user's edit of the `_divergence_kernel!` comment (says `lm==lmax`, the code uses `lmmax`), this doc.

Running when the session ended: `test/GPU/reactant.jl` with `env_reactant_main` on login02, detached;
result in `claude_reactant_gpu_2026-10-13/gpu_reactant_test_standard.log` (look for "Test Summary").
`initialize!` passes raised, the time step with the new RingGrids fix is untested in this script.

Open items:
1. Check `gpu_reactant_test_standard.log`. If `run!` fails, compare with `reactant_wet_compare_raise.jl`,
   which passes when called with the same settings.
2. Implicit kernel at 64 layers on CUDA: the split accumulators (committed in `e7f6942f`, pairs in
   `ntuple`s) are 1.25-1.5x slower than the old kernel at L64, equally fast up to L32. Variant b
   (`implicit_variants.jl`, separate real `ntuple`s) is 5-13% faster up to L32, but also 1.45-1.52x
   slower at L64. To decide: accept, keep the old kernel for non-Reactant architectures, or restructure.
3. CPU Reactant tests (`test/reactant`) with `raise = true` in `@maybe_jit` and `time_stepping!`: not run
   yet. They need `Reactant_jll` >= 0.0.419 (a Reactant release is expected in the coming days, then raise
   compat and update the CHANGELOG line of #1304, which says v0.2.292). `test_correctness.jl` and
   `benchmark/reactant_runner.jl` compile `time_step!` without `raise`; decide whether they should raise.
4. Report upstream: `reactant-issue-draft-complex-store-signext.md` (complex store sign extension), plus
   an issue for `ifelse` on `Complex` (aggregate select, can't be raised, worked around with
   `componentwise_ifelse`) and for raised `dot_general`s without precision config (TF32 on GPU,
   worked around with `xla_gpu_enable_triton_gemm = false`).
5. Shortwave/cloud diagnostics differ from the CPU by 1.7% on Reactant in general (not only raised or on
   the GPU), to investigate.

Environment: `claude_reactant_gpu_2026-10-13/env_reactant_main` (Reactant `main` d2285670 developed
from `reactant_main/` with jll compat widened to 0.0.419, `Reactant_jll` 0.0.419, CUDA). Scripts for
all the steps above are in `claude_reactant_gpu_2026-10-13/`.

> continue, i accept the currently fixed implicit version for now. its temporary anyway

Continued 2026-10-13 evening:
- Item 1 done: `test/GPU/reactant.jl` (standard script, `raise = true` defaults) passes on the A40,
  5/5 prognostic variables NaN-free after 10 steps (12m56s incl. compilation).
- Item 2 closed: the committed implicit kernel (split accumulators, pairs in `ntuple`s) stays, it's a
  temporary workaround for the upstream bug.
- Item 3 done: `test/reactant` on the Reactant CPU backend with the `raise = true` defaults passes,
  environment `env_reactant_tests` (test/reactant/Project.toml + Reactant `main` + `Reactant_jll`
  0.0.419 + Statistics), runner `run_reactant_tests.jl` (same as `runtests.jl`, CPU backend set
  explicitly as login02 has a GPU), log `reactant_tests_cpu.log`. Results: `@maybe_jit` nesting 8/8,
  Clock copy 14/14, views 4/4, MatrixSpectralTransform 2/2, Barotropic 8/8, ShallowWater 14/14,
  PrimitiveDry 17/17, PrimitiveWet 20/20, output writers 29/29, `make_zero` 4/4. Differentiability
  (BarotropicModel): 2 pass and 2 **unexpected passes**. Loss and gradient of the raised Reactant time
  step now agree with Enzyme on the CPU (rtol 1e-4 / 1e-3): the complex-constant fix unblocked it.
  TODO: change the two `@test_broken` in `test/reactant/differentation.jl` to `@test` (needs the
  Reactant release with `Reactant_jll` >= 0.0.419).
- Item 4: issue drafts with MWEs for the other two upstream problems:
  `reactant-issue-draft-complex-ifelse.md` (`ifelse` on `Complex` → aggregate `llvm.select`, can't be
  raised) and `reactant-issue-draft-raised-dot-precision.md` (raised `dot_general`s have no
  `precision_config`). The precision MWE (A40, n = 16, relative errors) corrects the earlier conclusion
  that turning off Triton GEMMs is enough. Triton GEMM: 3.1e-4, ignoring `HIGHEST` and
  `NVIDIA_TF32_OVERRIDE`. Triton off (cuBLAS): 2.7e-4, TF32 at default precision too. Only Triton off
  *and* `NVIDIA_TF32_OVERRIDE=0`: 8.3e-8. The implicit kernel was fine with Triton off alone,
  presumably because XLA didn't route its small dot to cuBLAS. In the MWE the loops are only raised
  with `inbounds = true` and static trip counts (`Val(n)`), and with n = 8 LLVM unrolls them so no
  `dot_general` appears.
- Item 5 in progress (shortwave/cloud differences), `cloud_top_compare.jl` (GPU, raised, HIGHEST,
  Triton off): after the first `parameterization_tendencies!`, `cloud_top` differs at 208 of 3168
  points, sometimes by several levels (e.g. 7 vs 2). Away from those points the shortwave agrees to
  ~1e-4 (`surface_shortwave_down` 1.2e-4, `outgoing_shortwave` 2.3e-4; max 1.7e-2 and 1.9e-1 overall).
  `cloud_cover` max 7e-4 relative everywhere, rain rates identical. `cloud_top` is the argmax of the
  relative humidity over layers 2 to nlayers-1, lowered by deep convection, so near ties can flip under
  rounding. Hypothesis: the initial humidity is close to a constant relative humidity with height.
  Not yet confirmed: the script still has to compare the relative humidity at the two chosen levels.
  It accesses the clouds component as `model.shortwave_radiation.clouds`, but PrimitiveWetModel has no
  field `shortwave_radiation`; find the right field name (`fieldnames(PrimitiveWetModel)`), then
  rerun. If ties are confirmed: harmless, maybe compare cloud diagnostics with a tolerance on the
  level or after the spin-up only.
  **Confirmed (2026-10-14), harmless:** the initial humidity has a relative humidity of ~0.70 at all
  levels (e.g. 0.700056, 0.700055, 0.699830, ... for levels 2, 3, 4, ...). At the 208 points the
  relative humidity of the levels chosen by CPU and Reactant differs by a median of 2.3e-6 (max
  1.1e-4), so the argmax is a near tie that ulp-level differences (e.g. in `exp` of the saturation
  humidity) flip. The shortwave differences follow from the different cloud top only. Nothing to fix;
  the CPU-vs-Reactant correctness tests compare after a 10-day spin-up and pass.
- 2026-10-14: `test/reactant/differentation.jl`: the two `@test_broken` are now `@test`, the
  `try_or_nothing` wrapper is removed and the header comment updated. Passes 4/4 on the Reactant CPU
  backend with `env_reactant_tests` (`run_reactant_tests_diff.jl`: only the differentiation tests, with
  `sync_variables!` from `test_correctness.jl` extracted to `sync_variables.jl`).

> (Should `test_correctness.jl` and `benchmark/reactant_runner.jl` raise too?) Both raise.
> (Next?) GPU benchmarks with Reactant, but on the H100 via submit_bench.sh

> other todos: differentiability tested on the PrimitiveWetModel instead of the Barotropic model in
> the unit tests; Reactant GPU tests also included in Metal tests (but not AMD)

- `@compile raise = true` in `test/reactant/test_correctness.jl` (3x) and
  `benchmark/reactant_runner.jl` (`build_runner` and `compiled_call`).
- Reactant GPU benchmarks on the H100: Slurm job 2483557 (`claude_reactant_gpu_2026-10-13/submit_bench_reactant.sh`,
  log `benchmarks-reactant-<jobid>.log` there). The driver activates its own directory's environment, so the job
  runs a copy of `SpeedyWeather/benchmark` (`benchmark_reactant/`) whose environment uses Reactant
  `main` d2285670 + `Reactant_jll` 0.0.419 (paths absolute, precompiled on login02). When it is done, copy
  `benchmark_reactant/README.md` and `benchmark_reactant/assets/benchmark_results.json` back into
  `SpeedyWeather/benchmark/` if the repository's copies haven't changed in the meantime (they include all
  other architectures plus the new `reactant-gpu` section).
- Job 2483557 done (H100, results in `benchmark_reactant/README.md`, section `reactant-gpu`). Failed: all
  PrimitiveWet configurations with >= 16 layers (T86/T128 with L16/L24, T32L16 in the layers suite; L12
  works), the vertical diffusion `while` loop, fixed since. Reactant (raise) vs CUDA on the same H100 (SYPD;
  CUDA results from 2026-10-06, `gpu-nvidia`): Barotropic T32 28405 vs 10929; ShallowWater T32 24420
  vs 11638, T64 8045 vs 3997, T86 2256 vs 701, T128 809 vs 292; PrimitiveWet L8 T32 1762 vs 5291, T43 2255
  vs 4342, T64 1361 vs 2920, T86 425 vs 652 (CUDA matrix transform 1075), T128 268 vs 262 (matrix 196);
  Float64 T32 2756 vs 4612. Reactant compile times 2-270 s per configuration. Memory isn't comparable
  (Reactant counts its device buffers, CUDA's `summarysize` doesn't). Rerun with the vertical diffusion fix:
  job 2485035.
- Todo: differentiability test with the PrimitiveWetModel instead of the BarotropicModel.
- Todo: include the Reactant GPU test in the Metal GPU tests (not AMDGPU). Not possible yet: Reactant
  has no Metal backend (jax-metal plugin disabled in Reactant.jl#1967, the MPSGraph backend #2489 is
  an open PR untouched since July), `set_default_backend("gpu")` finds no GPU client on a Mac. Decision:
  skip, comment in the Metal branch of `test/GPU/runtests.jl` to include `reactant.jl` once #2489 is in.
- Benchmark job 2483557 (H100, reactant-gpu): the PrimitiveWet resolution benchmarks with 16 layers failed
  to raise `column_parameterizations_kernel!`. `get_diffusion_coefficients!` (vertical diffusion)
  finds the boundary layer depth with `while kₕ > 0 && Ri[ij, kₕ] < Ri_c` (data-dependent early exit,
  `scf.while`). At 8 layers it raised (presumably unrolled), at 16 not.

> yes, try to make the vertical diffusion loop branch free. Check for identical results and performance
> on both CPU and GPU

- Fixed-trip loop `for k in nlayers:-1:1` with a flag `below_critical &= Ri[ij, k] < Ri_c` and
  `kₕ = ifelse(below_critical, k - 1, kₕ)`: same kₕ (lowest layer with Ri ≥ Ri_c scanning up from the
  surface, 0 if none, NaN stops both). A/B check in one process (`vertical_diffusion_ab.jl`:
  `parameterization_tendencies!` from the same state after a 2-day spin-up, new code vs the HEAD version
  redefined into SpeedyWeather, bitwise comparison of the grid tendencies and boundary layer heights,
  median time) on CPU and A40, T31L8, T63L16, T127L16, T127L32.
  A40 (`vertical_diffusion_aba.jl`, timed new → old → new as the first try showed an order effect):
  bitwise identical at all sizes. `parameterization_tendencies!` in ms (new / old / new): T31L8 0.215 / 0.207 / 0.213,
  T63L16 0.306 / 0.291 / 0.290, T127L16 0.422 / 0.419 / 0.424, T127L32 0.597 / 0.581 / 0.595, i.e.
  ratio (best new / old) 1.00-1.03, within noise. CPU (1 thread, `vertical_diffusion_ab.jl`): bitwise
  identical, ms (old / new): T31L8 4.67 / 4.65, T63L16 26.9 / 27.3, T127L16 103.8 / 103.0, T127L32 143.3 / 142.6,
  ratio 0.99-1.02, within noise. **Fix kept.** (Observation, the same in old and new:
  at T127L32 all columns have boundary layer height 33 = no boundary layer.)
- PrimitiveWet differentiability test (`reactant_tests_diff_wet.log`): stopped after 1.5 h without
  result (the Barotropic one took 3 min).

> about the differentiation tests, regular Enzyme won't run in Julia 1.12 and Julia 1.13 at the moment.
> We could only do the direct comparison in Julia 1.10 or maybe Julia 1.11. This test will also
> generally take too long to run in the unit tests. What we want: in the regular unit tests run the
> gradient comparison on the Barotropic model; provide an extra script based on the same code (adjust
> it so that it just calls functions defined in the script the CI runs) that makes it possible to test
> it on the PrimitiveWetModel as well.

- `test/reactant/differentation.jl` restructured: `differentiability_simulations(ModelType)` (CPU and
  Reactant simulation in the same state) and `test_differentiability(ModelType; rtol_loss, rtol_gradient)`
  (the testset). CI runs `make_zero` and `test_differentiability(BarotropicModel)`, skipped if the including script
  sets `DIFFERENTIATION_FUNCTIONS_ONLY = true`. New `test/reactant/differentiation_primitive_wet.jl`:
  activates the test environment, sets that flag, includes `setup.jl` and `differentation.jl`, calls
  `test_differentiability(PrimitiveWetModel)`; `REACTANT_BACKEND` selects the backend. `sync_variables!`
  moved from `test_correctness.jl` to `setup.jl` (also used by `test_output.jl` and the new script).
  Note: in Julia 1.12.2 the Barotropic Enzyme CPU reference did work (4/4); the PrimitiveWet run
  didn't finish in 1.5 h. CI path on Julia 1.12 (`run_reactant_tests_diff.jl`): `make_zero` 4/4,
  `test_differentiability(BarotropicModel)` 4/4. PrimitiveWet on Julia 1.11.7 (`env_reactant_tests_111`,
  Reactant `main` + `Reactant_jll` 0.0.419, wrapper `run_differentiation_primitive_wet.jl` = the script
  without `Pkg.activate`, log `differentiation_primitive_wet_111.log`): the first try failed compiling
  `implicit_primitive_leapfrog_kernel!`, "unsupported dynamic function invocation (call to zero)". The
  split-accumulator version (`e7f6942f`) uses `zero(NF)` inside the `ntuple` closure with `NF` captured,
  and Julia 1.11 doesn't infer the captured type (1.12 does). This would also break CUDA runs on Julia
  1.10/1.11. Fixed: `zero(real(eltype(temp_tend)))` inside the closure (committed by the user in `d85e908e`).
  Rerun: 3 of 4 pass (raised loss = CPU loss, Enzyme CPU gradient finite and non-zero, 16 min). The
  Reactant gradient fails to compile: "Unsupported operation in reduction rev autodiff", a
  `stablehlo.reduce` with select/or (a masked reduction) from the raised branches of
  `dpdy_kernel!` (SpeedyTransforms, meridional gradient in `∇!`), as earlier for the UV kernels.
- `dpdy_kernel!` branch-free (as `_divergence_kernel!`: neighbours via `previous_lm`/`next_lm`,
  missing terms zeroed with `componentwise_ifelse`). It was the last kernel with the `if lm == 1 ... elseif
  lm == lmmax` pattern. `dpdy_ab.jl` (old vs new, T31L8 to T255L64): results `==` everywhere (bitwise except
  signed zeros at the first/last element). A40: same speed (1.00-1.02). **CPU (1 thread): 1.9-2.4x
  slower** (e.g. T127L16 235 vs 512 μs), as the `rem`-based wrap is an integer division per element that
  dominates this cheap kernel. `∇!` runs once per time step on the 1-layer surface pressure, so about
  +17 μs per step at T127 on the CPU, negligible against the time step. Kept.
- Rerun on Julia 1.11: again 3 of 4, now "WhileOp does not have known iteration count for cache removal"
  in `longwave_radiation.jl` (OneBandLongwaveRadiativeTransfer, upward beam).

> (How to proceed with the PrimitiveWet differentiation?) Keep going

- MWEs (`mwe_reverse_steprange.jl`, `mwe_reverse_view_bounds.jl`, `mwe_reverse_dynamic_bounds.jl`,
  `mwe_reverse_beam.jl`): reversed `StepRange`s, views, runtime bounds and a beam updating rows k and k-1 of
  a tendency array all differentiate fine, so the loop form isn't the problem. In the dumped module
  the failing `stablehlo.while` carries an i1 exit flag and an i32 exit code computed from unsigned
  `index >= length` comparisons: **bounds checks**. The early exit makes the trip count unknown.
  Source: `pressure_thickness` (via `flux_to_tendency`), `pressure`, `pressure_half`, `sigma`,
  `pressure_below/above` in `vertical_coordinates.jl` were `@inline` but not `@propagate_inbounds`, so the
  kernel's `inbounds = true` didn't reach them. Fixed: all ten are `@propagate_inbounds` (removes checks
  only in an `@inbounds` caller, e.g. kernels).
- Faster iteration: `reverse_raise_check.jl` compiles only the Reactant reverse pass on Julia 1.12 (no
  Enzyme CPU reference), `MODE=physics` (loss of `parameterization_tendencies!`) or `MODE=step` (full
  time step as `test_differentiability`). With the `@propagate_inbounds` fix: `MODE=physics` compiles and
  runs (gradient finite, non-zero, 3.5 min), `MODE=step` too (3.7 min).
- **PrimitiveWet differentiation passes 4/4 on Julia 1.11** (`run_differentiation_primitive_wet.jl`,
  test configuration: constant longwave transmissivity, no convection): raised Reactant loss = CPU loss
  (rtol 1e-4), Enzyme CPU gradient finite and non-zero, Reactant gradient = Enzyme CPU gradient
  (rtol 1e-3), 16.6 min. Default PrimitiveWetModel (Frierson transmissivity, Betts-Miller convection):
  reverse check in progress (`reverse_raise_check.jl`, `MODEL=default`).

> Add a section to the differentiability section of `differentiability.md` about Reactant. Note that it
> is still work in progress, which Reactant options have to chosen and that this does also work with
> Julia 1.12 and Julia 1.13 in contrast to regular Enzyme without Reactant ("vanilla Enzyme")

- `docs/src/differentiability.md`: new section "Differentiating with Reactant" after the checkpointing
  section (work-in-progress warning, Reactant version, example with a `ReactantDevice` PrimitiveWetModel,
  loss = squared vorticity after one `time_step!`, `@jit raise = true raise_first = true autodiff(...)`,
  what the options do, note on TF32 matrix products on GPUs, the CI and PrimitiveWet scripts), and a
  pointer to it in the Julia version warning at the top. Julia 1.13 is stated as the user said, only
  1.12 (and 1.11) were run here.

> You reported that PrimitiveWet with reactant trails non-reactant at low resolution due to column
> physics. I want to investigate this further. What's the status on this?

- Not investigated before, "column physics" was inferred. The physics-only/dynamics-only suite of job
  2483557 confirms it. PrimitiveWet T32L8 on the H100, time per step from SYPD (Δt = 2400 s), Reactant
  vs CUDA: full 3.7 vs 1.2 ms, dynamics only 0.50 vs 0.90 ms (Reactant 1.8x faster), physics only 3.3 vs
  0.71 ms (Reactant 4.6x slower). Lead: in raised modules the layer loops of the column kernel become
  `stablehlo.while` loops whose iterations gather/scatter into the whole fused tendency array; on a
  GPU, each iteration means more kernel launches and possibly a host sync for the condition.
  A40, T31L8, `feedback = nothing`, 200 calls of the compiled `time_step!` (`physics_perf.jl`, raised modules
  in `hlo/physics_perf/`) and the same with CUDA.jl (`cuda_physics_perf.jl`), ms per step:

  | | Reactant raise = true | Reactant raise = false | CUDA (FFT) | CUDA (matrix) |
  |---|---|---|---|---|
  | physics only | 1.43 | 2.07 | 0.93 | 1.02 |
  | dynamics only | 0.64 | 1.59 | 1.09 | 1.06 |
  | full | 1.76 | 2.08 (NaN, layout bug) | 1.21 | 1.36 |

  So raising is not the cause: unraised is slower for both. On the A40 with this harness, Reactant's
  physics is 1.5x slower than CUDA (the H100 benchmark suggested 4.6x), dynamics 1.7x faster, and the full
  step 1.45x slower. The H100 benchmark had 3.3 ms for Reactant's physics alone, more than twice the A40
  value here. Suspect: the benchmark harness, i.e. `time_stepping!` calls `progress!` and `output!` on the
  host after every compiled step, possibly with a device sync (feedback). Next: Reactant with the
  default feedback vs `feedback = nothing` through `time_stepping!` (`harness_perf.jl`).
- `harness_perf.jl` (A40, PrimitiveWet T31L8, default precision): compiled `time_step!` called directly
  1.90 ms per step, through `time_stepping!` (as `run!`, `feedback = nothing`) **7.58 ms per step**.
  With the default feedback: **error "Scalar indexing is disallowed"**, so `run!` with the default
  progress feedback fails on Reactant GPU. Causes:
  - `do_output!` (Reactant method, SpeedyWeatherReactantExt) read `Int(clock.time_step_counter)` from the
    device before checking `output.active`: a device-to-host sync after every compiled step, also with
    output off, so the host can't queue the next step while the GPU computes. Fixed: `output.active ||
    return false` first. (With output on it still syncs every step; a host-side counter would avoid that.)
  - `progress!` → `max_speed`/`temperature_range` → `extrema(::Field)` iterates elementwise over a view of
    a concrete Reactant array: allowed (slowly) on Reactant's CPU backend, not on the GPU. First fix:
    `extrema` on the host in RingGridsReactantExt, replaced as follows.

> Just disable the extrema read on Reactant by default, only make it opt-in

  - `Feedback`: `show_umax` and `show_temperature_range` are `Union{Bool, Nothing} = nothing`, `nothing` =
    the architecture's default via `show_diagnostic(option, architecture)` in `feedback.jl`: shown for
    CPU/GPU, not for `ReactantDevice`. An explicit `true` opts in (on Reactant GPU that then hits the scalar
    indexing error, there is no host fallback). The `extrema` method in RingGridsReactantExt is removed.
    Test in `test/output/feedback.jl` ("Feedback diagnostics opt-in for Reactant", 8/8).
  Rerun with the `do_output!` fix (`harness_perf_fixed.log`): direct 1.90, `time_stepping!` without
  feedback **still 8.01 ms**, with feedback 13.2 ms per step. So `do_output!` isn't the cause (or not the
  only one). In progress: `harness_parts_perf.jl` times the parts of the loop separately (step + `output!`,
  step + `progress!`, `output!` alone, `time_stepping!`) and prints the dispatched methods.

> This note on precision of matrix products, should rather go into the general `architecture` docs.
> But also please explain to me and in the docs the consequences of that. Does the reactant model has
> to be run with those settings or not?

- Consequences measured (`reactant_wet_compare_raise.jl` with PRECISION=default, no Triton flag, log
  `wet_raise_gpu_tf32.log`, A40, PrimitiveWet T31L8, max rel. difference to the CPU, step 1 → step 10;
  in brackets HIGHEST + Triton off): vorticity 1.7e-4 → 1.4e-2 (3.4e-6 → 7.7e-5), divergence 4.8e-1 →
  3.9e-2 (9.3e-4 → 2.2e-4), temperature 2.1e-4 → 5.0e-4 (2.4e-6 → 1.8e-6), humidity 1.3e-2 → 8.4e-3
  (2.6e-5 → 1.6e-5), surface pressure 4.5e-6 → 4.7e-6 (9.6e-8 → 2.7e-8). Stable (no NaN), but effectively
  a lower precision model, needed: Float32 precision for comparisons, tests, gradients.
- Docs: the note moved from `differentiability.md` (now a short pointer) to `architectures_gpu.md`,
  Reactant Support, new subsection "Precision of matrix products on GPUs" (what TF32 is, which matrix
  products, the table above, when to use Float32, the three settings, example wrapping `initialize!` and
  `run!` in `with_config`, environment variables before starting Julia). Also updated there: Reactant
  version (newer than v0.2.292, `Reactant_jll` >= 0.0.419), `raise = true` default, GPU NaN test instead
  of "not yet tested".
- Speed cost (`precision_perf.jl`, A40, PrimitiveWet L8, direct compiled steps, ms per step): default (TF32)
  T31 1.81, T63 3.08; + HIGHEST T31 2.00 (+10%), T63 3.80 (+23%); + Triton off + `NVIDIA_TF32_OVERRIDE=0`
  T31 4.18 (+130%), T63 4.87 (+58%). The accuracy table (Float32 column) used HIGHEST + Triton off without
  the override, so the override may not be needed for the model. HIGHEST + Triton off without the override:
  T31 2.35 ms (+30%), so the override is the expensive part. In progress: T63 of that, accuracy of HIGHEST
  alone.
- Benchmark rerun job 2485035 (H100) **cancelled after 3.6 h**: the log stopped at 13:09 after compiling
  the first 16-layer configuration (T86L16, 174 s, ptxas: registers spilled, e.g. 27324 bytes spill stores in
  one fusion), then 3 h at 100% GPU without finishing it. So at 16 layers the raised time step compiles
  (with the branch-free vertical diffusion) but runs extremely slowly. The first run's L8/L12 results stand.
  Next: Reactant at 12 vs 16 layers on the A40 (`physics_perf_layers.jl`, physics only / dynamics only /
  full, 20 steps, raised modules in `hlo/physics_perf_L*`).
- `harness_parts_perf.jl` (A40, T31L8; may have overlapped with another GPU run): compiled step 2.15,
  + `output!` 1.92, + `progress!` 2.17, `output!` alone 0.001, but `time_stepping!` 7.52 ms per step. So neither
  call is the cause, it's how `time_stepping!` runs the loop. Next: `harness_variants_perf.jl`.
- GPU runs overlapped because the chained waits started at the first output line of the previous run.
  Now `gpu_run.sh` (flock on `gpu.lock`) and one sequential chain.

> also try to diagnose and improve the performance of the parameterizations / column physics with
> reactant. are they still slower with Reactant?

- Yes: A40 T31L8 physics only 1.43 ms (Reactant, raised) vs 0.93 ms (CUDA). The raised physics-only module
  (`hlo/physics_perf/physics only_raise_true`) has only 5 `stablehlo.while` loops, **all from the vertical
  diffusion**: the K loop and the four `_vertical_diffusion!` calls (u, v, humidity, temperature), all
  `for k in kₕ:nlayers` with the column-dependent boundary layer top `kₕ`, so they can't be unrolled. Each
  carries the whole fused tendency array (`1x73x3168`) through the iterations with gather/scatter. The
  other physics loops have static trip counts and become parallel array operations.
- Rewrite (bitwise identical): both loops over `1:nlayers` with the mask `k >= kₕ`, indices that need
  `kₕ` clamped to `min(kₕ, nlayers)`, Rifac of `Ri[kₕ]` computed once before the loop, the diffusion stores the
  unchanged tendency outside the boundary layer. CPU A/B (1 thread): bitwise identical, but the whole
  `parameterization_tendencies!` 13-29% slower (T31L8 5.36 vs 4.61 ms, T127L32 186 vs 144 ms), as the masked
  loops now cover all layers instead of the few boundary layers. Hence: the loop range depends on the
  architecture, `boundary_layers(kₕ, nlayers, ::AbstractArchitecture) = kₕ:nlayers` and
  `boundary_layers(kₕ, nlayers, ::ReactantDevice) = 1:nlayers`, the mask stays (a no-op for `kₕ:nlayers`),
  `model.architecture` is passed through `vertical_diffusion!`, `get_diffusion_coefficients!` and
  `_vertical_diffusion!` (only internal callers), the K reset loop over all layers is back.
  CPU A/B of the whole call chain (`vertical_diffusion_aba.jl`, old = HEAD from `# function barrier` to the
  end of `_vertical_diffusion!`): bitwise identical, `parameterization_tendencies!` new/old 0.99-1.00
  (T31L8 4.67 vs 4.68 ms, T127L32 144.9 vs 146.4 ms). GPU A/B and the Reactant speed-up in progress.
- `harness_variants_perf.jl` (A40, full PrimitiveWet T31L8, raise, **with the masked vertical diffusion**):
  loop of compiled steps 1.08 ms (1.20 after `initialize!`), loop + `progress!` + `output!` 3.39 ms,
  `time_stepping!` 6.92 ms the first time, **0.92 ms the second time**. So the `time_stepping!` overhead
  was Julia compiling the host code on the first call inside the timed loop, an artifact (the benchmark
  runner warms up first). And the full step went from 1.76 ms (before the vertical diffusion rewrite) to
  ~1.1-1.2 ms, on par with or faster than CUDA (1.21 ms).
- Nsight: XLA runs the whole time step as one CUDA graph (one `cuGraphLaunch` per step, plus 2
  `cuMemcpyDtoDAsync` and 5 `cuLaunchHostFunc`), kernels inside graphs need `--cuda-graph-trace=node`. The
  CUDA profile failed (`@compile` parsed in the non-Reactant branch, now `@eval`). Rerun queued.
- Correction: `physics_perf_layers.jl` L8 with the new code (A40, 50 steps): physics only 1.46 ms (before
  the vertical diffusion rewrite 1.43), dynamics only 0.75, full 1.24 ms. So the rewrite did not make the
  physics faster on its own. The full step's 1.76 → 1.24 ms has another cause, or the runs are noisy;
  to be checked with the profiles.

> we want to submit the benchmark jobs on the H100 next

- Job 2490777 (combined CUDA then Reactant) cancelled right away.

> don't do a combined job, do two separate jobs so that they'll finish more quickly

- Two jobs in parallel, each in its own directory so they don't overwrite each other's
  `benchmark_results.json`: **2490788** CUDA (`gpu`) with the user's `submit_bench.sh` in the repository's
  `SpeedyWeather/benchmark` (log `benchmarks-2490788.log` in the repository root), **2490789** Reactant
  (`reactant-gpu`) with `submit_bench_reactant.sh` in `benchmark_reactant/` (log there). Both use the current
  working tree. Afterwards merge the `reactant-gpu` entry of `benchmark_reactant/assets/benchmark_results.json`
  into the repository's and regenerate the README. Watch the 16/24-layer Reactant configurations.
- Precision: `HIGHEST` alone (no Triton flag, no override, `wet_raise_gpu_highest_only.log`) gives the
  same agreement as HIGHEST + Triton off after 10 steps: vorticity 8.0e-5, divergence 2.6e-4, temperature
  1.8e-6, humidity 1.9e-5, surface pressure 2.7e-8. So for the model `HIGHEST` alone suffices (+10% T31,
  +23% T63); the raised-kernel matrix products (implicit solver) don't matter at this level. Docs
  updated: recommend `with_config(dot_general_precision = HIGHEST)` around `initialize!` and `run!`,
  table column `HIGHEST`, the environment variables only as optional (1.6-2.3x slower, no noticeable
  effect on the model).
- Found on the way: `_vertical_diffusion!` uses `k₋ = max(k, 1)` and `k₊ = min(k, nlayers)`, i.e.
  `k₋ = k₊ = k`, so the vertical diffusion adds exactly zero (since 80926431, 2024-03-21, "Vertical
  diffusion operator streamlined"). **The user confirms, it's fixed independently in another PR**, so not
  touched here. Note for that PR: with the masked loop (k >= kₕ) the flux out of the top of layer kₕ into kₕ - 1
  is missing; a loop over all layers with K = 0 above the boundary layer would include it.
- Future: with Reactant.jl#3466 (`disable_dot_formation_passes`) no dots from raised loops, then only
  setting 1 would be needed, and SpeedyWeather could set both in its compile calls by default.

Next steps (state at 2026-10-13 23:00, last commit `e7f6942f`):
1. Commit the working tree: `raise = true` defaults (both Reactant extensions), RingGrids
   `_interpolate_kernel!` fix, `test/GPU/reactant.jl` + `runtests.jl`, this doc. The untracked issue
   drafts and MWEs in `docs/dev/2026-10/` are new too (`reactant-issue-*`).
2. Finish item 5 (cloud_top near-tie check, see above).
3. When the Reactant release with `Reactant_jll` >= 0.0.419 is out: compat bound for Reactant (SpeedyWeather,
   SpeedyWeatherInternals, RingGrids, LowerTriangularArrays extensions and the test environments),
   update the CHANGELOG line of #1304 ("requires Reactant v0.2.292"), `@test_broken` → `@test` in
   `test/reactant/differentation.jl`, run `test/reactant` and `test/GPU` in the repo's own environments.
4. Post the upstream issues (`reactant-issue-draft-complex-store-signext.md`,
   `-complex-ifelse.md`, `-raised-dot-precision.md`; the layout one has PR #3377).
5. Decide whether `test/reactant/test_correctness.jl` and `benchmark/reactant_runner.jl` should also
   compile `time_step!` with `raise = true` (they currently don't).
6. Precision of raised dots on GPU: document or default the workaround (Triton GEMM off *and*
   `NVIDIA_TF32_OVERRIDE=0`), or wait for an upstream fix.

2026-10-14 morning: Reactant `main` is at 0.2.293 with `Reactant_jll` 0.0.420 (not registered yet,
General ends at 0.2.292 / 0.0.419). Merged overnight: EnzymeAD/Reactant.jl#3466, a compile option
`disable_dot_formation_passes` (default `false`) that switches off `reduce_mul_to_dot_general` and
`split_reduce_add_mul_to_add_dot_general`, the passes that turn our raised reduction loops into
`dot_general`s. It needs a `Reactant_jll` built with Enzyme-JAX#3527 (Reactant.jl#3467). Once released,
`disable_dot_formation_passes = true` in our raised compilations should be the fix for item 6: the
raised contractions stay fused reductions in Float32 (also faster for small contractions on GPU according to
the PR), and traced dots (MatrixSpectralTransform) keep their `precision_config`.

2026-10-14 afternoon:

- Nsight profiles (A40, physics-only step T31L8, 50 steps, `run_profiles.sh`; the Reactant kernel CSV had
  to be regenerated by hand, `nsys/physics_reactant2_cuda_gpu_kern_sum.csv`): **Reactant 81 kernels and
  0.22 ms GPU time per step, CUDA.jl 232 kernels and 0.66 ms per step**. Wall time per step under the
  profiler: Reactant 0.94 ms, CUDA 1.52 ms (more launches → more profiler overhead). The two biggest
  Reactant kernels are `gemm_fusion_dot_2/3` (32-36 µs each, the MatrixSpectralTransform), no single
  physics kernel above 8 µs. The biggest CUDA.jl kernel is `gpu_column_parameterizations_kernel_` with
  60 µs. So the Reactant GPU work is ~3x smaller; its wall time is host-side overhead per call (one
  `cuGraphLaunch`, 5 `cuLaunchHostFunc`, 10 `cuEventRecord`, 2 DtoD copies per step, plus the Julia
  side of calling the compiled function). Unprofiled `physics_perf_layers.jl` L8 numbers vary between runs
  (physics only 1.46 ms > full 1.24 ms, which can't be compute), consistent with overhead dominating.
  Queued `call_overhead_perf.jl`: compiled no-op taking the whole `simulation` vs compiled `time_step!`.

> Do we really need a separate `initialize_traceable!`? If it's identical is it really that much slower
> on CPU and regular GPU?

- `implicit_init_compare.jl`, both from the same inputs, max |ΔS⁻¹| / max |S⁻¹| over the used degrees and
  the identity error ‖S⁻¹S − I‖ (Float64 S) of each, median of 15 calls, CPU 1 thread:

  | T, L     | LU ms  | traceable ms | ratio | rel. ΔS⁻¹ | ‖S⁻¹S − I‖ LU / traceable |
  | -------- | ------ | ------------ | ----- | --------- | ------------------------- |
  | T31 L8   | 0.18   | 0.07         | 0.40  | 5.3e-7    | 2.6e-6 / 2.9e-6           |
  | T63 L16  | 0.79   | 1.31         | 1.67  | 5.1e-7    | 4.7e-6 / 3.7e-6           |
  | T127 L32 | 4.6    | 28.1         | 6.1   | 8.9e-7    | 5.9e-6 / 5.8e-6           |
  | T255 L64 | 23.7   | 709          | 30    | 2.5e-6    | 7.8e-6 / 1.2e-5           |

  Same results to Float32 precision, equally good inverses. The traceable kernel has one thread per degree
  doing O(nlayers³) serially (R*L recomputed per l + Gauss-Jordan), and `S⁻¹[l, k, r]` is strided by
  `truncation + 1` in k and r (coalesced on GPU, cache misses on CPU); LU calls LAPACK on contiguous
  matrices. On CPU/GPU `reinitialize!` only calls it when Δt changes (~3 times per run), so the 30x
  would cost ~2 s startup at T255L64, but LU also pivots.

> okay, implement your suggestion for the implicit initialization and continue to analyse performance
> differences between CUDA.jl and Reactant

- Implemented (the user committed it as `3c1dd133`): one `initialize!` that builds the operators via
  broadcasting (the former traceable code, works on all architectures) and then calls
  `invert_S!(implicit, ξ, architecture)`: `invert_S_lu!` (transfer R, L, U, W, S, S⁻¹ to the CPU, LU per l,
  `R*L + U*W'` hoisted out of the l loop, NaN guard kept, S⁻¹ transferred back; a scratch matrix instead of
  reusing `L1`) or, for `ReactantDevice`, `invert_S_gauss_jordan!` (the kernel). `initialize_traceable!`
  is gone, -77 lines. CPU A/B against the old `initialize!` (`scratchpad/implicit_ab.jl`): **bitwise
  identical** operators and S⁻¹ for Float32/Float64 at T31L8, T63L16, T127L32. Unit test
  `test/dynamics/implicit.jl` rewritten: S⁻¹S ≈ I for every l and Gauss-Jordan vs LU, 66/66.
  Queued: CUDA (device-built operators == CPU, `implicit_init_compare.jl`), CPU `test/reactant`.
- H100 benchmarks: CUDA job 2490788 finished (committed by the user as `f9468d64`), no regressions vs
  the previous `gpu-nvidia` results (all within -8%..+4%, T86 +67% (L8) / +32% (L16) where the old default
  transform values were outliers next to the matrix transform). Reactant job 2490789: T86L16 now **fails**
  (instead of hanging for hours) in XLA's GPU codegen: `INTERNAL: RET_CHECK failure
  (xla/backends/gpu/codegen/emitters/mlir_kernel_emitter.cc:705) llvm_module != nullptr Failed to
  translate module to LLVM IR`, recorded as failed, the job continues.
- The login node (128 cores, 21 users, load averages up to 65) makes host-bound timings meaningless:
  `physics_perf_layers.jl` L12 physics only 1.00 ms vs L8 1.46 ms. Clean timings now in Slurm job
  **2491189** (`submit_perf_compare.sh`, `perf_compare.jl`): CUDA.jl vs Reactant in separate processes from
  the same environment, T31L8/T63L8/T127L8 × full/physics/dynamics, 7 × 100 steps (min/median), host
  dispatch time before the final sync, and for Reactant compiled no-ops on the whole `simulation` vs on
  the clock only.
- The compiled physics-only step (T31L8) has **183 arguments (91 donated), 97 results**, 33 MB, 50
  arguments ≤ 1 KB. The A40 profile's CUDA API time adds up to ~0.13 ms per step (graph launch 28 µs, 5 host
  functions, 2 DtoD copies 35 µs, 24 graph node updates), wall time 0.94 ms, so most of the host time is
  in the Julia wrapper or XLA's runtime (per-buffer handling). `profile_host.jl` queued: allocations/GC
  per call + Julia sampling profile with C frames.
- GPU checks of the unified implicit initialization (A40, `implicit_init_compare.jl`, `implicit_R_check.jl`):
  works on CuArrays; `U, L, W, L1, L3` bitwise equal to the CPU, `R` (hence S⁻¹) not, because its inputs
  `Δp_geopot_full/half` are computed on the GPU and differ by 1e-8/7e-8 relative; R built on the GPU from
  the CPU's inputs is bitwise equal to the CPU's. A40 timings: `invert_S_gauss_jordan!` 0.07 / 0.7 / 11 /
  169 ms vs `invert_S_lu!` (incl. transfers) 0.3 / 0.9 / 4.0 / 21.5 ms at T31L8 / T63L16 / T127L32 /
  T255L64. CPU `test/reactant` suite on the unified code: all pass (PrimitiveWet 20/20, differentiability 4/4).

> We need to further work on "Reactant does less GPU work but spends its time on the host".

- **H100 job 2491189** (`perf_compare.jl`, clean node, min ≈ median): ms per step, CUDA.jl / Reactant

  | step     | T31L8       | T63L8       | T127L8      |
  | -------- | ----------- | ----------- | ----------- |
  | full     | 1.09 / 1.21 | 2.15 / 0.97 | 4.80 / 4.09 |
  | physics  | 0.77 / 1.53 | 1.53 / 0.69 | 3.51 / 2.73 |
  | dynamics | 1.00 / 0.44 | 2.05 / 0.75 | 4.70 / 3.57 |

  CUDA.jl's host dispatch is 86-91% of the wall time (host bound too), Reactant's 99% at all resolutions:
  the calls are effectively synchronous. Reactant is faster except T31 full/physics; T31 physics only is
  slower than T31 full and than T63 physics, so something specific is slow there.
- **H100 job 2491225** (`call_overhead_perf.jl`): compiled no-op on the whole `simulation` (183 buffers)
  0.10 ms per call, on the clock only 0.05 ms. So flattening/passing the buffers is cheap; the time step's
  1.0-1.5 ms is lost inside the execution.
- Calls are synchronous because the generated argument code calls `XLA.synced_buffer` on every input,
  which waits on the buffer's future, i.e. for the previous step to finish. PJRT orders executions on the
  stream itself. `no_wait_perf.jl` patches `synced_buffer` to not wait (queued on the A40).
- Per-step timeline from the A40 Nsight sqlite (physics only T31L8, 928 µs launch to launch): 24
  `cuGraphExecNodeSetParams` spread over ~330 µs before every launch, `cuGraphLaunch` itself 160-290 µs,
  then events/host functions/2 DtoD copies on other threads ~150 µs, ~100 µs Julia; the GPU runs its 81
  kernels (224 µs busy) starting 118 µs after the launch. So XLA **updates its CUDA graph every call**.
  Graph updates are needed when buffer addresses change: in the compiled physics step 6 of 97 results are
  not aliased to inputs, all scalars (5 × i64, 1 × f32, the clock fields and probably a scalar `Ref`),
  i.e. new buffers every call. Reactant aliases a result to an argument only if the result is the same
  traced object (`Compiler.jl`, `tf.aliasing_output`): true for arrays updated in place, false for
  `clock.time_step_counter += 1`. Test queued: `XLA_FLAGS=--xla_gpu_enable_command_buffer=` (no CUDA
  graphs).
- Multi-step calls (`multistep_perf.jl`, CHUNK steps in a `@trace for`): (1) `@trace for` traces all
  numbers of the captured `simulation` by default and fails on our views (`Base.OneTo{TracedRNumber}`),
  fixed with `@trace track_numbers = false`; (2) then XLA export fails on `tensor.empty`: a zero-length
  array carried through a traced loop (Reactant bug, MWE `reactant_issue_zero_size_loop_mwe.jl`, 10
  lines). Our only zero-length array is `geometry.lond` (unused, `Float64[]` for reduced grids); the
  prototype overrides `get_lond` for reduced grids to a 1-element vector. CPU backend (T31L8, CHUNK = 5):
  compiles (39 s vs 12 s for one step), agrees with single steps to 2e-7 after 5 steps (3.7e-4 after
  64, not bitwise, a different program), no speed-up on CPU (compute bound, ~9.2 ms per step). GPU queued.
- No-wait patch (A40, `no_wait_perf.jl`): **bitwise identical** results after 629 steps; timing noisy on
  the login node (0.61-1.23 vs 0.92-1.53 ms), host dispatch still 100%: the host work (~0.7 ms) is longer
  than the GPU work, so not waiting only hides the GPU time.
- **H100 host profile** (`profile_host.jl`, T31L8 full, 1.28 ms per step, 135 KB and 1036 allocations per
  step, GC 2.6%), main thread split: XLA execute (`RunGpuAsync`: CUDA graph update + launch) 39%;
  Reactant result write-back 21% (`traced_setfield_buffer_at_parent!` → `traced_getfield`, which is
  `@nospecialize` → dynamic `getfield` that allocates per result); `donate_argument!` 7.5%
  (`Tuple(d.buffer for d in carg.data)`, dynamic dispatch per argument); waiting for the previous step
  (`synced_buffer` → `FutureAwait`) 7%; PJRT `ExecutePrepare` 6%; output buffer creation 3%;
  `CopyBufferToDevice` 3% (the 2 DtoD copies per step). So ~0.4 ms per call is type-unstable Julia code in
  Reactant's generated call wrapper. `reactant_call_patches.jl` patches the three Reactant functions
  stage by stage.
- Multi-step on the A40 (`multistep_perf.jl`, CHUNK = 10, T31L8 full): compiles (46 s vs 18 s), **0.40 ms
  per step vs 1.05-1.09 ms** for single-step calls (2.6-2.7x, 2 of 3 repetitions; the 3rd 1.17 ms, login
  node). Agreement 3.4e-4 after 10 steps, 2.3e-2 after 609 (CPU backend: 2e-7 after 5): presumably TF32
  rounding in a differently fused program; to check with `HIGHEST`.
- Upstream issue draft for the zero-length loop-carried array: `reactant-issue-draft-zero-size-loop.md`
  + `reactant-issue-zero-size-loop-mwe.jl`.
- H100, `reactant_call_patches.jl` T31L8 full (patches cumulative, median ms per step / allocations per
  step): baseline 1.21 / 1036; specialized `traced_getfield` 1.23 / 1039 (no effect, the allocations come
  from the runtime field index, not the missing specialization); type-stable `donate_argument!` 1.13 / 932
  (-7%); **no wait for inputs 0.77** (-32%), i.e. 30% faster than CUDA.jl's 1.09 ms. The wait dates back
  to Reactant's original async design (2024, "Refactor compile and test buffer donation"); nothing in
  the history says it is needed for correctness on GPUs. Upstream issue, not a patch in SpeedyWeather
  (type piracy on Reactant internals).
- **H100 job 2492144 results** (median ms per step, PrimitiveWet full):
  - `reactant_call_patches.jl` T63L8: baseline 1.16, specialized `traced_getfield` 1.10, type-stable
    `donate_argument!` 1.04, **no wait 0.67** (-42%; CUDA.jl 2.15).
  - `alias_clock_perf.jl` (T31L8): bitwise identical, but no measurable effect (1.0-1.4 ms per step in
    both variants, the repetitions vary by more than any difference). Dropped.
  - Without command buffers (`XLA_FLAGS=--xla_gpu_enable_command_buffer=`): T31L8 1.01, T63L8 0.91,
    T127L8 3.56 vs 1.21 / 0.97 / 4.09 with (job 2491189), i.e. 7-17% faster on the H100, but on the
    A40 12-30% slower (T31 full 1.19 vs 1.06, physics 0.92 vs 0.70). Not robust, leave the default.
  - **`multistep_perf.jl` CHUNK = 10 (T31L8): 0.41 ms per step vs 1.14-1.69 ms single-step calls, 2.8x**
    (2.7x faster than CUDA.jl's 1.09 ms). Compile 36 s vs 16 s. Agreement 2.1e-4 after 10 steps, 3.1e-2
    after 609. Same 0.41 ms on the A40, so the loop body is now launch bound. With `HIGHEST` (A40):
    **1.2e-6 after 10 steps** (TF32: 2e-4), 5e-3 after 609 (TF32: 3e-2), so the differences are TF32
    rounding in two differently fused programs; 0.44 ms per step chunked vs 1.0-1.7 single-step.
  - Implementation plan for chunked time stepping: `2026-10-14-reactant-multistep-time-stepping.md`
    (status planned, awaiting sign-off).
- H100 job **2492144** (`submit_host_variants.sh`): `reactant_call_patches.jl` (T31, T63),
  `alias_clock_perf.jl` (clock fields and `implicit.Δt` updated in place via `TracedUtils.set_mlir_data!`
  so that Reactant aliases them; first attempt failed on a method ambiguity, `Δt::Dates.Period` now),
  `perf_compare.jl` without command buffers, `multistep_perf.jl` CHUNK = 10.
- Reactant benchmark job 2490789 finished: everything ran except **T86L24** (not T86L16: the Reactant run
  skips T171/T256, so its "Benchmark 8" is T86L24) with the XLA `RET_CHECK` in `mlir_kernel_emitter.cc`.
  Compile times up to 173 s (T86L16), 154 s (T128L24). `reactant-gpu` merged into the repository's
  `benchmark_results.json` (README regenerated with `scratchpad/merge_reactant_benchmarks.jl`: the README
  part of `manual_benchmarking.jl` with `benchmark_suite.jl` included). Note: the overview table shows the
  failed T86L24 as `0.0` SYPD (the per-architecture table says "failed"). Reactant vs CUDA.jl (both with the
  matrix transform, H100, SYPD): T32 5256/5343, T43 5974/3986, T64 3067/2804, T86L8 1025/1085, T86L16
  869/746, T128L8 399/199, T128L16 327/138, T128L24 303/118; vs CUDA.jl's default FFT+Legendre transform
  at T128: 399/271, 327/243, 303/227.

## Handoff 2026-10-14

State at the end of 2026-10-14 (last commit `f5082066`). Read this first, details in the log above
("2026-10-14 afternoon").

Done today (committed by the user): unified implicit initialization (`invert_S_lu!` /
`invert_S_gauss_jordan!`, bitwise identical on CPU, CPU `test/reactant` passes, CUDA checked), its unit
test; CUDA H100 benchmarks.

Uncommitted in the working tree:
- `SpeedyWeather/benchmark/assets/benchmark_results.json` + `README.md`: `reactant-gpu` (H100) merged.
  Overview table shows the failed T86L24 as `0.0` SYPD (maybe make `write_overview` print "failed").
- This doc; new untracked docs: `2026-10-14-reactant-multistep-time-stepping.md` (implementation plan),
  `reactant-issue-draft-call-overhead.md`, `reactant-issue-draft-zero-size-loop.md` +
  `reactant-issue-zero-size-loop-mwe.jl` (and the older issue drafts/MWEs).

Performance findings (H100, PrimitiveWet T31L8, ms per step): Reactant 1.21 vs CUDA.jl 1.09, but the GPU
is busy only ~0.2-0.3 ms; the rest is host time per call (Reactant's call wrapper ~0.4 ms, XLA graph
update + launch ~0.5 ms, plus a wait for the previous step's outputs). Reactant patched to not wait:
0.77 (T63: 0.67, CUDA.jl 2.15). **10 steps per compiled call (`@trace` loop): 0.41**, matches single
steps to 1.2e-6 after 10 steps with `HIGHEST` (2e-4 with TF32).

Open items:
1. **Sign-off on `2026-10-14-reactant-multistep-time-stepping.md`**, three decisions: default
   `steps_per_call`; compile chunk + single step (2.3x compile time) or chunk only; workaround for the
   empty `geometry.lond` on reduced grids. Then implement (prototype: `claude_reactant_gpu_2026-10-13/
   multistep_perf.jl`).
2. Post the upstream issues (call overhead: no wait for inputs -32..42%, type-stable `donate_argument!`
   -7%; zero-length array in traced loops; the earlier drafts).
3. Reactant T86L24 on the H100 fails in XLA's GPU codegen (`RET_CHECK ... mlir_kernel_emitter.cc:705
   Failed to translate module to LLVM IR`); not investigated (module dump was on the compute node's /tmp).
4. Earlier open items still open: Reactant release with `Reactant_jll` >= 0.0.419 (compat bounds,
   `@test_broken` → `@test`), `disable_dot_formation_passes` once released.

Scripts and logs: `/p/tmp/maxgelbr/code/Speedy2/claude_reactant_gpu_2026-10-13/` (`perf_compare.jl`,
`call_overhead_perf.jl`, `profile_host.jl`, `reactant_call_patches.jl`, `no_wait_perf.jl`,
`alias_clock_perf.jl`, `multistep_perf.jl`, `submit_*.sh`, `*-24xxxxx.log`). Clean timings need a
compute node (login node load averages up to 65 make host-bound timings meaningless).

## Handoff 2026-10-16

State on 2026-10-16. Everything since the last commit (`f5082066`) is uncommitted. The details of the work
since 2026-10-14 are in `2026-10-14-reactant-multistep-time-stepping.md` (implementation plan with the
logs of 2026-10-15/16), read that next.

Done since 2026-10-14:
- **Reactant 0.2.293** (Reactant_jll 0.0.421): compat `Reactant = "0.2.293"` in LowerTriangularArrays,
  RingGrids, SpeedyWeatherInternals and SpeedyWeather. The gitignored benchmark manifest was regenerated
  (it pinned 0.2.292 / jll 0.0.418 and couldn't be updated in place).
- **Several time steps per compiled program**: `time_steps!(simulation, n; checkpointing = false,
  mincut = true)` runs `n` steps in one `@trace` loop; `time_stepping!` compiles one program for the
  steps between two outputs and one program for the whole run without output; feedback (progress, NaN
  check) per program. (An `AutomaticDifferentiation` options struct was added and removed again: the
  keywords go straight to `@trace`.)
- **Central compile options**: `SpeedyWeather.reactant_compile_options(; kwargs...)` (=
  `Reactant.CompileOptions(; raise = true, excluded_passes = REACTANT_EXCLUDED_PASSES, kwargs...)`) with
  `SpeedyWeather.REACTANT_EXCLUDED_PASSES = ["concat_view_cse"]`, used by `@maybe_jit`, `time_stepping!`,
  the tests and the benchmark runner. `concat_view_cse` (new in 0.2.293) doesn't terminate on several of
  our programs (`run!` of every ShallowWaterModel hung in `initialize!`).
- **Output with Reactant**: views of concrete arrays are transferred as their parent and copied out on
  the host, the mean sea level pressure is computed on the output's architecture
  (`on_output_architecture`), ~0.65 s → a few ms per snapshot.
- **Differentiation through the time loop** works with `excluded_passes = [REACTANT_EXCLUDED_PASSES;
  "pad_concat_to_concat_pad"; "slice_pad"]` and `raise_first = true` (two Enzyme-JAX bugs, see below):
  BarotropicModel gradients of 4 and 6 steps agree with vanilla Enzyme on the CPU to 2e-7 to 6e-7
  without checkpointing and with `Reactant.Periodic(2)` / `Reactant.Binomial(2)`, 2e-4 on the A40 (TF32).
  `mincut = false` doesn't finish compiling when the loop is a `stablehlo.while`; Reactant unrolls loops
  of at most 4 iterations, so the tests use `n = 6`.
- **Tests**: `test_time_steps.jl` (new), `differentation.jl` with the loop gradients and a
  `compare_with_enzyme` keyword; CI additionally checks that the PrimitiveWetModel loop gradient compiles
  (`compare_with_enzyme = false`, vanilla Enzyme of the primitive equations takes too long for CI).
- **Verification on 0.2.293** (CPU backend, before the last API cleanups): the whole `test/reactant`
  suite passed in two parts (109/109 and 96/96); `test/GPU/reactant.jl` (PrimitiveWet standard script,
  A40) 5/5. After the `on_output_architecture` cleanup: `utilities.jl` + `test_output.jl` 49/49.

Upstream issues, drafted in this folder, not posted:
1. `reactant-issue-draft-pad-concat.md` (+ MWE): `pad_concat_to_concat_pad` hoists the padding out of a
   concatenation also along the concatenation dimension and drops elements. Also in 0.2.292.
2. `reactant-issue-draft-strided-slice-adjoint.md` (+ MWE): the reverse rule of a strided
   `stablehlo.slice` pads the adjoint as `dim - limit`, one short when the limit isn't tight; the loose
   limits come from `GatherOpCanon` (`limit = start + count * stride`, a one-line fix). Also in 0.2.292.
3. `reactant-issue-draft-concat-view-cse-hang.md` (+ MWE): a 0.2.293 regression; the cycle is
   `concat_view_cse` + `broadcast_to_reshape` + `reshape_concat` (minimized by delta debugging), a
   7-op StableHLO module reproduces it.
The earlier drafts (call overhead, zero-size loop, complex `ifelse`, complex store, kernel operand
layout, raised dot precision) are still unposted too.

The PrimitiveWet CI check passes (7/7) but takes 31 min on the loaded login node (three gradient
compiles); the whole `differentation.jl` in its final form passes (Barotropic 6/6 and 12/12,
PrimitiveWet 7/7). The `reactant-gpu` benchmarks on an H100 were rerun (Slurm job 2516570, results
in `SpeedyWeather/benchmark/README.md`, the 2026-10-14 ones backed up in
`claude_reactant_gpu_2026-10-13/bench_0293/`): every SYPD row is faster, 1.2-4x per suite (T32L8
16563 SYPD, 3.1x CUDA.jl on the same GPU), T86L24 now works, details in the multistep plan.

Open items:
1. The PrimitiveWet CI check costs ~30 min (three gradient compiles); reduce to one checkpointing
   option if that is too much for CI.
2. Post the upstream issues; remove the pattern exclusions once they are fixed.
3. Commit the work (everything since `f5082066`).
4. `disable_dot_formation_passes` (released in 0.2.293) is not used yet: evaluate it against the
   TF32 precision of raised reductions on GPUs.
5. **GPU CI times out** (Buildkite CUDA job, `timeout_in_minutes: 40` incl. instantiate): the whole
   GPU suite takes 39.2 min on the A40 here (`claude_reactant_gpu_2026-10-13/gpu_ci_2026-10-16/`,
   timestamped log), of which `reactant.jl` 17.7 min (~3 min precompiling Reactant, 12.9 min testset:
   two XLA programs, the 10-step loop program ~9 min to compile) and `output.jl` 5.3 min. The `@trace`
   loop compiles as fast as the old single step (12.9-13.3 min before, 12.9-13.7 min now; H100
   benchmark 192 vs 196 s); `XLA_FLAGS=--xla_gpu_autotune_level=0` doesn't help (13.8 min). Options:
   own Buildkite step for `reactant.jl`, all output variables with one writer only in `output.jl`,
   kernel-level checks instead of whole models in `vertical_integration.jl`, the second
   `primitive_wet.jl` testset and the graph-cache test.

Scripts and logs: `/p/tmp/maxgelbr/code/Speedy2/claude_reactant_gpu_2026-10-13/` (see the multistep plan
for which script does what; `run_reactant_tests_*.jl` run parts of `test/reactant` in
`env_reactant_0293`).

## Reactant 0.2.294 (2026-10-16 evening)

> A new Reactant version 0.2.294 was just released that should include fixes for several of the upstream
> problems we had. Please check and remove the workarounds. We want to keep the option to add
> compiler_options easily later though, so don't completely remove the option for that

Reactant 0.2.294 / Reactant_jll 0.0.423 closes #3476 (zero-length array in a traced loop), #3478
(strided-slice adjoint) and #3479 (`pad_concat_to_concat_pad`). Every MWE of our issue drafts rerun on it
(CPU backend, `claude_reactant_gpu_2026-10-13/reactant_0294/mwe/`):

| Draft / workaround | 0.2.294 | Action |
| --- | --- | --- |
| pad-concat (AD exclusion `pad_concat_to_concat_pad`) | fixed (1e-7 vs host) | exclusion removed |
| strided-slice adjoint (AD exclusion `slice_pad`) | fixed (gradient error 0) | exclusion removed |
| zero-size loop (`geometry.lond = nothing`, `unused_length` = 1) | fixed | zero-length arrays again |
| `concat_view_cse` loop (`REACTANT_EXCLUDED_PASSES`) | MWE terminates (0.11 s); ShallowWater `initialize!` 32 s and `time_steps!` compile 68 s without exclusions | list emptied, mechanism kept |
| complex `ifelse` in kernels (`componentwise_ifelse`) | not fixed (fails to raise) | kept |
| `ComplexF32` accumulator store (separate real/imag accumulators) | not fixed (NaN imaginary part) | kept |
| `donate_argument!` type stability (#3477/#3484) | fixed | nothing to remove |

The while-loop unroll threshold on the Julia side is still 4 (`Reactant.Compiler.WHILE_UNROLL_THRESHOLD`).
`reactant_compile_options(; kwargs...)` and `REACTANT_EXCLUDED_PASSES` (now `String[]`) stay, so options or
pattern exclusions can be added in one place again; the AD tests' `REACTANT_AD_COMPILE_OPTIONS` is
`reactant_compile_options(; raise_first = true)`. Compat `Reactant = "0.2.294"` in all five packages,
CHANGELOG updated. The CUDA test and benchmark manifests were re-resolved from scratch (old ones in
`reactant_0294/`).

Verification on 0.2.294 (CPU backend unless noted, `reactant_0294/test_*.log`):
- `test_maybe_jit`, `utilities`, `test_transform`, `test_correctness`: 165/165; `test_output`,
  `test_time_steps`: 76/76; differentiation without the PrimitiveWet loop check (make_zero, one step with
  both transforms, Barotropic 6-step loop with `(;)`, `Periodic(2)`, `Binomial(2)`): 48/48, no exclusions.
- A40: `test/GPU/reactant.jl` 5/5, PrimitiveWet NaN check with the `SpectralTransform` 5/5, transform test
  with `HIGHEST` 2/2 + 48/48.
- **New with 0.2.294: the PrimitiveWet loop gradient with `checkpointing = Reactant.Periodic(2)` doesn't
  finish compiling** (killed after 60 min, with and without the old exclusions; the CI runner sat in it for
  2.5 h in the greedy pattern rewriter, e.g. `AddTermsToReduce`). Without checkpointing and with
  `Binomial(2)` it compiles in 3.5 / 3.6 min, the BarotropicModel with `Periodic(2)` is fine. On 0.2.293 the
  whole CI check (three options) took ~30 min. `test_differentiability_time_steps` got an `options` keyword,
  the PrimitiveWet CI check and `differentiation_primitive_wet.jl` run without `Periodic(2)` (comment in the
  test); the PrimitiveWet CI check then passes (`test_C3.log`, 21 min incl. setup). Not minimized or reported
  upstream yet.

> the concat_view_cse should fixed in Reactant 0.2.295. test that.

Reactant 0.2.295 (2026-10-17, `Reactant_jll` 0.0.424; release notes: only a macOS symbol fix, #3487), no
exclusions (`claude_reactant_gpu_2026-10-13/reactant_0295/`): the `concat_view_cse` MWE terminates in 0.1 s
(530 patterns), the ShallowWater `initialize!` takes 31.8 s and the compile of `time_steps!(simulation, 5)`
67.5 s, as with 0.2.294. So `concat_view_cse` works with both; `REACTANT_EXCLUDED_PASSES` stays empty, the
compat (`0.2.294`) allows 0.2.295. The PrimitiveWet loop gradient with `Periodic(2)` is being rechecked on
0.2.295.

Benchmark fix (2026-10-17): the merge of main with #1313 (`68e4088f`) dropped `pressure_gradient_flux!` from
the `functions_and_arguments` of the "Individual dynamics functions" suite (its signature became
`(vars, time_stepping)`), so the name assertion failed; entry restored, the suite runs again on the CPU
(12/12 functions).

## Known limitations

- Particles are not supported on Reactant.
- Differentiation with Reactant needs `raise_first = true` (and Reactant ≥ 0.2.294 for the time loop). It is verified against vanilla Enzyme
  on the CPU for the BarotropicModel (one step, and several steps with each checkpointing option) and
  for one step of the PrimitiveWetModel (`differentiation_primitive_wet.jl`, 2026-10-14); the
  PrimitiveWet loop gradient is only checked to compile and be finite. ShallowWater and PrimitiveDry
  are untested. `mincut = false` doesn't compile for loops of more than 4 steps.
- On NVIDIA GPUs with TF32 the raised matrix products run in TF32 unless disabled (see the docs,
  "Precision of matrix products on GPUs").
- The prescribed seasonal climatologies (e.g. `SeasonalOceanClimatology`) are untested inside
  the compiled time step.

## Future work

- Move the `# TODO: move upstream to Reactant` fixes into Reactant, remove the `inv` and `extrema`
  workarounds once Reactant.jl#2560 and #2387 are fixed.
- On Reactant the implicit operators are recomputed every time step (`reinitialize!`): ~0.5 ms per
  step and ~8 s compile time, although Δt only changes in the spin-up steps. With the spin-up steps
  already outside the compiled loop (`time_stepping!`), initializing them once for the constant Δt of
  the loop is the next step (multistep plan, "Future work").
- Tighten the tolerances of the correctness tests, differences are orders of magnitude smaller
  than the current tolerances.
