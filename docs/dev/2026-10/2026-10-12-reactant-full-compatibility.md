# Full Reactant compatibility incl. the primitive equation models

> Status: **in progress**. All models run on Reactant with `raise = true` (now the default): the CPU-backend test suite passes, and on CUDA the PrimitiveWetModel runs a standard script without NaNs, matching the CPU to ~1e-4 for 10 steps. Needs the upcoming Reactant release with `Reactant_jll` >= 0.0.419. See "Handoff 2026-10-13" and the evening update below it for open items.

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
  Step 1-2 in progress (`physics_perf.jl`, A40, T31L8: physics-only / dynamics-only / full time step,
  raise = true vs false, raised modules dumped to `hlo/physics_perf/`).

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

## Known limitations

- Reactant on GPU (CUDA) is not tested in this PR, `SpeedyWeather/test/GPU/reactant.jl` is still
  not included in the GPU test suite.
- Particles are not supported on Reactant.
- Differentiation with Reactant (Enzyme reverse mode, kernels raised to StableHLO) needs
  `Reactant_jll` >= 0.0.419 (raising kernels that store complex constants, Reactant.jl#3450). With it,
  the loss and gradient of a `BarotropicModel` time step agree with Enzyme on the CPU
  (`test/reactant/differentation.jl`, now `@test` instead of `@test_broken`). ShallowWater and
  primitive equation models are untested.
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
