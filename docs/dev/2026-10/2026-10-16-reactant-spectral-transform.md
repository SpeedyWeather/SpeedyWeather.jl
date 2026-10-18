# FFT-based `SpectralTransform` with Reactant

> Status: **in progress**. See "Handoff 2026-10-16 (to the HPC node with GPU)" first. The FFT-based `SpectralTransform` runs, differentiates and passes the CPU-vs-Reactant correctness tests for all models on Reactant's CPU backend; GPU tests/benchmarks and the default transform on Reactant are open.

Date of initial draft: 2026-10-16

Base revision: c0471b16 (branch `mg/reactant-fft`)

## Originating prompt

> In this branch we made most of SpeedyWeather compatible with Reactant. However this is so far
> limitied to `MatrixSpectralTransform`. Now we want to tackle the regular `SpectralTransform` that
> uses FFTs with FFTW on CPU and CUFFT on GPU. It also uses CUDA Graphs on GPU that might need to be
> treated in a special way for Reactant. Maybe it works with them, maybe we have to turn them off.
> Try to make the `SpectralTransform` run with Reactant, keep protocol in a plan and report to me
> and in the plan

## Revision log

> Why is the dense copy of the legendre polynomials needed for reactant?

Answered in the protocol (2026-10-16): needed for the batched matrix multiply, precomputed rather than
gathered in every compiled program; the GPU kernels with the transposed copy don't raise.

> If they are not needed on Reactant, the regular legendre_polynomials and legendre_polynomials_transposed
> should be set to zero-d arrays

On Reactant `legendre_polynomials` (built, used for the dense copy, then replaced) and
`legendre_polynomials_transposed` are unused. Zero-length arrays break `run!` on Reactant (a traced loop
can't carry them, see the protocol), so they have one latitude / one element there (`unused_length`),
zero-length on CPU/GPU as before.

> one of the differentiation tests, should also test the regular SpectralTransform

`test_differentiability(BarotropicModel; Transform = SpectralTransform)` runs in CI.

> abort the benchmark here. we'll run it on a HPC node

> in the markdown file write a handoff to the HPC node with GPU, so that we can continue the tests there.

Benchmark stopped after T31 and T63 (results in the protocol), "Handoff 2026-10-16 (to the HPC node with GPU)".

> Okay, we now also have the FFTs ready for Reactant. Inspect the plan, we also need to run the tests on GPU
> first, and then if everything passes we want to rerun the main benchmarks

On the HPC node (merged into `mg/new-reactant` as #1312): GPU checks in "Protocol on the GPU node" below.

> oh, then we have adjust the benchmark suite to also run FFTs with Reactant

Open decision 2 resolved: the benchmark suite's default rows run with the `SpectralTransform` on Reactant,
the matrix rows with the `MatrixSpectralTransform` (see the protocol).

## Problem description

On `ReactantDevice` all models use the `MatrixSpectralTransform` (`WhichTransform` in
`SpeedyWeather/src/dynamics/spectral_grid.jl`): dense forward/backward matrices, one matrix
multiply per direction. Its cost grows with the square of the number of grid points times the
number of harmonics, O(N⁴) in the truncation N, and its memory too, so on Reactant runs above
~T150 are infeasible (the benchmarks drop them). The FFT + Legendre `SpectralTransform` (O(N³))
does not run on Reactant yet.

## Background

How the `SpectralTransform` works (`SpeedyTransforms/src/`):
- Fourier (`fourier.jl`): one real FFT per latitude ring with FFTW (CPU) or cuFFT/rocFFT/Metal
  (GPU), plans per ring, batched over the layers (`K`), stored in the transform (serial `K = 1`
  plans and batched plans in a `Dict{Int, Vector{Plan}}` keyed by `K`). On GPU, plans for an
  unplanned `K` are created on first use.
- Legendre: on CPU scalar loops (`legendre.jl`), on GPU KernelAbstractions kernels (`legendre_ka.jl`),
  which need the transposed polynomials and thread index tables built only for `GPU`.
- CUDA/HIP graphs (`gpu_graphs_common.jl`, `SpeedyTransformsCUDAExt.jl`): the batched Fourier
  transform on `CuArray`/`ROCArray` scratch memory is captured once per field buffer and replayed.

How Reactant handles FFTs (Reactant 0.2.293, `ext/ReactantAbstractFFTsExt.jl`, `ReactantFFTWExt.jl`):
- `rfft`, `fft`, `ifft`, `irfft` and `plan_*` on traced arrays lower to `stablehlo.fft` (along the
  trailing dimensions, other dimensions are permuted there first). XLA runs these with DUCC on the
  CPU and cuFFT on CUDA.
- FFTW plans captured in traced structs are converted to Reactant plans (which only store the
  dimensions) by `make_tracer`.

## Protocol

Environment: Reactant 0.2.293 / `Reactant_jll` 0.0.421, Julia 1.11.9, Apple M3 (Reactant CPU backend
only, no GPU here), scratch environment with the five monorepo packages as path sources (the
gitignored `test/reactant` manifest pinned Reactant 0.2.264). Scripts are in the session scratchpad
(not part of the repository).

2026-10-16
- **As is**, `SpectralTransform(spectral_grid)` with `ReactantDevice` fails in the constructor:
  `plan_FFTs` calls FFTW's `plan_rfft` on a view of a `ConcretePJRTArray`, "conversion to pointer
  not defined for ConcretePJRTArray". The rest of the transform would not trace either: the CPU
  Legendre transform is scalar loops (every scalar index a traced slice), the GPU Legendre kernels
  need the transposed polynomials only built for `GPU`, and the Fourier transform falls back to the
  serial (one FFT per ring *and* layer) path for every unplanned batch size `K` on non-GPU
  architectures.
- **CUDA graphs**: the graph-accelerated Fourier transform is selected by dispatch on `CuArray`
  (`ROCArray`) scratch memory, which Reactant arrays never are, so it can't be hit with Reactant.
  Reactant/XLA already runs a whole compiled program as one CUDA graph (Nsight, see
  `2026-10-12-reactant-full-compatibility.md`), so there is nothing to gain from our graphs.
  `gpu_graphs` defaults to `true` for any non-CPU backend though, to be set `false` for
  `ReactantDevice` for clarity.
- **Prototype** (`proto_dense_legendre.jl`, `proto_reactant_transform.jl`), traced Julia, no package
  changes:
  - Fourier: rings of equal length in one `rfft`/`brfft` each, i.e. every north/south pair (octahedral
    grids: `nlat_half` FFTs instead of `2nlat_half`), all rings of a full grid in one FFT; the per-group
    spectra are zero-padded to `nfreq_max` and concatenated, then split into the north/south scratch
    arrays; the inverse gathers the rings back into grid order with one constant permutation. (Later
    replaced by runs of latitudes without gathers, see "Differentiation" below.)
  - Legendre as one `dot_general` per direction: coefficients gathered from the packed lower-triangular
    layout into a dense `(parity, l'', m, k)` array (column m shifted to start at its first degree,
    zero-padded), dense polynomials `(parity, l'', m, j)` with the Legendre shortcut applied as zeros;
    batched over parity (0-based even/odd degree, which gives the north/south symmetry) and order m,
    contracting the degree (inverse) or the latitude (forward). Complex coefficients as stacked
    real/imaginary parts (real `dot_general`).
  - Plain-array version vs `_legendre!` on the CPU, T31, K = 3, all six grids: Float32 ≤ 1.0e-7,
    Float64 ≤ 3.1e-16 relative (max abs), both directions.
  - Traced version, `@compile raise = true` (Reactant CPU backend) vs the CPU `SpectralTransform`, T31, K = 4,
    Float32, relative max abs error: grid→spectral 4e-9 to 7e-8, spectral→grid (with `unscale_coslat`)
    1.5e-7 to 3.3e-7 on OctahedralGaussian, FullGaussian, HEALPix, OctaHEALPix, FullClenshaw and
    OctahedralClenshaw grids (24, 1, 16, 32, 1, 32 FFT groups).
- Timing (`proto_bench.jl`), OctahedralGaussianGrid, K = 8, Float32, median ms. FFTW on the CPU with 1 Julia
  thread, Reactant on its CPU backend (XLA uses its own thread pool, so not a like-for-like comparison with FFTW):

  | T    | FFTW CPU s→g / g→s | prototype s→g / g→s | MatrixSpectralTransform (Reactant) s→g / g→s |
  | ---- | ------------------ | ------------------- | -------------------------------------------- |
  | T31  | 0.086 / 0.077      | 0.12 / 0.14         | 0.69 / 0.51                                  |
  | T63  | 0.50 / 0.46        | 0.36 / 0.35         | 10.8 / 10.5                                  |
  | T127 | 3.3 / 3.2          | 1.5 / 1.3           | 87 / 130                                     |

  So on the Reactant CPU backend the FFT transform is 4-5x faster than the matrix transform at T31 and 60-100x at
  T127. Compile times for both directions 10-45 s (prototype) vs 0.2-3 s (matrix).
- **Reusing the GPU Legendre kernels** (KernelAbstractions, transposed polynomials) instead of the dense layout
  (`proto_ka_kernels.jl`; Reactant's kernel support needs `using CUDA` also on its CPU backend): neither raises.
  The forward kernel aborts the process (`LLVM ERROR: AffineApplyNormalizer: cannot move operand to become a valid
  affine symbol`, the loop bounds are loaded from `lm_indices`), the inverse kernel fails with
  `'affine.load' op operand cannot be used as a dimension id` (T31, T63, T127). Both have data-dependent trip
  counts (the column length per order, the rings per order). The dense Legendre alone: 0.043 / 0.054 ms
  (inverse / forward) at T31, 0.11 / 0.14 ms at T63 (K = 8).
- **Why a dense copy of the polynomials** (question by the user): a batch dimension of `dot_general` must be
  rectangular, the packed lower-triangular columns have `lmax - m + 1` degrees for order m. The dense copy is the
  same polynomials with each column shifted to start at its first degree and zero-padded, split by parity (one
  multiply gives both hemispheres) and with the Legendre shortcut as zeros. It is precomputed for speed, not
  strictly needed: the same array could be gathered from `legendre_polynomials` with constant indices in every
  compiled program instead. Memory about twice the packed polynomials: 0.1 MB (T31), 6 MB (T127), 50 MB (T255),
  400 MB (T511).
- **Package implementation**: `SpeedyTransformsReactantExt` (new, Reactant weak dependency of SpeedyTransforms):
  `plan_FFTs(::ReactantDevice, ...)` returns no plans, `dense_legendre_polynomials(::ReactantDevice, ...)` builds the
  dense copy (new field `legendre_polynomials_dense`, an empty vector on CPU/GPU, no new type parameter),
  `_fourier!`/`_legendre!` for both directions as in the prototype, and `_transform_grid!`/`_transform_spec!`
  wrap the Fourier + Legendre transform in `@maybe_jit` (compiled on the fly outside of compilation).
  `default_gpu_graphs(::ReactantDevice) = false`. Pitfalls on the way: traced arrays have element type
  `TracedRNumber{T}`, so signatures `AbstractArray{<:Complex, 3}` don't match them (the generic method is picked
  instead) and `eltype(field) == eltype(S)` fails (`Reactant.unwrapped_eltype`). `step2_package.jl` vs the CPU
  `SpectralTransform` (T31, K = 1 and 3, six grids): eager and compiled (on views into a fused parent array),
  grid→spectral ≤ 6.7e-8, spectral→grid (`unscale_coslat`) ≤ 5.9e-7 relative, 60/60 tests pass. Construction
  0.3 s (3.2 s for the first, compilation).
- The packed `legendre_polynomials` are empty (zero latitude rings) on Reactant after the dense copy is built
  from them, the transposed copy was already empty (only built for `GPU`). Pitfall: `isempty` of a concrete
  Reactant array checks whether its buffer was deleted (e.g. donated), not whether its length is zero, so
  `isempty(ConcretePJRTArray(Float32[]))` is `false`; the test checks `length(...) == 0`.
- Repository tests: `test/reactant/test_transform.jl` "SpectralTransform on Reactant" (OctahedralGaussian,
  FullGaussian, HEALPix; eager and compiled vs the CPU transform), `create_model`/`test_model` in
  `test/reactant/setup.jl`/`test_correctness.jl` take a `Transform` keyword (default `MatrixSpectralTransform`).
- First model run (`test_model(BarotropicModel; Transform = SpectralTransform)`): "Scalar indexing is
  disallowed" in `gather_rings` during `initialize!`. Single-layer variables are vectors `(npoints,)` and
  `data[range, :]` on a 1-dimensional traced array falls back to scalar indexing; the unit test only had
  `(npoints, 1)` matrices. Fixed: fields and coefficients are read as `(n, columns)` matrices and written back
  in their own shape; the unit test now also covers vectors.

- Barotropic correctness run (`run!` spin-up, several time steps per compiled program): XLA export fails,
  `'tensor.empty' op unsupported op for export to XLA`. The known Reactant bug with a zero-length array
  carried through a traced loop (`reactant-issue-draft-zero-size-loop.md`, worked around for
  `geometry.lond` with `nothing`): here the zero-latitude `legendre_polynomials` and the empty
  `legendre_polynomials_transposed` (that one was empty on Reactant before, too). `nothing` would change the
  field types (a new type parameter, but the type name has only ~10 characters margin to the 1024 limit
  of `test/type_name_length.jl`), so unused arrays have length `unused_length(architecture)`: 0, but 1 on
  Reactant (TODO in the docstring: 0 once fixed upstream). A recursive scan of a Reactant PrimitiveWetModel
  simulation with the `SpectralTransform` finds no zero-length arrays after that (`find_empty.jl`).

- **All four models pass the CPU-vs-Reactant correctness tests with the `SpectralTransform`**
  (`run_correctness_fft.jl`: the repository's `test_model(ModelType; Transform = SpectralTransform)`,
  both models with the `SpectralTransform`, T31, 10-day spin-up, tendencies after one step and variables
  after 100 steps): Barotropic 8/8, ShallowWater 14/14, PrimitiveDry 17/17, PrimitiveWet 20/20, 125 / 88 /
  202 / 277 s each (incl. compilation). CI (`test_correctness.jl`) runs ShallowWater and PrimitiveWet with
  the `SpectralTransform` (+~6 min), the BarotropicModel is covered by the differentiability test.
- `test/reactant/test_transform.jl` passes 48/48, now also with vector-shaped (single-layer) data.
- SpeedyTransforms CPU test suite (Julia 1.11, `--check-bounds=yes`, scratch copy of the test environment
  with the monorepo packages as path sources): all pass (incl. JET dispatch, type name length 1008-1014
  of 1024, AD rules 8151/8151), i.e. the struct and `plan_FFTs` changes leave CPU unaffected.

- `test_differentiability(BarotropicModel; Transform = SpectralTransform)` in the CI part of
  `differentation.jl` (`Transform` keyword for `differentiability_simulations`/`test_differentiability`).

- **Differentiation** (`run_diff_fft.jl`, `test_differentiability(BarotropicModel; Transform = SpectralTransform)`
  with the AD compile options of `differentation.jl`): Julia aborts compiling the gradient,
  `LLVM ERROR: Failed to infer result type(s): "stablehlo.add"(...) : (tensor<57x2x2xcomplex<f32>>,
  tensor<57x2x1xcomplex<f32>>)` (57 = `nfreq_max` at T31). Bisected (`bisect_ad.jl`, gradient of each part
  separately, one process per part): only the forward Fourier transform. Its north/south pairs made picking
  the northern rings out of the concatenated spectra a stride-2 gather, which `GatherOpCanon` turns into a
  strided slice, whose reverse rule pads the adjoint one element short: the Enzyme-JAX bug of
  `reactant-issue-draft-strided-slice-adjoint.md`, not avoided by the `slice_pad` exclusion here.
  Rewritten without any gather (`latitude_runs`): one FFT per run of latitudes j1:j2 with rings of equal length,
  northern and southern rings stacked; the northern rings of a run are contiguous grid points and so are its
  southern ones, so only contiguous slices, `reverse` and concatenations (also no permutation gather in the
  inverse: northern runs, then southern runs reversed, are the grid order). Same number of FFTs as before
  for octahedral (one per latitude) and full grids (one).
  Then all six parts differentiate, directional derivative vs a central finite difference (Float32, T31, K = 2,
  OctahedralGaussian): Fourier forward/inverse 4.0e-5 / 1.6e-4, Legendre forward/inverse 9.5e-5 / 2.2e-5,
  transform grid→spectral / spectral→grid 6.6e-4 / 5.8e-4 relative.
- The first gather-free version stacked northern and southern rings along a new dimension, `(nlon, n, 2, K)`.
  The 3-layer case of the transform test then failed in Reactant's optimization ("MLIR pass pipeline "all"
  failed: 'stablehlo.multiply' op requires compatible types", `tensor<2x3x1x20xf32>` times
  `tensor<3x1x2x20xf32>`, from the `rmul!` scaling of `brfft`). `brfft`/`irfft` on such arrays alone are
  correct for all shapes tried (`mwe_brfft.jl`), so a pass rewrites the larger program wrongly; not minimized
  further. Northern and southern rings concatenated along the ring dimension instead, `(nlon, 2n, K)`, and
  split by contiguous slices: the transform test passes 48/48 again.
- AD of the final version (`bisect_ad.jl`, K = 3, Float32): Fourier forward/inverse and transform
  spectral→grid agree with central finite differences to 6e-6 to 4e-4 on OctahedralGaussian, HEALPix and
  FullGaussian grids; grid→spectral 1.2e-3 (full grid), 3% (HEALPix) and 8% (octahedral) in Float32, where
  the directional derivative is tiny compared to the loss and the Float32 difference quotient is rounding
  noise: in Float64 (the loss is quadratic, so the central difference is exact up to rounding) AD and FD agree
  to 1.3e-12 (octahedral) and 1.4e-11 (HEALPix).

- PrimitiveWetModel time step (`bench_model.jl`, 8 layers, OctahedralGaussianGrid, `feedback = nothing`,
  compiled `time_step!`, median of 10, ms; CPU = SpeedyWeather on the CPU with FFTW, 1 Julia thread; Reactant
  on its CPU backend, which uses all cores, so not comparable with the CPU column). Stopped by the user before
  T127 (to be run on an HPC node instead):

  | T   | CPU (FFTW) | Reactant matrix | Reactant FFT | compile matrix / FFT [s] |
  | --- | ---------- | --------------- | ------------ | ------------------------ |
  | T31 | 6.03       | 7.39            | 3.53         | 27.8 / 26.8              |
  | T63 | 24.0       | 68.5            | 13.0         | 10.7 / 16.6              |

  So for the whole model on Reactant's CPU backend the `SpectralTransform` is 2.1x (T31) and 5.3x (T63) faster
  than the `MatrixSpectralTransform`, at similar compile times.

## Protocol on the GPU node (2026-10-16, login02 A40 and Slurm H100)

Branch `mg/new-reactant` at `296bb05d` (#1312 merged). Environments `SpeedyWeather/test/GPU/CUDA` and
`SpeedyWeather/benchmark` re-resolved (both record `SpeedyTransformsReactantExt`; Reactant 0.2.293, jll
0.0.421), backups of the old manifests and all scripts/logs in
`/p/tmp/maxgelbr/code/Speedy2/claude_reactant_gpu_2026-10-13/gpu_fft_2026-10-16/`.

- The committed `test/GPU/reactant.jl` runs only the `MatrixSpectralTransform`
  (`for Transform in (MatrixSpectralTransform, ) # SpectralTransform`), so the NaN check with the
  `SpectralTransform` runs from a scratch copy (`nan_check_spectral.jl`).
- **First H100 run (job 2526764): the `SpectralTransform` constructor fails on Reactant's GPU backend**,
  "Scalar indexing is disallowed", `setindex!(::ConcretePJRTArray, ...)` from
  `legendre_polynomials[:, j] = LowerTriangularArray(...)` (`spectral_transform.jl`, filling the
  polynomials latitude by latitude). On Reactant's CPU backend a scalar `setindex!` writes straight into
  the host buffer (no `assertscalar`), which is why the M3 runs never saw it; on the GPU backend it is an
  error (with `allowscalar` it would compile a program per element). Fix: on Reactant the packed
  polynomials are computed on a CPU spectrum (they are only used there to build the dense copy, which
  transferred them to the host anyway, and are then replaced by the minimal array); CPU/GPU unchanged.
- **With the fix (H100, job 2526880): the GPU checks pass.**
  - PrimitiveWet NaN check with the `SpectralTransform` (default resolution, 10 steps, `run!`): 5/5, 10.9 min
    (as long as with the matrix transform).
  - `test/reactant/test_transform.jl` on the GPU backend with `dot_general_precision = HIGHEST`: all pass
    (matrix testset 2/2, `SpectralTransform` 48/48). With the default precision (TF32) the matrix testset
    fails 1 of 2 (`≈`), which aborts the file before the `SpectralTransform` testset, as expected.
  - Relative errors vs the CPU `SpectralTransform` (T32, 8 layers), grid→spectral / spectral→grid:

    | Grid | precision | `MatrixSpectralTransform` | `SpectralTransform` |
    | --- | --- | --- | --- |
    | OctahedralGaussian | TF32 | 1.5e-4 / 9.6e-4 | 4.2e-4 / 6.3e-4 |
    | FullGaussian | TF32 | 7.8e-5 / 6.2e-4 | 4.4e-4 / 9.0e-4 |
    | HEALPix | TF32 | 7.3e-5 / 5.6e-4 | 4.5e-4 / 1.1e-3 |
    | OctahedralGaussian | HIGHEST | 2.7e-7 / 1.7e-6 | 6.6e-8 / 4.1e-7 |
    | FullGaussian | HIGHEST | 2.0e-7 / 2.3e-6 | 5.1e-9 / 2.3e-7 |
    | HEALPix | HIGHEST | 2.7e-7 / 1.5e-6 | 3.3e-8 / 2.2e-7 |

    So with `HIGHEST` the `SpectralTransform` on the GPU is 4-40x more accurate than the matrix transform;
    with TF32 both are at the TF32 level (1e-4 to 1e-3).
- CUDA.jl GPU suite on the A40 with the merged code (committed `test/GPU/reactant.jl`, matrix transform):
  **44/44 testsets pass, 34.1 min** (spectral transforms, CUDA graphs, models, output, Reactant NaN check
  with the matrix transform 12.8 min).
- Benchmarks (code `4c0257ff`, H100, submitted after all GPU checks passed): the main suite
  `manual_benchmarking.jl reactant-gpu` (Slurm 2527277, `benchmarks-reactant-fft-2527277.log` in the claude
  dir, writes `SpeedyWeather/benchmark/README.md` and `assets/benchmark_results.json`, previous versions in
  `gpu_fft_2026-10-16/*_before_fft.*`) and the two-transform comparison of the handoff below, T31-T255,
  TF32 and HIGHEST (Slurm 2527278, `bench_transforms.jl`, `gpu_fft_2026-10-16/bench-transforms-2527278.log`).
- **Benchmark suite** (`SpeedyWeather/benchmark`, open decision 2): on Reactant `build_model` now honours
  the row's transform, `:default` builds the FFT-based `SpectralTransform` explicitly (the models' default on
  Reactant, `WhichTransform`, is the matrix transform at every resolution), `:matrix` the
  `MatrixSpectralTransform`. The Reactant method of `adapt_to_architecture!` (forced every row to the matrix
  transform and dropped T > 150) is removed: `define_benchmarks.jl` creates matrix rows only up to
  `MAX_MATRIX_TRUNCATION` anyway, so nothing is left to adapt (the generic no-op hook stays). So the Reactant
  resolution suite has both transforms up to T128 and the `SpectralTransform` up to T255, and every other
  suite (default rows) runs with the `SpectralTransform` on Reactant. Driver docs and README text adjusted.
- Finding, not changed: on CUDA.jl the resolution suite's "default" rows at T ≤ 100 are the
  `MatrixSpectralTransform` too (the model default on GPUs below T100), their SYPD equals the matrix rows
  (T32L8 5439 vs 5343), while the overview labels them "LT+FFT".

### GPU performance analysis (2026-10-16 evening)

> I am surprised the FFT based transform is faster than matrix transform even for small problems with
> Reactant. Can you try to analyse and profile if there's anything we can still optimize?

On the GPU it isn't: the matrix transform is faster wherever both run. PrimitiveWet, 8 layers, H100, TF32,
ms per step (`bench_transforms.jl`, job 2527486; the "FFT faster" numbers above are Reactant's CPU backend):

| T | Matrix, 10 steps per call | FFT, 10 steps per call | matrix / FFT |
| --- | --- | --- | --- |
| T31 | 0.451 | 1.12 | 2.5x faster |
| T63 | 0.687 | 2.25 | 3.3x faster |
| T127 | 5.27 | 5.67 | 1.08x faster |
| T255 | (too large) | 13.0 | |

Profiles (A40, Nsight with `--cuda-graph-trace=node`, scripts and reports in
`claude_reactant_gpu_2026-10-13/profile_fft_2026-10-16/`):
- **Per time step both transforms make the same four calls**: one grid→spectral (73 columns, all tendencies
  batched) and three spectral→grid (2, 16 and 33 columns). Matrix: four matrix products reading 53 MB (T31) /
  729 MB (T63) of dense matrices per step, memory-bound. FFT: 3-4 batched Legendre `dot_general`s (1 MB / 8 MB)
  and **96 (T31) / 192 (T63) FFT ops per step**, one per run of equal-length rings and call (octahedral grid:
  one per latitude).
- **XLA can't capture FFTs in its CUDA graphs** (command buffer command types in this build: FUSION, CUBLAS,
  CUBLASLT, CUDNN, CUSTOM_CALL, COLLECTIVES, CONDITIONAL, WHILE, DYNAMIC_SLICE_FUSION, no FFT). With the matrix
  transform the whole step is one graph launch; with FFTs the program is split and most kernels are launched
  one by one from the host.
- Full step, A40, T31: matrix 0.418 ms per step (350 µs GPU time, 106 kernels, 1.2 graph launches per step);
  FFT 1.466 ms (835 µs GPU time, 385 kernels, 4.2 graph launches + 269 separate kernel launches per step).
  FFT step GPU time: FFT kernels 302 µs (103 kernels), concatenate/slice/pad 172 µs (88), transposes 119 µs
  (64), **cuBLAS `scal` 91 µs (72 kernels: XLA normalizes every inverse FFT, `brfft` then multiplies back)**,
  matrix products 10 µs, other fusions 141 µs. Matrix at T63: 1.48 ms per step, its matrix products 983 µs of
  1386 µs GPU time (bandwidth bound).
- Transform alone (forward + inverse, 73 columns, traced loop of 20, A40): T31 matrix 0.084 ms vs FFT 0.534 ms
  (6 vs 135 kernels per roundtrip, FFT GPU time 282 µs: FFTs 145, concat/slice/pad 92); T63 matrix 0.81 ms vs
  FFT 1.32 ms.

### Re-audit of the Reactant `SpectralTransform` (2026-10-16 night)

> We don't see much speed ups with Reactant LT+FFT versus CUDA.jl. Is this expected because we already use
> CUDA Graphs, or is there still something we can do? Re-audit the SpectralTransform for Reactant

Main benchmark suite on the H100 (`reactant-gpu` 2026-10-16 vs `gpu-nvidia` 2026-10-14, PrimitiveWet, SYPD):
Reactant LT+FFT vs CUDA.jl LT+FFT at T128-T171 0.93-0.98x, T256 1.09-1.21x; at T ≤ 86 CUDA.jl's "default"
rows are its matrix transform (`WhichTransform`), there Reactant LT+FFT is 0.51-1.09x; Reactant's matrix
transform is 1.2-3.1x CUDA.jl. Reactant only pulls ahead where the whole step runs as one CUDA graph, which
the FFT path never does (XLA can't capture FFTs). At high resolution both spend most of the transform time in
the same cuFFT kernels.

Inefficiencies of `SpeedyTransformsReactantExt.jl`:
1. One FFT op per latitude run and transform call (octahedral grids: per latitude, 96-384 per step); each
   breaks XLA's CUDA graph, adds a cuBLAS `scal` (XLA normalizes `IRFFT`, `brfft` multiplies back) and pre/post
   kernels.
2. All `nfreq_max` frequencies computed and zero-padded, only `mmax` used (e.g. 32 of 57 at T31).
3. Fourier and Legendre communicate through the scratch arrays `f_north`/`f_south`/`g_north`/`g_south` with a
   `permutedims` on each side, real device writes/reads of loop-carried buffers inside the traced loop.
4. Per run slices, `cat`, zero padding and `reverse`, then concatenation over runs (T31 step: 172 µs
   concatenate/slice/pad + 119 µs transposes of 835 µs GPU time; T63 transform alone: 354 µs of 834 µs).
5. `dot_general_complex` stacks and unstacks real/imaginary parts per product.

**Prototype: Fourier as one batched matrix product** (`claude_reactant_gpu_2026-10-13/dft_proto_2026-10-16/
dft_transform_proto.jl`): since only `mmax` frequencies are needed, the Fourier transform of every ring is a
truncated DFT, precomputed per latitude: forward `W[2mmax, nlon_max, nlat_half]` (real/imaginary rows, quadrature
weights and longitude offsets folded in), inverse `B[nlon_max, 2mmax, nlat_half]` (brfft of the first `mmax`
frequencies, factor 2 and Nyquist handling); one constant-index gather of the grid into padded rings
`(nlon_max, nlat_half, north/south, K)` and one back; Legendre unchanged. Two gathers and four `dot_general`s
per roundtrip, no FFT ops, so XLA keeps the program in one CUDA graph. Memory O(N³): 2 × 0.7 / 5 / 39 MB at
T31 / T63 / T127 (≈ 2 × 300 MB at T255).

Forward + inverse, 73 columns, traced loop of 20, A40, TF32 (relative error vs the CPU transform 4-7e-4 for the
prototype, the TF32 level of the other two):

| T | DFT prototype | current LT+FFT | MatrixSpectralTransform |
| --- | --- | --- | --- |
| T31 | 0.045 ms | 0.452 ms | 0.083 ms |
| T63 | 0.149 ms | 1.250 ms | 0.797 ms |
| T127 | 0.851 ms | (lost in the log) | 10.24 ms |
| T255 | 5.42 ms | 9.76 ms | (too large) |

With `dot_general_precision = HIGHEST`: errors vs the CPU transform 1.3e-7 / 5.1e-7 (T31) and 2.0e-7 / 1.1e-6
(T127) forward / inverse; times T31 0.051 (DFT) / 0.449 (FFT) / 0.111 ms (matrix), T127 0.957 / 3.231 /
17.67 ms. So the prototype is 1.8x (T255) to 10x (T31) faster than the current FFT path and 2-18x faster than
the matrix transform in the transform alone, at the same accuracy. Open: the full time step, odd `nlat` (equator ring: its southern slot →
the zero row), gradients (constant-index gathers, scatter-add adjoint). Implementing it needs an
implementation plan (it replaces `_fourier!`/the scratch round trip on Reactant), independent of it: merging
the three spectral→grid calls per step (see above).

## Handoff 2026-10-16 (to the HPC node with GPU)

State: branch `mg/reactant-fft`, code committed by the user as `14a458a1` ("FFTs working with Reactant"),
uncommitted: the CHANGELOG line and the last edits of this document. **The branch isn't pushed yet**:
commit and push here, then `git fetch && git checkout mg/reactant-fft` in the HPC clone
(`/p/tmp/maxgelbr/code/Speedy2/`).

Done (Reactant CPU backend, Apple M3, see the protocol): the FFT-based `SpectralTransform` constructs, runs
eagerly and compiled, passes the CPU-vs-Reactant correctness tests of all four models, and its gradients agree
with vanilla Enzyme on the CPU (BarotropicModel, one step). Never run on a GPU.

### Environment

- Reactant 0.2.293 / `Reactant_jll` 0.0.421 (the earlier `claude_reactant_gpu_2026-10-13/env_reactant_0293`
  should do; `SpeedyWeather/test/GPU/CUDA` for the GPU tests). Julia 1.11 or 1.12: Julia 1.10 ignores
  `[sources]` (would load the registered SpeedyWeather, check `pathof(SpeedyWeather)`).
- **Re-resolve every environment that has the monorepo packages as path sources**
  (`julia --project=<env> -e 'using Pkg; Pkg.resolve()'`): SpeedyTransforms has a new weak dependency on
  Reactant and a new extension `SpeedyTransformsReactantExt`, and the manifest records a path package's
  extensions. Without it the extension doesn't load and `SpectralTransform(spectral_grid)` on
  `ReactantDevice` fails in FFTW: "conversion to pointer not defined for ConcretePJRTArray".
  Check: `Base.get_extension(SpeedyWeather.SpeedyTransforms, :SpeedyTransformsReactantExt) !== nothing`.
- `using CUDA` is needed for Reactant to compile KernelAbstractions kernels, also on its CPU backend.
- Timings only on a compute node (Slurm, e.g. as `submit_perf_compare.sh`), the login node's load makes
  them meaningless.

### To do on the GPU node

1. **NaN check on CUDA**: `SpeedyWeather/test/GPU/reactant.jl` runs the PrimitiveWetModel (default
   resolution) for 10 steps with both transforms now:
   `julia --project=SpeedyWeather/test/GPU/CUDA SpeedyWeather/test/GPU/runtests.jl` (or only that file with
   `using Test, SpeedyWeather` first). Expected: 10/10.
2. **Transform on CUDA vs the CPU**, both directions, with and without
   `Reactant.with_config(dot_general_precision = Reactant.PrecisionConfig.HIGHEST)`: the Legendre transform
   is a traced `dot_general`, so TF32 by default on the A40/H100 (expect ~1e-4 relative with TF32, ~1e-6
   with HIGHEST). The quickest way is `test/reactant/test_transform.jl` with
   `Reactant.set_default_backend("gpu")` (it compares against the CPU with rtol 1e-5, so it should fail
   with TF32 and pass with HIGHEST).
3. **Benchmark** on a compute node (H100, script below): PrimitiveWet time step with the
   `MatrixSpectralTransform` vs the `SpectralTransform` on Reactant (GPU), T31-T255 (matrix only up
   to T127), with TF32 and with HIGHEST. CUDA.jl reference: the `gpu` results (H100) in
   `SpeedyWeather/benchmark/README.md`. Note: a single compiled `time_step!` per call is host-bound at low
   resolution on GPUs (see `2026-10-14-reactant-multistep-time-stepping.md`), so also time
   `time_steps!(simulation, 10)` per call (compiled loop), which is how `run!` runs.
4. **Decide the default transform on Reactant** (`WhichTransform` in
   `SpeedyWeather/src/dynamics/spectral_grid.jl`, currently always the matrix transform for
   `ReactantDevice`) from 3., see "Open decisions".
5. Optional: the whole `test/reactant` suite in one go on the CPU backend (only run in parts so far):
   `Reactant.set_default_backend("cpu")` before including `runtests.jl` (the GPU node would pick the GPU);
   the new parts are the ShallowWater/PrimitiveWet runs with the `SpectralTransform` in `test_correctness.jl`
   (~3.5 min) and `test_differentiability(BarotropicModel; Transform = SpectralTransform)` (~3.5 min).

Benchmark script (the one used for the table in the protocol, for the GPU; `TRUNCS=31,63,127,255`,
`PRECISION=highest` for Float32 matrix products):

```julia
using CUDA, Reactant, SpeedyWeather
Reactant.set_default_backend(get(ENV, "REACTANT_BACKEND", "gpu"))
using SpeedyWeather: ReactantDevice
const options = SpeedyWeather.reactant_compile_options()
const precision = get(ENV, "PRECISION", "default") == "highest" ? Reactant.PrecisionConfig.HIGHEST : Reactant.PrecisionConfig.DEFAULT

function median_ms(f, sync; n = 20)
    f(); sync()
    ts = map(1:n) do _
        t0 = time_ns(); f(); sync(); (time_ns() - t0) / 1.0e6
    end
    return sort(ts)[cld(n, 2)]
end

# a Reactant model at truncation `trunc` with `Transform`, spun up, clock set for one step
function reactant_simulation(Transform, trunc)
    spectral_grid = SpectralGrid(; truncation = trunc, nlayers = 8, architecture = ReactantDevice())
    model = PrimitiveWetModel(spectral_grid; spectral_transform = Transform(spectral_grid), feedback = nothing)
    simulation = initialize!(model)
    run!(simulation, steps = 2)
    initialize!(simulation, steps = 1)
    return simulation
end

# parent: synchronize(::SubArray) is a no-op
sync(simulation) = () -> Reactant.synchronize(parent(simulation.variables.prognostic.vorticity.data))

Reactant.with_config(; dot_general_precision = precision) do
    for trunc in parse.(Int, split(get(ENV, "TRUNCS", "31,63,127,255"), ','))
        for Transform in (MatrixSpectralTransform, SpectralTransform)
            Transform == MatrixSpectralTransform && trunc > 127 && continue
            simulation = reactant_simulation(Transform, trunc)
            tc = @elapsed step! = @compile compile_options = options SpeedyWeather.time_step!(simulation)
            t1 = median_ms(() -> step!(simulation), sync(simulation))
            tc10 = @elapsed steps! = @compile compile_options = options SpeedyWeather.time_steps!(simulation, 10)
            t10 = median_ms(() -> steps!(simulation, 10), sync(simulation)) / 10
            nans = any(isnan, Array(parent(simulation.variables.prognostic.vorticity.data)))
            println("T$trunc $(nameof(Transform)): $(round(t1, sigdigits = 3)) ms per step (single steps), ",
                "$(round(t10, sigdigits = 3)) ms per step (10 per call), compile $(round(tc, digits = 1)) / ",
                "$(round(tc10, digits = 1)) s", nans ? ", NaNs!" : "")
        end
    end
end
```

## Summary of changes

SpeedyTransforms (version stays `1.0.0+DEV`)
- `Project.toml`: Reactant as weak dependency (compat 0.2.293), new extension `SpeedyTransformsReactantExt`.
- `SpectralTransform`: new field `legendre_polynomials_dense::VectorType` (no new type parameter, an empty
  vector on CPU/GPU), on Reactant `legendre_polynomials` and `legendre_polynomials_transposed` minimal
  (`unused_length`: one latitude / one element, a traced loop can't carry zero-length arrays); the FFT planning goes through
  `plan_FFTs(architecture, planned_K, NF, grid, nlayers, scratch, rings, nlons)` (allocates the planning
  array itself, so Reactant allocates nothing for it); `default_gpu_graphs(::ReactantDevice) = false`;
  `dense_legendre_polynomials(architecture, ...)` hook (empty fallback in `legendre.jl`); `show` counts the
  dense copy.
- `SpeedyTransformsReactantExt`: no FFT plans, the dense polynomials, Fourier and Legendre transforms for
  `SpectralTransform{NF, <:ReactantDevice}` (Legendre as in the prototype, Fourier per run of latitudes
  without gathers, see "Differentiation" in the protocol), `_transform_grid!`/`_transform_spec!`
  compiled with `@maybe_jit` outside of compilation.

SpeedyWeather
- No source changes. `WhichTransform` still returns the `MatrixSpectralTransform` for `ReactantDevice`
  (see "Open decisions").
- Tests: `test/reactant/test_transform.jl` (new testset), `Transform` keyword for `create_model`,
  `test_model`, `differentiability_simulations` and `test_differentiability`; `test/GPU/reactant.jl` runs
  the PrimitiveWet NaN check with both transforms.
- Docs: `architectures_gpu.md` (Reactant support: `SpectralTransform` works, how, CUDA graphs off; the TF32
  note also covers its Legendre transform). CHANGELOG line (references #1304, as this branch continues
  `mg/new-reactant`; to change if it becomes its own PR).
- Docs: also `differentiability.md` (Reactant differentiation works with either transform).

## Testing and verification

All on Reactant's CPU backend (Julia 1.11.9, Reactant 0.2.293), final code:
- `test/reactant/test_transform.jl`, "SpectralTransform on Reactant": 48/48 (three grids, no layer dimension /
  one / three layers, eager and compiled vs the CPU transform, no FFT plans, `gpu_graphs == false`, minimal
  unused polynomial arrays).
- `test/reactant/test_correctness.jl` functions with `Transform = SpectralTransform`: Barotropic 8/8,
  ShallowWater 14/14, PrimitiveDry 17/17, PrimitiveWet 20/20 (CI runs ShallowWater and PrimitiveWet).
- `test/reactant/differentation.jl`, `test_differentiability(BarotropicModel; Transform = SpectralTransform)`:
  6/6 (loss vs CPU, Reactant gradient vs vanilla Enzyme on the CPU), 3.5 min.
- Gradients of each part of the transform vs finite differences on three grids (see the protocol).
- SpeedyTransforms CPU test suite: all pass.
- Not run here: the whole `test/reactant/runtests.jl` in one go (its parts above were run separately with
  scratch runners, the MatrixSpectralTransform parts unchanged), anything on a GPU.

**To run on a GPU** (CUDA, login02/H100): `test/GPU/reactant.jl` (PrimitiveWet NaN check, now with both
transforms) and a benchmark of both transforms with Reactant on the GPU, which decides the default (below).

## Documentation changes

`docs/src/architectures_gpu.md` (Reactant support: the `SpectralTransform` works, how it is computed,
CUDA graphs off, CPU-backend speed, the TF32 note covers its Legendre transform), `docs/src/differentiability.md`,
CHANGELOG.

## Known limitations

- Not tested on a GPU (no GPU here). The batched Legendre transform is a `dot_general` that Reactant traces
  itself, so on TF32 GPUs it runs in TF32 unless `dot_general_precision = HIGHEST` (documented).
- Compile time of the `SpectralTransform` is longer than of the `MatrixSpectralTransform` (both directions
  10-45 s vs 0.2-3 s in the prototype; one FFT per run of latitudes, i.e. per latitude on octahedral grids).
- The Legendre transform does about twice the floating point operations of the CPU loops (the dense
  layout's zero padding of the triangle) and the dense polynomials take about twice the memory of the packed
  ones (400 MB at T511).
- Workarounds for Reactant/Enzyme-JAX bugs: no gathers in the Fourier transform (strided slice adjoint, fixed
  in 0.2.294, but the gather-free version is kept); no 4-d arrays around the FFTs (an optimization pass fails,
  not minimized). The minimal unused arrays (`unused_length`) are gone with Reactant 0.2.294: zero-length again.

## Open decisions

1. **Default transform on Reactant**: `WhichTransform` still picks the `MatrixSpectralTransform` for every
   `ReactantDevice` model. On Reactant's CPU backend the `SpectralTransform` is faster at all resolutions
   tested, on GPUs untested. Proposal: the GPU rule (matrix up to T100, `SpectralTransform` above) once
   a GPU benchmark confirms it, or the `SpectralTransform` always if it also wins on GPUs at low resolution.
2. The benchmark runner (`benchmark_suite.jl`) uses the `MatrixSpectralTransform` for Reactant and skips
   T > 150; with the `SpectralTransform` the Reactant benchmarks could cover all resolutions.
3. CHANGELOG references #1304 (this branch continues `mg/new-reactant`), to change if it becomes its own PR.

## Future work

- Post upstream: the strided-slice adjoint issue now also hit by gathers of north/south ring pairs (draft
  exists); minimize the 4-d FFT pass failure; Reactant's `isempty` meaning "buffer deleted" for concrete
  arrays is surprising (maybe by design).
- Pack the dense Legendre layout tighter (e.g. pair column m with column mmax - m + 1 to fill the triangle).
