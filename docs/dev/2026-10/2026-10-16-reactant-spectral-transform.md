# FFT-based `SpectralTransform` with Reactant

> Status: **in progress**. See "Handoff 2026-10-16 (to the HPC node with GPU)" first. The FFT-based `SpectralTransform` runs, differentiates and passes the CPU-vs-Reactant correctness tests for all models on Reactant's CPU backend; GPU tests/benchmarks and the default transform on Reactant are open. A Fourier stage as truncated DFT (no FFTs) is planned for a follow-up PR, see the last section.

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

> For the SpectralTransform, why do you estimate a DFT approach to be faster than the pure MatrixSpectralTransform
> approach? MatrixSpectralTransform just does a single transform

Measured (transform alone), not estimated. One product, but O(N⁴) work and memory, bandwidth bound; the DFT +
Legendre factorisation is O(N³). Numbers in "Follow-up plan: Fourier transform as a truncated DFT on Reactant".

> Would you also estimate the DFT to be faster without Reactant? With regular CUDA.jl

Estimate only: no gain at T ≤ 64 (CUDA.jl is host-bound there), ~1.5x at T86 against the matrix transform, ~1.1x
above T100 against LT+FFT. Out of scope of the follow-up plan, see there.

> Would it be possible to somehow run Reactant FFTs with CUDA Graphs though?

Not without changing XLA or Reactant: `FftThunk` can't be recorded into XLA's command buffers, neither in the
Reactant_jll build nor upstream. Capturing them would remove only the launch overhead. See the follow-up plan.

> okay make a more detailed plan of the DFT approach to implmeent later in the already existing markdown document.
> We defer this to a follow up PR

Added "Follow-up plan: Fourier transform as a truncated DFT on Reactant" (deferred, not part of this PR).

> okay we swtichted to a new branch. In this branch, to implement the new transform for Reactant. Thsis should be
> implemented for reactant only. If you implement a new FFT logic, it should be tested against standard ffts in
> correctness and precision. Both while implmenting but also in the unit tests. this is very important.

Implementation started on branch `mg/reactant-dft`, see "Protocol (implementation)" in the follow-up plan.

> Why do we now have two versions of legendre and fourier in the extension?

Two Fourier transforms by design (the FFTs stay for `fourier_dft = false`, above `DFT_MAX_TRUNCATION`), but the
Legendre transform only needed one: the first implementation kept the FFT path's `_fourier!`/`_legendre!` (complex
scratch arrays, `dot_general_complex`) next to the new all-real DFT functions. Proposed to share the all-real
Legendre transform.

> yes

Refactored: one Legendre transform (`legendre_forward`/`legendre_inverse`) for both Fourier transforms
(`dft_fourier_*`, `fft_fourier_*`, same interface), no scratch memory on Reactant at all, see the protocol.

> do the tests really need to implement yet another fft variant? Can't we just reuse the regular CPU spectraltransform
> of the model?

No: the ring-by-ring references of `test_transform.jl` now call the CPU `SpectralTransform`'s `_fourier!` (FFTW) and
only apply the solid angles and longitude offsets as its Legendre transform does. Same values as the earlier
FFTW-per-ring helpers (forward ≤ 2.8e-16, inverse identical, six grids).

> continue with the spectral transform work

Steps 6-8 of the follow-up plan (GPU benchmarks, profiles, defaults, docs), see the protocol (2026-10-18 evening).

> re-audit the code of the extension. We don't to keep code duplication at a minimum but also not change much in the
> main non-extension code.

Audit of `SpeedyTransformsReactantExt` (protocol): one dense (parity, order, degree) layout and index function for
coefficients and polynomials (three copies of the same loop before), `stack_pair` for five `cat`/`reshape` pairs,
`S.nlon_max`/`S.nfreq_max` instead of recomputing them; main code: one comment.

> I also want you to prepare a markdown for a github PR with an overview of both the performacne across T32, T64,
> T128, T256 but also the memory usage. Also write down the basic math of the DFT you solve

`2026-10-18-reactant-dft-pr-description.md` (next to this plan): math, H100 performance, memory of the precomputed
arrays and of the model run, accuracy, changes.

> add to that the precision of the transform, how accurate it is compared to a F32 baseline. Also include numbers
> for doing the transform both in TF32 and F32

`precision_table.jl` (H100): matrix, FFT path and DFT on Reactant with TF32 and with Float32 (`HIGHEST`) vs the
Float32 CPU `SpectralTransform` (FFTW) baseline, all against a Float64 reference, T32 to T256, both directions; a
"Precision" section in the PR description.

> we also need to test if the model is table for long roll outs with TF32. Test if T256L8 is stable for 10 years
> with a job

`long_run_stability.jl` (Slurm): PrimitiveWetModel T256L8, DFT, TF32, 10 years of 365 days, a compiled program per
30 days with a NaN/Inf check and global diagnostics after each (CSV), see the protocol.

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
the zero row), gradients (constant-index gathers, scatter-add adjoint). Implementation plan: "Follow-up
plan: Fourier transform as a truncated DFT on Reactant" at the end (deferred to a follow-up PR); independent of it: merging
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
- Fourier stage as a truncated DFT on Reactant (follow-up PR): see the next section.

## Follow-up plan: Fourier transform as a truncated DFT on Reactant

> Status: **in progress** on branch `mg/reactant-dft` (2026-10-18). Implemented for Reactant only and verified
> against FFTW on Reactant's CPU backend, see "Protocol (implementation)" at the end; GPU checks and benchmarks open.

Date: 2026-10-18. Base revision: `e22793e3` (branch `mg/new-reactant`). Prototype:
`/p/tmp/maxgelbr/code/Speedy2/claude_reactant_gpu_2026-10-13/dft_proto_2026-10-16/dft_transform_proto.jl`.

### Problem

On Reactant's GPU backend the FFT-based Fourier stage is why the `SpectralTransform` loses to the matrix transform
at T ≤ 127 (H100, PrimitiveWet L8, 10 steps per call: T31 1.12 vs 0.451 ms, T63 2.25 vs 0.687 ms per step):
- every `stablehlo.fft` becomes an XLA `FftThunk`, which can't be recorded into a CUDA graph (below), so the step
  splits into several graphs plus ~270 separately launched kernels per step (T31, A40);
- XLA normalizes every inverse FFT with a cuBLAS `scal` that `brfft` multiplies back. The per-run slices, zero
  padding, concatenations and the scratch round trip with `permutedims` add ~380 µs. Of the 835 µs GPU time per
  step at T31 (A40), 91 µs are `scal` and 291 µs concatenate/slice/pad and transposes.

Only the first `mmax` of a ring's `nlon/2 + 1` frequencies are used. So the Fourier transform of a ring is a
truncated DFT, a small dense matrix. All rings together are one batched `dot_general`, which XLA records into its
CUDA graph like any matrix product.

**Why not FFTs inside CUDA graphs** (checked 2026-10-18). In the depot's `libReactantExtra.so` (Reactant_jll
build of 2026-10-17), the thunks that XLA can record into a command buffer (CUDA graph) have a `Record` method:
kernel, custom kernel, custom call, cuDNN, while, conditional, memzero, device-to-device copy, collectives.
`FftThunk` has none. Upstream XLA `main` is the same: `FftThunk` derives only from `Thunk`, and
`GetCommandBufferCmdType` (`command_buffer_conversion_pass.cc`) has no case for `kFft`. So there is no `FFT`
command type to enable with `XLA_FLAGS`. Two ways would be possible, both out of scope:
- an XLA patch making `FftThunk` a `TracedCommand` with stream-captured cuFFT calls (plans created before capture);
- a command-buffer-compatible cuFFT FFI custom call in Reactant.

Either removes only the launch overhead. The FFT path's GPU time would stay: 835 µs per step at T31, about twice
the matrix transform's whole step (418 µs, A40).

**Why the DFT is cheaper than the matrix transform, although that is a single product.** A spherical harmonic is
separable, Y_lm(λ, φ) = P_lm(sin φ) e^{imλ}. The `MatrixSpectralTransform` stores the product for every harmonic
and every grid point, `(2nharmonics, npoints)`, O(N⁴) in work and memory. Applied in two stages (a DFT per ring,
then a Legendre transform per order) the factors are O(N³). The matrix transform is memory-bandwidth bound: it
reads the whole matrix on every call, whatever K. The DFT matrices fit, or nearly fit, into L2.

| T | matrix, per direction | DFT matrix `W`, per direction | FLOPs per column: matrix / DFT + Legendre |
| --- | --- | --- | --- |
| T31 | 14 MB | 0.7 MB | 7.1 M / ≈ 0.8 M |
| T63 | 188 MB | 5 MB | 94 M / ≈ 6 M |
| T127 | 2.7 GB | 39 MB | 1.4 G / ≈ 45 M |
| T255 | ≈ 41 GB | 307 MB | 20 G / ≈ 0.36 G |

### Design

**Selection.** A new `Bool` field `fourier_dft`, a keyword of the `SpectralTransform` constructor like `gpu_graphs`.
It adds no type parameter (the type name is at 1008-1014 of 1024 characters, `test/type_name_length.jl`). Defaults
from `default_fourier_dft(architecture, truncation)`:
- `false` on CPU/GPU: the CUDA.jl path is unchanged.
- on `ReactantDevice`, `true` up to `dft_max_truncation` (proposal T255, decided by the benchmark).

The FFT path stays for larger truncations, memory-limited runs and possibly Reactant's CPU backend. `if
S.fourier_dft` is resolved while tracing.

**Precomputed matrices.** `dense_fourier_matrices(architecture, ...)` has an empty fallback in `src/` and its
method in `SpeedyTransformsReactantExt`. The matrices are computed on the host in Float64, converted to `NF` and
stored flattened as `VectorType` like `legendre_polynomials_dense`. They are empty on CPU/GPU and with
`fourier_dft = false`.
- Forward `W[2mmax, nlon_max, nlat_half]`: row m + 1 holds the real part and row mmax + m + 1 the imaginary part
  of `solid_angles[j] * conj(lon_offsets[m + 1, j]) * cis(-2π m (i - 1) / nlon_j)`. It is zero for longitudes
  i > nlon_j (padding) and for orders beyond the ring's Nyquist frequency (2m > nlon_j, as the zero-padded `rfft`
  is now). This folds in the quadrature weights and longitude offsets that the forward `_legendre!` applies today.
- Inverse `B[nlon_max, 2mmax, nlat_half]`: `brfft` of the first mmax frequencies with the longitude offset folded
  in. Columns are `c_m cos(θ + m λ_j)` (real part) and `-c_m sin(θ + m λ_j)` (imaginary part), where:
  - θ = 2π m (i - 1) / nlon_j;
  - λ_j is the first longitude of ring j;
  - c_m = 1 for m = 0 and the Nyquist frequency, 2 otherwise;
  - entries beyond the ring are zero.

  The prototype multiplied `lon_offsets` separately as complex numbers. Folding them in makes the whole traced
  transform real.
- Memory 2 × 0.7 / 5 / 39 / 307 MB at T31 / T63 / T127 / T255, ≈ 2 × 2.4 GB at T511. Reported by `show`.

**Index tables.** Built on the host while tracing (like `dense_indices`), constants in the program:
- `ring_gather_indices(S)`, `(nlon_max, nlat_half, 2)`: the grid point at longitude i of northern/southern ring j.
  Padding points to npoints + 1, an appended zero row. For odd nlat the southern copy of the equator ring does too
  (it is transformed once, as now).
- `ring_positions(S)`, `(npoints,)`: the position of every grid point in the inverse product's output
  `(nlat_half, nlon_max, 2)`; `dot_general` puts the batch dimension first. Equator points come from the northern
  slot only.

**Traced transform.** New methods of `_transform_nonchunked!` for `ReactantSpectralTransform`, both directions
(`_transform_grid!`/`_transform_spec!` already wrap them in `@maybe_jit`). They branch on `S.fourier_dft`. The FFT
branch calls `_fourier!`/`_legendre!` as now. The DFT branch fuses the Fourier and Legendre transforms and skips
the scratch memory.

Grid → spectral:
1. Read the field as an (npoints, K) matrix and append a zero row.
2. Gather into padded rings `(nlon_max, nlat_half, 2, K)`.
3. `dot_general` with `W`, batched over j, contracting the longitude: `(nlat_half, 2mmax, 2, K)`.
4. Even/odd = north ± south.
5. Legendre `dot_general` with the dense polynomials, batched over parity and order, contracting j. The real and
   imaginary parts are a free dimension of the right-hand side, so neither `dot_general_complex` nor complex arrays
   are needed.
6. `complex(re, im)`, pack with `dense_of_lm`, `write_matrix!`.

Spectral → grid:
1. Unpack to the dense layout with `lm_of_dense`, real and imaginary parts as a free dimension.
2. Legendre `dot_general`, contracting the degree.
3. North/south = even ± odd, times `coslat⁻¹[j]` with `unscale_coslat`.
4. `dot_general` with `B`, batched over j, contracting the order and the real/imaginary part.
5. Gather with `ring_positions`, `write_matrix!`.

Per direction that is one gather, two `dot_general`s and a few fused elementwise operations. There are no FFTs, no
`scal`, no slices or concatenations per latitude run and no loop-carried scratch arrays. Layouts are chosen so one
product's output feeds the next without transposes; check the optimized HLO for leftover `transpose`s. Arrays are
real throughout, except the packed coefficients at the boundary. That also avoids the complex-number Reactant
issues (`reactant-issue-draft-complex-*.md`).

**Scratch memory.** Unused with `fourier_dft = true`. Allocate it with zero layers in that case (zero-length arrays
are fine in traced loops since Reactant 0.2.294). Do the same for SpeedyWeather's `vars.scratch.transform_memory`
if it is sized from the transform (check).

**Default transform on Reactant.** If the full time step with the DFT path beats the matrix transform at all
resolutions (prototype, transform alone: 1.8x at T31), `WhichTransform` returns the `SpectralTransform` for every
`ReactantDevice` model. That closes open decision 1.

### Implementation steps

0. **Go/no-go: AD through the gathers.** Differentiate the prototype's `dft_forward`/`dft_inverse` on Reactant's
   CPU backend and compare with central finite differences in Float64. These are constant-index gathers with
   scatter-add adjoints. The strided-slice adjoint bug that forced the gather-free FFT path is fixed in 0.2.294,
   the current compat. If AD fails, minimize and report it, then fall back to padded rings from contiguous slices
   per ring (slower, as in the FFT path).
1. `dense_fourier_matrices` and the index tables, with host-only unit tests:
   - `W` times a ring equals the weighted, rotated `rfft` (FFTW) of that ring for its first mmax frequencies;
   - `B` matches `brfft` in the same way.
2. The `fourier_dft` field and keyword, `default_fourier_dft`, `show`.
3. The DFT branch of `_transform_nonchunked!`, both directions, all real.
4. Odd nlat (Clenshaw grids) and the HEALPix/OctaHEALPix ring offsets. These live only in the tables and matrices
   and are covered by the tests of step 5.
5. Tests (below) on Reactant's CPU backend.
6. On the GPU node: GPU checks, Nsight, benchmarks (below).
7. Defaults from the benchmarks: `dft_max_truncation`, the CPU-backend rule, `WhichTransform`.
8. Docs: the Reactant section of `architectures_gpu.md`, i.e. how the Fourier transform is computed and when.
   CHANGELOG line. The SpeedyTransforms version stays `1.0.0+DEV` (additive change).

### Testing and verification

- `test/reactant/test_transform.jl`: the DFT path vs the CPU `SpectralTransform`.
  - Grids: OctahedralGaussian, FullGaussian, HEALPix, OctaHEALPix, FullClenshaw and OctahedralClenshaw (the last
    two with odd nlat).
  - No layer dimension, one and three layers; eager and compiled; with and without `unscale_coslat`.
  - Also the DFT vs the FFT path on Reactant in Float64 (≤ 1e-12).
- `test/reactant/test_correctness.jl`: ShallowWater and PrimitiveWet with the DFT path in CI; Barotropic and
  PrimitiveDry once by hand.
- `test/reactant/differentation.jl`: `test_differentiability(BarotropicModel; Transform = SpectralTransform)` with
  the DFT path. Also finite-difference checks of both directions in Float64, by hand as for the FFT path.
- SpeedyTransforms CPU suite (type name length, JET, AD rules): CPU must be unaffected.
- GPU (H100/A40):
  - `test/GPU/reactant.jl` NaN check.
  - Errors vs the CPU transform with TF32 and `HIGHEST`; expected as in the prototype, 4-7e-4 and ~1e-7.
  - Nsight with `--cuda-graph-trace=node`: one graph launch per compiled call, no cuFFT kernels, no `scal`.
  - `bench_transforms.jl`: matrix vs FFT path vs DFT path, T31-T255, TF32 and `HIGHEST`, single steps and
    `time_steps!(simulation, 10)`. Then the main suite, `manual_benchmarking.jl reactant-gpu`.
- Reactant CPU backend: the same benchmark, to decide the default there.

### Expected results (estimates, to be measured)

- Transform alone: as in the prototype. 10x (T31) to 1.8x (T255) faster than the FFT path, 2-18x faster than the
  matrix transform.
- Full step (H100, 10 steps per call):
  - T31: about the matrix transform's 0.45 ms/step (both run as one CUDA graph, and the transform is a small part).
  - From T63: faster than the matrix transform, whose products take 983 of 1386 µs GPU time per step at T63 (A40).
  - T127/T255: faster than the FFT path's 5.67 / 13.0 ms, by less than the transform-alone ratio, because the
    Legendre products and the rest of the step are unchanged.

### Risks and open questions

- AD of the gathers (step 0).
- TF32: the DFT adds one TF32 product per direction. The prototype's errors (4-7e-4) are at the level of the
  current path with TF32. `HIGHEST` cost 10-15% in the prototype.
- Memory is O(N³), and every call reads all of `W`/`B` whatever K. A K = 2 call at T255 reads 307 MB (~0.1 ms on
  an H100). See `2026-10-16-fuse-inverse-transforms.md` for merging small calls.
- Padding: on octahedral grids rings average about half of `nlon_max`, so about half of `W`/`B` is zeros. Optional:
  split the latitudes into 2-4 bands, each with its own `nlon_max`. That means less memory and work for a few more
  operations.
- Full grids (all rings the same length) could share one `W`/`B` across latitudes (O(N²) memory) if the weights and
  offsets are applied separately. Optional.
- Reactant's CPU backend has no tensor cores, so the O(N³) DFT may be slower than XLA's FFTs (DUCC) from about
  T127. The default there may differ; check how to query the backend when the transform is constructed.
- Compile time: fewer operations than the FFT path (which has one per latitude run), but constant index arrays of
  npoints entries (1.2 MB at T255). Measure.

### Out of scope

- The DFT Fourier stage for CUDA.jl. Rough estimate (H100, from the benchmark README and
  `docs/dev/2026-08/gpu-idle-time-reduction.md`):
  - T ≤ 64: no gain. A CUDA.jl step is bound by host dispatch at ~1.1 ms whatever the transform.
  - T86: ~1.5x against the bandwidth-bound matrix transform (585 MB per direction).
  - Above T100: ~1.1x against LT+FFT. Its FFTs already run in CUDA graphs: cuFFT and Legendre are ~2 ms each of
    10.5 ms GPU time per step at T255, and Legendre would be unchanged.

  CUDA.jl's default math mode probably runs Float32 gemms without TF32 (to check). This needs a separate analysis
  and benchmark.
- FFT command-buffer support in XLA or Reactant (see above). An upstream issue could be drafted.

### Protocol (implementation)

Branch `mg/reactant-dft` (base `c32de0ee`), login02 (A40) and Slurm (H100). Scripts and logs in
`/p/tmp/maxgelbr/code/Speedy2/claude_reactant_gpu_2026-10-13/dft_impl_2026-10-18/`. Environments: `env_reactant_0295`
(Reactant 0.2.295, Julia 1.12.2, Reactant's CPU backend) and `env_dft_gpu` (the same, GPU, see below).

2026-10-18
- **Step 0, go: AD through the gathers** (`ad_gather_check.jl`, Reactant CPU backend, Float64, quadratic losses):
  forward gather into padded rings + batched `dot_general`, and inverse `dot_general` + gather of the grid points, on
  OctahedralGaussian, HEALPix and FullClenshaw (odd nlat) grids at T15. Values vs host loops ≤ 1.4e-14, Reactant
  gradients vs central finite differences ≤ 4e-11 relative, the forward gradient vs the analytic one 3e-16.
- **Implementation** (Reactant only):
  - `src/` (SpeedyTransforms): fields `fourier_matrix_forward`, `fourier_matrix_inverse` (flattened `VectorType`,
    empty elsewhere) and `fourier_dft::Bool` (no new type parameter), keyword `fourier_dft` with
    `default_fourier_dft(architecture, spectrum)`: `false` on CPU/GPU, on `ReactantDevice` `true` up to
    `DFT_MAX_TRUNCATION = 256` (1-based `mmax`, 2 × 308 MB in Float32 at 256), `ArgumentError` for `true` elsewhere;
    the hook `dense_fourier_matrices` (empty fallback in `fourier.jl`); `solid_angles` computed before it;
    `_scratch_nlayers(::ReactantDevice, ...) = 0` (no scratch memory on Reactant, see the refactor); `show` adds a
    "Fourier: truncated DFT, matrices ..." line only for `fourier_dft = true` (CPU/GPU output unchanged).
  - `SpeedyTransformsReactantExt`: `fourier_matrices` (host, Float64: phases `cispi(2 mod(m(i-1), nlon)/nlon)`
    with the argument reduced exactly, solid angle and conjugate longitude offset folded into the forward, the
    longitude offset and the `brfft` factors 1/2/1 into the inverse, zero beyond a ring's Nyquist frequency and
    beyond the ring), `ring_indices` (gather into padded rings with a zero row for the padding and the southern
    copy of the equator ring; positions of the grid points in the inverse result), `dft_fourier_forward`,
    `dft_fourier_inverse`, the all-real Legendre transform (`dense_of_lm_order_first` maps the forward result
    without a transpose).
- **Found and fixed in the FFT path** (pre-existing): the inverse `_fourier!` failed on grids whose equator ring is
  a run on its own (OctaHEALPix, OctahedralClenshaw: dropping its southern copy left an empty slice, which
  Reactant's indexing rejects, "BoundsError: attempt to access 0-element Vector{Union{}}"). The tests only covered
  OctahedralGaussian, FullGaussian and HEALPix. Fixed (the run is dropped), both grids are tested now.
- **Verification against FFTW** (`dft_verify.jl`, all six grids), relative max errors:
  - A, host matrices vs FFTW `rfft`/`brfft` ring by ring (with the weights/rotations of the FFT path), Float64:
    ≤ 1.0e-15 (T32), ≤ 1.5e-15 (T64).
  - B, traced Fourier stages vs FFTW ring by ring: Float64 ≤ 1e-15, Float32 ≤ 4.8e-7 (T32), ≤ 6.3e-7 (T64).
  - C, whole transform, DFT vs CPU FFTW and vs the Reactant FFT path: Float64 ≤ 1.0e-15 both directions; Float32
    error vs the Float64 FFTW transform for DFT / Reactant FFT / FFTW: grid→spectral 0.6-1.6e-7 / 0.6-1.3e-7 /
    0.5-2.2e-7, spectral→grid 1.3-3.7e-7 / 1.3-4.2e-7 / 1.2-5.3e-7 (T32, T64): the same order, the DFT's error at
    most 1.9x the smaller of the two FFT errors.
  - The same on the A40 with `dot_general_precision = HIGHEST`: Float64 ≤ 1e-15, Float32 stages ≤ 6.3e-7, Float32
    transform errors ≤ 6.6e-7, at most 2.4x the smaller of the two FFT errors.
- **Tests**: `test/reactant/test_transform.jl`: the existing testset for both Fourier transforms and five grids
  (OctaHEALPix and OctahedralClenshaw new), new testset "Truncated DFT against FFTW" (defaults; the matrices vs
  the CPU SpectralTransform's Fourier transform (FFTW) ring by ring on six grids, Float64 rtol 1e-13; the traced
  Fourier stages (DFT and FFTs) vs the same, Float64 1e-13, Float32 2e-6; the whole transform, Float64 DFT vs the
  CPU SpectralTransform and vs the FFT path 1e-13, Float32 error vs the Float64 CPU SpectralTransform
  transform < 1e-6 and < 4x that of FFTW or the FFT path). Reactant CPU backend: 2 + 79 + 79 + 1461 pass (8.6 min,
  of which the new testset 3.1 min). SpeedyTransforms CPU suite with the new testset "Fourier transform as truncated
  DFT only on Reactant": all pass (type name length, JET, AD rules 8175/8175).
- **GPU environment**: the shared depot's precompiled `Reactant_jll` (and `CUDA_Runtime_jll`) caches had been built
  where no GPU driver was found (`gpu = "none"`, "No functional client was found" on the GPU backend, CUDA.jl "could
  not find an appropriate CUDA runtime"). Not rebuilt (shared with other jobs); instead `env_dft_gpu`, a copy of
  `env_reactant_0295` with `LocalPreferences.toml`: `Reactant_jll` `gpu = "cuda"`, `gpu_version = "13.1"`,
  `mode = "opt"` (`gpu_version` isn't a recorded compile-time preference, `mode` forces the separate cache;
  `Reactant_jll` under `[extras]` so that its preferences are found), `CUDA_Runtime_jll` `version = "13.0"`
  (`CUDA.set_runtime_version!`).
- **Refactor** (user: one Legendre transform): `legendre_forward`/`legendre_inverse` (all real) for both Fourier
  transforms; `fft_fourier_forward`/`fft_fourier_inverse` with the interface of the DFT (weights and rotations in
  the Fourier stage); `_transform_reactant!` for both directions; `_fourier!`, `_legendre!` and
  `dot_general_complex` of the extension removed; no scratch memory on Reactant at all. The traced-stage testset
  covers both Fourier transforms.
- **After the refactor**, Reactant CPU backend: `test_transform.jl` 2 + 79 + 79 + 1479 pass (7.4 min); `dft_verify.jl`
  T32/T64, six grids: Float64 DFT vs FFTW and vs the FFT path ≤ 1.0e-15, Float32 errors ≤ 3.7e-7 (DFT), ≤ 4.2e-7
  (FFT path).
- **Model tests**, Reactant CPU backend, CPU-vs-Reactant correctness (`test_correctness.jl`'s `test_model`, the CPU
  model with FFTW) and differentiability (`test_differentiability`, gradients vs Enzyme on the CPU):
  - before the refactor, DFT: Barotropic 8/8, ShallowWater 14/14, PrimitiveDry 17/17, PrimitiveWet 20/20, AD
    Barotropic 6/6;
  - after the refactor: PrimitiveWet 20/20 with the DFT and 20/20 with the FFTs (`fourier_dft = false`), AD
    Barotropic 6/6 with the DFT and 6/6 with the FFTs.
- **A40** (`env_dft_gpu`): `test_transform.jl` with `HIGHEST` (refactored code) 2 + 79 + 79 + 1479 pass (8.3 min).
  With the default TF32: Float64 unaffected (≤ 1e-15), the Float32 Fourier stages 2-4.4e-4 (the TF32 level), the whole
  transform 4.1-7.8e-4 (DFT) vs 2.6-6.9e-4 (FFT path, its TF32 product is only the Legendre transform), i.e. up to
  ~1.8x with one more TF32 product per direction. PrimitiveWet NaN check with the DFT (default resolution, 10 steps,
  `run!`): 6/6, 12.7 min (incl. compilation).
- **`show`** of the transform reported the size of the Reactant array wrappers (`sizeof` of a `ConcretePJRTArray`),
  e.g. "matrices 64 bytes", also for the dense Legendre polynomials before: now `nbytes` = length × size of the element
  type (unchanged for CPU/GPU arrays). T32 on Reactant: DFT matrices 1.38 MB, polynomials 105.78 KB, scratch 0 bytes.

2026-10-18 evening
- **H100 benchmark** (Slurm 2566174, `bench_dft.jl`, PrimitiveWet L8, OctahedralGaussian, TF32, ms per step with 10
  steps per call; matrix / FFT path / DFT): T32 0.391 / 0.852 / 0.284, T64 0.527 / 1.76 / 0.401, T128 2.45 / 4.07 /
  0.745. Transform alone (forward + inverse, 73 columns): T32 0.0406 / 0.359 / 0.0312, T64 0.218 / 0.809 / 0.066,
  T128 1.84 / 2.10 / 0.233 ms. So the DFT path is the fastest at every resolution so far (T256, `HIGHEST` and the CPU
  backend still running).
- **Audit of `SpeedyTransformsReactantExt`** (user: minimal duplication, little change in the main code), 519 → 499
  lines:
  - one dense layout (parity, order, degree) for the coefficients and the polynomials, so one `dense_indices` for
    both maps; it builds the dense polynomials too. Before: three copies of the same (l, m) loop and two layouts, the
    forward result was (parity, order, degree) by `dot_general`'s batch-first output, the inverse and the
    polynomials (parity, degree, order), hence the extra `dense_of_lm_order_first`;
  - `stack_pair(a, b, dims)` for the five `cat(reshape(a, ..., 1, ...), reshape(b, ..., 1, ...); dims)` (real and
    imaginary parts, north and south, even and odd);
  - `S.nlon_max` and `S.nfreq_max` instead of recomputing them from `nlons` (`fourier_matrices` keeps
    `maximum(nlons)`, a host function of `nlons` that the tests call);
  - the unused `_jit` import removed (`@maybe_jit` resolves it in its own module);
  - kept: the duplicated boundschecks of the two `_transform_reactant!` methods (as in the main
    `_transform_nonchunked!`), the empty-vector branch of `dense_fourier_matrices` (one line);
  - main code: only the layout in the comment of `legendre_polynomials_dense`.
  
  Transform alone before / after the audit (A40, TF32, `bench_transform_alone.jl`, DFT and FFT path, ms): T31 0.0416
  / 0.0416 and 0.495 / 0.491, T63 0.137 / 0.136 and 1.24 / 1.24, T127 0.766 / 0.766 and 2.95 / 2.95: unchanged.
- **Memory and FLOPs** of the precomputed arrays (`memory_flops_table.jl`, formulas from the array sizes, Float32,
  OctahedralGaussian; T32 checked against `show`): matrix transform 28 MB / 375 MB / 5.4 GB / 82 GB at T32 / T64 /
  T128 / T256; DFT matrices 1.4 / 10 / 79 / 617 MB (54 to 44% nonzero, the rest padding), dense Legendre 0.10 / 0.81
  / 6.4 / 51 MB. FLOPs per column and direction, matrix vs DFT + Legendre: 7.1 vs 0.79, 94 vs 5.9, 1352 vs 46,
  20470 vs 359 MFLOP (9 to 57× fewer).
- **Tests after the audit**, Reactant CPU backend (`env_reactant_0295`): `test_transform.jl` 2 + 89 + 89 + 1479 pass
  (7.6 min; 89 instead of 79 since the `show` tests); PrimitiveWet CPU vs Reactant 20/20 with the DFT and 20/20 with
  the FFT path; Barotropic gradients vs Enzyme on the CPU 6/6 and 6/6.
- **Nsight, DFT path** (A40, `dft_impl_2026-10-18/profile/`, `run_nsys.sh`; matrix and FFT references from
  `profile_fft_2026-10-16/`): the optimized StableHLO of the 10-step loop has no `fft` (FFT path: 96), 29
  `dot_general` (matrix 27), 33 transposes (matrix 13, FFT path 123), 6380 lines (matrix 6171, FFT 7686). Per time
  step (20 calls of 10 steps): 1.2 graph launches and 0.5 kernel launches outside graphs (matrix 1.2 / 0.5, FFT path
  4.2 / 270), 107.5 kernels, no cuFFT, no `scal`; GPU time 304 µs at T31 (matrix 350, FFT 835) and 513 µs at T63
  (matrix 1386, FFT 1833); wall time 0.367 / 0.585 ms per step at T31 / T63 (matrix 0.418 / 1.482, FFT 1.466 /
  3.139). Transform alone (73 columns, traced loop): GPU time 38 / 134 µs (matrix 78 / 791, FFT 282 / 834), 11
  kernels per forward + inverse, one graph launch.
- **Precision** (`precision_table.jl`, Slurm 2567378, H100, OctahedralGaussian, 8 columns of random coefficients and
  the band-limited fields from them, relative max error vs a Float64 CPU reference): CPU Float32 baseline 1.4 to
  5.8e-7; DFT with Float32 (`HIGHEST`) 1.8 to 6.7e-7 grid → spectral (0.8 to 1.3× the baseline) and 2.9 to 8.7e-7
  spectral → grid (1.7 to 3.3×), the FFT path 0.7 to 1.6× the baseline, the matrix transform 4.3e-7 to 2.2e-6 (up to
  9×, growing with resolution); with TF32 all at 2 to 4.4e-4, the DFT 1.0 to 1.5× the FFT path. The T128 matrix TF32
  line was lost in the job's filtered log (a ptxas warning printed onto it, deterministic: also lost in a filtered
  rerun); unfiltered rerun on the A40: 2.0e-4 grid → spectral, but 1.6e-6 spectral → grid (XLA apparently chose a
  GEMM without TF32 for that shape on the A40).
- **`HIGHEST` time step** (Slurm 2566174, H100, ms per step, 10 per call; matrix / FFT path / DFT): T32 0.415 /
  0.823 / 0.280, T64 0.763 / 1.76 / 0.424, T128 5.67 / 4.07 / 0.862. Float32 instead of TF32: DFT +0 to 16%, matrix
  +6% / +45% / 2.3×, FFT path unchanged.
- **GPU NaN check** `test/GPU/reactant.jl`, BarotropicModel now with both transforms: 1/1 and 1/1 (A40).
- **Long-run script** (`long_run_stability.jl`), smoke test T31L8 for 1 year on the A40: stable, 2975 SYPD, global
  mean surface temperature 285 → 275 K in the spin-up, then 274-278 K, surface pressure 962.6 to 963.6 hPa, humidity
  ≥ 0, |u| ≤ 87 m/s. The 10-year T256L8 TF32 run is Slurm 2568499 (queued).
- Pitfall: waiting with `while pgrep -f "julia.*precision_table.jl"` matches the waiting shell's own command line and
  never ends (two A40 runs didn't start for ~40 min).
- **`HIGHEST` T256** (Slurm 2566174): FFT path 9.8, DFT 3.1 ms per step (TF32: 9.68 / 2.29), i.e. Float32 costs the
  DFT path 35% at T256 (0 to 16% up to T128).
- **Reactant CPU backend** (16 cores of the H100 node, same job), FFT path / DFT, ms per step with 10 per call
  (transform alone): T32 4.14 / 4.82 (1.84 / 1.58), T64 12.7 / 12.7 (5.84 / 5.06), T128 44.9 / 43.8 (21.1 / 20.7).
  Within 15%, no backend-specific rule: `default_fourier_dft` stays `true` on every Reactant backend (querying the
  backend at construction isn't worth ±15%).
- **Step 7, defaults**: `WhichTransform` returns the `SpectralTransform` for `ReactantDevice` (faster than the matrix
  transform at every resolution: H100 TF32 1.3 to 3.3×, `HIGHEST` 1.5 to 6.6×, Reactant CPU backend 2.1 / 5.3× at
  T31 / T63 for the FFT path, M3); `DFT_MAX_TRUNCATION = 256` kept (617 MB of DFT matrices there, the H100 runs
  it 4.2× faster than the FFT path). New testset "Default transform on Reactant" in `test_transform.jl`: 4/4 (CPU
  backend). Open decision 1 closed.
- **Step 8, docs**: `architectures_gpu.md` Reactant section rewritten (default transform, truncated DFT and FFT path,
  memory, H100 performance, CPU backend; the TF32 model table marked as measured with the matrix transform, the cost
  of `HIGHEST` for the default transform; the GPU NaN test is the BarotropicModel with both transforms);
  `differentiability.md` (default transform in the example). The PrimitiveWet gradient in CI uses the matrix
  transform; a one-off check of the one-step PrimitiveWet gradient with the DFT is running.
