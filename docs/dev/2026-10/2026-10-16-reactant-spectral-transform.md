# FFT-based `SpectralTransform` with Reactant

> Status: **in progress**. Investigating what breaks when the FFT + Legendre `SpectralTransform` is used on `ReactantDevice`.

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
  `SpectralTransform{NF, <:ReactantDevice}` (see the prototype above), `_transform_grid!`/`_transform_spec!`
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
