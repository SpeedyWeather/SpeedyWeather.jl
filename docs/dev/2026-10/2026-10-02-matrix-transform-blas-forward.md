# MatrixSpectralTransform: real stacked matrices so both multiplies hit BLAS on every backend

> Status: **completed**. Implemented on `mg/matrix-transform-regression` (based on `origin/main`),
> all CPU and GPU tests pass, before/after transform benchmarks on CPU (EPYC 9554, 16 cores) and GPU
> (H100) below: the grid → spectral multiply is 84–108× faster on CPU and 2.6–13× faster on GPU,
> the `PrimitiveWetModel` with the matrix transform runs 15× (T32) to 41× (T86) faster on CPU and
> 2.4× (T64) to 3.5× (T86) faster on GPU, where it now beats the FFT transform up to at least T86.

Date of initial draft: 2026-10-02

Base revision: `8e99099f` (`origin/main`, after rebase; drafted against `8a29efdd` v0.23.0)

## Originating prompt

> Ok in a new branch `mg/matrix-transform-regression` we need to investigate problem 2. First
> investigate the change we have done in PR #1201, then we need to find a version that restores old
> CPU performance (or improves it) but also keeps GPU compatability and keeps or improves GPU
> performance

## Revision log

- 2026-10-02: initial draft.
- 2026-10-02 (review 2): on request the FFT transform's scratch/`nlayers` uses the same expression
  `max(maximum(transform_batch), max_transform_batch(nlayers))` as the matrix transform (one line in
  `spectral_grid.jl`); on CPU this only raises the `ismatching` limit (the FFT scratch is the largest
  planned K, wider calls are chunked), on GPU it closes a corner case where a custom `transform_batch`
  without the 9L+1 batch left the scratch too narrow for the tendency batch (the 4L+1 floor dated from
  #1099, before the tendency batch existed).
- 2026-10-02 (review): column blocking dropped on request, the matrix transform's scratch is sized
  to the widest model batch instead (`SpectralGrid` constructor, `max_transform_batch`).
- 2026-10-02: branch rebased onto `origin/main` (`8e99099f`) on request; the release-branch base only
  differed in versions/changelog/benchmark results. Implemented; SpeedyTransforms unit tests pass
  (209 tests, `--check-bounds=yes`). Validation extended to a before/after transform benchmark
  across resolutions on CPU and GPU (`claude_bisect/bench_matrix_transform.jl`), results below.

## Problem description

The [benchmark regression analysis](2026-10-02-x86-benchmark-regression-analysis.md) found that the
`MatrixSpectralTransform` rows of the x86 benchmark halved (T32 L8: 107 → 56 SYPD) with PR #1201
(`b1ad362b`, 2026-08-20). The grid→spectral `transform!` multiplies the **complex** forward matrix
`F` (nharmonics × npoints) with the **real** field. `LinearAlgebra.mul!` cannot hand a complex × real
product to BLAS `gemm!` (equal element types required), so it uses Julia's generic triple loop. On
Julia 1.12 the dedicated complex × real path in `LinearAlgebra` (which reinterprets the complex
matrices as real) is unreachable: its `generic_matmatmul_wrapper!` methods take `::Val{true}` while
`_mul!` now passes `Val(BlasFlag.GEMM)`. #1201 replaced

```julia
coeffs_matrix = reshape(coeffs.data, size(coeffs.data, 1), :)
```

by `_as_matrix`, which on CPU returns a 2-D array unchanged. The model's fused variables are 2-D
`SubArray`s of a 3-D parent, and the generic loop's element-wise `setindex!` is twice as slow through
the bare view (123 ms) as through the `ReshapedArray` the old code wrapped it in (62 ms, T32 L8
tendency batch, see `mt_micro.jl`). A real `sgemm` of the same size takes 1 ms.

On GPU the same multiply is in the same situation: CUDA.jl calls CUBLAS only when `eltype(A) ==
eltype(B) == eltype(C)` and otherwise falls back to `GPUArrays.generic_matmatmul!` (a naive kernel).
That is why T128 L8 runs at 33 SYPD with the matrix transform against 261 with the FFT transform on
the H100, although the matrix transform is the default for `truncation <= 64` on GPU.

## Background

- `#1201` also introduced the `_as_matrix(x, ::GPU)` materialization of non-contiguous views
  (nested `SubArray`s of a `CuArray` are not `StridedCuArray`s and crash CUBLAS). That part must
  stay. Contiguous views of GPU arrays are plain `CuArray`/`ROCArray`s, so the common fused-variable
  case never copies.
- The spectral→grid transform already avoids the complex × real problem by splitting `B` into
  `backward_real`/`backward_imag` and doing two real `gemm!`s plus two broadcasts into a real
  scratch `(nharmonics × K)`.
- `reinterpret(NF, complex_matrix)` would give the forward product in one real `gemm!` with no
  scratch (measured 0.98 ms vs 122 ms on CPU), but only CUDA.jl defines `reinterpret` for its array
  type; AMDGPU.jl, Metal.jl and GPUArrays do not, so a reinterpreted `ROCArray` would not reach
  rocBLAS. A backend-agnostic formulation is preferred.
- The transform's `scratch_memory` is shared with the model as `vars.scratch.transform_memory`; it
  is sized by `SpectralGrid` to `max(maximum(transform_batch), 4 nlayers + 1)` columns, which on CPU
  (`transform_batch = [1, nlayers]`) is smaller than the grid→spectral tendency batch (`9 nlayers +
  1` for `PrimitiveWet`), so a forward path through scratch must block over columns.

## Summary of changes

All in `SpeedyTransforms/src/matrix_transform.jl` (plus `show.jl`, tests, docs, changelog):

1. **Real stacked matrices instead of complex ones.** Store
   `forward_stacked = [Re F; Im F]` (`2 nharmonics × npoints`, real) and
   `backward_stacked = [Re B  -Im B]` (`npoints × 2 nharmonics`, real). The complex `forward`,
   `backward` and the `backward_real`/`backward_imag` copies are dropped, which reduces the matrix
   memory from three to two complex-equivalents (T128 L8: 8.2 GB → ≈5.5 GB). The
   `MatrixComplexType` type parameter goes away.
2. **Scratch** becomes `(2 nharmonics × nlayers)` real, used by both directions. The matrix
   multiply has no batch restriction (unlike FFT plans), so the `SpectralGrid` constructor sizes the
   matrix transform's scratch to the widest batch any model emits (`max_transform_batch(nlayers) =
   9 nlayers + 1`, the `PrimitiveWet` tendency batch, a few MB) instead of chunking; the FFT
   transform keeps its 4L+1 scratch and chunking. A batch wider than the scratch throws a
   `DimensionMismatch` saying which `nlayers` to construct the transform with.
3. **Grid→spectral**: `mul!(scratch, forward_stacked, field_matrix)` (one real `gemm!` / CUBLAS /
   rocBLAS) followed by one broadcast `coeffs .= complex.(scratch[1:n, :], scratch[n+1:2n, :])`.
4. **Spectral→grid**: two broadcasts fill the scratch halves with `real`/`imag` of the
   coefficients and one `mul!(field_matrix, backward_stacked, scratch)` replaces the two `gemm!`s
   (measured 1.03 → 0.74 ms on CPU). `unscale_coslat` unchanged.
5. **`_as_matrix`** on CPU unchanged (BLAS and broadcasts handle strided views; no copies). On GPU
   materialize whenever `x` is not an `AbstractGPUArray` (also for 2-D nested views, which the
   `ndims == 2` shortcut used to skip), otherwise reshape for free.
6. `show` reports `sizeof(forward_stacked) + sizeof(backward_stacked)`.
7. Reactant keeps working through the same code: `@maybe_jit` wraps the `mul!` calls as before, and
   the broadcasts/views are of the kind the current spectral→grid path already executes eagerly.

Expected effect: CPU forward multiply ≈100× faster than today (and than July), backward ≈1.4×;
on GPU the forward multiply moves from the generic GPUArrays kernel to CUBLAS, so the default
low-resolution GPU configurations (`truncation <= 64`) and all `matrix` rows should speed up.

## Testing and verification

- `SpeedyTransforms/test/matrix_transform.jl`: adapt the field/size checks; add a short test that
  (a) a 9L+1-column batch is transformed when the transform is constructed for it and a too narrow
  scratch throws, (b) transforms into/out of views (2-D view of a 3-D parent, 1-D slot view) agree
  with the plain-array result, on CPU. JET `@test_opt` guard for both directions in `dispatch.jl`.
- `SpeedyWeather/test/dynamics/matrix_transform.jl` (20-day `PrimitiveWetModel` run, no NaN).
- GPU: `SpeedyWeather/test/GPU/runtests.jl` (CUDA env) on an H100 node, in particular the existing
  "transform! on views into a larger backing array" regression test and the model tests that use
  the matrix transform.
- Benchmarks: the reduced benchmark from the regression analysis (`claude_bisect/bench_subset.jl`,
  T32 L8 MT) on the same node as before (baseline 107, HEAD 57 SYPD), and the GPU suite's matrix
  rows (`manual_benchmarking.jl gpu`) against the committed README.

## Results

Before = `origin/main` (`8e99099f`), after = this branch. `claude_bisect/bench_matrix_transform.jl`:
`transform!` wall time (minimum of 10 calls, 3 at T128, device-synchronised) with the transform as the
model constructs it (`MatrixSpectralTransform(spectral_grid)`, 8 layers) for the batch widths the
`PrimitiveWetModel` uses (K = 8 layers, 33 = prognostic batch, 73 = tendency batch), plus the model's
SYPD with the matrix transform measured like the benchmark suite. CPU: `csp14c04` (AMD EPYC 9554,
`--cpus-per-task=16`, Julia 1.12.2); "BLAS 16" = `BLAS.set_num_threads(16)` instead of Julia's default
of 64 threads on the 16 allocated cores. GPU: NVIDIA H100 80GB (`csl14c246`), CUDA.jl 6.2.2.

### CPU, grid → spectral (ms)

| T | K | before | after | after, BLAS 16 |
| --- | --- | --- | --- | --- |
| 32 | 8 / 33 / 73 | 9.3 / 38.6 / 84.6 | 0.44 / 0.59 / 1.00 | 0.29 / 0.40 / 0.75 |
| 64 | 8 / 33 / 73 | 128 / 523 / 1162 | 6.6 / 7.8 / 12.5 | 3.7 / 4.2 / 7.9 |
| 86 | 8 / 33 / 73 | 393 / 1624 / 3600 | 22.2 / 27.0 / 39.1 | 12.3 / 14.4 / 24.3 |
| 128 | 8 / 33 / 73 | – / 7724 / 17038 | 79.9 / 94.0 / 157 | 52.1 / 57.0 / 104 |

### CPU, spectral → grid (ms; before could not take K = 73, scratch too narrow)

| T | K | before | after | after, BLAS 16 |
| --- | --- | --- | --- | --- |
| 32 | 8 / 33 / 73 | 0.42 / 0.55 / – | 0.40 / 0.54 / 0.95 | 0.25 / 0.34 / 0.67 |
| 64 | 8 / 33 / 73 | 6.8 / 7.9 / – | 6.5 / 7.7 / 12.3 | 3.7 / 4.5 / 7.7 |
| 86 | 8 / 33 / 73 | 17.0 / 20.3 / – | 17.1 / 20.3 / 34.0 | 10.4 / 11.9 / 21.5 |
| 128 | 8 / 33 / 73 | – / 93.7 / – | 80.5 / 94.3 / 158 | 52.1 / 57.8 / 104 |

### CPU, `PrimitiveWetModel` L8 with the matrix transform (SYPD)

| T | July README | before | after | after, BLAS 16 | FFT transform (Oct README) |
| --- | --- | --- | --- | --- | --- |
| 32 | 107 | 56.9 | 882 | 980 | 830 |
| 64 | 3.7 | 2.2 | 64.8 | 82.9 | 104 |
| 86 | 0.9 | 0.5 | 20.7 | 27.2 | 39 |

### GPU (H100), `transform!` (ms) and model SYPD

| T | K | grid → spectral before → after | spectral → grid before → after |
| --- | --- | --- | --- |
| 32 | 8 / 33 / 73 | 0.189 / 0.196 / 0.200 → 0.056 / 0.055 / 0.078 | 0.062 / 0.075 / 0.084 → 0.055 / 0.062 / 0.068 |
| 64 | 8 / 33 / 73 | 1.23 / 1.23 / 1.78 → 0.180 / 0.179 / 0.294 | 0.140 / 0.189 / 0.307 → 0.134 / 0.177 / 0.285 |
| 86 | 8 / 33 / 73 | 2.13 / 2.14 / 6.26 → 0.301 / 0.444 / 0.782 | 0.344 / 0.522 / 0.853 → 0.321 / 0.492 / 0.782 |
| 128 | 33 / 73 | 13.9 / 42.9 → 1.70 / 3.26 | 1.92 / 3.36 → 1.84 / 3.29 |

| T | before | after | FFT transform (Oct README) |
| --- | --- | --- | --- |
| 32 | 5653 | 5836 | 5621 |
| 64 | 1206 | 2836 | 1181 |
| 86 | 310 | 1095 | 653 |

The backward transform was already BLAS-bound and is unchanged (slightly faster: one `gemm!` instead of
two); the forward transform now costs the same as the backward one, as it should. On CPU the matrix
transform is now within 1.3× of the FFT transform up to T64 and beats it on GPU at T64 and T86, so the
`WhichTransform` threshold (`truncation <= 64` on GPU) could be re-evaluated with the full GPU suite.

### Tests

- `SpeedyTransforms/test/matrix_transform.jl` (209 tests incl. the new wide-batch/view tests) and
  `dispatch.jl` (32 incl. the new JET guard), `--check-bounds=yes`, against the local package
  (`SpeedyTransforms/test/Project.toml` now has `[sources]` for the monorepo packages; without them
  `--project=SpeedyTransforms/test` silently tested the registry release).
- `SpeedyWeather/test/dynamics/matrix_transform.jl` (20-day `PrimitiveWetModel`, no NaN) and a one-day
  run of all four models with the FFT transform after the scratch-width change.
- `SpeedyWeather/test/GPU/runtests.jl` on an H100 (CUDA env): all 36 testsets pass, including the
  matrix-transform round trip, the nested-view regression test from #1201, Barotropic with the matrix
  transform and the default GPU `PrimitiveWetModel` (`WhichTransform` → matrix transform).
- Not run: Reactant (`test/reactant` environment has no Reactant installed here) and AMDGPU/Metal;
  the code path is backend-agnostic (`mul!` + broadcasts over views, as the previous backward path).

## Documentation changes

- `docs/src/speedytransforms.md`: update the two formulas (both directions are now one real
  matrix multiply with stacked matrices) and the memory remark.
- `CHANGELOG.md` entry; `SpeedyTransforms` version `1.0.0` → `1.0.0+DEV`.

## Known limitations

- `M.forward`, `M.backward`, `M.backward_real`, `M.backward_imag` no longer exist (internal
  fields; only `show` and the SpeedyTransforms tests used them).
- BLAS thread count inside a Slurm allocation is Julia's default `Sys.CPU_THREADS ÷ 2` (64 threads on
  the 16 allocated cores). Now that both multiplies are BLAS-bound this oversubscription costs 1.3–1.5×
  (tables above); the benchmark driver could set `BLAS.set_num_threads` from `SLURM_CPUS_PER_TASK`.
  Not changed here.

## Future work

- Instability of the default `PrimitiveWetModel` with 24 layers (all T) and 16 layers at T171+
  (problem 1 of the regression analysis).
- Make the benchmark suite record NaN for crashed runs.
