# MatrixSpectralTransform: real stacked matrices so both multiplies hit BLAS on every backend

> Status: **in progress**. Plan drafted after the regression analysis; implementation, CPU tests,
> GPU tests and CPU/GPU benchmarks to follow on branch `mg/matrix-transform-regression`.

Date of initial draft: 2026-10-02

Base revision: `8e99099f` (`origin/main`, after rebase; drafted against `8a29efdd` v0.23.0)

## Originating prompt

> Ok in a new branch `mg/matrix-transform-regression` we need to investigate problem 2. First
> investigate the change we have done in PR #1201, then we need to find a version that restores old
> CPU performance (or improves it) but also keeps GPU compatability and keeps or improves GPU
> performance

## Revision log

- 2026-10-02: initial draft.
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
2. **Scratch** becomes `(2 nharmonics × nlayers)` real (`nlayers` = the scratch width the
   `SpectralGrid` passes), used by both directions.
3. **Grid→spectral**: for each block of at most `size(scratch, 2)` columns,
   `mul!(scratch, forward_stacked, field_block)` (one real `gemm!` / CUBLAS / rocBLAS) followed by
   one broadcast `coeffs_block .= complex.(scratch[1:n, :], scratch[n+1:2n, :])`.
4. **Spectral→grid**: per block, two broadcasts fill the scratch halves with `real`/`imag` of the
   coefficients and one `mul!(field_block, backward_stacked, scratch)` replaces the two `gemm!`s
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
  (a) a batch wider than the scratch (`nlayers = 3` transform, 8-layer data) still round-trips, i.e.
  the column blocking works, (b) transforms into/out of views (2-D view of a 3-D parent, 1-D slot
  view) agree with the plain-array result, on CPU.
- `SpeedyWeather/test/dynamics/matrix_transform.jl` (20-day `PrimitiveWetModel` run, no NaN).
- GPU: `SpeedyWeather/test/GPU/runtests.jl` (CUDA env) on an H100 node, in particular the existing
  "transform! on views into a larger backing array" regression test and the model tests that use
  the matrix transform.
- Benchmarks: the reduced benchmark from the regression analysis (`claude_bisect/bench_subset.jl`,
  T32 L8 MT) on the same node as before (baseline 107, HEAD 57 SYPD), and the GPU suite's matrix
  rows (`manual_benchmarking.jl gpu`) against the committed README.

## Documentation changes

- `docs/src/speedytransforms.md`: update the two formulas (both directions are now one real
  matrix multiply with stacked matrices) and the memory remark.
- `CHANGELOG.md` entry; `SpeedyTransforms` version `1.0.0` → `1.0.0+DEV`.

## Known limitations

- `M.forward`, `M.backward`, `M.backward_real`, `M.backward_imag` no longer exist (internal
  fields; only `show` and the SpeedyTransforms tests used them).
- BLAS thread count inside a Slurm allocation is still `Sys.CPU_THREADS ÷ 2` (64 threads on 16
  cores); now that the forward multiply is BLAS-bound this oversubscription matters more for the
  benchmark numbers. Measured as part of verification; not changed here.

## Future work

- Instability of the default `PrimitiveWetModel` with 24 layers (all T) and 16 layers at T171+
  (problem 1 of the regression analysis).
- Make the benchmark suite record NaN for crashed runs.
