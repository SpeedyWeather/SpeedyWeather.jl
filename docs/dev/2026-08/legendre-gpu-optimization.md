# GPU Legendre transform optimization

> Status: **completed**. Steps 0–2 are implemented, unit tested and verified against the CPU reference: the Legendre transforms are 5–37× faster at T127–T511 and are no longer the bottleneck of `transform!`. Step 3 was dropped in favour of `MatrixSpectralTransform` at low resolution.

Date of initial draft: 2026-08-07

Base revision: `12fbce2e` (drafted on `mg/fix-filter-and-diffusion`); implementation is being carried on top of `8adf1052` (`main`)

## Originating prompt

> We want to optimize the legendre! for better GPU performance across all truncations, analyse the
> functions for potential performance bottlenecks, profile the function and make a plan to further
> optimize it.

## Revision log

- **2026-08-07, initial draft.** Analysis, benchmark sweep (T31–T511) and profiling of the two GPU
  Legendre kernels; five root causes identified and a four-step plan proposed.
- **2026-08-07, "Start implementing the plan, after each step make a quick benchmark to make sure
  it's really useful".** Implementation started.
  - *Step 0 implemented and measured.* Inverse transform improved 1.3×–3.2× across truncations.
    Step 0 also introduced a **forward-transform regression at T511** (16.0 ms → 68.0 ms): moving
    from a layer-major 1-D launch to a 2-D `(row, layer)` launch puts all layers in flight
    simultaneously, so the atomically-accumulated coefficient array (8.4 MB at T511) no longer fits
    the A40's 6 MB L2, whereas the old launch kept a single 1.05 MB layer column resident. This is
    moot once Step 1 removes the atomics, and it is the reason Step 1 was not deferred.
  - *Step 1 design changed.* The plan assumed the rings contributing to a given order `m` form a
    band `j_start(m):nlat_half` reaching the equator, which would follow from `mmax_truncation`
    being monotonic pole→equator. **Measured and disproved**: `LegendreShortcutLinCubCoslat` and
    `LegendreShortcutLinQuadCoslat²` are non-monotonic, because their denominator grows toward the
    equator (`nlon/(2+2cos(lat))`, `nlon/(2+cos²(lat))`) faster than `nlon` does once `nlon`
    saturates. Of 120 grid × shortcut × truncation combinations, **16 have a contributing band that
    stops short of the equator** — always a contiguous interval, never with an interior hole
    (0/120), but summing to `nlat_half` would still wrongly include the equatorward rings. The
    kernel therefore walks an explicit, precomputed list of contributing rings instead of a range,
    which costs one warp-broadcast `Int32` load per iteration and assumes nothing about the
    shortcut.
- **2026-08-07, "Inspect the CLAUDE.md of the project again and archive the plan in accordance to it
  in a docs/dev folder".** Plan moved from the repository root to
  `docs/dev/2026-08/legendre-gpu-optimization.md` and restructured to the template in `CLAUDE.md`.
  This is the first document under `docs/dev/`.
- **2026-08-07, Step 2 implemented, with a changed approach.** The plan called for a
  workgroup-cooperative kernel staging tiles of the polynomials and coefficients in `@localmem`.
  Implemented instead as a **second, transposed copy of the Legendre polynomials** read by the
  inverse kernel, whose thread table is now ordered by order `m` first so that neighbouring threads
  hold neighbouring latitude rings. This needs no shared memory, no `@synchronize` and no
  compile-time-constant tile sizes (KernelAbstractions requires static `@localmem` extents, which
  the runtime layer count would have forced through a `Val` parameter), and it landed the expected
  speedup: the inverse transform reaches 37–41 % of peak bandwidth, up from 1.4–8 %. Its cost is
  memory, see *Memory consumption* below.
  - *Two type parameters had to be removed again.* Adding `MatrixType` and `ArrayTypeIntVector` to
    `SpectralTransform` pushed its type name to 1046 characters and **broke the regression guard in
    `SpeedyTransforms/test/type_name_length.jl`** (cap 1024, margin was only ~10 characters before
    this work). Exceeding the cap silently disables nested Enzyme AD. Fixed by storing the
    transposed polynomials flattened into the existing `VectorType`, and by dropping the separate
    `m_rings` array: the rings inside the Legendre truncation are now ordered first within each
    order in `jm_indices`, so they form the contiguous row range the forward transform sums over.
  - *Decision on memory (user, 2026-08-07):* "The factor 2 in memory consumption is still fine for
    now" — both layouts are kept; the numbers are recorded below so the trade can be revisited.
- **2026-08-07, `Atomix` dependency removed.** The forward transform's atomics were its only user
  anywhere in the repository, so it was dropped from `SpeedyTransforms` and from `SpeedyWeather`
  (which declared but never used it).
- **2026-08-07, "step 3 is skipped, for very low resolutions we have the matrix transform".** Step 3
  dropped; `MatrixSpectralTransform` is the answer at low resolution, not shaving launch overhead
  off the FFT + Legendre path.
- **2026-08-07, loop form and fused multiply-adds.** The forward kernel's ring loop was written as a
  `while` to match the comment in the inverse kernel claiming that is faster inside a kernel
  (from #815). Re-measured: **the claim does not hold on this stack.** The apparent 1–7 % advantage
  reversed when the measurement order was swapped (`for` ~244 µs consistently, `while` 246–254 µs at
  T255), i.e. it was inside the ~4 % run-to-run spread. The generated PTX in fact favours `for`
  (161 instructions vs 205, 34 loads vs 43), because `r += 1` adds an `Int64` literal to an `Int32`
  counter and drags 64-bit conversions through the loop, while iterating a `UnitRange{Int32}` stays
  32-bit. The forward loop was switched to `for`; the inverse keeps its `while` (two induction
  variables advancing together plus the odd-degree tail) with the performance claim replaced by the
  real reason.
- **2026-08-07, "check if using FMA in the kernel further speeds it up".** Julia does not contract
  `a + x*z` into fused multiply-adds by itself. Fusing the two real components explicitly gives a
  consistent **3.6–7.4 % on the forward transform** and nothing measurable on the inverse (1–2 %,
  inside noise — that one is bandwidth-bound at 41 % of peak, with little arithmetic to fuse).
  Applying `muladd` to the *complex* number instead promotes the real factor and does a full complex
  multiply-add, measuring 12–16 % *slower* than plain `*`/`+`.
- **2026-08-07, "test `muladd` within `fma_complex` to make sure it's really performant across all
  architectures and vendors".** The helper originally used `fma` on the components, which is a
  portability hazard: `fma` requires a genuine single-rounding fused instruction and degrades to
  software emulation where the hardware lacks one, while this code has to run on CUDA, AMDGPU,
  Metal and the KernelAbstractions CPU backend. Switched to `muladd` (renaming the helper to
  `muladd_complex`) and **verified that this costs nothing on CUDA**: the generated PTX has an
  identical instruction mix for both spellings (inverse 189 instructions, 6 `fma.rn.f32`, 34 loads;
  forward 155 instructions, 4 `fma.rn.f32`, 34 loads), differing only in load scheduling. So the
  contraction happens either way where hardware supports it, and the fallback stays safe where it
  does not.
  - *A measurement lesson worth recording:* a standalone rewrite of the inverse kernel suggested a
    17–22 % FMA win, but its own no-FMA baseline was ~20 % slower than the production kernel
    (loads hoisted into locals changed scheduling). Kernel micro-experiments must be validated
    against the real kernel through the normal benchmark before their numbers are believed.

## Problem description

The Legendre transform is the dominant cost of `transform!` on GPU at production resolutions: it is
40–50 % of a spectral↔grid transform at T255 and 75–81 % at T511. Both GPU kernels sustain only
1–7 % of the useful memory bandwidth the hardware can deliver, i.e. they are 6× (T127) to 66× (T511)
away from a memory-traffic roofline. The goal is to close most of that gap without giving up
portability (the transform must keep running through KernelAbstractions on CUDA, AMDGPU and Metal).

## Background

**Hardware and setup for all numbers below:** NVIDIA A40 (Ampere, 84 SMs, ~696 GB/s peak DRAM,
6 MB L2), Julia 1.12.2, CUDA.jl 6.2.1, KernelAbstractions 0.9.42, `Float32`, `nlayers = 8`, default
`SpectralGrid` (`OctahedralGaussianGrid`, quadratic dealiasing, linear Legendre shortcut). Timings
are min-of-30 samples of 10-call averages taken with CUDA events. Nsight Compute hardware counters
are admin-restricted on this cluster (`ERR_NVGPUCTRPERM`), so the claims below were confirmed with
targeted experiments instead of counter data.

### Measured baseline (before any change)

| trunc | nlat_half | nlm | inverse `_legendre!` | of which fills | forward `_legendre!` | Legendre share of `transform!` (inv / fwd) |
|------:|----------:|--------:|------------:|-----------:|------------:|:----------------:|
| T31   | 24  | 560     | 20.3 µs   | 13.7 µs (68 %) | 14.7 µs   | 10 % / 11 % |
| T63   | 48  | 2 144   | 27.9 µs   | 14.0 µs (50 %) | 52.6 µs   | 6 % / 14 %  |
| T127  | 96  | 8 384   | 110.4 µs  | 44.1 µs (40 %) | 292.9 µs  | 9 % / 27 %  |
| T255  | 192 | 33 152  | 1 809.7 µs | 159.8 µs (9 %) | 2 054.7 µs | 40 % / 50 % |
| T511  | 384 | 131 840 | 26 968.7 µs | 811.2 µs (3 %) | 16 003.5 µs | 81 % / 76 % |

Two regimes: below T255 the Legendre step is launch- and fill-bound but only ~10 % of `transform!`
(the Fourier side dominates, already CUDA-graph accelerated); from T255 up it dominates and is
limited by memory access inefficiency rather than DRAM capacity.

### Roofline

Irreducible traffic per call assuming perfect on-chip reuse (read Λ once, read/write the
coefficients once, write the Fourier scratch once), at the ~520 GB/s this GPU sustains for coalesced
streaming:

| trunc | floor bytes | floor time | measured kernel | headroom |
|------:|------------:|-----------:|----------------:|---------:|
| T63   | 1.1 MB  | ~2 µs   | 14 µs    | ~6×  |
| T127  | 5.7 MB  | ~11 µs  | 66 µs    | ~6×  |
| T255  | 32 MB   | ~62 µs  | 1 670 µs | ~27× |
| T511  | 206 MB  | ~400 µs | 26 160 µs | ~66× (inverse), ~40× (forward) |

### Root causes

The original design ran one thread per `(latitude ring j, order m, layer k)` triple, each doing a
serial dot product over degree `l` (19–513 iterations, mean ~150 at T255):

1. **Uncoalesced loads of both operands.** Warp lanes were consecutive orders `m`, and adjacent
   lower-triangular columns start ~`lmax` elements apart, so each warp load requested 32 scattered
   addresses. A synthetic probe replicating exactly this access pattern measured **4.4× (T255-like)
   to 5.2× (T511-like)** bandwidth loss versus the same data walked coalesced (99–119 GB/s vs
   513–526 GB/s).
2. **No on-chip reuse by construction.** The coefficients are needed by every ring and the Legendre
   polynomials by every layer, but nothing was staged in shared memory, so both were re-fetched.
   Layer-count scaling at T255 confirmed it: the inverse transform grew **super-linearly**, 13.2×
   the time for 8× the layers (134.9 µs at K=1 → 1 775.8 µs at K=8). Even K=1 — every 2-D surface
   field transform — was ~20× off its floor.
3. **Contended atomics in the forward kernel.** Every coefficient was accumulated by up to
   `nlat_half` (=192 at T255) threads through `Atomix.@atomic` on a `reinterpret`ed real array. The
   forward kernel was 2.3× slower than the inverse at T31, where loads are cheap — that gap is
   purely atomic cost.
4. **Whole-array zero-fills plus a read-modify-write.** Each output element of the Fourier scratch
   is written by exactly one thread, so `+=` (and the two `g .= 0` broadcasts it required) were
   unnecessary except for the frequencies past the Legendre truncation, which the transform never
   touches but the FFT does read. The broadcast fills ran at ~62 GB/s (`fill!` on the same array
   does ~640 GB/s) and were 40–68 % of the call below T255.
5. **Oversized index metadata.** `kjm_indices` stored three `Int64` per thread (24 B), replicated
   over layers — ~25 MB per launch at T511, comparable to the entire Legendre polynomial array,
   despite the layer index being derivable from the thread index.

### Target design

Per order `m`, both directions are a batched GEMM over the contiguous degree index `l`:

- inverse: `g[j,k] = Σ_l Λ[l,j] · specs[l,k]`, with the even/odd accumulator pair giving the
  northern and southern hemisphere;
- forward: `specs[l,k] = Σ_j Λ[l,j] · f̃[j,k]`, the transposed reduction.

cuBLAS grouped-batched GEMM was considered and rejected: it is CUDA-only, handles the variable
per-`m` sizes awkwardly, and cannot fuse the even/odd recombination, the longitude-offset rotation
or the north/south pair. Portable KernelAbstractions kernels (with `@localmem` tiling for Step 2)
are used instead.

## Summary of changes

### Step 0 — index and fill rework (implemented)

`SpeedyTransforms/src/spectral_transform.jl`, `src/legendre_ka.jl`, `src/legendre.jl`:

- Replaced `kjm_indices` (`Int64`, `(rows × layers, 3)`) with `jm_indices` (`Int32`,
  `(rows, 4)` = ring, order, column offset, column length), not replicated over layers; the layer
  became the second dimension of a 2-D `ArrayWorkOrder` launch. The column offset and length are
  precomputed, so the kernels no longer recompute `get_lm_range`.
- The index table now also covers the Fourier frequencies **past** the Legendre truncation, with
  column length 0. Those threads write the zeros the FFT reads, which removes both `g .= 0`
  broadcasts; the inverse kernel writes with `=` instead of `+=`.
- `unscale_coslat!` now scales only the layers the transform actually wrote.
- Fixed vacuous bounds checks: `length(nlayers)` on an `Int` is always 1, so the scratch-size
  checks in both GPU `_legendre!` methods never fired. They now compare against `nlayers`.
- Precomputation is allocation-free (the old `kjm_indices[i, :] .= [k, j, m]` allocated a small
  vector per row, ~1 M allocations at T511).

### Step 1 — atomic-free forward transform (implemented)

- The forward kernel is now parallelised over its **output**: one thread per (coefficient `lm`,
  layer `k`), summing over the contributing rings. Every coefficient is written exactly once, so
  there are no atomics, no `fill!(specs, 0)` and no `reinterpret` to a real array. Warp lanes are
  consecutive `lm`, which are contiguous in memory, so the Legendre polynomial loads and the
  coefficient stores are coalesced.
- Two new precomputed fields on `SpectralTransform`: `lm_indices` (`Int32`, `(nlm, 4)` = order,
  hemisphere sign, first and last index into `m_rings`) and `m_rings` (`Int32`, the contributing
  rings grouped by order). The explicit ring list rather than a range is what makes this correct for
  non-monotonic Legendre shortcuts (see revision log).
- Summation now proceeds in ascending ring order, matching the CPU reference exactly, instead of the
  nondeterministic atomic order.

### Step 2 — coalesced inverse transform (implemented)

- `SpectralTransform` gained `legendre_polynomials_transposed`: the polynomials in
  (latitude ring, harmonic) order, flattened to a vector (see the type-name note in the revision
  log) and only built on GPU.
- `jm_indices` is now ordered by order `m` first, latitude ring `j` second, so neighbouring threads
  hold neighbouring rings at one order. The inverse kernel then reads the transposed polynomials
  fully coalesced, and the coefficients — identical across the rings of one order — become a warp
  broadcast. Within an order, the rings inside the Legendre truncation come first, followed by the
  zero-only rows, which is what gives the forward transform its contiguous row range.
- The cost is the write to the Fourier scratch, which is now strided (neighbouring threads hold
  neighbouring rings, and the scratch is indexed frequency-first). That is 2 writes against `L_m`
  reads per thread, and the measurements show it does not dominate.

### Fused multiply-adds

Both kernels accumulate through `muladd_complex(x::Real, z::Complex, a::Complex)`, a small helper
computing `a + x*z` with the real and imaginary parts fused separately. Julia does not contract that
expression on its own. Two spellings matter and are easy to confuse:

- **`muladd` applied to the complex number** promotes the real factor to complex and performs a full
  complex multiply-add — 4 multiplies where 2 suffice. Measured **12–16 % slower** than not fusing
  at all. This is the trap.
- **`muladd` applied to the two real components** is what the helper does, and is correct.

`muladd` rather than `fma` on the components is deliberate and portability-driven: `fma` demands a
true single-rounding fused instruction and falls back to slow software emulation where the hardware
has none, whereas `muladd` contracts where a fused instruction exists and degrades to a plain
multiply and add otherwise. On CUDA the two are equivalent — the generated PTX has an identical
instruction mix (inverse 189 instructions with 6 `fma.rn.f32`, forward 155 with 4; both spellings
agree exactly, differing only in load scheduling) — so `muladd` costs nothing here and keeps
AMDGPU, Metal and the CPU path safe.

Worth 3.6–7.4 % on the forward transform, nothing measurable on the inverse (1–2 %, inside noise:
it is bandwidth-bound at 41 % of peak with little arithmetic to hide). Accuracy improves marginally
(~1e-7 relative, fewer roundings) and now matches the CPU reference, which already used `muladd`.

### Dependency cleanup

The forward transform's atomics were the only user of `Atomix` anywhere in the repository, so the
dependency was dropped from `SpeedyTransforms` (`import Atomix` plus its `[deps]`/`[compat]`
entries). `SpeedyWeather` declared it too but never used it — dead weight predating this work —
and was cleaned up in the same pass. `Atomix` remains in the manifests as a transitive dependency
of `KernelAbstractions`; only the two direct dependencies were removed, so the manifest diff is
minimal rather than a full re-resolve.

### Step 3 — launch-count reduction for small truncations (dropped)

The plan proposed folding the Legendre kernels into the existing CUDA-graph capture in
`SpeedyTransformsCUDAExt` to recover the last ~5–10 µs of launch latency per call, which only
matters at small truncations. **Dropped**: at very low resolution the answer is not to shave launch
overhead off the FFT + Legendre path but to use `MatrixSpectralTransform` instead, which replaces
the whole ring-by-ring transform with one dense matrix-matrix multiply and is already the
recommended choice there (see the "When to use it" section of `docs/src/speedytransforms.md`; it is
also the transform the Reactant/XLA path requires). Its dense matrices scale as
`O(n_grid × n_harmonics)`, which is exactly why it suits low resolution and not high — the
complementary regime to the kernels optimised here.

## Testing and verification

- `SpeedyWeather/test/GPU/legendrebench/verify_legendre.jl` compares the GPU `_legendre!` in both
  directions against the CPU reference, plus full `transform!` round trips and the `unscale_coslat`
  path, over 6 grids × 2 truncations × 2 layer counts. It deliberately pre-fills the scratch memory
  with `NaN` and with large garbage values so that any element the kernels fail to write is caught
  rather than masked by a leftover zero — which is what makes it a real test of the removed fills.
  **All checks passed after every step** (121 after Step 0, 122 after Steps 1 and 2; max relative
  error 6.4e-7, 1.2e-5 on the 10× round trip).
- `SpeedyTransforms/test/type_name_length.jl` is the guard that caught the type-parameter problem;
  it passes again.
- Unit test `"legendre kernels: compare to CPU"` added to
  `SpeedyWeather/test/GPU/spectral_transform.jl`: both directions against the CPU reference over
  2 grids × 2 layer counts, with the scratch memory and the output coefficients dirtied first so a
  missing write cannot pass, plus a round trip on `HEALPixGrid` with `LegendreShortcutLinCubCoslat`
  that asserts (via `!issorted(mmax_truncation)`) it really exercises the non-monotonic case.
- Benchmarks use `SpeedyWeather/test/GPU/legendrebench/bench_legendre.jl` (sweep T31–T511 with
  bandwidth models and `CUDA.@profile` per-kernel tables). Supporting experiments:
  `micro_experiments.jl` (layer-count scaling, coalescing probe) and `ncu_target.jl` (minimal Nsight
  Compute target, unusable on this cluster).
- The Enzyme AD rules only wrap `_fourier!` and the `transform!` entry points, whose signatures are
  unchanged, but `test/differentiability/speedy_transforms.jl` still needs a re-run.

### Measured results

Kernel times per call (`nlayers = 8`, `Float32`, default grid), original code → after each step:

| trunc | inverse: before | Step 0 | Step 2 | +FMA | **total** | forward: before | Step 0 | Step 1 | +FMA | **total** |
|------:|----------:|---------:|--------:|------:|------:|-----------:|-----------:|---------:|------:|------:|
| T31   | 20.3 µs   | 6.6 µs   | 7.1 µs  | 6.9 µs  | **2.9×** | 14.7 µs    | 16.6 µs    | 9.7 µs   | 9.8 µs   | **1.5×** |
| T63   | 27.9 µs   | 16.1 µs  | 9.6 µs  | 9.6 µs  | **2.9×** | 52.6 µs    | 65.9 µs    | 16.5 µs  | 17.0 µs  | **3.1×** |
| T127  | 110.4 µs  | 73.5 µs  | 21.7 µs | 22.2 µs | **5.0×** | 292.9 µs   | 301.9 µs   | 35.1 µs  | 35.4 µs  | **8.3×** |
| T255  | 1 809.7 µs | 566.4 µs | 120.7 µs | 119.2 µs | **15.2×** | 2 054.7 µs | 1 834.1 µs | 246.6 µs | 235.2 µs | **8.7×** |
| T511  | 26 968.7 µs | 20 726.5 µs | 722.8 µs | 713.2 µs | **37.8×** | 16 003.5 µs | 67 982.8 µs | 1 689.8 µs | 1 543.7 µs | **10.4×** |

The Step 0 column for the forward transform shows the L2 regression described in the revision log,
which Step 1 removes. Useful bandwidth of the inverse transform went from 1.4–8 % of peak to
37–41 %, of the forward transform from 2–3 % to 18–23 % — both are now ordinary memory-bound
kernels rather than pathological ones.

Whole-transform effect, and the share the Legendre step still takes:

| trunc | spec→grid before → after | grid→spec before → after | Legendre share now (inv / fwd) |
|------:|------------------------:|------------------------:|:------------------:|
| T31   | 206.7 → 200.0 µs (1.03×) | 140.3 → 135.6 µs (1.03×) | 3.5 % / 7.1 % |
| T63   | 508.3 → 451.7 µs (1.13×) | 365.1 → 328.6 µs (1.11×) | 2.1 % / 5.0 % |
| T127  | 1 242.7 → 1 150.7 µs (1.08×) | 1 092.9 → 858.8 µs (1.27×) | 1.9 % / 4.1 % |
| T255  | 4 507.0 → 2 808.9 µs (1.60×) | 4 113.6 → 2 307.6 µs (1.78×) | 4.3 % / 10.7 % |
| T511  | 33 350.7 → 7 128.1 µs (4.68×) | 21 179.9 → 6 826.9 µs (3.10×) | 10.1 % / 24.8 % |

The Legendre transform was 76–81 % of `transform!` at T511 and is now 10–25 %; the Fourier stage
dominates at every truncation. Further work on `transform!` belongs there, not here.

### Memory consumption

The transposed copy costs exactly as much as the Legendre polynomials themselves, which are already
the largest precomputed array in a `SpectralTransform` (at T511 they are 5× the Fourier scratch).
They scale as **O(T³)** — `nlm ~ T²/2` times `nlat_half ~ T` — so every doubling of the truncation
multiplies them by ~8:

| trunc | Λ (Float32) | with transposed copy | Λ (Float64) | with copy | Fourier scratch (K=8, F32) |
|------:|------------:|---------------------:|------------:|----------:|---------------:|
| T31   | 52.5 KB  | 105 KB   | 105 KB   | 210 KB   | 171 KB |
| T63   | 402 KB   | 804 KB   | 804 KB   | 1.6 MB   | 630 KB |
| T127  | 3.1 MB   | 6.1 MB   | 6.1 MB   | 12.3 MB  | 2.4 MB |
| T255  | 24.3 MB  | 48.6 MB  | 48.6 MB  | 97.1 MB  | 9.2 MB |
| T511  | 193 MB   | 386 MB   | 386 MB   | 772 MB   | 36.4 MB |
| T1023 | 1.50 GB  | 3.01 GB  | 3.01 GB  | 6.02 GB  | 145 MB |

Only the GPU builds the copy; on CPU the field is empty. Up to T255 the doubling is immaterial
(≤ 49 MB). At T511 it is 386 MB, which is nothing on a 46 GB A40 but noticeable on an 8–16 GB card,
and at T1023 (or Float64 at T511) it becomes a real constraint.

**The single-copy alternative was measured** rather than assumed. Keeping only the transposed layout
and letting the forward transform read it strided costs (T63/T127/T255/T511):

| trunc | forward, both layouts | forward, transposed only | ratio | memory saved |
|------:|----------------------:|-------------------------:|------:|-------------:|
| T63   | 16.6 µs   | 17.0 µs   | 1.02× | 0.4 MB |
| T127  | 36.0 µs   | 69.6 µs   | 1.93× | 3.1 MB |
| T255  | 245.9 µs  | 414.1 µs  | 1.68× | 24.3 MB |
| T511  | 1 675.4 µs | 3 045.5 µs | 1.82× | 193 MB |

That is far cheaper than the ~4–5× the coalescing probe suggested, because each thread still walks
contiguous memory in the transposed layout and only the 32 lanes of a warp diverge, so sectors are
reused within a lane. In whole-transform terms it would cost nothing on spec→grid and about +7 %
(T255) to +20 % (T511) on grid→spec. **Decision: keep both copies for now** (factor 2 accepted); if
memory becomes the binding constraint at T1023 or in Float64, dropping to the transposed copy alone
is a small, well-understood change that retains most of the speedup.

## Documentation changes

None so far. The `SpectralTransform` struct gained three precomputed fields (`jm_indices`,
`lm_indices`, `m_rings`) which are internal and documented inline via `TYPEDFIELDS`. If the
non-monotonicity of the coslat-dependent Legendre shortcuts is worth surfacing to users, the place
would be the shortcut docstrings in `SpeedyTransforms/src/legendre_shortcuts.jl`.

## Known limitations

- **GPU memory doubles for the Legendre polynomials** (see *Memory consumption*): +24 MB at T255,
  +193 MB at T511, +1.5 GB at T1023 in Float32, and twice that in Float64. Accepted for now.
- All measurements are from a single A40. The balance between the two regimes (launch-bound vs
  memory-bound) will shift on other hardware; nothing in the changes is A40-specific, but the
  crossover truncation is.
- Non-CUDA backends (AMDGPU, Metal) are untested — the kernels stay pure KernelAbstractions with no
  warp intrinsics or shared memory, so they should work, but nobody has run them.
- ~~`import Atomix` is now unused~~ — removed, see *Summary of changes*.
- Step 0 changed the meaning of `SpectralTransform.jm_index_size` (it now counts the extra
  zero-writing rows), and Step 2 changed the row *order* of `jm_indices` from ring-major to
  order-major. Any external code reading those fields would need updating; nothing in the
  repository does.
- The `SpectralTransform` type name is close to the LLVM value-name cap that nested Enzyme AD
  depends on. This work consumed the remaining margin twice over before being reworked to fit; the
  next person adding a field should keep it within an existing type parameter and re-run
  `SpeedyTransforms/test/type_name_length.jl`.

## Future work

- **The Fourier stage is now the bottleneck of `transform!` at every truncation** (75–98 % of it).
  It runs ~2 FFTs plus a gather/scatter per latitude ring; batching multiple rings into one
  execution is the obvious next project, and it is where the remaining time is. This matters most
  at mid and high resolution, where `MatrixSpectralTransform` is not an option because its dense
  matrices grow as `O(n_grid × n_harmonics)`.
- Dropping to a single (transposed) copy of the polynomials if memory becomes binding — the cost is
  measured and recorded above.
- The layer-count scaling measurement suggested the K=1 case (2-D surface fields) was ~20× off the
  roofline; Steps 0–2 should have fixed that too, but it was measured only at K=8 and deserves a
  re-check at K=1.
