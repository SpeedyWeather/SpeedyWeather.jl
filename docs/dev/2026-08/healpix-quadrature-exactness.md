# Exact HEALPix spectral transforms via per-order quadrature weights

> Status: **completed**. Implemented, unit tested and verified. The round-trip error on all four
> HEALPix grids drops from ~5·10⁻³ to roundoff (10⁻¹⁵ in Float64, 10⁻⁷ in Float32) with every
> quadrature weight strictly positive and within ±18% of equal area. Grids with an exact quadrature
> rule are bit-identical to before and the forward transform is no slower.

Date of initial draft: 2026-08-31

Base revision: `20467269` (drafted on `mg/healpix-exactness`)

## Originating prompt

> I want to work on optimizing the quadrature weights of the Healpix spectral transform to ensure
> higher transform exactness. Currenlty all rings are equally weighted and the Healpix transform has
> a relatively high inaccurary that also leads to problems in model stability (espacially for low
> resolutions and low dealiasing). Can we optimize the quadrature weights in a way that reduces the
> transform error? See also https://arxiv.org/pdf/2510.01785 whether they already implemented
> something like this

> We have four target resolutions: T32, T64, T128, T256

> What's different between T31 and T32? If that is preferable, we might be able to change the target
> to one degree less

## Revision log

- **2026-08-31, initial draft.** Root-cause analysis, exact operator model validated to machine
  precision against `transform!`, weight optimisation prototyped, end-to-end verification at the
  four target resolutions.
- **2026-08-31, reduced formulation.** The full `O(n_m²)`-row least-squares system was replaced by
  an `n_m`-row one after noticing the off-diagonal exactness conditions are implied by the diagonal
  ones (see "Per-order quadrature weights"). This removed the regularisation knob entirely, cut the
  condition number from ~10¹⁰ to ~10¹–10², and made every weight strictly positive. All tables below
  are from the reduced formulation.
- **2026-08-31, "Okay I reviewed the plan, execute it".** Option B selected (`dealiasing = 3.5`,
  `nlat_half = 36/72/144/288` for T32/T64/T128/T256). Implementation started.
- **2026-08-31, "keep changes to the non-HEALPix transform absolutely minimal".** Three changes
  followed. (i) The sub-Nyquist truncation was moved out of `LegendreShortcutLinear` and into the
  `SpectralTransform` constructor, gated on `optimize_quadrature(grid)`. Applied globally it
  degraded `OctahedralGaussianGrid` at dealiasing 2 (T32: 1.5·10⁻¹⁰ → 2.4·10⁻¹⁰; T128:
  5.7·10⁻⁷ → 8.5·10⁻⁷), because there it removes ring contributions without a refit to compensate.
  (ii) The weights are fused into
  `conj(lon_offsets)` instead of living in their own array — listed as future work in the first
  draft, but pulled forward because adding a type parameter broke the Enzyme value-name cap (see
  below). (iii) The fuse is done after rounding to `NF`, reproducing the kernel's original
  arithmetic exactly, so non-HEALPix output is bit-identical rather than merely equivalent.
- **2026-08-31, restriction to the solvable case.** The reduced formulation is equivalent to the
  full system only where a solution exists. Applying it as a least-squares fit below that made the
  `dealiasing = 2` HEALPix transforms *worse* (caught by the existing "inexact transforms" tests),
  because satisfying the diagonal conditions in a least-squares sense leaves the off-diagonal ones
  unconstrained. The fit now declines unless at least as many rings contribute as there are
  degrees, leaving those orders exactly as accurate as before.

## Problem description

The HEALPix and OctaHEALPix grids use `equal_area_weights`: every grid point is assigned the same
solid angle `4π/npoints`, so every ring carries the same weight per point. The resulting spherical
harmonic transform is not exact. Measured relative L2 error of the spectral → grid → spectral
round-trip (Float64, random band-limited field, default `dealiasing = 3`):

| grid | T32 | T64 | T128 | T256 |
|---|---|---|---|---|
| `HEALPixGrid` | 7.5·10⁻³ | 3.7·10⁻³ | 9.8·10⁻⁴ | 5.5·10⁻⁴ |
| `OctaHEALPixGrid` | 5.8·10⁻³ | 2.8·10⁻³ | 8.2·10⁻⁴ | 4.7·10⁻⁴ |

The worst-case amplification is the more relevant number for model stability: the round-trip
operator deviates from the identity by `‖T − I‖₂ = 0.12` at *every* resolution — it does not
converge away — and `‖T‖₂ = 1.005 > 1`, i.e. the transform pair injects energy.

## Background

### What exactness requires

Write the analysis (grid → spectral) and synthesis (spectral → grid) operators as `A` and `S`.
Exactness of the transform pair means `A∘S = I` on the retained spectral space. Because the Fourier
step is exact per ring, the condition decouples by order `m` and reads

```
Σ_j  g_j · λ_lm(μ_j) · λ_l'm(μ_j)  =  δ_ll'        for all l, l' ≥ m up to lmax
```

where `λ_lm` are the normalised associated Legendre polynomials, `μ_j = sin(lat_j)` and `g_j` is the
total solid angle of ring pair `j` (north + south), with `Σ_j g_j = 4π`. In the code
`g_j = nlon_j · solid_angles[j] · (2 unless j is the equator)`.

Two consequences matter:

1. **The condition is linear in `g`.** Minimising `‖A∘S − I‖_F` over the weights is an ordinary
   linear least-squares problem, not a nonlinear optimisation.
2. **It has far fewer independent conditions than it looks.** Substituting `x = μ²`, the products
   `λ_lm λ_l'm` with `l + l'` even span `(1−x)^m ×` {polynomials of degree ≤ `lmax₀ − m` in `x`}.
   So order `m` imposes only `n_m = lmax₀ − m + 1` independent conditions on the weights, and the
   exactness problem is *solvable* whenever the number of contributing rings satisfies
   `|J_m| ≥ n_m`. For `m = 0` that is `nlat_half ≥ truncation`.

This is why the problem is worth attacking at all: HEALPix is usually described as a grid that
"does not support exact spherical harmonic transforms" (cuHPX, §V), but that statement is about
*equal-area* weights, not about the grid.

### Two error sources, not one

The equal-area weighting is only the larger of two independent defects.

**(a) The Nyquist bin of each ring.** A ring with `nlon_j` longitudes has its Fourier bin
`m = nlon_j÷2` at the Nyquist frequency, where the real-FFT pair `rfft ∘ brfft` keeps only *one* of
the two components. `LegendreShortcutLinear` returns `nlon ÷ 2`, so that bin *is* retained. On
HEALPix polar-cap rings (`nlon_j = 4j`, retained up to `m = 2j`) this is hit on every single ring,
and because those rings carry a half-pixel longitude offset the surviving component is
`o = exp(i·m·lon₁) = i` — i.e. the component *orthogonal* to the one the quadrature needs.

The practical consequence is that at even orders `m` the real and imaginary parts of the same
coefficient are analysed over **different sets of rings**. No choice of a single weight per (ring,
order) can be consistent for both, so exactness is unreachable while this bin is retained. This was
confirmed by probing the true operator with unit basis vectors: at `T16`, `‖T_re − T_im‖` is
2.1·10⁻² at `m = 2`, 1.5·10⁻² at `m = 4`, 9.2·10⁻³ at `m = 6`, and exactly zero at every odd `m`.

**(b) The equal-area weights themselves**, which are simply not a quadrature rule.

### What the literature does

`arXiv:2510.01785` (cuHPX, Cheng et al., *GPU-Accelerated Differentiable Spherical Harmonic
Transforms on HEALPix Grids*) evaluates exactly this question in §IV-A. It compares equal, ring and
pixel weights and finds ring weights best: *"reducing the round-trip error by nearly an order of
magnitude compared to equal weighting"*, and *"adopts ring weights as the default"*. The weights
are order-independent by construction — the paper's own definition is *"ring (all pixels in a ring
share the same weight)"*.

The paper does **not** discuss making the weights depend on the order `m`, anywhere. Every mention
of weights in it uses `w_j`, indexed by ring alone — *"W is a diagonal matrix of quadrature
weights"* — and §V ("Limitation & Future Work") proposes only three extensions: multi-GPU
distributed transforms, using cuHPX as the core of a PDE solver, and iterative refinement.

Their own comparison points at why. The three schemes they try vary the weights **spatially**, with
strictly increasing degrees of freedom: equal (one value everywhere), ring (one per ring), pixel
(one per pixel). The most spatially-resolved of the three is not the best — the reported ordering is
equal < pixel < ring, i.e. *"pixel weights provide moderate improvement. Ring weights, however,
provide the best results"*. Spatial freedom stops paying off. The axis never varied is the spectral
one, and that is where the freedom is needed: the exactness condition is a separate constraint
system per order `m`, so per-order weights are not *finer* weights than pixel weights, they are
weights along a different axis.

Two things to keep straight about what the paper does and does not say:

- It never states where its ring weights come from. The presumption here — that they are the
  classical HEALPix ones (Górski et al. 2005, shipped as FITS tables, used by healpy/ducc) — is an
  inference from the standard meaning of the term and from the fact that cuHPX benchmarks against
  healpy and ducc. It is not a quote. The comparison in the table below does not depend on it: the
  row labelled "classical ring weights" is that scheme reimplemented and measured here, and stands
  as the best-known order-independent ring weighting regardless of which variant cuHPX ships.
- **Iterative refinement is not merely future work in that literature.** §II-A ("Standard Methods
  in HEALPix Software") presents it with explicit equations (4) and (5) as what traditional CPU
  libraries already do — *"This allows the system to be solved iteratively or via gradient descent,
  with accuracy improved by increasing the number of iterations"* — and rejects it on cost:
  *"its computational cost scales as O(ℓ³max), which quickly becomes prohibitive at high resolutions
  (e.g., ℓmax ≳ 2000)"*. §V then revisits it as an extension cuHPX could explore. So the established
  fallback for HEALPix inexactness is an O(ℓ³) iteration paid at **every transform**; per-order
  weights buy exactness in a one-off precomputation and cost nothing per transform. That is the
  sharper contrast, and the first draft of this document understated it.

The classical ring weights solve a strictly weaker problem than the one above: they only require
`Σ_j g_j λ_l0(μ_j) = 2√π δ_l0`, i.e. that the quadrature integrate every band-limited function
exactly. That is the `m = 0, l' = 0` column of the condition — one row per even degree instead of
the full orthogonality system.

Reimplemented here and measured on the same footing:

| | T32 | T64 | T128 | T256 |
|---|---|---|---|---|
| equal area (current) | 5.99·10⁻³ | 2.99·10⁻³ | 1.49·10⁻³ | 7.43·10⁻⁴ |
| classical ring weights (healpy/ducc) | 4.56·10⁻³ | 2.26·10⁻³ | 1.13·10⁻³ | 5.64·10⁻⁴ |
| per-order weights (this plan) | 7.20·10⁻⁴ | 3.36·10⁻⁴ | 1.62·10⁻⁴ | 7.99·10⁻⁵ |

(RMS of `T − I`; the per-order row here is the conservatively regularised variant, *not* the exact
one.) The classical ring weights buy only 1.3× in this metric. That is not in conflict with cuHPX's
"nearly an order of magnitude": the baseline is the same (equal weighting) but the quantity is not —
they measure `‖f − iSHT(SHT(f))‖` on a map, in RMS and L∞, where this table measures the deviation
of the round-trip *operator* from the identity over spectral coefficients. Per-order weights are
**7× better than the classical ring weights** in this metric even before exploiting exactness.

So: no, the paper has not implemented anything equivalent, and does not suggest it. It uses ring
weights that are one weight per ring, order-independent, and fitted to a weaker condition. The
per-order idea here comes from two properties of this codebase rather than from the literature: the
Legendre shortcut already makes the set of contributing rings depend on `m`, so one weight per ring
cannot be consistent across orders; and the exactness condition decouples by `m`, so each order is
an independent, small, and exactly solvable problem.

## Summary of changes

### 1. Drop the Nyquist bin from the Legendre shortcut

`LegendreShortcutLinear(nlon, latd) = nlon ÷ 2` → `(nlon - 1) ÷ 2`.

On full and octahedral grids this changes nothing that matters: there `nlon ≫ 2·mmax` on the rings
that carry weight, and where the shortcut does bite, `λ_lm ≈ 0` near the poles anyway (measured:
`OctahedralGaussianGrid` stays at 2.5·10⁻¹⁴). On HEALPix it is a prerequisite for everything below.

Note this change is **not** an improvement on its own — with equal-area weights it makes the HEALPix
error marginally *worse* (5.70·10⁻³ → 5.99·10⁻³ at T32), because it removes contributions without
re-fitting the weights. It only pays off together with step 2, so the two must land together.

### 2. Per-order quadrature weights

Replace the single `solid_angles[j]` with `solid_angles[m, j]`.

The exactness system looks like it has `O(n_m²)` rows, but **the off-diagonal conditions are implied
by the diagonal ones**. In `x = μ²` the products `λ_lm λ_l'm` (`l+l'` even) span
`(1−x)^m ×` {polynomials of degree ≤ `n_m−1`}, and the diagonal products `λ_lm²` are already a
triangular basis of that space. So the whole system collapses to `n_m` rows — *"integrate every
`|Y_lm|²` exactly"*:

```
Σ_j  g_j · λ_lm(μ_j)²  =  1        for l = m … lmax
```

Solve that `n_m × |J_m|` system for the minimum **relative** change from equal area, i.e. for
`u = (g − g₀)/g₀` rather than for `g` (the weights span three orders of magnitude between polar and
equatorial rings, so minimising `‖g − g₀‖` lets a polar ring swing past zero while barely moving the
objective). This needs no regularisation parameter: the reduced system is well conditioned
(`cond` 12–660 for `HEALPixGrid` at the target resolutions) and its minimum-norm solution is
already the one we want.

Only the forward (grid → spectral) kernel reads these weights, in `legendre.jl:169` and
`legendre_ka.jl:154`. The change is `ΔΩ = solid_angles[j]` (hoisted out of the `m` loop) →
`ΔΩ = solid_angles[m, j]` (inside it). On GPU `m` is the fast dimension, matching `lon_offsets`
which is already indexed `[m, j]`, so access stays coalesced.

The weights are stored **pre-multiplied by `conj(lon_offsets)`**, in a field named
`solid_angles_rotated`, which is exactly the factor the forward kernel needs — so it reads one array
instead of two and drops a multiply from the inner loop. This also keeps the storage in the
`MatrixComplexType` type parameter the struct already has, which turned out to be load-bearing:
adding a parameter for a separate real `[m, j]` array pushed the `SpectralTransform` type name from
1008 to 1025 characters, past the 1024-character LLVM value-name cap that
`test/type_name_length.jl` guards, which would have blocked nested Enzyme AD entirely.

Memory: `mmax × nlat_half × Complex{NF}` replaces an `nlat`-element vector — 512 KB at T256 in
Float32, against the ~250 MB of Legendre polynomials at that resolution. Negligible.

Construction cost (Float64, laptop CPU, prototype), against the current `SpectralTransform`
constructor: 0.0015 s at T32 (+181%), 0.011 s at T64 (+405%), 0.085 s at T128 (+663%), 1.08 s at
T256 (+1120%). Absolutely small below T128 but a large *relative* increase, and ~1 s at T256 is
noticeable. It scales as roughly `nlat_half³`. The prototype uses an SVD per order for the
minimum-norm solve; since the reduced system is well conditioned a cheaper factorisation should
recover most of that. T512+ may still want caching.

### 3. Keep `get_solid_angles` as-is

`RingGrids.get_solid_angles` is the *geometric* area of a grid cell and is used for area-weighted
means elsewhere. The optimised weights are a property of the *quadrature*, not of the grid, and must
not silently redefine what "the area of a HEALPix pixel" means. They belong in `SpeedyTransforms`,
computed from the grid + spectrum + shortcut, which is where all three are known.

### Choosing the truncation

Exactness is reachable only when there is spectral slack, `slack = nlat_half − truncation`. With
zero slack the `m = 0` system is square, its unique solution is an interpolatory quadrature at
maximal degree, and the weights oscillate the way Newton–Cotes does. This is a property of the
problem, not of the solver — it persists with the well-conditioned reduced formulation:

| grid | `nlat_half` | truncation | slack | cond | g/g₀ range | negative |
|---|---|---|---|---|---|---|
| `HEALPixGrid` | 32 | T32 | 0 | 1.4·10³ | [−7.8, 9.5] | 5 |
| `HEALPixGrid` | 32 | T30 | 2 | 1.2·10¹ | [0.905, 1.178] | 0 |
| `HEALPixGrid` | 64 | T64 | 0 | 1.3·10⁶ | [−3791, 3745] | 14 |
| `HEALPixGrid` | 64 | T60 | 4 | 2.0·10¹ | [0.904, 1.178] | 0 |
| `HEALPixGrid` | 128 | T128 | 0 | 1.5·10⁹ | [−1.1·10⁶, 1.1·10⁶] | 112 |
| `HEALPixGrid` | 128 | T120 | 8 | 5.6·10¹ | [0.903, 1.179] | 0 |
| `HEALPixGrid` | 256 | T256 | 0 | 9.3·10⁹ | [−9.2·10⁴, 9.2·10⁴] | 1696 |
| `HEALPixGrid` | 256 | T240 | 16 | 6.6·10² | [0.903, 1.179] | 0 |

`HEALPixGrid` needs `slack ≥ nlat_half/16`. `OctaHEALPixGrid` needs more — at `nlat_half = 256`,
T240 (`nlat_half/16`) is still unusable for it (cond 9.8·10⁹, 4926 negative weights) — but less than
`nlat_half/8`: `nlat_half/9` is verified sufficient at all four target resolutions.

**Decision (2026-08-31): Option B — keep the truncation, grow the grid.**

`dealiasing = 3.5` maps all four target resolutions onto valid HEALPix grids exactly, and gives
`slack = nlat_half/9`, which clears both grids' requirements:

| truncation | `nlat_half` | even | `nlon_max` | HEALPix `npoints` | slack |
|---|---|---|---|---|---|
| T32 | 36 | ✓ | 72 = 2³·3² | 3888 | 4 |
| T64 | 72 | ✓ | 144 = 2⁴·3² | 15552 | 8 |
| T128 | 144 | ✓ | 288 = 2⁵·3² | 62208 | 16 |
| T256 | 288 | ✓ | 576 = 2⁶·3² | 248832 | 32 |

Cost is ~27% more grid points than `dealiasing = 3`, paid in every grid-space parameterisation.
The rejected alternative (Option A) was to keep the grid and run T30/T60/T120/T240 at
`dealiasing = 3.267`, which is free but changes what `truncation` means on HEALPix grids.

Implementation: `default_dealiasing(::Type{<:Union{HEALPixGrid, OctaHEALPixGrid}}) = 3.5`, plus the
full/`Full…` HEALPix variants. Note `get_nlat_half` already returns an even `nlat_half` for every
input (`roundup_fft` starts even and steps by 2), which `HEALPixGrid` requires.

## Testing and verification

The prototype validates an explicit matrix model of `A∘S` against the real `transform!` to
**machine precision** (max relative deviation 3·10⁻¹⁶ … 1·10⁻¹⁵) for `HEALPixGrid`,
`OctaHEALPixGrid`, `FullHEALPixGrid`, `OctahedralGaussianGrid` and `FullClenshawGrid` at T16/T32/T64
— including the Nyquist real/imaginary asymmetry described above. All numbers in this document come
from that validated model or from `transform!` directly.

End-to-end round-trip with the real `transform!`, random band-limited field, Float64, at the
Option B resolutions (reduced formulation, no regularisation):

| grid | `nlat_half` | truncation | equal-area L2 | optimised L2 | g/g₀ range | negative | cond |
|---|---|---|---|---|---|---|---|
| `HEALPixGrid` | 36 | T32 | 5.64·10⁻³ | 7.90·10⁻¹⁶ | [0.913, 1.174] | 0 | 1.2·10¹ |
| `HEALPixGrid` | 72 | T64 | 2.83·10⁻³ | 1.37·10⁻¹⁵ | [0.913, 1.174] | 0 | 1.8·10¹ |
| `HEALPixGrid` | 144 | T128 | 6.55·10⁻⁴ | 2.58·10⁻¹⁵ | [0.913, 1.174] | 0 | 2.7·10¹ |
| `HEALPixGrid` | 288 | T256 | 3.75·10⁻⁴ | 4.86·10⁻¹⁵ | [0.913, 1.174] | 0 | 3.9·10¹ |
| `OctaHEALPixGrid` | 36 | T32 | 4.12·10⁻³ | 8.22·10⁻¹⁶ | [0.961, 1.147] | 0 | 1.3·10¹ |
| `OctaHEALPixGrid` | 72 | T64 | 2.07·10⁻³ | 1.46·10⁻¹⁵ | [0.960, 1.148] | 0 | 1.9·10¹ |
| `OctaHEALPixGrid` | 144 | T128 | 4.55·10⁻⁴ | 2.38·10⁻¹⁵ | [0.959, 1.148] | 0 | 1.9·10² |
| `OctaHEALPixGrid` | 288 | T256 | 2.63·10⁻⁴ | 4.56·10⁻¹⁵ | [0.929, 1.160] | 0 | 1.0·10⁵ |

In **Float32**, the model's default number format, with weights computed in Float64 and stored as
Float32, the prototype reached 1.2·10⁻⁷ … 2.8·10⁻⁷ — Float32 roundoff, a 2,000–60,000×
improvement. `‖T‖₂` becomes exactly 1, so the transform pair no longer injects energy. To be
re-confirmed against the implementation.

Unit tests added in `SpeedyTransforms/test/quadrature.jl` (23489 assertions, 3 s):

- round-trip exactness for all four HEALPix grids at T32/T64 in Float32 and Float64, held to
  `1e-10`/`1e-4` — tolerances the pre-existing equal-area weights (~5·10⁻³) miss by orders of
  magnitude, so the test cannot pass by accident;
- every weight strictly positive and within `[QUADRATURE_MIN_RATIO, QUADRATURE_MAX_RATIO]` of the
  grid's geometric solid angle;
- `Σ_j g_j = 4π` at `m = 0`, i.e. the global mean is exactly conserved;
- no ring is kept at its own Nyquist order on any HEALPix grid;
- grids with an exact quadrature rule keep their geometric solid angles, identical at every order;
- an under-resolved grid (`dealiasing = 2`) warns and falls back rather than returning unusable
  weights;
- `default_dealiasing == 3.5` for the HEALPix grids, giving an even `nlat_half` with `nlat_half/9`
  of slack at all four target resolutions, and every order fitted without a warning at T128.

Regression checks beyond the unit tests:

- **Non-HEALPix output is bit-identical.** 40 configurations (`OctahedralGaussianGrid`,
  `FullGaussianGrid`, `OctaminimalGaussianGrid`, `FullClenshawGrid`, `OctahedralClenshawGrid` ×
  T32/T64 × dealiasing 2/3 × Float32/Float64) hash-compared against `20467269`: all identical.
- **Full `SpeedyTransforms` suite passes**, including the Enzyme AD rules (8143 assertions) and the
  type-name guard (1008 characters, 16 of margin — unchanged).
- **No performance change.** Forward transform at T64/T128/T256 on three non-HEALPix grids: within
  ±2% of baseline, mostly marginally faster.
- **Models run.** `ShallowWaterModel` on `HEALPixGrid` and `OctaHEALPixGrid` at truncation 31
  (`nlat_half = 36`) integrates 5 days with finite, sane vorticity.

Measured round-trip error after the change, at each grid's default dealiasing:

| grid | T32 | T64 | T128 | T256 |
|---|---|---|---|---|
| `HEALPixGrid` before | 5.60·10⁻³ | 2.77·10⁻³ | 6.38·10⁻⁴ | 3.75·10⁻⁴ |
| `HEALPixGrid` after | 8.1·10⁻¹⁶ | 1.4·10⁻¹⁵ | 2.6·10⁻¹⁵ | 4.9·10⁻¹⁵ |
| `OctaHEALPixGrid` after | 8.3·10⁻¹⁶ | 1.5·10⁻¹⁵ | 2.4·10⁻¹⁵ | 4.6·10⁻¹⁵ |

In Float32: 1.4·10⁻⁷ … 3.2·10⁻⁷, i.e. Float32 roundoff. `‖T‖₂` becomes 1, so the transform pair no
longer injects energy.

`healpix_quadrature/transform_error.jl` reproduces these numbers against the current code.

Still worth doing separately: confirm the premise that transform inexactness was behind the HEALPix
stability problems, with a long model run at the configuration that previously failed soonest.

## Documentation changes

- `equal_area_weights` docstring currently states the transform "is not exact with these grids but
  errors reduce for higher resolution". Both halves are wrong: it does not converge (`‖T − I‖₂`
  is 0.12 at every resolution), and it can be made exact.
- Document the `slack ≥ nlat_half/16` rule and whichever of Option A / B is chosen.
- `CHANGELOG.md` entry under `## Unreleased`.
- Version bumps: `SpeedyTransforms` `0.2.1` → `0.3.0-DEV` (the `solid_angles` field is replaced, so
  this is breaking). `RingGrids` and `SpeedyWeather` are already at `+DEV` and need no further bump —
  the only change to either is a docstring and a compat bound.
- `SpeedyWeather/Project.toml` pins `SpeedyTransforms = "0.2"`, which no longer resolves against
  `0.3.0-DEV`; raised to `"0.3"`. No `SpeedyWeather` source uses the changed API (`solid_angles`
  appears nowhere outside `SpeedyTransforms`), so the bound is the only coupling.

## Known limitations

- **Dealiasing below 3 cannot be fixed by weights.** On `HEALPixGrid` the equatorial belt has only
  `2·nlat_half` longitudes, so `LegendreShortcutLinear` caps the representable order at
  `m ≤ nlat_half`. At `dealiasing = 2`, T32 on `nlat_half = 24` retains at most `m = 24`, so orders
  `m = 25…31` (7 of 32) are covered by **no ring at all** and are annihilated outright:
  `‖T − I‖₂ = 1.0`, RMS 0.27. At T64 on `nlat_half = 48` it is 15 orders of 64. That is a resolution
  deficit, not a quadrature deficit, and no weighting can repair it. This is worth stating plainly
  in the docs, since it is the regime the original report singles out as unstable.
- Weight computation adds ~2–11× to `SpectralTransform` construction time (0.0015 s at T32, 0.085 s
  at T128, 1.08 s at T256) and scales as roughly `nlat_half³`. The transform itself is unaffected.
  T512+ likely wants caching or a cheaper factorisation than the per-order SVD.
- The exact weights are specific to the `(grid, truncation, shortcut)` triple. Changing the Legendre
  shortcut after construction would invalidate them.
- Only the forward transform changes. Synthesis was always exact and stays untouched.
- At the old default `dealiasing = 3` the gain is small (T64: 3.75·10⁻³ → 3.51·10⁻³): every order
  but `m = 0` becomes exact, and `m = 0` — the one order with no slack — then dominates the error.
  The benefit is tied to the dealiasing change, not separable from it.

## Future work

- The same machinery applies unchanged to any ring grid; it could replace `clenshaw_curtis_weights`
  and reduce the ring count needed for exactness on Clenshaw grids.
- Iterative refinement (cuHPX §II-A for the method, §V for the proposal) becomes unnecessary at the
  sanctioned truncations, but would be the fallback if the dealiasing change is rejected and zero
  slack has to be supported. It costs an extra analysis **and** synthesis per iteration, i.e. it is
  paid on every transform, where the weights are a one-off precomputation.
