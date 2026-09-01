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
- **2026-09-01, "I also want to be able to reproduce the ring-weight from cuHPX as a baseline".**
  The weights are now built by a selectable scheme, `SpectralTransform(...; Quadrature = ...)`:
  `EqualAreaQuadrature` (the status quo), `RingQuadrature` (the classical HEALPix ring weights that
  healpy, ducc and cuHPX use) and `PerOrderQuadrature` (the default on the HEALPix grids).
  `healpix_quadrature/plot_exactness.jl` plots round-trip error against truncation for dealiasing
  1.5 … 4.0, one panel per grid × scheme. Dropping a ring's Nyquist bin stays tied to the *grid*
  rather than the scheme, so all three see the same rings and the comparison isolates the weights.
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

throughout this document `l` and `m` are the 0-based degree and order, so `lmax = truncation − 1`,
where `λ_lm` are the normalised associated Legendre polynomials, `μ_j = sin(lat_j)`, `δ_ll'` is the
Kronecker delta and `g_j` is the total solid angle of ring pair `j` (north + south), with
`Σ_j g_j = 4π`. The full derivation, and why this system is much smaller than its `O(n²)` rows
suggest, is in "Per-order quadrature weights" below; what follows here is only what is needed to
see that the problem is tractable.

Two consequences matter:

1. **The condition is linear in `g`.** Minimising `‖A∘S − I‖_F` over the weights is an ordinary
   linear least-squares problem, not a nonlinear optimisation.
2. **It has far fewer independent conditions than it looks.** Substituting `x = μ²`, the products
   `λ_lm λ_l'm` with `l + l'` even span `(1−x)^m ×` {polynomials of degree ≤ `lmax − m` in `x`}.
   So order `m` imposes only `n_m = lmax − m + 1` independent conditions on the weights, and the
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

`RingQuadrature` reproduces that baseline, so the comparison below is measured rather than argued.
Two caveats on how faithful it is: the HEALPix package's tabulated weights depend on `nside` alone
where these are fitted at the transform's own band limit, and cuHPX has no Legendre shortcut at all
so its ring sets differ from SpeedyWeather's regardless of the weights.

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

`mmax_truncation[j]` is capped at `(nlons[j] - 1) ÷ 2` rather than `nlons[j] ÷ 2`, in the
`SpectralTransform` constructor and only for grids whose weights are refitted
(`optimize_quadrature(grid)`).

Why that bin has to go. Take one ring `j` with `n = nlons[j]` longitudes and look at its Nyquist
order `m = n/2`.

1. **The real FFT stores half a spectrum.** `rfft` of a real ring of length `n` returns
   `nfreq = n÷2 + 1` coefficients `c_0 … c_{n/2}`; the rest are recovered from Hermitian symmetry,
   `C_{n−m} = conj(C_m)`.

2. **Two bins are their own conjugate partner.** For even `n`, both `m = 0` and `m = n/2` satisfy
   `n − m ≡ m (mod n)`. Hermitian symmetry therefore *forces* `C_0` and `C_{n/2}` to be real — a
   complex value in those bins does not describe a real ring at all.

3. **So the inverse FFT discards their imaginary parts.** `brfft` computes

   ```
   x_i = Re(c_0) + 2 Σ_{m=1}^{n/2−1} [Re(c_m) cos(2πmi/n) − Im(c_m) sin(2πmi/n)] + (−1)^i Re(c_{n/2})
   ```

   in which `Im(c_{n/2})` simply never appears. Transforming back, `rfft(x)` returns `n·Re(c_{n/2})`
   there. The round trip on the Nyquist bin is `c ↦ n·Re(c)` — a projection, not the identity.

4. **But the transform does not work in the ring's own frame.** Ring `j` starts at longitude
   `lon1_j`, and the transform rotates by `o = exp(i·m·lon1_j)` (`lon_offsets`) on the way out and
   by `conj(o)` on the way back. The value actually sitting in the FFT bin is `o·F_m`, where
   `F_m = Σ_l a_lm λ_lm(μ_j)` is the ring's Fourier coefficient referred to the prime meridian, and
   step 3 keeps `Re(o·F_m)`. Rotating that back:

   ```
   conj(o) · Re(o F)  =  conj(o) · (o F + conj(o) conj(F)) / 2  =  (F + conj(o)² conj(F)) / 2
   ```

5. **On a HEALPix polar-cap ring the projection lands on the wrong component.** Those rings carry a
   half-pixel offset, `lon1_j = π/n`, so at `m = n/2`

   ```
   o = exp(i · (n/2) · (π/n)) = exp(iπ/2) = i ,      conj(o)² = (−i)² = −1
   ```

   and the bracket collapses to `(F − conj F)/2 = i·Im(F)`. The ring contributes through `Im(F)`
   only — precisely the component *orthogonal* to the `Re(F)` an unoffset ring would contribute.

6. **The analysis therefore stops being complex-linear.** Split `a_lm = x_lm + i·y_lm`, so
   `Re F = Σ_l λ_lm x_lm` and `Im F = Σ_l λ_lm y_lm`. An ordinary ring feeds both; an offset
   Nyquist ring feeds only the `y` equations and an unoffset one only the `x` equations. The
   round-trip operator at that order is thus two *different* real matrices — `T_re` acting on
   `Re(a)`, `T_im` on `Im(a)` — differing by exactly that ring's rank-one term:

   ```
   T_re − T_im  =  ± g_j · λ_·m(μ_j) λ_·m(μ_j)ᵀ
   ```

7. **No choice of weight can reconcile them.** Exactness needs `T_re = T_im = I`, but their
   difference is proportional to `g_j` and vanishes only for `g_j = 0` — that is, only if the ring
   does not participate at that order. Dropping the bin *is* the fix; there is no weighting
   alternative.

8. **On HEALPix this fires on every polar-cap ring.** They have `n = 4j` longitudes while
   `LegendreShortcutLinear` retains up to `m = n÷2 = 2j`, landing exactly on Nyquist for every
   `j ≤ nside`. The equatorial belt has `n = 2·nlat_half` and Nyquist at `m = nlat_half`, which sits
   above `mmax−1` at the sanctioned dealiasing, so belt rings are never affected.

Measured directly by probing the true operator with unit basis vectors at T16: `‖T_re − T_im‖_∞` is
2.1·10⁻² at `m = 2` (ring `j = 1`), 1.5·10⁻² at `m = 4` (`j = 2`), 9.2·10⁻³ at `m = 6` (`j = 3`) —
and *exactly zero at every odd `m`*, since `n` is always even and no ring has an odd Nyquist order.
That odd/even signature is what identified the mechanism.

**Why this is gated per grid.** Dropping the bin removes a contribution, which only pays off if the
weights are refitted afterwards. Applied globally it made `OctahedralGaussianGrid` at
`dealiasing = 2` worse (T32: 1.5·10⁻¹⁰ → 2.4·10⁻¹⁰; T128: 5.7·10⁻⁷ → 8.5·10⁻⁷), so grids that keep
their geometric weights keep the bin too. Even on HEALPix it is not an improvement by itself: with
equal-area weights it makes the error marginally *worse* (5.70·10⁻³ → 5.99·10⁻³ at T32). It only
pays off together with step 2, so the two must land together.

### 2. Per-order quadrature weights

Replace the single `solid_angles[j]` with `solid_angles[m, j]`.

#### What exactness actually asks of the weights

1. **Compose the two halves of the transform.** Synthesis puts `F_m(θ_j) = Σ_l a_lm λ_lm(μ_j)` into
   ring `j`'s Fourier bin `m`; analysis reads that bin back and accumulates
   `Σ_j ΔΩ_j conj(o_{m,j}) λ_lm(μ_j) ·` (bin `m` of ring `j`). With no aliasing and no Nyquist
   pathology `rfft ∘ brfft = nlon_j · Id` on the retained bins, and the rotations cancel because
   `conj(o)·o = 1`.

2. **Collect the per-ring factor.** The `nlon_j` from the FFT pair cancels the `1/nlon_j` inside
   `ΔΩ_j = w_j · 2π/nlon_j`. What survives is a single number per ring — its total solid angle,
   north and south together:

   ```
   g_j = nlon_j · ΔΩ_j · (2 for a north/south pair, 1 for the equator) ,      Σ_j g_j = 4π
   ```

3. **Read off the round-trip operator.** Order `m` round-trips through

   ```
   â_lm = Σ_l' T^m_{l,l'} a_l'm ,        T^m_{l,l'} = Σ_j g_j λ_lm(μ_j) λ_l'm(μ_j)
   ```

   which is **linear in `g`** — so fitting the weights is ordinary linear algebra, not a nonlinear
   optimisation. Exactness is `T^m = I` at every order:

   ```
   Σ_j g_j λ_lm(μ_j) λ_l'm(μ_j) = δ_ll'          for all l, l' ≥ m
   ```

   which is just the discrete form of the continuous orthonormality `2π ∫ λ_lm λ_l'm dμ = δ_ll'`.
   It asks the quadrature to integrate products of Legendre polynomials exactly — a statement about
   the integration rule alone, not about the round trip. That is what makes this a legitimate fix
   rather than a correction tuned to band-limited input: better weights improve the analysis of
   *any* field, where a post-hoc correction of `T` would only help fields that came from synthesis.

#### Why the system is far smaller than it looks

4. **Half the conditions are automatic.** `λ_lm(−μ) = (−1)^{l+m} λ_lm(μ)`, so for `l+l'` odd the
   northern and southern contributions cancel and the condition holds for *any* symmetric weights.
   Only `l+l'` even is a real constraint.

5. **The surviving products live in a small space.** Write `λ_lm(μ) = (1−μ²)^{m/2} p_{l−m}(μ)` with
   `p_k` of degree `k` and parity `(−1)^k`. For `l+l'` even,

   ```
   λ_lm λ_l'm = (1−μ²)^m · p_{l−m}(μ) p_{l'−m}(μ)
   ```

   and `p_{l−m} p_{l'−m}` is an *even* polynomial of degree `(l−m)+(l'−m)`. Substituting `x = μ²`
   it becomes a polynomial in `x` of degree `d = ((l−m)+(l'−m))/2`, with `0 ≤ d ≤ n_m − 1` where
   `n_m = lmax − m + 1` is the number of degrees at order `m`. So every product lies in

   ```
   (1−x)^m · P_{n_m−1}(x) ,        dimension n_m
   ```

   The map `g ↦ T^m` factors through those `n_m` numbers: the system has at most `n_m` independent
   rows, not `O(n_m²)`.

6. **The diagonal products already span that space.** Take `l = l' = m + d`. Then
   `λ_{m+d,m}² = (1−x)^m p_d(μ)²`, and `p_d²` has exact degree `d` in `x`. So
   `{λ_{m+d,m}² : d = 0 … n_m−1}` is triangular in `x`-degree and therefore a basis. Satisfying the
   diagonal conditions consequently *implies* all the off-diagonal ones, and the system collapses
   to `n_m` rows — *"integrate every `|Y_lm|²` exactly"*:

   ```
   Σ_j  g_j · λ_lm(μ_j)²  =  1          for l = m … lmax
   ```

   The equivalence holds only where a solution exists. Below that — `|J_m| < n_m`, too few
   contributing rings — satisfying the diagonal rows in a least-squares sense leaves the
   off-diagonal ones unconstrained and can make the transform *worse*, which is why the code
   declines to fit rather than fitting anyway in that case.

#### Solving it

7. **Count unknowns against equations.** `n_m` equations against `|J_m|` unknowns, the rings the
   Legendre shortcut keeps at order `m`. A solution exists when `|J_m| ≥ n_m`; at `m = 0` that reads
   `nlat_half ≥ truncation`, which is the binding case and drives the choice of dealiasing below.

8. **Pick the solution closest to equal area — relatively.** With `|J_m| > n_m` the system is
   underdetermined, so solve for the relative change `u = (g − g₀)/g₀` instead of for `g`:

   ```
   (Λ² · diag(g₀)) u = 1 − Λ² g₀ ,        Λ²[l, j] = λ_lm(μ_j)²
   ```

   taking the minimum-norm `‖u‖₂` via a thin SVD. Relative rather than absolute, because `g₀` spans
   orders of magnitude between polar and equatorial rings: minimising `‖g − g₀‖` lets a polar ring
   with tiny `g₀` swing far past zero while barely moving the objective.

9. **No regularisation parameter is needed.** The reduced `n_m × |J_m|` system is well conditioned:
   `cond` 12 to 39 for `HEALPixGrid` across the four target resolutions (up to 1·10⁵ for
   `OctaHEALPixGrid` at T256). Forming the full `O(n_m²)`-row system and solving it through its
   Gram matrix instead squares the conditioning — at `nlat_half = 288`, T256 that route returned
   weights spanning `[−5.6, 6.9]` with 383 negative entries, where the reduced system lands inside
   `[0.913, 1.174]` with none.

10. **Two guards, because "solved" is not "exact".** A solution within
    `[QUADRATURE_MIN_RATIO, QUADRATURE_MAX_RATIO]` of equal area is accepted; anything outside is
    discarded in favour of the geometric weights, since an oscillating quadrature amplifies
    grid-space noise and a negative weight breaks the reading of the analysis as an integral. And
    because weights can sit *inside* those bounds while the solve still leaves a residual, the
    residual is checked separately — without that check `OctaHEALPixGrid` at T512 returned weights
    in `[0.119, 1.882]`, comfortably within bounds, for a transform that was not exact. Both paths
    report through one warning naming the `nlat_half` that would be required.

#### Implementation

Only the forward (grid → spectral) kernel reads these weights — `_legendre!` in `legendre.jl` on
CPU and `forward_legendre_kernel!` in `legendre_ka.jl` on GPU. The change is
`ΔΩ = solid_angles[j]`, hoisted out of the `m` loop, becoming `ΔΩ = solid_angles[m, j]` inside it. On GPU `m` is the fast dimension, matching `lon_offsets`
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

Construction cost (Float64, laptop CPU), against the current `SpectralTransform` constructor:
0.0015 s at T32 (+181%), 0.011 s at T64 (+405%), 0.085 s at T128 (+663%), 1.08 s at T256 (+1120%).
Absolutely small below T128 but a large *relative* increase, and ~1 s at T256 is noticeable. It
scales as roughly `nlat_half³`, dominated by one thin SVD per order. Since the reduced system is
well conditioned, a cheaper factorisation should recover most of that; T512+ may still want
caching. The transform itself is unaffected.

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

Round-trip accuracy is asserted in the pre-existing `Transform roundtrip accuracy` testset in
`SpeedyTransforms/test/spectral_transform.jl`, which this change extends rather than duplicates:

- the hardcoded `dealiasing = 3` becomes `max(3, default_dealiasing(Grid))`. Dropping the override
  outright would have *loosened* the Gaussian rows badly — at their own default of 2 the Legendre
  shortcut drops polar rings and `OctahedralGaussianGrid` reaches only 2.6·10⁻⁸ against its 10⁻¹³
  tolerance — so the Gaussian and Clenshaw rows keep exactly the footing they had and only the
  HEALPix grids move to their 3.5;
- the HEALPix tolerances drop from **10⁻² to 10⁻¹³**, the bound the Gaussian grids are held to;
- `FullHEALPixGrid` and `FullOctaHEALPixGrid` join the grid list;
- a `Float32` sweep is added, where every exact grid bottoms out at the same arithmetic noise floor
  (~2·10⁻⁷ in this metric) and so shares one tolerance.

The testset grows from 14 to 36 assertions. Its T43 rows were *already* exact at `dealiasing = 3`
(`roundup_fft` bumps 43 to `nlat_half = 48`, giving slack 5) but were being held to 10⁻², 13 orders
looser than the measurement; only the T64 rows ever needed that tolerance.

Unit tests for the weights themselves in `SpeedyTransforms/test/quadrature.jl` (23475 assertions,
~15 s), deliberately not re-testing the round trip:

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

`healpix_quadrature/transform_error.jl` reproduces these numbers against the current code, and
`healpix_quadrature/plot_exactness.jl` sweeps dealiasing 1.5 … 3.5 across all three schemes. That
sweep shows the cliff plainly — `HEALPixGrid` with per-order weights, relative L2:

| dealiasing | 1.5 | 2.0 | 2.5 | 3.0 | 3.5 |
|---|---|---|---|---|---|
| T32 | 4.1·10⁻¹ | 2.7·10⁻¹ | 8.4·10⁻² | 7.1·10⁻³ | **8.1·10⁻¹⁶** |
| T64 | 4.2·10⁻¹ | 2.5·10⁻¹ | 7.5·10⁻² | 3.5·10⁻³ | **1.4·10⁻¹⁵** |
| T128 | 4.1·10⁻¹ | 2.5·10⁻¹ | 7.0·10⁻² | 9.6·10⁻⁴ | **2.6·10⁻¹⁵** |
| T256 | 4.2·10⁻¹ | 2.6·10⁻¹ | 6.6·10⁻² | 5.0·10⁻⁴ | **4.9·10⁻¹⁵** |

Equal-area and per-ring weights show no such cliff: both sit at 10⁻⁴–10⁻³ at every dealiasing and
every truncation, and are nearly indistinguishable from each other.

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
