# Fitting the optimal contractive quadrature: half the `EqualArea` transform error at `‖A‖ = 1`

> Status: **completed** (implemented and unit tested; **not** yet run in a model). Replaces the
> scalar-shrink construction inside `ContractiveQuadrature` with a
> per-order fit to the *provably optimal* contractive analysis operator. Keeps `‖A‖ = 1.000000`
> exactly — the property the scheme exists for — while cutting the round-trip error to **1.75–1.92×
> below `EqualAreaQuadrature`**, where the shipped scalar shrink is 1.07–1.45× *above* it. The
> systematic round-trip damping that the contractive plan lists as an untested risk shrinks by 8.5×
> to 26× at the same time. Measured better than the scalar shrink at all 24 grid × truncation ×
> dealiasing combinations tried, and never worse.

Date of initial draft: 2026-09-01

Base revision: `5aa60fc3` (`mg/healpix-exactness`)

## Originating prompt

> First, of all, I am confused how the transform can be supposdely exact yet still inject energy
> into the system

> How is that different to how the transform with e.g. gaussian grids work. Why don't they have that
> problem?

> Is there really no way we can optimize both for contractive behaviour and exactness at the same
> time?

> the transform doesn't need to be really exact, we just want the transform error to be lower than
> with the `EqualArea` quadrature while still being contractive enough so that the model also
> doesn't loose stability because of that

## Revision log

- 2026-09-01: initial draft.
- 2026-09-01, implementation. Three things the draft did not anticipate, all found by testing and
  all changing the design (see "Guards, as built"): the ratio band and the residual check are both
  blind to the failure mode here, so acceptance is now decided by the *achieved* round trip; the
  `σmax` used for the shrink must come from an SVD of the block rather than from the Gram's
  eigenvalues; and the shortcut that keeps exact grids bit-identical has to be gated on
  non-expansiveness at roundoff, not on a fitting tolerance. Construction cost came in at the high
  end of the estimate.

## Problem description

[`ContractiveQuadrature`](2026-09-01-contractive-quadrature.md) achieves `‖A‖ = 1` by shrinking the
geometric weights with one scalar per order, `g = g⁰/σmax`. That hits the norm target but its own
"Known limitations" concede two costs:

- *"The scheme is a rescaled equal-area rule. Everything wrong with equal-area accuracy is still
  wrong; only the sign of the error is fixed."* Measured, it is slightly **worse** than equal area:
  6.73·10⁻⁴ against 5.76·10⁻⁴ at `HEALPixGrid` T128 d3.5.
- *"The gain deficit is a new, small, systematic damping … its effect on the simulated climate is
  untested."* A round trip returns 0.03% (T128) to 0.14% (T32) less than it was given.

Both are artefacts of the *construction*, not of contractivity. `ContractiveQuadrature` is not the
best contractive scheme; it is the simplest one. The stated requirement — transform error below
`EqualAreaQuadrature` while remaining contractive — is reachable, and there is a closed-form bound
saying exactly how far it can be pushed.

## Background

### The dual bound

The contractive plan's §3b establishes one half of a pair. Both are governed by `σmin(S)`, which
depends only on the grid's geometric ring areas:

| | statement | `HEALPixGrid` T64 d3.5 |
|---|---|---|
| §3b (already documented) | `A S = I` ⟹ `‖A‖ ≥ 1/σmin` | `‖A‖ ≥ 1.050914` |
| **dual (this plan)** | `‖A‖ ≤ 1` ⟹ `‖A S − I‖₂ ≥ 1 − σmin` | `‖A S − I‖₂ ≥ 0.048447` |

Proof of the dual. With `B = Λ diag(g⁰)^(1/2) = Sᵀ` and its SVD `B = U Σ Vᵀ`, put `Y = Uᵀ A V`, so
`‖Y‖ ≤ ‖A‖ ≤ 1` and `A S − I = U (Y Σ − I) Uᵀ`. Then

```
‖A S − I‖₂ = ‖Y Σ − I‖₂ ≥ max_i |Y_ii σ_i − 1| ≥ max_i (1 − σ_i) = 1 − σmin
```

using `|Y_ii| ≤ ‖Y‖ ≤ 1`. Equality holds at `Y = I`, i.e. `A = U Vᵀ` — the **orthogonal polar
factor** of `Sᵀ`, whose singular values are all exactly 1. Its round-trip operator is `A S = G^(1/2)`
with `G = Λ diag(g⁰) Λᵀ` the Gram matrix.

So the polar factor is the *unique optimum* of the problem the user posed: among all analysis
operators with `‖A‖ ≤ 1` — dense ones included, not just weight schemes — none has a smaller
round-trip error. That is the target this plan fits to.

### The frontier is a straight line

Solving the same problem at a budget `‖A‖ ≤ c` gives `A = φ(G) B` with `φ(λ) = min(1/λ, c/√λ)`, and

```
‖A S − I‖₂  =  max(0, 1 − c·σmin)
```

a straight line from `c = 1` to `c = 1/σmin`, where it reaches zero and `PerOrderQuadrature` sits.
Verified with a weights-only fit at `HEALPixGrid` T64 d3.5:

| budget `c` | achieved `‖A‖` | `‖A S − I‖₂` | frontier `1 − c·σmin` |
|---|---|---|---|
| 1.0 | 1.0 | 0.049819 | 0.048447 |
| 1.02 | 1.02 | 0.031596 | 0.029416 |
| 1.05 | 1.05 | 0.005211 | 0.000869 |
| 1.0509 | 1.0509 | 0.004397 | 1.3·10⁻⁵ |

This plan takes `c = 1`. The knob is recorded because it makes the exactness/contractivity trade a
continuum rather than the binary the two existing schemes suggest, but a budget above 1 gives back
the energy-injection channel both post-mortems diagnosed and is not proposed.

## Summary of changes

### 1. Fit the weights to `G^(1/2)` instead of scaling them by `1/σmax`

Per order `m` and parity of `l − m`, form the Gram matrix `G = Λ diag(g⁰) Λᵀ` and its square root
`G^(1/2)` from a symmetric eigendecomposition. Then solve the **reduced** `n_m`-row system — the same
one `PerOrderQuadrature` uses, with the target `1` replaced by the polar round-trip's diagonal:

```
Σ_j  g_j · λ_lm(μ_j)²  =  (G^(1/2))_{ll}          for l = m … lmax
```

taking the minimum-norm relative change `u`, `g = g⁰(1 + u)`, through the existing
`minimum_relative_change`. Finally shrink, `g ← g / max(1, σmax(A))`, which makes `‖A‖ ≤ 1` true
**by construction** rather than by the fit succeeding.

The reduction is legitimate here for the reason the exactness plan gives: the map `g ↦ T^m` factors
through `n_m` numbers and the diagonal products are triangular in `x`-degree, so matching the `n_m`
diagonal entries determines the whole achievable `T^m`. Measured, it agrees with the full
`O(n_m²)`-row system to 5–6 digits and is **14× cheaper** (0.32 s versus 4.5 s at T128).

### 2. Guards, as built

Unlike the scalar shrink, the fit has a resolution requirement — the same one `PerOrderQuadrature`
has. Below `dealiasing = 3.5` it returns weights that are negative or wildly oscillating, with
round-trip errors far worse than equal area (`HEALPixGrid` T128 d3 unguarded: weights in
`[−1.02, 2.26]`, error 0.227 against equal area's 9.4·10⁻⁴). On failure it drops to the scalar
shrink rather than to the geometric weights, so `‖A‖ ≤ 1` stays unconditional at every resolution
and dealiasing and only the accuracy gain depends on spectral slack.

**The guards `PerOrderQuadrature` uses are not sufficient here**, which the draft assumed they would
be. Three findings from implementing it:

1. **The ratio band and the residual check are both blind to this failure mode.** `[0.1, 2.0]` is
   far too wide, and the reduced system solves to machine precision while still landing far from the
   target — because the diagonal only determines the whole round trip where the target is
   *reachable*, and below the sanctioned dealiasing it is not. `OctaHEALPixGrid` T64 d3 passed both
   guards with weights in `[0.24, 1.40]` and came out 28× worse than equal area. A tuned bound on
   the weight excursion suppressed the worst of it but still let 4.5× regressions through.
   Acceptance is therefore decided on the quantity that matters: form the achieved round trip
   `Λ diag(g) Λᵀ` per parity block and keep the fit only if it beats the fallback it would
   otherwise get. That is one `n_m² × |J_m|` product per order, affordable next to the SVD, and it
   makes "never worse than the scalar shrink" a property of the construction rather than of tuning.
   The fallback's own error is free — a scalar rescaling takes `G`'s eigenvalues with it, so
   `‖G/σmax − I‖_F² = Σ_k (λ_k/σmax − 1)²`.
2. **`σmax` must come from an SVD of the block, not from the Gram's eigenvalues.** Taking
   `√λmax(B Bᵀ)` squares the condition number, which put `σmax` wrong in the 9th digit and silently
   skipped a shrink that was needed — `‖A‖` came out at `1 + 5.9·10⁻⁹`, past the `1 + 10⁻¹²` the
   testset asserts. `G^(1/2) = U Σ Uᵀ` for `B = U Σ Vᵀ`, so one SVD yields the target diagonal and
   `σmax` together, more accurately and at the same cost.
3. **The shortcut that keeps exact grids bit-identical must be gated on non-expansiveness at
   roundoff.** Gated instead on the fitted diagonal being within `QUADRATURE_RTOL = 10⁻⁸` of the
   target, it also swallowed HEALPix orders that are merely *nearly* exact — `HEALPixGrid` T64 d3.5
   at `m = 55` sits 6·10⁻⁹ off — and skipped the shrink they needed. Gated on
   `σmax ≤ 1 + QUADRATURE_NONEXPANSIVE_TOL` with the tolerance at `10⁻¹²`, the guarantee is
   manifest: an order is left alone only when it already satisfies it. Gating on the block being an
   *isometry* instead is too strict — `OctahedralGaussianGrid`'s reduced rings put it 9·10⁻⁸ from
   one — and would have changed Gaussian-grid output.

The warning also had to change. Declining an order is graceful here, so it fires only where the grid
is genuinely under-resolved; a few top orders declining on a well-resolved grid is normal operation,
and the inherited `PerOrderQuadrature` wording claimed "the grid is too coarse" while reporting an
`nlat_half` that exceeded its own stated requirement.

### 3. Where it lives

Recommended: **replace the body of `ContractiveQuadrature`**, keeping the scalar shrink as its
fallback path. One type with a fitted primary path and a guaranteed fallback is the structure
`PerOrderQuadrature` already has, and `ContractiveQuadrature` is unreleased (added on this branch),
so nothing external depends on its current numerics. The alternative — a fifth `AbstractQuadrature`
subtype — leaves a strictly dominated scheme in the API for users to pick by mistake.

Two of the current `ContractiveQuadrature` test assertions describe the scalar shrink rather than
contractivity and would have to change: weights no longer *only* shrink (they land in
`[0.956, 1.083]` of equal area on `HEALPixGrid`, `[0.979, 1.072]` on `OctaHEALPixGrid`) and are no
longer within 1%. The assertion that matters, `σmax(A) ≤ 1 + 10⁻¹²`, is unaffected.

## Testing and verification

Implemented, unit tested, and re-measured through `plot_exactness.jl` — so the numbers below are in
the **same metric as the other two plan docs** (relative L2 of a spectral → grid → spectral round
trip of a random band-limited field, `dealiasing = 3.5`). The first draft of this document quoted a
grid-space metric instead, which was not comparable with them; those numbers are superseded.

`HEALPixGrid`, Float64:

| scheme | T32 | T64 | T128 | T256 |
|---|---|---|---|---|
| equal area | 5.644·10⁻³ | 2.827·10⁻³ | 6.553·10⁻⁴ | 3.748·10⁻⁴ |
| per ring | 5.070·10⁻³ | 2.373·10⁻³ | 1.227·10⁻³ | 5.998·10⁻⁴ |
| per order | 7.02·10⁻¹⁶ | 1.14·10⁻¹⁵ | 1.75·10⁻¹⁵ | 3.04·10⁻¹⁵ |
| contractive, scalar shrink (was) | 5.83·10⁻³ | 2.91·10⁻³ | 7.39·10⁻⁴ | 4.21·10⁻⁴ |
| **contractive, fitted (now)** | **2.976·10⁻³** | **1.498·10⁻³** | **3.639·10⁻⁴** | **2.142·10⁻⁴** |

Round-trip gain `‖a₂‖/‖a‖` on the same fields:

| scheme | T32 | T64 | T128 | T256 |
|---|---|---|---|---|
| equal area | 0.99978512 | 0.99993966 | **1.00000150** | 0.99999967 |
| contractive, scalar shrink (was) | 0.99864410 | 0.99935732 | 0.99970338 | 0.99985111 |
| **contractive, fitted (now)** | **0.99984019** | **0.99994185** | **0.99998876** | **0.99999408** |

`OctaHEALPixGrid`, Float64:

| scheme | T32 | T64 | T128 | T256 |
|---|---|---|---|---|
| equal area | 4.121·10⁻³ | 2.066·10⁻³ | 4.551·10⁻⁴ | 2.634·10⁻⁴ |
| **contractive, fitted** | **2.150·10⁻³** | **1.081·10⁻³** | **2.479·10⁻⁴** | **1.910·10⁻⁴** |
| equal area, gain | 0.99983897 | 0.99995288 | **1.00000502** | **1.00000041** |
| **contractive, gain** | 0.99988899 | 0.99996046 | 0.99999566 | 0.99995297 |

Float32 agrees with Float64 to four digits at every entry, confirming the error is quadrature- and
not roundoff-dominated at these resolutions.

Against the stated requirement:

- **Error below `EqualAreaQuadrature`:** yes at every configuration — **1.75–1.90×** on
  `HEALPixGrid`, **1.38–1.92×** on `OctaHEALPixGrid`. The scalar shrink it replaces was 1.03–1.13×
  *above* equal area. Against that predecessor the gain is a clean ~2× (1.94–2.03× on `HEALPixGrid`).
- **Contractive:** `‖A‖ = 1` to 5·10⁻¹³, at every grid and truncation, and by construction rather
  than by the fit succeeding. Equal area is 1.0008–1.0018 and `PerOrderQuadrature` 1.0398–1.0563.
- **Optimal:** `‖A S − I‖₂` lands 2.3–2.9% above the `1 − σmin` floor, so at most ~3% remains
  available to any construction whatsoever, dense operators included.
- **Damping:** the round-trip gain deficit on `HEALPixGrid` falls **8.5× (T32) to 26× (T128)**
  against the scalar shrink. Equal area's gain exceeds 1 at `HEALPixGrid` T128 and at
  `OctaHEALPixGrid` T128 *and* T256 — visible injection in the shipped path — while the fitted
  scheme is strictly below 1 everywhere.
- **Never worse than the scheme it replaces**, checked directly against the staged implementation at
  all 24 {`HEALPixGrid`, `OctaHEALPixGrid`} × {T32, T64, T128} × {d2, d2.5, d3, d3.5} combinations:
  better at all 24, worse at none. Below `dealiasing = 3.5` most orders decline and the result
  converges on the scalar shrink's, which is worse than equal area — that is the pre-existing trade
  of this scheme, not a regression.

Construction cost, `HEALPixGrid`, Float64, second call (compilation excluded), on a shared node so
±30% between repeats:

| | T64 | T128 | T256 |
|---|---|---|---|
| `EqualAreaQuadrature` | 0.005 s | 0.02 s | 0.4 s |
| `PerOrderQuadrature` | 0.025 s | 0.21–0.52 s | 2.2–2.5 s |
| `ContractiveQuadrature` | 0.045 s | 0.39–0.69 s | 4.3–5.2 s |

Roughly 2× `PerOrderQuadrature`, as the draft predicted: one SVD per parity block for the target and
`σmax`, plus the acceptance check's round-trip product. The transform itself is unaffected.

Still to do:

- [x] Implemented in `SpeedyTransforms/src/quadrature.jl`, replacing `ContractiveQuadrature`'s body.
- [x] `SpeedyTransforms/test/quadrature.jl` — 46623 assertions pass with `--check-bounds=yes`,
      including a new "more accurate than equal area" testset, which is the requirement this plan
      exists to meet and which no previous test covered.
- [x] Clean construction-cost measurement (table above).
- [x] Fallback verified contractive at `dealiasing` 2 / 2.5 / 3 on both HEALPix grids.
- [x] Gaussian and Clenshaw-Curtis grids bit-identical to `EqualAreaQuadrature`.
- [x] Float32 — agrees with Float64 to four digits, no interaction.
- [x] `plot_exactness.jl` regenerated (`exactness_Float64.png`, `exactness_Float32.png`,
      `analysis_norm.png`). Its header pointed at `healpix_quadrature/`, which no longer exists;
      corrected to the script's actual location and invocation.
- [x] `CHANGELOG.md` entry under `## Unreleased`. `SpeedyTransforms` stays at `0.3.0-DEV`.
- [ ] The low-`m` declines on `OctaHEALPixGrid` T256, if T256 on that grid matters.
- [ ] **A model run.** As with the contractive plan, all evidence here is at the operator and
      transform level.

## Documentation changes

- `ContractiveQuadrature`'s docstring currently describes the scalar shrink and claims *"no failure
  mode to guard against and no resolution requirement"*. After this change that is true of the
  fallback only; the docstring must say the primary path is a fit that can decline.
- The contractive plan's "Known limitations" bullet *"The scheme is a rescaled equal-area rule"* no
  longer holds and should be struck, along with the note that it is dominated in its own class.
- The dual bound belongs in the contractive plan's §3b alongside the exactness bound — it is what
  turns "contractive costs you accuracy" into a quantified floor, and it retires the open question
  the same way §3b retired the `‖A‖`-optimisation question.

## Known limitations

- **The transform is still not exact.** ~3·10⁻⁴ at T128 against `PerOrderQuadrature`'s 10⁻¹⁵. This
  plan improves the constant, not the order. Users needing an accurate HEALPix transform for
  diagnostics, interpolation or output regridding — where analysis only sees near-band-limited
  input — should still choose `PerOrderQuadrature`.
- **A ~1.9× accuracy gain may not be worth a constructor cost increase** on its own. The stronger
  argument is the 8.5–26× smaller damping, which removes an untested systematic forcing bias from the
  scheme, rather than the error ratio.
- **`‖A‖ ≤ 1` still does not imply `‖A S‖ ≤ 1`**, unchanged from the contractive plan.
- **`OctaHEALPixGrid` at T256 is the weak point, and it is a low-`m` problem.** 1.38× rather than
  ~1.9×, with a gain deficit of 4.7·10⁻⁵ against T128's 4.3·10⁻⁶. Diagnosed: 14 of 256 orders
  decline there, and they are `m = 3 … 16` — the orders carrying most of the energy, so a handful of
  them dominates both metrics. Contrast `HEALPixGrid` T256, which declines *more* orders (49 of 256)
  with barely any effect, because they are all `m = 207 … 255`. `OctaHEALPixGrid` T128 declines
  none. This tracks the conditioning the exactness plan already reported for this grid — cond 1·10⁵
  at T256 against `HEALPixGrid`'s 3.9·10¹ — so the fit is failing where the system is worst
  conditioned, not at random. The result is still better than equal area and still contractive; the
  fix, if it is wanted, is a better-conditioned solve at low `m` rather than a change of target.
- **Operator-level evidence only.** Both schemes have `‖A‖ = 1`, so this changes nothing about the
  stability guarantee — it is an accuracy and bias improvement within a class already believed safe.
  It does not add evidence that the class *is* safe; only a model run does that.
- **The default is unchanged.** `default_quadrature` still returns `PerOrderQuadrature` for the
  HEALPix grids. Flipping it remains the human decision the contractive plan flagged.

## Future work

- Decide the default, unchanged as an open question from the contractive plan — but note this
  strengthens the case, since the contractive arm now beats equal area on accuracy as well as on
  `‖A‖` and no longer trades one for the other.
- Re-run `healpix-run.jl` on this scheme. Superrotation onset was ~day 1400, so 5 years
  discriminates; the controlled comparison wants equal-area and `PerOrderQuadrature` arms alongside.
- The budget `c` from §"The frontier is a straight line" could be exposed as a scheme parameter if
  a configuration ever wants to buy accuracy back at a bounded, *chosen* injection rate, rather than
  the unbounded one exactness implies.
