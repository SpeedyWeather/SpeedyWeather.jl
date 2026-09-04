# The HEALPix instability is polar ring truncation, not quadrature

> Status: **in progress**. Three experiments, two of them decisive at the operator level and one
> already decisive in the model.
>
> The five-document search for a better HEALPix quadrature was looking in the wrong place. A
> per-ring decomposition of the error an analysed *product* carries shows that **the entire
> zonal-mean error sits on the innermost 8 rings of each hemisphere, and it is not a quadrature
> error at all** — it is content lost because a cap ring with `4j` longitudes cannot represent the
> zonal orders the field actually carries there. Longitudinal folding is zero to machine precision
> on every ring; latitude quadrature is at `1e-14` once the dense operator is used. What is left is
> plain truncation loss, and it is one-signed for squared quantities.
>
> **Experiment 1 (done, decisive):** `OctaminimalGaussianGrid` — *exact* Gaussian latitude
> quadrature, but the same `4j` cap rings — diverges at **4.89 years**, inside the 4.5–6.1 year
> band of every HEALPix arm, with an identical angular-momentum trajectory. Latitude quadrature is
> exonerated as the cause.
>
> **Experiment 2 (operator level done, model runs in flight):** `HEALPixPaddedGrid`, HEALPix
> latitudes and belt and ring areas with `4j + 16` longitudes on the cap rings, drives the
> zonal-mean product error from `3.2e-6` to `1.1e-15` — a factor of 3 million, and seven orders
> below the octahedral Gaussian control — for **+3.6 %** grid points.
>
> **Experiment 3 (done):** a field integrated on the padded grid is delivered onto true HEALPix
> pixels by per-ring Fourier resampling with a relative difference of `5.0e-4`, all of it the
> discarded polar content, and the re-analysed spectral state is identical to `3.2e-16`. That is
> **58× more accurate** than the interpolator `NetCDFOutput` uses today.

Date of initial draft: 2026-09-04

Base revision: `c2a9f97f6b9d1558cceb0db8b58acc2b1ac217b8` (`mg/healpix-exactness`)

## Originating prompt

> I need new ideas and approach to get simulations stable with HEALPix grids. So far I focused on
> finding better quadrature rules based on optimizing for exactness at different truncations.
> Exactness isn't actually that central, the transforms just need to be exact enough. Stability is
> what really matters. Do you have any ideas?

> okay do one proposed experiment after the other and present me a summary in the end

> Also research if there are any other spectral models that can run with a native healpix resolution

## Revision log

- 2026-09-04: initial draft, after the per-ring diagnostic and the three experiments.

## Problem description

Five documents on this branch have tried to fix a HEALPix superrotation runaway by changing the
Legendre quadrature weights, and all five failed in the model:

| scheme | what it optimised | model outcome at T128, 4 h diffusion |
|---|---|---|
| `EqualAreaQuadrature` (status quo) | nothing | fails 4.9–6.1 y |
| `RingQuadrature` | classical HEALPix ring weights | fails 5.1 y |
| `PerOrderQuadrature` | `A∘S = I` to `1e-14` | fails 4.5–5.6 y |
| `ContractiveQuadrature` | `‖A‖ ≤ 1` exactly | fails 5.2 y (since removed) |
| `DenseQuadrature` | alias-free in latitude to `2T` | **fails 4.6 y** (confirmed this session) |

`DenseQuadrature` is the sharpest disconfirmation. It drives the latitude quadrature error on a
product to `1.2e-14`, matching the Gaussian control and beating the shipped weights by ten orders
of magnitude, and the model fails at 4.61 years — indistinguishable from everything else. The
[alias-free analysis](2026-09-03-healpix-alias-free-analysis.md) predicted exactly this and pointed
at longitudinal folding in the equatorial belt as the remaining suspect.

That suspect is wrong too. This document identifies the actual mechanism.

## Diagnosis: split the product error per ring, into truncation and folding

The dynamical core does not analyse band-limited fields, it analyses products. For a product of two
fields on one ring, write `F_m` and `G_m` for the exact per-ring Fourier coefficients of the two
factors. Three quantities are then computable per ring:

```
exact    = Σ_m F_m conj(G_m)                       the true ring mean of the product
truncated = Σ_{m ≤ M} F_m conj(G_m)                 what the ring can represent, M = ring Nyquist
sampled  = mean over the ring's nlon points of the truncated fields multiplied pointwise
```

- `truncated - exact` is **truncation loss**: orders the ring cannot hold, discarded at synthesis.
- `sampled - truncated` is **folding**: classical aliasing from sampling the product.

`healpix_quadrature/polar_truncation.jl` computes both, validated against the shipped `transform`
to `1.9e-14` on a belt ring. T128, dealiasing 3.5, atmosphere-like `(1+l)^-1` spectra, Float64.

### 1. Folding is zero. Truncation loss is everything, and it is polar.

| grid | fold, all rings | trunc loss ring 2 (8 pts) | ring 8 (32) | ring 24 | ring 48+ |
|---|---|---|---|---|---|
| `HEALPixGrid` 144 | `1e-16` | **2.8e-3** | 2.8e-7 | 8e-12 | `1e-16` |
| `OctaHEALPixGrid` 144 | `1e-16` | **1.1e-3** | 7.5e-10 | `1e-15` | `1e-16` |
| `OctahedralGaussianGrid` 96 | ≤ `1e-7` | `1e-15` | 2e-8 | 7e-8 | `1e-15` |

Two readings, both against the previous documents:

- **The belt is innocent.** The fold column is at roundoff on *every* ring of both HEALPix grids,
  including the 288-point belt rings that
  [the alias-free plan](2026-09-03-healpix-alias-free-analysis.md) nominated as the load-bearing
  problem. Its §2 proposal — more longitudes in the belt, `nside ≥ 0.75T` — addresses a channel
  that measures zero here.
- **The caps carry all of it, and the profile is brutal.** Ring 2 has 8 longitudes and an error of
  `2.8e-3`; by ring 24 (96 longitudes) it is `8e-12`. Removing the innermost 8 rings of each
  hemisphere leaves **0.4 %** of the total zonal-mean error.

### 2. Through the real transform, and why hyperdiffusion cannot touch it

Pushing the per-ring errors through the grid's own analysis, relative to the exactly analysed
product's zonal-mean coefficients:

| grid | all degrees | degree ≤ 10 |
|---|---|---|
| `HEALPixGrid` d3.5 | 2.1e-5 | **3.2e-6** |
| `OctaHEALPixGrid` d3.5 | 7.1e-6 | 1.2e-6 |
| `OctaminimalGaussianGrid` d3.5 | 8.3e-6 | 1.8e-6 |
| `OctahedralGaussianGrid` d2 | 8.8e-8 | **2.4e-8** |

A defect localised on a handful of polar rings is *not* a high-degree defect once projected onto
spherical harmonics: a polar spike has power at every degree, including `l ≤ 10` where `power = 4`
hyperdiffusion is by design near-zero. That is the property both blow-up post-mortems demanded of
any candidate mechanism, and it is why 14 diffusion settings could not fix it.

### 3. The error is one-signed for squared quantities

Repeating with the same field twice (`a·a`, the structure of a kinetic-energy or momentum-flux
term), the truncation loss is **negative on every cap ring**: ring 1 `-1.4e-3`, ring 2 `-1.7e-3`,
ring 3 `-1.5e-4`, monotonically to roundoff by ring 24. That is not noise averaging out over a long
integration; it is a systematic polar sink for a positive-definite quantity, present at every
timestep, with the same sign for ten years.

Measured directly on `u²` through the shipped `transform`, as the relative loss of each ring's
mean of `u²` (the ring's own sampled mean against the exact `Σ_m |F_m|²` from the spectral state):

| ring | nlon | lat | `HEALPixGrid` | `HEALPixPaddedGrid` (nlon) |
|---|---|---|---|---|
| 1 | 4 | 89.35° | **−1.38e−2** | +6.7e−16 (20) |
| 2 | 8 | 88.70° | **−1.36e−2** | −2.7e−15 (24) |
| 3 | 12 | 88.05° | −8.7e−5 | −1.0e−14 (28) |
| 5 | 20 | 86.75° | −1.2e−5 | −6.8e−15 (36) |
| 16 | 64 | 79.59° | −1.9e−9 | −1.1e−16 (80) |
| 32 | 128 | 69.09° | −3.3e−16 | −4.4e−16 (144) |

The two innermost HEALPix rings lose **1.4 % of their kinetic energy**, every timestep, always
negative. The padded grid is at roundoff on every ring. This is the mechanism in its most direct
form, and the sharpest single contrast in this document.

Chaining the pieces: a persistent one-signed polar defect projects onto all zonal degrees →
hyperdiffusion removes only the high ones → what survives is a smooth zonal-mean forcing that
includes the tropics → a spurious angular-momentum source. The last link is still an inference, as
it was in the previous documents; what is new is a mechanism that survives every test the earlier
candidates failed.

### 3b. Two distinct errors, and only one of them is the quadrature's

A global `<u²>` test initially looked like it contradicted the above — the padded grid's total
error came out *larger* than HEALPix's, and positive. Decomposing it settles the matter; the two
channels partly cancel on HEALPix, which is why the total is a poor diagnostic:

| grid | (A) sampling loss | (B) analysis quadrature | total |
|---|---|---|---|
| `HEALPixGrid` d3.5 | +4.81e−6 | −5.95e−6 | −1.14e−6 |
| `OctaHEALPixGrid` d3.5 | +4.03e−6 | −4.45e−6 | −4.17e−7 |
| `OctaminimalGaussianGrid` d3.5 | −1.74e−7 | **+2.2e−16** | −1.74e−7 |
| `HEALPixPaddedGrid` d3.5 | +5.85e−6 | **−1.8e−16** | +5.85e−6 |
| `OctahedralGaussianGrid` d2 | −7.61e−9 | **+1.8e−16** | −7.61e−9 |

(A) is the area-weighted mean of the sampled `u²` against the exact value; (B) is what the analysis
adds on top. Two things follow. **The padded grid's `l = 0` quadrature is exact to roundoff**,
joining the Gaussian grids, because its per-point solid angles are now consistent with its ring
counts — HEALPix's `−5.9e−6` is gone. And **(A) is dominated by the equal-area weights, not by the
polar sampling**: `OctahedralGaussianGrid`, which has both exact weights and adequate polar rings,
is the only row small in both columns.

So the padded grid fixes the *quadrature* channel and the *per-ring polar* channel (§3), while
inheriting HEALPix's equal-area weighting error in the area-weighted sum. That last one is the
channel `PerOrderQuadrature` was built for, and it is orthogonal to this work — the two are
composable, not competing.

### 3c. The padded grid and `PerOrderQuadrature` compose

Because the two channels are orthogonal, the grid fix and the weight fit can be used together. The
per-order fit works on the padded grid without modification (it is a reduced ring grid like any
other), and the combination is the best row measured on every metric:

| configuration | round trip | product L2 | zonal, l ≤ 10 |
|---|---|---|---|
| `HEALPixGrid` d3.5 equal area | 2.2e−4 | 8.06e−4 | 4.3e−5 |
| `HEALPixGrid` d3.5 per order | 2.1e−15 | 7.99e−4 | 2.8e−6 |
| `HEALPixPaddedGrid` d3.5 equal area | 1.9e−4 | 1.43e−4 | 4.6e−5 |
| **`HEALPixPaddedGrid` d3.5 per order** | **2.6e−15** | **5.01e−5** | 4.1e−6 |

The padded grid alone cuts the product error 5.6× against HEALPix; with the per-order weights it is
16× better, at an exact round trip. Note the `OctahedralGaussianGrid` control has a *worse* product
L2 (1.7e−2) than any HEALPix row while surviving ten years, which is a reminder that the total
product norm is not the quantity that predicts stability — the zonal low-degree column is, and there
the Gaussian grid leads by three orders of magnitude.

### 3d. The defect barely improves with resolution — which explains the dealiasing-7 result

The mechanism predicts that raising `dealiasing` helps only weakly: the cap ring's `nlon` grows in
step with the truncation, so the ratio of the field's order content to the ring's Nyquist limit is
roughly preserved. Measured, as the relative loss of ring-mean `u²`:

| grid | T | dealiasing | `nlat_half` | ring 1 | ring 2 |
|---|---|---|---|---|---|
| `HEALPixGrid` | 64 | 3.5 | 72 | −1.9e−3 | −3.9e−2 |
| `HEALPixGrid` | 64 | 5.5 | 108 | −2.1e−4 | −3.6e−3 |
| `HEALPixGrid` | 128 | 3.5 | 144 | −1.4e−2 | −1.4e−2 |
| `HEALPixGrid` | 128 | 5.5 | 216 | −3.2e−3 | −1.1e−3 |
| `HEALPixGrid` | 256 | 3.5 | 288 | −1.2e−2 | −4.1e−3 |
| `HEALPixGrid` | 256 | 5.5 | 432 | −2.6e−3 | −2.0e−4 |
| **`HEALPixPaddedGrid`** | 64 / 128 / 256 | 3.5 | 72 / 144 / 288 | **≤2e−15** | **≤3e−15** |

Going from T64 to T256 at fixed dealiasing does *not* reduce the loss — it is still ~1 % at T256 —
and going from dealiasing 3.5 to 5.5 buys only a factor of 4–5 while costing 2.25× the grid points.
This retro-explains the
[2026-08-31 analysis](../2026-08/2026-08-31-healpix-superrotation-blowup.md), where the composite
eigenvalue only dropped below 1 at `dealiasing = 7`: throwing resolution at a fixed-shape polar
defect is the expensive way to fix it. The padded grid removes it outright at every truncation, for
3.6 %.

### 4. Why no quadrature scheme could ever have fixed this

The loss happens at **synthesis**, when the spectral state is evaluated on a ring that cannot
represent its own orders, before the analysis quadrature runs at all. `A` cannot recover information
that `S` never put on the grid. Every scheme on this branch optimises `A`. This is the structural
reason the five-way matrix showed no ranking.

## The experiments

### Experiment 1 — negative control: exact latitudes, HEALPix cap rings

`OctaminimalGaussianGrid` is already in the codebase and is the perfect discriminator: **exact
Gaussian latitude quadrature**, and cap rings of 4, 8, 12, ... longitudes, exactly HEALPix's. If
latitude quadrature were the cause it must survive; if polar truncation is the cause it must fail.

At T128, dealiasing 3.5, `nlat_half = 144`, 4 h diffusion, matched to the failing HEALPix
configuration in every other respect:

| run | latitude quadrature | cap ring `nlon` | diverges |
|---|---|---|---|
| `OctahedralGaussianGrid` d2 (control, earlier) | exact | 20, 24, 28, … | **survives 10 y** |
| **`OctaminimalGaussianGrid` d3.5 (this experiment)** | **exact** | **4, 8, 12, …** | **4.89 y** |
| `HEALPixGrid` d3.5 equal area | inexact | 4, 8, 12, … | 5.23 y (4.89–6.07, 4 members) |
| `HEALPixGrid` d3.5 per order | exact round trip | 4, 8, 12, … | 4.57 y (4.53–5.63, 4 members) |
| `HEALPixGrid` d3.5 dense | alias-free in latitude | 4, 8, 12, … | 4.61 y |

The angular-momentum trajectory is indistinguishable from the HEALPix arms:

| relative angular momentum | y1 | y2 | y3 | y4 | y5 |
|---|---|---|---|---|---|
| `OctaminimalGaussianGrid` | 0.44 | 0.55 | 0.56 | 0.57 | NaN |
| `HEALPixGrid` equal area | 0.42 | 0.54 | 0.53 | 0.60 | 1.16 |
| `HEALPixGrid` per order | 0.45 | 0.56 | 0.56 | 0.58 | NaN |

**A grid with a provably exact latitude quadrature reproduces the HEALPix failure, at the same
time, in the same way.** Two further members are queued to bound the scatter, but a single member
landing mid-band already settles the direction: the shared property is the `4j` cap ring, and the
distinguishing property is the only one that changed.

### Experiment 2 — `HEALPixPaddedGrid`

`RingGrids/src/grids/healpix_padded.jl`, a new reduced grid. HEALPix latitudes, HEALPix belt,
HEALPix ring *areas*; the only change is that cap ring `j` carries `min(4j + 16, 2·nlat_half)`
longitudes. 16 matches the `OctahedralGaussianGrid`'s pole offset, which exists for precisely this
reason. It is deliberately **not** equal area — a padded cap ring shares HEALPix's ring area among
more points — so `get_solid_angles` falls through to the generic per-ring method rather than
HEALPix's `4π/npoints` shortcut.

Cost, at `nlat_half = 144`: 64432 points against HEALPix's 62208, **+3.6 %**.

Product error, same diagnostic as above:

| grid | zonal-mean error, all degrees | degree ≤ 10 |
|---|---|---|
| `HEALPixGrid` d3.5 | 2.1e-5 | 3.2e-6 |
| `OctahedralGaussianGrid` d2 (10-y survivor) | 8.8e-8 | 2.4e-8 |
| **`HEALPixPaddedGrid` d3.5** | **2.6e-15** | **1.1e-15** |

A factor of **3 million** against HEALPix and seven orders below the grid that survives ten years,
for 3.6 % more points. The per-ring table is flat at roundoff from ring 1 onward. A `+8` padding
was also measured and gives `~1e-9`, so the effect is graded and 16 is comfortably past the knee.

Model runs at the failing configuration (3 members) are queued; the model itself runs on the grid
(verified 10 days at T31, finite, `max|vor|` within 2 % of HEALPix's).

### Experiment 3 — delivering true HEALPix pixels from the padded grid

The padded grid is only useful if native HEALPix output survives it. Two rings at the same latitude
with different equidistant samplings are related by a Fourier resampling along the ring, exact for
everything below the shorter ring's Nyquist limit. `healpix_quadrature/padded_to_healpix.jl`
implements it (with a self-test asserting exactness to `1e-12` up- and down-sampling) and compares
three ways of producing the identical `HEALPixGrid` pixel set from one spectral state:

| delivery | relative rms difference from direct spectral synthesis onto `HEALPixGrid` |
|---|---|
| **Fourier ring resampling** | **4.96e-4** |
| `AnvilInterpolator` (what `NetCDFOutput` uses today) | 2.90e-2 |

and, re-analysing the delivered map back to spectral on the HEALPix grid:

| | all degrees | degree ≤ 20 |
|---|---|---|
| Fourier ring resampling | **3.2e-16** | **1.8e-16** |
| | | |

The `5e-4` grid-space difference is **entirely** the polar content the HEALPix pixels cannot
represent — the per-ring table shows it falls from 8.2e-2 on the 4-point ring to 1.7e-13 by ring 48
— and the spectral state is recovered to roundoff. So delivery costs nothing a HEALPix-native run
would have had, and is **58× more accurate** than today's interpolator.

## Are there other spectral models running natively on HEALPix?

Researched this session. **No.** The result is consistent across the literature:

- **cuHPX** (arXiv 2510.01785, GPU-accelerated differentiable SHTs on HEALPix) uses ring weights and
  explicitly lists solving PDEs on HEALPix as *future* work: "a promising avenue is to leverage
  cuHPX as the computational core for solving physical PDEs on the sphere, such as the shallow water
  equations on HEALPix grids." No nonlinear integration is attempted, and no stability result is
  reported.
- **Reinecke & Seljebotn's libsharp** and its successor **ducc0** support HEALPix, but ducc0's
  accurate-analysis machinery targets Clenshaw-Curtis, Fejer-1 and McEwen-Wiaux grids — the
  sub-classes where exactness is reachable. HEALPix is supported, not made exact.
- **The 2019 fast/accurate HEALPix algorithm** (arXiv 1904.10514, non-uniform FFT + double Fourier
  sphere + Slevinsky's fast SHT) improves complexity and convergence for **CMB analysis**. It does
  not modify the grid and does not consider fluid dynamics.
- **HEALPix's own software** offers ring weights, pixel weights (3.40+) and an iterative `map2alm`.
  All are analysis-side corrections — the side that experiment 1 exonerates.
- **nextGEMS / DYAMOND** use HEALPix heavily, but as an **output and analysis** grid. ICON-Sapphire
  integrates on its icosahedral grid and writes HEALPix, for chunking, hierarchical zoom levels and
  cross-model comparison. That is exactly the split experiment 3 supports.

So SpeedyWeather is, as far as this search goes, the only spectral model integrating natively on
HEALPix — which is why this failure mode has no prior art to borrow from. Notably the project's own
EGU 2026 abstract (Klöwer, Gelbrecht, Leland, Groenke, Hotta, "Variants of HEALPix grids for global
climate modelling") states that "the inexact transform with the HEALPix grids does not pose any
problems in simulations where other sources of error dominate" — the multi-year runs on this branch
contradict that, and this document supplies the specific reason.

## Summary of changes

- `RingGrids/src/grids/healpix_padded.jl`: new `HEALPixPaddedGrid`, registered and exported in
  `RingGrids.jl`, re-exported from `SpeedyWeather.jl`, default dealiasing 3.5 in
  `SpeedyTransforms/src/aliasing.jl`.
- `RingGrids/test/healpix_padded.jl`: 2352 assertions covering ring structure, the HEALPix
  latitudes, `npoints` and its inverse, ring-index partitioning, that ring weights are HEALPix's and
  the sphere integrates to `4π`, that the grid is deliberately *not* equal area, and the longitude
  offset convention. Added to `RingGrids/test/runtests.jl`.
- `healpix_quadrature/polar_truncation.jl` and `polar_truncation_T128.txt`: the per-ring
  truncation/folding decomposition, validated against the shipped `transform`.
- `healpix_quadrature/padded_to_healpix.jl`: the ring-resampling delivery path and its self-test.
- `healpix_quadrature/run_case.jl`: `OctaminimalGaussianGrid` and `HEALPixPaddedGrid` added to the
  grid table, with comments recording what each one tests.

## Testing and verification

- [x] `RingGrids` full suite passes with `--check-bounds=yes`, including the 2352 new assertions.
- [x] Per-ring diagnostic validated against the shipped `transform` to `1.9e-14`.
- [x] Ring resampler self-tested exact to `1e-12` in both directions (it had a phase-origin sign
      error on first writing, caught by that test — the uncorrected version reported a `3.4e-2`
      delivery error, which would have wrongly killed experiment 3).
- [x] `HEALPixPaddedGrid` runs `PrimitiveWetModel` for 10 days at T31, finite.
- [x] Experiment 1 member 1: `OctaminimalGaussianGrid` diverges at 4.89 y.
- [x] Per-ring `u²` loss measured through the shipped `transform`: −1.4 % on HEALPix's two
      innermost rings, roundoff on the padded grid's (§3).
- [x] The padded grid composes with `PerOrderQuadrature` (§3c).
- [x] `SpeedyTransforms` full suite passes with `--check-bounds=yes`.
- [ ] Experiment 1 members 2–3 (queued, SLURM array 2052347).
- [ ] Experiment 2, `HEALPixPaddedGrid` T128 × 3 members (queued, SLURM array 2052428). **This is
      the test that decides the whole line of argument.**
- [x] GPU: forward and inverse transform on `HEALPixPaddedGrid` agree with the CPU path to
      `1.2e-7` / `1.7e-7` in Float32 at T42 — identical to `HEALPixGrid`'s `1.2e-7` / `1.7e-7` on
      the same test, i.e. Float32 roundoff. The per-ring offset vector and the FFT planning handle
      the new ring lengths without changes.

## Known limitations

- **Experiment 2's model outcome is not in yet.** Everything about the padded grid is operator level
  so far. The operator improvement is 3 million-fold and the mechanism is now identified by a
  positive *and* a negative control, but the previous five schemes also looked good at the operator
  level. Only the 10-year runs settle it.
- **One member for experiment 1.** 4.89 y sits mid-band rather than at an extreme, which is the most
  informative single draw, but two more are queued.
- **The causal chain from a one-signed polar sink to a tropical angular-momentum source is still an
  inference**, exactly as in the two blow-up post-mortems. It is now a much better-supported one:
  the mechanism is present on every failing grid, absent on the surviving one, unreachable by
  hyperdiffusion, and unaffected by every analysis-side scheme — the four properties the runs demand.
- **The padded grid is not equal area**, so it gives up the property HEALPix users may have chosen
  HEALPix for. Experiment 3 is the answer: integrate on the padded grid, deliver equal-area HEALPix
  pixels. But a user wanting equal-area *pixels during integration* is not served by this.
- **`OctaminimalGaussianGrid`'s divergence is a new finding in its own right** and is not a HEALPix
  problem: a shipped grid fails a 10-year primitive-equation integration at the default diffusion.
  That deserves its own note regardless of what happens to HEALPix.

## Future work

- If experiment 2 survives, decide whether `HEALPixPaddedGrid` becomes the recommended grid for long
  HEALPix integrations, and wire the ring resampling into `NetCDFOutput` as a HEALPix delivery path
  (it is 58× better than the interpolator there today, independent of this whole argument).
- The same `min(4j + pad, belt)` treatment applies to `OctaHEALPixGrid` and to
  `OctaminimalGaussianGrid`, which experiment 1 shows needs it too.
- Revisit `OctahedralClenshawGrid`, whose 20-point pole rings suggest it should be fine, against the
  `λmax = 1.106` reported for it in the 2026-08-17 analysis — that number now looks like a different
  problem.
- The EGU 2026 abstract's claim about the inexact transform posing no problems should be revised in
  light of the branch's run record.
