# Recompute the Legendre polynomials on the fly

> Status: **in progress**, opened as a draft for discussion. The design is validated numerically
> (see the accuracy table below); the `ScaledLegendre` module and its tests are done, the
> `SpectralTransform` integration and the CPU recompute path are drafted, and the GPU recursion
> kernels (section 4 and 5) are **not implemented yet** — `recompute_legendre = true` on a GPU
> architecture currently throws. The full `SpeedyTransforms` test suite has not been run against
> this branch yet.

Date of initial draft: 2026-09-01

Base revision: c3f443af8c22aca1506ea2ec0ff553278a03c4d5

## Originating prompt

> We currently precompute the Legendre polynomials in SpeedyWeather but I would like you to
> explore options to recompute them on the fly. The Legendre polynomials are big at higher
> resolution and scale badly in terms of required memory. We currently use
> AssociatedLegendrePolynomials.jl to precompute them, have a look at that code whether
> computing them on the fly on the GPU would be feasible. There's several problems we need to be
> aware of. On the GPU we'd ideally compute the polynomials in float32 but lower precision can
> easily be a problem here because of the recursion relation to compute them, see Justin
> Willmert's blogposts that are linked in the docs. In the end the recursion relation creates
> intermediate values that are really small before growing again. I'm wondering whether parts of
> it could be computed in log-space, or the sectoral modes (on the diagonal) could be precomputed
> in higher precision but all others in float32 ... I want you to be creative. If you find
> changes to AssociatedLegendrePolynomials necessary, then create a module in SpeedyTransforms.
> Don't do any SpeedyWeather unit testing or docs, just write some code that could efficiently
> compute the polynomials on the GPU, and draft how that can be integrated into
> SpeedyTransforms.jl's interface so that we have a simple precompute/recompute flag. Then create
> a draft pull request, plan with opus, execute with Sonnet

## Revision log

- 2026-09-01: initial draft, after a numerical prototype on CPU established that plain `Float32`
  recursion is unusable, that extended-exponent scaling alone is not enough, and that
  double-single (`Float32` pair) arithmetic reaches the accuracy of the current precomputation.

## Problem description

`SpectralTransform` precomputes the associated Legendre polynomials
``\lambda_\ell^m(\cos\theta_j)`` for every spherical harmonic ``(\ell, m)`` and every latitude
ring ``j`` of the northern hemisphere. That array has

```
nonzeros(spectrum) * nlat_half  ≈  (T²/2) * T  =  O(T³)
```

entries, so its memory grows with the *cube* of the truncation while every other array in the
transform grows with the square. Measured for a square spectrum with `nlat_half = T+1` in
`Float32`:

| truncation | precomputed polynomials |
|-----------:|------------------------:|
| T31        | 0.07 MB                 |
| T127       | 4.2 MB                  |
| T255       | 34 MB                   |
| T511       | 269 MB                  |
| T1023      | ~2.1 GB                 |

On GPU there is a second copy (`legendre_polynomials_transposed`, needed because the forward and
inverse kernels want opposite layouts), so the real figure is twice that. At T1023 that is more
than 4 GB of a GPU's memory spent on a table that is pure function of ``(\ell, m, j)``, before a
single prognostic variable has been allocated. This is what caps the resolution we can reach on
a single GPU.

Recomputing the polynomials inside the transform kernels would replace that `O(T³)` table with
`O(T²)` of recursion coefficients. The obstacle is numerical: the recursion is run in the number
format of the transform, which on GPU we want to be `Float32`.

## Background

### The recursion

SpeedyWeather uses the spherical-harmonic normalisation ``\lambda_\ell^m`` provided by
[AssociatedLegendrePolynomials.jl](https://github.com/jmert/AssociatedLegendrePolynomials.jl)
(`Legendre.λlm!`), documented in Justin Willmert's blog series linked from
`docs/src/spectral_transform.md`. Three recurrences are used:

```math
\lambda_0^0 = \frac{1}{\sqrt{4\pi}}, \qquad
\lambda_m^m = -\mu_m \sqrt{1-x^2}\, \lambda_{m-1}^{m-1}, \qquad
\lambda_\ell^m = \alpha_\ell^m x \lambda_{\ell-1}^m - \beta_\ell^m \lambda_{\ell-2}^m
```

with ``x = \cos\theta``. The coefficients (Willmert's well-conditioned forms, reproduced in the
new module) are

```math
\mu_\ell = 1 + \left[2\ell + \sqrt{2\ell(2\ell+1)}\right]^{-1}, \quad
\alpha_\ell^m = \sqrt{\frac{4\ell^2-1}{\ell^2-m^2}}, \quad
\beta_\ell^m = \sqrt{\frac{(2\ell+1)\left[(\ell-1)^2-m^2\right]}{(2\ell-3)(\ell^2-m^2)}}
```

Two facts make this a good fit for our GPU kernels:

1. The one-term recurrence for ``\ell = m+1``, ``\lambda_{m+1}^m = \nu_{m+1} x \lambda_m^m`` with
   ``\nu_\ell = \sqrt{2\ell+1}``, is the *same* two-term recurrence with
   ``\alpha_{m+1}^m = \nu_{m+1}`` and ``\beta_{m+1}^m = 0`` (verified analytically and
   numerically). So a single branch-free recursion with ``\lambda_{m-1}^m \equiv 0`` covers the
   whole column.
2. The recursion marches *down a column of fixed order ``m``*, which is exactly the direction the
   inverse Legendre kernel already walks: one thread per (ring ``j``, order ``m``, layer ``k``)
   dotting a lower-triangular column against the coefficients. The polynomial can be produced in
   registers as the dot product accumulates.

### Why `Float32` fails, in two separate ways

**(a) Range.** ``\lambda_m^m \sim (\sin\theta)^m``. On the ring closest to the pole of a T511
full grid that is ``\sim 10^{-1500}``, far below the `Float32` (and `Float64`) minimum. Starting
the recursion from a flushed zero gives zero for the whole column — but ``\lambda_\ell^m`` grows
back to ``O(1)`` once ``\ell \gtrsim m/\sin\theta``, so real coefficients are silently lost. In
`Float32` this bites early: at ``\theta = 30°`` the sectoral mode underflows at ``m \approx 126``
while the column becomes ``O(1)`` again at ``\ell \approx 252``, i.e. already at T255. This is
why the current code gets away with `Float32` *storage*: the values are produced by
AssociatedLegendrePolynomials.jl in `Float64` (its work arrays promote against the `Float64`
colatitudes) and only rounded on write.

**(b) Accumulated round-off.** Near the poles at low ``m`` the two solutions of the
three-term recurrence are nearly degenerate (the characteristic roots coalesce as ``x \to 1``),
so injected round-off grows like ``\ell^2`` rather than ``\sqrt{\ell}``. With ``\epsilon_{32}``
and ``\ell = 512`` that is ``\epsilon_{32}\ell^2/2 \approx 8\cdot10^{-3}`` — which is exactly
what the prototype measures.

### Prototype measurements

Max ``|\lambda_\text{recomputed} - \lambda_\text{Float64}|`` over the whole triangle and all
rings, square spectrum, `nlat_half = T+1`:

| variant                                          | T127     | T255     | T511      |
|--------------------------------------------------|---------:|---------:|----------:|
| `Float32` recursion, no scaling (naive)          | 5.3e-4   | 1.9e-2   | **5.5e35** |
| `Float32` recursion, extended-exponent scaling   | 5.3e-4   | 2.5e-3   | 9.5e-3    |
| `Float64` recursion, `Float32` coefficients      | 4.4e-5   | 5.7e-4   | 5.0e-3    |
| **double-single `Float32` (hi/lo pairs)**        | **2.0e-7** | **2.8e-7** | **4.0e-7** |
| `Float64` recursion, `Float64` coefficients      | 2.0e-7   | 2.8e-7   | 4.0e-7    |
| *current: `Float64` precompute rounded to `Float32`* | 1.2e-7 | 2.4e-7 | 2.4e-7 |

Reading the table:

- Scaling alone fixes the blow-up (5.5e35 → 9.5e-3) but not the accuracy.
- Rounding the *coefficients* to `Float32` costs as much as `Float32` arithmetic does, so both
  have to be carried in higher precision.
- Double-single arithmetic is indistinguishable from a full `Float64` recursion and lands on the
  `Float32` output-rounding floor, i.e. **as accurate as what we store today**. The residual
  3.99e-7 at T511 is the `Float32` rounding of values of magnitude ``\sqrt{2\ell+1}/\sqrt{4\pi}
  \approx 9``, plus the `Float32` sectoral start; storing the sectoral start as a hi/lo pair too
  should bring it to ~1.2e-7.

Double-single ("two-float") arithmetic represents a number as an unevaluated sum of two
`Float32`s, giving ~48 bits of mantissa using only `Float32` FMAs. It costs roughly 30 `Float32`
operations per recursion step against the 2 FMAs of a naive step — but it buys back a global
memory load per step, and the Legendre kernels are strongly memory bound (currently ~1.5 flop
per byte against a machine balance of ~13 on an A100), so the arithmetic is expected to be
close to free. On the L4s available here `Float64` runs at 1/64 rate, so double-single is the
only viable route to `Float64`-grade accuracy in the kernel.

### Why not log space

Log space is the right tool for the *sectoral* start, where the whole difficulty is dynamic
range and the value is a product, and it is wrong for the ``\ell``-recursion, which subtracts.
But an ``\exp`` of an argument of magnitude ``\sim m\ln\sin\theta \sim -10^3`` needs the argument
to ~1e-11 absolute to give `Float32`-accurate output, which is itself a double-single
computation. The extended-exponent representation below is the same idea in base ``2^{32}``,
computed exactly once on the host in `Float64`, and needs no transcendentals in the kernel.

## Summary of changes

### 1. New module `SpeedyTransforms.ScaledLegendre` (`src/scaled_legendre.jl`)

Self-contained, no GPU or SpeedyWeather dependencies, so it can be tested and reasoned about on
its own. It replaces our use of AssociatedLegendrePolynomials.jl inside the transform (that
package stays a dependency for now, as the reference in tests).

*Recursion coefficients.* `coeff_μ`, `coeff_α`, `coeff_β`, `coeff_ν` for the spherical-harmonic
normalisation, in the well-conditioned forms above, always evaluated in `Float64`.

*Double-single arithmetic.* `two_prod`, `two_sum`, `quick_two_sum`, `df_mul`, `df_add` on
`Float32` pairs, all `@inline` and FMA-based so they are GPU-safe (no branches, no allocations).

*Extended-exponent scaling.* Constants and helpers for the representation

```
true value  =  p · 2^(-SHIFT·s),    SHIFT = 32,    s ≥ 0 integer
invariant:  |p| ≤ 1 whenever s > 0
```

Rescaling rule inside the recursion: `if s > 0 && |p| > 1` then `p *= 2^-SHIFT` (exact, a power
of two) for *both* running values and `s -= 1`. The invariant guarantees that anything with
`s > 0` has a true magnitude below `2^-32 ≈ 2.3e-10`, i.e. below `Float32` eps relative to the
`O(1)` values in the same column, so it can be treated as zero. Once `s == 0` the value is
exact-as-stored and no further rescaling happens. `SHIFT = 32` leaves 94 bits of headroom above
the `Float32` minimum normal, so the running values can never underflow.

*Host-side builders.*

- `recursion_coefficients(NF, spectrum) -> (αhi, αlo, βhi, βlo)`: four vectors of length
  `nonzeros(spectrum) + 2`, indexed by the same `lm` running index as the polynomials
  themselves, zero-padded at the end so a kernel may run one recursion step past the end of a
  column without a bounds check. Computed in `Float64` and split into hi/lo `Float32` pairs.
  Diagonal entries are zero (unused).
- `sectoral_modes(NF, cos_colat, mmax) -> (hi, lo, scale)`: `(nlat_half, mmax)` arrays, `j`
  fastest so that neighbouring GPU threads (neighbouring rings at the same order) read
  contiguous memory. The ``\mu`` recursion is run once per ring in `Float64` with the scaling
  above, then split into a hi/lo pair; `scale` is `Int16` (the largest scale at T1023 is ~330).
- `legendre_column!(out, x, coefficients, lm_offset, ncol, sectoral, scale)`: reference
  implementation of the column recursion, used by the CPU path and by tests.

### 2. `AbstractLegendrePolynomials` and the `recompute` flag

```julia
abstract type AbstractLegendrePolynomials{NF} end

struct PrecomputedLegendre{NF, LTA, V} <: AbstractLegendrePolynomials{NF}
    polynomials::LTA                # (lm, j), as today
    polynomials_transposed::V       # flattened (j, lm), GPU only, empty on CPU
end

struct RecomputedLegendre{NF, V, M, MI, S} <: AbstractLegendrePolynomials{NF}
    αhi::V; αlo::V; βhi::V; βlo::V  # length nonzeros+2
    sectoral_hi::M                  # (nlat_half, mmax)
    sectoral_lo::M
    sectoral_scale::MI              # (nlat_half, mmax), Int16
    x::V                            # cos(colat), length nlat_half
    tile::S                         # scratch buffer for the forward transform, see below
    tile_orders::Vector{UnitRange{Int}}   # order blocks, host side
end
```

`SpectralTransform` loses its `legendre_polynomials` and `legendre_polynomials_transposed`
fields and gains a single `legendre::LegendrePolynomialsType` field (one type parameter replaces
one, since `LowerTriangularArrayType` was only there for the polynomials). The constructor gains

```julia
SpectralTransform(spectrum, grid; recompute_legendre::Bool = false, kwargs...)
```

`false` keeps today's behaviour bit-for-bit. Everything downstream reaches the polynomials
through the two `_legendre!` methods only, so the blast radius is `legendre.jl`,
`legendre_ka.jl`, `spectral_transform.jl` and `show.jl` (the only other reader, for the memory
report — which should now print the recompute footprint and the mode).

### 3. CPU path (`legendre.jl`)

Both CPU `_legendre!` methods already loop over rings `j` on the outside and take a
`view(legendre_polynomials.data, lm:lm_end, j)` on the inside. Introduce

```julia
legendre_ring!(scratch, legendre::AbstractLegendrePolynomials, j) -> AbstractVector
```

which returns a `view` into the precomputed array for `PrecomputedLegendre`, and for
`RecomputedLegendre` fills a `nonzeros(spectrum)`-long scratch vector with all columns for ring
`j` and returns it. Hoisted to the top of the `j` loop it is amortised over every order `m` and
every layer `k`, so the CPU cost is one recursion pass per ring per transform. The inner loops
are unchanged apart from indexing the returned vector.

On CPU the recursion runs in `Float64` when `NF == Float64` and in double-single otherwise, so
CPU and GPU agree.

### 4. GPU inverse transform — true on the fly, no buffer at all

`inverse_legendre_kernel!` already has one thread per (`j`, `m`, `k`) walking a column. Replace
the polynomial load with the recursion carried in registers:

```
p1 = (sectoral_hi[j,m], sectoral_lo[j,m]);  p2 = 0;  s = sectoral_scale[j,m]

phase 1 — climb out of the underflow region, no accumulation:
    while s > 0 and q < ncol:  advance one step, rescale

phase 2 — hot loop, no branches, two steps per iteration to keep the
          existing even/odd (north/south symmetry) split:
    a = Σ over even positions,  b = Σ over odd positions
    (even, odd) = isodd(q_first) ? (a, b) : (b, a)
```

If phase 1 never reaches `s == 0` the whole column is negligible and the thread writes zeros —
which is the correct answer, and is what the current `Float64` precompute effectively produces
too.

`α`, `β` are read by `lm`, which is *identical for every thread in a warp* (a warp holds
consecutive rows of `jm_indices`, i.e. neighbouring rings at the same order), so those loads
broadcast and hit cache; the per-thread memory traffic of the inverse Legendre transform drops
to essentially zero. `legendre_polynomials_transposed` disappears entirely in this mode.

### 5. GPU forward transform — recompute in tiles over the order `m`

The forward kernel is parallelised over the output coefficients (one thread per `lm`, summing
over rings), which is orthogonal to the recursion direction, so it cannot recompute in
registers. Instead, tile over **blocks of orders `m`**:

- A block of orders is a *contiguous* `lm` range in the lower-triangular layout, so no index
  table needs restructuring: `lm_indices` and `jm_indices` are used exactly as they are.
- Each coefficient belongs to exactly one block, so each is written exactly once — no
  accumulation across tiles, no `fill!`, and the result is bit-identical regardless of tile size.
- `tile::Matrix` of size `(nnz_tile, nlat_half)` where `nnz_tile` is chosen from a memory budget
  (default 32 MB, so ~20 orders per tile at T1023). A recompute kernel with one thread per
  (`m`, `j`) fills it; the existing `forward_legendre_kernel!` then runs over the block's `lm`
  range with a new `lm_offset` argument (`0` for the precomputed path, so no cost there).
- The recompute cost is per *transform call*, not per layer, so it is amortised over `nlayers`.

## Testing and verification

Per the prompt, no SpeedyWeather-level unit tests or docs in this PR. What goes in:

- `SpeedyTransforms/test/scaled_legendre.jl`: the recursion against `Legendre.λlm!` in `Float64`
  over several truncations and both a full and a reduced grid; the extreme-underflow rings
  explicitly; and `transform!` round-trips with `recompute_legendre = true` against
  `recompute_legendre = false` for both `Float32` and `Float64`.
- A standalone script under `SpeedyTransforms/benchmark/` (not part of the test suite) that
  reproduces the accuracy table above and times forward/inverse transforms in both modes on CPU
  and, if a GPU backend is loaded, on GPU.

## Documentation changes

None in this PR beyond docstrings, per the prompt.

## Known limitations

- The inverse kernel repeats the recursion for every layer `k`, since `k` is a thread dimension.
  For large batches it would be better to give each thread a block of layers and amortise the
  recursion; deferred (see Future work).
- The sectoral table is `O(nlat_half · mmax)`, i.e. still `O(T²)` — small next to the `O(T³)` it
  replaces, but it could be dropped entirely by recursing in `m` inside the kernel at roughly
  2× the arithmetic.
- Reduced grids with aggressive Legendre shortcuts get less benefit from the tiled forward
  transform, since the tile is sized for the full ring set.
- `MatrixSpectralTransform` is untouched and always precomputes.

## Future work

- Layer-blocked inverse kernel (amortise the recursion over several `k` per thread).
- `Int16`/packed sectoral storage, and dropping `β` in favour of
  ``\beta_\ell^m = \alpha_\ell^m / \alpha_{\ell-1}^m`` (verified identity) if coefficient memory
  ever matters.
- Making `recompute_legendre` the default above a resolution threshold once benchmarked.
