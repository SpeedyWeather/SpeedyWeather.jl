"""Self-contained building blocks for recomputing associated Legendre polynomials on the fly
instead of precomputing and storing them.

The precomputed table `legendre_polynomials[lm, j]` used by `SpectralTransform` has
`O(nonzeros(spectrum) * nlat_half) = O(T^3)` entries, so its memory grows with the cube of the
truncation `T` while every other array in the transform grows with the square. Recomputing the
polynomials inside the transform kernels replaces that table with `O(T^2)` of recursion
coefficients, at the cost of running the recursion itself in the number format of the transform
(`Float32` on GPU), which is numerically delicate. This module provides exactly the arithmetic
needed to do that safely: recursion coefficients in the spherical-harmonic normalisation,
double-single ("two-float") arithmetic to recover close to `Float64` accuracy from `Float32`
storage, and an extended-exponent scaling scheme so the recursion does not underflow near the
poles.

The coefficient forms (the "well-conditioned" `α`, `β`, `μ` below) are Justin Willmert's, from his
Legendre.jl / AssociatedLegendrePolynomials.jl blog series on numerically stable recursions for
associated Legendre polynomials (linked from `docs/src/spectral_transform.md`), and are ported
from `AssociatedLegendrePolynomials.jl`'s `src/norm_sphere.jl`
(https://github.com/jmert/AssociatedLegendrePolynomials.jl).

This module is deliberately self-contained: no GPU imports and no dependency on RingGrids or
LowerTriangularArrays for the core math, so it can be reasoned about and tested on its own. See
`docs/dev/2026-09/recompute-legendre-polynomials.md` for the full design and the numerical
prototype that validated this approach.
"""
module ScaledLegendre

using DocStringExtensions

export coeff_μ, coeff_α, coeff_β, coeff_ν,
    two_prod, two_sum, quick_two_sum, df_mul, df_add, df_mul_f,
    scale_shift, rescale_shift,
    recursion_coefficients, sectoral_modes,
    legendre_column!, legendre_recursion_step

# ============================================================================================
# RECURSION COEFFICIENTS (spherical-harmonic normalisation), 0-based degree l and order m.
# Ported from AssociatedLegendrePolynomials.jl's src/norm_sphere.jl (Justin Willmert), which
# derives these "well-conditioned" forms in the accompanying Legendre.jl blog series to avoid
# the catastrophic cancellation of the textbook formulas when l and m are both large. Always
# evaluated in Float64, even when the recursion itself later runs in Float32: rounding these to
# Float32 before use costs as much accuracy as running the whole recursion in Float32, since the
# coefficients multiply an O(1) polynomial value every step (see the design doc's accuracy table).
# ============================================================================================

"""$(TYPEDSIGNATURES)
Coefficient of the sectoral (diagonal, l == m) recursion
`λ_m^m = -μ_m sin(θ) λ_{m-1}^{m-1}`, 0-based `l` denoting the *target* order `m` of that step.
Written as `1 + 1/(2l + sqrt(2l(2l+1)))` rather than the textbook `sqrt((2l+1)/(2l))` so that no
cancellation occurs when the two terms under the square root are added; both forms agree to
machine precision but this one is stable to evaluate directly in `Float32` as well as `Float64`."""
coeff_μ(::Type{T}, l::Integer) where {T} = (xT = convert(T, 2l); one(T) + inv(xT + sqrt(muladd(xT, xT, xT))))

"""$(TYPEDSIGNATURES)
Coefficient of `x λ_{l-1}^m` in the two-term recursion `λ_l^m = α_l^m x λ_{l-1}^m - β_l^m λ_{l-2}^m`
(0-based degree `l`, order `m`). Willmert's well-conditioned form,
`sqrt[(2l+1)(4(l-1)^2-1) / ((2l-3)(l^2-m^2))]`, factored to avoid overflow/cancellation for large
`l`. Satisfies `coeff_α(l, m) == coeff_ν(l)` exactly for `l == m + 1` (the one-term sectoral-to-
first-off-diagonal step), which is what lets a single branch-free recursion cover a whole column
with `λ_{m-1}^m ≡ 0`, see `legendre_column!`."""
function coeff_α(::Type{T}, l::Integer, m::Integer) where {T}
    lT, mT = convert(T, l), convert(T, m)
    fac1 = (2lT + 1) / ((2lT - 3) * (lT^2 - mT^2))
    fac2 = 4 * (lT - 1)^2 - 1
    return sqrt(fac1 * fac2)
end

"""$(TYPEDSIGNATURES)
Coefficient of `λ_{l-2}^m` in the two-term recursion `λ_l^m = α_l^m x λ_{l-1}^m - β_l^m λ_{l-2}^m`
(0-based degree `l`, order `m`). Willmert's well-conditioned form,
`sqrt[(2l+1)((l-1)^2-m^2) / ((2l-3)(l^2-m^2))]`. Satisfies `coeff_β(l, m) == 0` exactly for
`l == m + 1`, so the two-term recursion degenerates to the correct one-term sectoral-to-first-off-
diagonal step without a branch. Also satisfies `coeff_β(l, m) ≈ coeff_α(l, m) / coeff_α(l - 1, m)`
for `l > m + 1` (verified numerically in the tests); that identity is not used here but is noted
in the design doc as a possible way to drop `β` from the stored coefficients in the future."""
function coeff_β(::Type{T}, l::Integer, m::Integer) where {T}
    lT, mT = convert(T, l), convert(T, m)
    fac1 = (2lT + 1) / ((2lT - 3) * (lT^2 - mT^2))
    fac2 = (lT - 1)^2 - mT^2
    return sqrt(fac1 * fac2)
end

"""$(TYPEDSIGNATURES)
Coefficient of the one-term recursion `λ_{m+1}^m = ν_{m+1} x λ_m^m` from the sectoral mode to the
first off-diagonal, 0-based degree `l`. Equal to `coeff_α(l, l - 1)`; provided separately because
that identity is exactly what makes the general two-term recursion below valid at the top of every
column, and it is cheaper to state directly than to evaluate `coeff_α` at the boundary."""
coeff_ν(::Type{T}, l::Integer) where {T} = sqrt(convert(T, 2l + 1))

# ============================================================================================
# LOWER-TRIANGULAR INDEXING
# Local copy of LowerTriangularArrays.lm2i's convention (1-based l, m; running index increases
# down a column first) so this module has no dependency on LowerTriangularArrays. Must be kept
# in sync with that definition by construction/tests, not by sharing code, since the whole point
# of this module is to be usable without pulling in the rest of the package.
# ============================================================================================

"""$(TYPEDSIGNATURES)
Running index into the lower-triangular storage layout for 1-based degree `l`, order `m`, and
number of degrees `lmax`. Matches `LowerTriangularArrays.lm2i`'s convention exactly (same formula,
kept as a local copy so this module stays free of a LowerTriangularArrays dependency)."""
@inline lm2i(l::Integer, m::Integer, lmax::Integer) = l + (m - 1) * lmax - m * (m - 1) ÷ 2

"""$(TYPEDSIGNATURES)
Number of nonzero (on-or-below-diagonal) entries of a lower-triangular spectrum with `lmax`
degrees and `mmax` orders. Matches `LowerTriangularArrays.nonzeros(lmax, mmax)`."""
@inline count_nonzeros(lmax::Integer, mmax::Integer) = lmax * mmax - (mmax - 1) * mmax ÷ 2

# ============================================================================================
# DOUBLE-SINGLE ("TWO-FLOAT") ARITHMETIC ON Float32 PAIRS
# Represents a real number as an unevaluated sum hi + lo of two Float32s (|lo| much-less-than
# |hi|), giving close to Float64 accuracy (~48 bits of mantissa) from Float32 hardware. All
# routines are @inline, FMA-based and branch-free (no allocation, no error paths), so they are
# safe to call from GPU kernels. They take and return plain Float32 tuples/values rather than a
# wrapping struct for the same reason: this is meant to inline into kernel registers, not to
# allocate a boxed number. Standard algorithms (Dekker 1971 / Shewchuk 1997), specialised to use
# `fma` instead of a Veltkamp split since target hardware has a fused multiply-add.
# ============================================================================================

"""$(TYPEDSIGNATURES)
Error-free product `a*b = hi + lo` in `Float32`: `hi` is the correctly rounded product, `lo` is
the exact rounding error, computed as `fma(a, b, -hi)` (exact because `hi` is the Float32 nearest
to the true product, so `a*b - hi` is representable and the fused multiply-add computes it without
an intermediate rounding step)."""
@inline function two_prod(a::Float32, b::Float32)
    p = a * b
    return p, fma(a, b, -p)
end

"""$(TYPEDSIGNATURES)
Error-free sum `a+b = hi + lo` in `Float32` for general `a`, `b` (unlike `quick_two_sum`, does not
require `|a| >= |b|`). Knuth/Møller's two-sum algorithm: 6 flops, no comparisons needed to decide
which operand is larger, so it is branch-free on top of being exact."""
@inline function two_sum(a::Float32, b::Float32)
    s = a + b
    bb = s - a
    return s, (a - (s - bb)) + (b - bb)
end

"""$(TYPEDSIGNATURES)
Error-free sum `a+b = hi + lo` in `Float32`, assuming `|a| >= |b|` (Dekker's quick-two-sum, 3
flops instead of `two_sum`'s 6). Used where that ordering is already guaranteed by construction,
e.g. combining a double-single leading term with a small correction."""
@inline function quick_two_sum(a::Float32, b::Float32)
    s = a + b
    return s, b - (s - a)
end

"""$(TYPEDSIGNATURES)
Double-single times plain `Float32`: `(ah, al) * b`, returned as a new double-single `(hi, lo)`."""
@inline function df_mul_f(ah::Float32, al::Float32, b::Float32)
    p, e = two_prod(ah, b)
    e = fma(al, b, e)
    return quick_two_sum(p, e)
end

"""$(TYPEDSIGNATURES)
Double-single times double-single: `(ah, al) * (bh, bl)`, returned as a new double-single
`(hi, lo)`. Drops the `al*bl` cross term (it is below `Float32` eps relative to the other three
products), which is the standard double-single multiplication truncation."""
@inline function df_mul(ah::Float32, al::Float32, bh::Float32, bl::Float32)
    p, e = two_prod(ah, bh)
    e = fma(ah, bl, e)
    e = fma(al, bh, e)
    return quick_two_sum(p, e)
end

"""$(TYPEDSIGNATURES)
Double-single plus double-single: `(ah, al) + (bh, bl)`, returned as a new double-single
`(hi, lo)`. The low-order correction `al + bl` is added as plain `Float32` (its own rounding error
is below the precision the result can represent), then folded back in with a final
`quick_two_sum`."""
@inline function df_add(ah::Float32, al::Float32, bh::Float32, bl::Float32)
    s, e = two_sum(ah, bh)
    e += al + bl
    return quick_two_sum(s, e)
end

# ============================================================================================
# EXTENDED-EXPONENT SCALING
# λ_m^m ~ sin(θ)^m underflows Float32 (and Float64) range near the poles at high m, long before
# λ_l^m has decayed to a value that is actually negligible (it grows back to O(1) once
# l ≳ m/sin(θ)). Representing the running recursion state as p * 2^(-scale_shift(NF)*s) with an
# integer s tracked alongside it lets the recursion carry values through that underflow region
# without ever flushing to zero, using only exact power-of-two rescalings (no transcendentals, no
# rounding of the mantissa).
# ============================================================================================

"""$(TYPEDSIGNATURES)
Exponent step, in bits, of one unit of the extended-exponent scale factor `s` in number format
`NF`: the true value of a scaled running state is `p * 2^(-scale_shift(NF) * s)`, with the
invariant `|p| <= 1` whenever `s > 0`. That invariant means anything still carrying `s > 0` has
true magnitude below `2^-scale_shift(NF)`, so it is safe to treat as exactly zero for output
purposes. Two constraints fix the value, and the shift has to sit between them:

- it must be **large enough** that `2^-scale_shift(NF)` is below `eps(NF)` relative to the `O(1)`
  values that appear later in the same column, so that flushing those values to zero loses
  nothing `NF` could represent anyway. Hence at least `precision(NF)` bits.
- it must be **small enough** that a running mantissa just under the rescale threshold of `1`
  still has headroom above `floatmin(NF)` after one more step down in `s`, so that the mantissa
  itself can never underflow between rescalings however many recursion steps pass before the next
  `|p| > 1` trigger.

`precision(NF) + 8` satisfies the first with margin and is capped at a third of the exponent range
for the second. This gives 32 bits for `Float32` (2.3e-10 neglected, 94 bits of headroom above
`2^-126`) and 61 for `Float64` (4.3e-19 neglected, far under `eps(Float64)`). Formats with a very
narrow exponent range relative to their mantissa — `Float16` in particular — cannot satisfy both
at once and are not supported by the recomputed transform."""
scale_shift(::Type{NF}) where {NF <: AbstractFloat} =
    min(precision(NF) + 8, -exponent(floatmin(NF)) ÷ 3)

# spell the two supported formats out so the shift is a literal inside GPU kernels
scale_shift(::Type{Float32}) = 32
scale_shift(::Type{Float64}) = 61

"""$(TYPEDSIGNATURES)
The rescaling factor `2^(-scale_shift(NF))` in number format `NF`, i.e. what a running recursion
state `p` is multiplied by (together with decrementing its scale `s`) whenever `|p|` grows back
past `1` while `s > 0`. Exact in any binary floating-point format since it is a power of two."""
rescale_shift(::Type{NF}) where {NF} = NF(ldexp(1.0, -scale_shift(NF)))

# ============================================================================================
# SPLITTING A Float64 COEFFICIENT INTO A NUMBER-FORMAT hi/lo PAIR
# ============================================================================================

"""$(TYPEDSIGNATURES)
Split a `Float64` value `v` into an `NF` double-single pair `(hi, lo)` with `hi = NF(v)` and
`lo = NF(v - Float64(hi))` the exactly-representable remainder. For `NF == Float64` this makes
`hi == v` exactly and `lo == 0`: the pair degenerates to a single value but the interface stays
uniform across `NF`, so callers do not need to special-case `Float64`."""
@inline function split_two_float(::Type{NF}, v::Float64) where {NF}
    hi = NF(v)
    lo = NF(v - Float64(hi))
    return hi, lo
end

# ============================================================================================
# HOST-SIDE BUILDERS
# ============================================================================================

"""$(TYPEDSIGNATURES)
Build the `α`, `β` recursion coefficient tables for a lower-triangular spectrum with `lmax`
degrees and `mmax` orders, split into number-format `NF` double-single pairs. Returns
`(αhi, αlo, βhi, βlo)`, four `Vector{NF}` of length `count_nonzeros(lmax, mmax) + 2`: the running
index `lm` (`LowerTriangularArrays.lm2i` convention, 1-based `l`, `m`) addresses them directly,
and the `+ 2` zero-padding at the end lets a column recursion take one step past the end of the
last column (`m == mmax`, `l == lmax`) without a bounds check, since the write of that
one-past-the-end value is always discarded. Diagonal entries (`l == m`, the sectoral modes
themselves) are left zero and are never read: the recursion starts from the precomputed sectoral
value instead, see `sectoral_modes`. Coefficients are evaluated in `Float64` (`coeff_α`,
`coeff_β`), so accuracy does not depend on `NF`; only the final split does."""
function recursion_coefficients(::Type{NF}, lmax::Integer, mmax::Integer) where {NF}
    n = count_nonzeros(lmax, mmax)
    αhi, αlo = zeros(NF, n + 2), zeros(NF, n + 2)
    βhi, βlo = zeros(NF, n + 2), zeros(NF, n + 2)
    for m in 1:mmax, l in m:lmax
        l == m && continue                              # diagonal: unused, keep zero
        i = lm2i(l, m, lmax)
        a = coeff_α(Float64, l - 1, m - 1)               # 0-based degree/order in the coefficient
        b = coeff_β(Float64, l - 1, m - 1)               # formulas, 1-based in the lm2i index
        αhi[i], αlo[i] = split_two_float(NF, a)
        βhi[i], βlo[i] = split_two_float(NF, b)
    end
    return αhi, αlo, βhi, βlo
end

"""$(TYPEDSIGNATURES)
Build the sectoral (diagonal, `l == m`) starting values `λ_m^m(cos_colat[j])` for every ring `j`
and order `m`, with the extended-exponent scaling of `scale_shift(NF)` applied. Returns
`(hi, lo, scale)`: `hi`, `lo` are `(nlat_half, mmax)` `Matrix{NF}` double-single pairs and `scale`
is the matching `Matrix{Int16}` of scale factors `s` (`Int16` because even at T1023 the largest
scale needed is of order a few hundred, far under the `Int16` range). All three are laid out with
`j` (ring) as the fastest (first) dimension rather than `m` (order), even though the recursion
below runs one ring at a time: on the GPU, neighbouring threads process neighbouring rings at the
same order, so laying rings out contiguously in memory turns those reads into coalesced accesses
instead of strided ones. The `μ` recursion that produces the sectoral values is run once per ring
in `Float64` (so its own accuracy does not depend on `NF`) and only split into an `NF`
double-single pair at the end."""
function sectoral_modes(::Type{NF}, cos_colat::AbstractVector, mmax::Integer) where {NF}
    nlat_half = length(cos_colat)
    hi = zeros(NF, nlat_half, mmax)
    lo = zeros(NF, nlat_half, mmax)
    scale = zeros(Int16, nlat_half, mmax)
    SMALL = ldexp(1.0, -scale_shift(NF))
    BIG = ldexp(1.0, scale_shift(NF))
    for j in 1:nlat_half
        x = Float64(cos_colat[j])
        y = sqrt(max(0.0, 1 - x^2))                       # sin(colat), clamped against roundoff
        p = inv(sqrt(4 * Float64(π)))                      # λ_0^0
        s = 0
        hi[j, 1], lo[j, 1] = split_two_float(NF, p)
        scale[j, 1] = 0
        for m in 1:(mmax - 1)                              # 0-based target order m, μ_m: (m-1,m-1) -> (m,m)
            p = -coeff_μ(Float64, m) * y * p
            while p != 0 && abs(p) < SMALL
                p *= BIG
                s += 1
            end
            hi[j, m + 1], lo[j, m + 1] = split_two_float(NF, p)
            scale[j, m + 1] = Int16(s)
        end
    end
    return hi, lo, scale
end

# ============================================================================================
# COLUMN RECURSION
# ============================================================================================

"""$(TYPEDSIGNATURES)
One step of the double-single column recursion: given the recursion coefficient at this step as
double-single pairs `(αh, αl)`, `(βh, βl)`, the evaluation point as a double-single pair
`(xh, xl)`, and the running state `(p1h, p1l)` (the current value, `λ_{l-1}^m`) and
`(p2h, p2l)` (the previous value, `λ_{l-2}^m`), returns the new value `λ_l^m = α x λ_{l-1}^m -
β λ_{l-2}^m` as a double-single pair. `@inline`d and allocation-free so it inlines into registers;
this is the piece meant to be reused verbatim inside the GPU recursion kernels (see the design
doc), where the surrounding loop becomes a `while`/unrolled loop over threads instead of the plain
`for` loop `legendre_column!` uses on the CPU."""
@inline function legendre_recursion_step(
        αh::Float32, αl::Float32, βh::Float32, βl::Float32,
        xh::Float32, xl::Float32,
        p1h::Float32, p1l::Float32, p2h::Float32, p2l::Float32,
    )
    ah, al = df_mul(αh, αl, xh, xl)          # α*x
    th, tl = df_mul(ah, al, p1h, p1l)        # α*x*λ_{l-1}^m
    uh, ul = df_mul(βh, βl, p2h, p2l)        # β*λ_{l-2}^m
    return df_add(th, tl, -uh, -ul)
end

"""$(TYPEDSIGNATURES)
Reference/CPU implementation of one Legendre column recursion (fixed order `m`, running over
degree `l`), writing the *true* (unscaled) polynomial values into `out[1:ncolumn]`. Dispatches its
arithmetic on `eltype(out)`: double-single (`Float32` hi/lo pairs, via `legendre_recursion_step`)
for `Float32` output, plain `Float64` arithmetic for `Float64` output — carrying the coefficients
and the running state in `Float32` alone is what produces the ~1e-2 errors documented in the
design doc, so this is not an optional refinement.

Arguments: `x` the evaluation point `cos(colat)` (as `Float64`, split internally); `αhi, αlo, βhi,
βlo` the coefficient tables from `recursion_coefficients`; `lm_offset` the running index just
before this column's first entry (so degree `l = m` is written to `out[1]`, and coefficient lookup
for step `q` is `lm_offset + q + 1`, i.e. the coefficient belonging to the *next* value produced);
`ncolumn` the number of degrees in this column (`lmax - m + 1`); `sectoral_hi, sectoral_lo, scale`
the starting `λ_m^m` double-single pair and its extended-exponent scale from `sectoral_modes`.

Writes zero wherever `scale` has not yet reached `0`: per the `scale_shift` invariant, that means
the true value is below `Float32` epsilon relative to the `O(1)` values that appear later in the
same column, so it is indistinguishable from zero at this working precision anyway."""
function legendre_column!(
        out::AbstractVector{Float32}, x::Real,
        αhi::AbstractVector{Float32}, αlo::AbstractVector{Float32},
        βhi::AbstractVector{Float32}, βlo::AbstractVector{Float32},
        lm_offset::Integer, ncolumn::Integer,
        sectoral_hi::Float32, sectoral_lo::Float32, scale::Integer,
    )
    x64 = Float64(x)
    xh = Float32(x64)
    xl = Float32(x64 - Float64(xh))
    p1h, p1l = sectoral_hi, sectoral_lo               # λ_{l-1}^m, starting at the sectoral mode
    p2h, p2l = 0.0f0, 0.0f0                            # λ_{l-2}^m ≡ 0 above the sectoral mode
    s = Int(scale)
    SMALL = rescale_shift(Float32)
    @inbounds for q in 1:ncolumn
        out[q] = s == 0 ? p1h + p1l : 0.0f0
        i = lm_offset + q + 1
        ph, pl = legendre_recursion_step(αhi[i], αlo[i], βhi[i], βlo[i], xh, xl, p1h, p1l, p2h, p2l)
        p2h, p2l = p1h, p1l
        p1h, p1l = ph, pl
        if s > 0 && abs(p1h) > 1.0f0                   # rescale: exact, both running values and s
            p1h *= SMALL
            p1l *= SMALL
            p2h *= SMALL
            p2l *= SMALL
            s -= 1
        end
    end
    return out
end

function legendre_column!(
        out::AbstractVector{Float64}, x::Real,
        αhi::AbstractVector{Float64}, αlo::AbstractVector{Float64},
        βhi::AbstractVector{Float64}, βlo::AbstractVector{Float64},
        lm_offset::Integer, ncolumn::Integer,
        sectoral_hi::Float64, sectoral_lo::Float64, scale::Integer,
    )
    x64 = Float64(x)
    p1 = sectoral_hi + sectoral_lo                     # lo is exactly 0 for Float64, kept for a
    p2 = 0.0                                            # uniform interface with the Float32 method
    s = Int(scale)
    SMALL = rescale_shift(Float64)
    @inbounds for q in 1:ncolumn
        out[q] = s == 0 ? p1 : 0.0
        i = lm_offset + q + 1
        a = αhi[i] + αlo[i]
        b = βhi[i] + βlo[i]
        p = a * x64 * p1 - b * p2
        p2, p1 = p1, p
        if s > 0 && abs(p1) > 1.0
            p1 *= SMALL
            p2 *= SMALL
            s -= 1
        end
    end
    return out
end

end # module
