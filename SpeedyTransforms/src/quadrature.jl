# Relative singular value below which a direction of the exactness system is treated as noise.
const QUADRATURE_RTOL = 1.0e-8

# Bounds on the optimised weights relative to the equal-area weights. A well-resolved
# configuration stays within [0.9, 1.2]; an under-resolved one overshoots by orders of magnitude
const QUADRATURE_MIN_RATIO = 0.1
const QUADRATURE_MAX_RATIO = 2.0

"""Quadrature scheme used by the forward (grid to spectral) transform to weight each latitude
ring. 

Implementations are `EqualAreaQuadrature`, `RingQuadrature` and `PerOrderQuadrature`."""
abstract type AbstractQuadrature end

"""$(TYPEDSIGNATURES)
The grid's geometric solid angles, used unchanged at every order `m`. Exact for grids whose
latitudes carry a quadrature rule (Gaussian, Clenshaw-Curtis) and a plain Riemann sum for the
equal-area HEALPix grids, where it is what SpeedyWeather used before per-order weights."""
struct EqualAreaQuadrature <: AbstractQuadrature end

"""$(TYPEDSIGNATURES)
Classical HEALPix ring weights: one weight per ring, the same at every order `m`, fixed by
requiring that the quadrature integrate every band-limited function exactly,

```
Σ_j g_j λ_l0(μ_j) = ∫ λ_l0 dΩ = 2√π δ_l0        for l = 0 … lmax
```

with `g_j` the total solid angle of ring pair `j`. Odd `l` vanish by north/south symmetry, so this
is one condition per even degree — the `m = 0`, `l' = 0` column of the full orthogonality system
that `PerOrderQuadrature` solves, and strictly weaker than it.

This is the scheme shipped as FITS tables by the HEALPix package and used by healpy, ducc and
cuHPX (arXiv:2510.01785, §IV-A), provided here to reproduce that baseline. The defining conditions
above are the classical ones; the system is underdetermined and the particular solution taken here
(minimum relative change from equal area) is this implementation's choice."""
struct RingQuadrature <: AbstractQuadrature end

"""$(TYPEDSIGNATURES)
Weights fitted separately for every order `m`, which makes the transform exact on the HEALPix
grids given enough spectral slack. See `quadrature_weights` for the derivation."""
struct PerOrderQuadrature <: AbstractQuadrature end

"""$(TYPEDSIGNATURES)
Default quadrature scheme for a grid. Gaussian and Clenshaw-Curtis grids already carry an exact
quadrature rule in their latitudes and keep their geometric solid angles; the HEALPix grids are
equal-area discretisations whose weights are a Riemann sum rather than a quadrature rule, and get
per-order weights."""
default_quadrature(Grid::Type{<:AbstractGrid}) = EqualAreaQuadrature
default_quadrature(
    ::Type{
        <:Union{
            HEALPixGrid, OctaHEALPixGrid,
            FullHEALPixGrid, FullOctaHEALPixGrid,
        },
    }
) = PerOrderQuadrature

default_quadrature(grid::AbstractGrid) = default_quadrature(typeof(grid))

"""$(TYPEDSIGNATURES)
Whether to stop the Legendre loop one order below a ring's Nyquist frequency `m = nlon_j/2`.
There the real FFT carries only one of the two components, so the ring contributes to only one of
the real/imaginary parts of the coefficients and no single weight is consistent for both, which
makes an exact transform unreachable while the bin is retained.

This follows the grid, not the quadrature scheme: only the HEALPix grids have rings short enough
for the Legendre shortcut to reach their Nyquist order, and they are exactly the grids whose
weights are fitted. Tying it to the grid also keeps every scheme on the same set of rings, so
comparing `RingQuadrature` against `PerOrderQuadrature` isolates the weights instead of changing
the ring sets at the same time."""
drop_nyquist(grid::AbstractGrid) = default_quadrature(grid) === PerOrderQuadrature

"""$(TYPEDSIGNATURES)
Total solid angle `g_j` covered by ring pair `j` (north and south together), so `Σ_j g_j = 4π`.
The equator, if the grid has one, is a single ring and enters once."""
function ring_solid_angles(nlons, nlat::Integer, nlat_half::Integer, solid_angles)
    return [((nlat - j + 1 == j) ? 1 : 2) * nlons[j] * Float64(solid_angles[j]) for j in 1:nlat_half]
end

"""$(TYPEDSIGNATURES)
Minimum-norm solution `u` of `A diag(g₀) u = b - A g₀`, the smallest *relative* change from the
weights `g₀` that satisfies the conditions `A g = b`. Relative rather than absolute because `g₀`
spans orders of magnitude between polar and equatorial rings, and minimising `‖g - g₀‖` would let
a polar ring swing past zero while barely moving the objective. Returns `u` and the residual."""
function minimum_relative_change(A::AbstractMatrix, b::AbstractVector, g₀::AbstractVector)
    Ã = A * LinearAlgebra.Diagonal(g₀)
    b̃ = b - A * g₀
    F = LinearAlgebra.svd(Ã)
    keep = 1:count(>(QUADRATURE_RTOL * F.S[1]), F.S)
    u = F.V[:, keep] * ((F.U[:, keep]' * b̃) ./ F.S[keep])
    return u, LinearAlgebra.norm(Ã * u - b̃) / max(LinearAlgebra.norm(b̃), eps())
end

"""$(TYPEDSIGNATURES)
Quadrature weights `ΔΩ[m, j]` = sinθ Δθ Δϕ for order `m` (1-based) and northern latitude ring `j`,
as read by the forward (grid to spectral) Legendre transform. Dispatches on `Quadrature`, one of
`EqualAreaQuadrature`, `RingQuadrature` or `PerOrderQuadrature`."""
function quadrature_weights end

# EQUAL AREA: the grid's geometric solid angles, order-independent
function quadrature_weights(
        ::Type{EqualAreaQuadrature},
        NF::Type{<:Real},
        grid::AbstractGrid,
        spectrum::AbstractSpectrum,
        legendre_polynomials::LowerTriangularArray,
        mmax_truncation::AbstractVector{<:Integer},
        nlons::AbstractVector{<:Integer},
        nlat::Integer,
        solid_angles::AbstractVector,
    )
    return NF[solid_angles[j] for _ in 1:spectrum.mmax, j in 1:grid.nlat_half]
end

"""$(TYPEDSIGNATURES)
Classical HEALPix ring weights, see `RingQuadrature`. Solves

```
Σ_j g_j λ_l0(μ_j) = 2√π δ_l0          for even l = 0 … lmax
```

for the minimum relative change from the equal-area weights. Odd degrees are omitted because they
vanish for any north/south symmetric weighting, and the `l = 0` row is already satisfied by the
equal-area weights (it is just `Σ_j g_j = 4π`), so the global mean stays exactly conserved.

Unlike `PerOrderQuadrature` this is not guarded against straying far from equal area: it is a fixed
reference scheme rather than a fit we are free to decline, so out-of-range weights are reported but
still returned."""
function quadrature_weights(
        ::Type{RingQuadrature},
        NF::Type{<:Real},
        grid::AbstractGrid,
        spectrum::AbstractSpectrum,
        legendre_polynomials::LowerTriangularArray,
        mmax_truncation::AbstractVector{<:Integer},
        nlons::AbstractVector{<:Integer},
        nlat::Integer,
        solid_angles::AbstractVector,
    )
    (; lmax, mmax) = spectrum                   # 1-based
    (; nlat_half) = grid

    g₀ = ring_solid_angles(nlons, nlat, nlat_half, solid_angles)
    λ = Array(legendre_polynomials.data)        # (nonzeros, nlat_half), on the CPU for the fit

    # λ_l0 at every ring, even degrees only (odd ones vanish by symmetry)
    λ₀ = Float64.(λ[LowerTriangularArrays.get_lm_range(1, lmax - 1), :])
    A = λ₀[1:2:end, :]
    b = zeros(size(A, 1)); b[1] = 2sqrt(π)      # ∫λ_00 dΩ = 2√π, all higher degrees integrate to 0

    u, _ = minimum_relative_change(A, b, g₀)

    if !all(isfinite, u) || !all(r -> QUADRATURE_MIN_RATIO <= r <= QUADRATURE_MAX_RATIO, 1 .+ u)
        @warn "Ring quadrature weights stray outside " *
            "[$QUADRATURE_MIN_RATIO, $QUADRATURE_MAX_RATIO] of the equal-area weights " *
            "(range $(round(minimum(1 .+ u), digits = 3)) … $(round(maximum(1 .+ u), digits = 3))) " *
            "on $(nonparametric_type(typeof(grid))) at truncation $lmax, nlat_half $nlat_half. " *
            "They are returned unchanged; RingQuadrature is a fixed reference scheme." maxlog = 3
    end

    return NF[solid_angles[j] * (1 + u[j]) for _ in 1:mmax, j in 1:nlat_half]
end

"""$(TYPEDSIGNATURES)
Per-order quadrature weights, see `PerOrderQuadrature`. The transform is exact when, for every
order `m`,

```
Σ_j g_j λ_lm(μ_j) λ_l'm(μ_j) = δ_ll'      for all l, l' ≥ m
```

with `λ_lm` the normalised associated Legendre polynomials, `μ_j = sin(lat_j)`, `g_j` the total
solid angle of ring pair `j` (north + south, so `Σ_j g_j = 4π`) and `δ_ll'` the Kronecker delta,
1 for `l == l'` and 0 otherwise. That is, the quadrature has to reproduce the exact orthonormality
of the spherical harmonics: 1 on the diagonal (right normalisation, no gain error) and exactly 0
off it (no leakage between degrees at the same order `m`) — which is what makes analysis ∘
synthesis the identity.

That system looks quadratic in the number of degrees, but it is not: substituting `x = μ²`, the
products `λ_lm λ_l'm` with `l + l'` even span `(1-x)^m` times the polynomials of degree
`≤ lmax - m` in `x`, and the diagonal products `λ_lm²` are already a triangular basis of that
space. So the off-diagonal conditions follow from the diagonal ones and the whole system collapses
to one row per degree — *integrate every `|Y_lm|²` exactly*:

```
Σ_j g_j λ_lm(μ_j)² = 1                    for l = m … lmax
```

which is solved here for the minimum relative change from the equal-area weights.

A solution exists only if enough rings contribute at order `m`. At `m = 0` that means
`nlat_half ≥ truncation`, and with equality the system is square: its unique solution is an
interpolatory quadrature at maximal degree, which oscillates the way Newton-Cotes does. A margin of
`nlat_half/16` is enough for `HEALPixGrid`, `OctaHEALPixGrid` needs a little more, and the
`nlat_half/9` that `default_dealiasing` provides keeps every weight positive and within 20% of
equal area on both. Orders where no acceptable solution exists keep the geometric weights, leaving
them exactly as accurate as before, and are reported in a warning."""
function quadrature_weights(
        ::Type{PerOrderQuadrature},
        NF::Type{<:Real},
        grid::AbstractGrid,
        spectrum::AbstractSpectrum,
        legendre_polynomials::LowerTriangularArray,
        mmax_truncation::AbstractVector{<:Integer},
        nlons::AbstractVector{<:Integer},
        nlat::Integer,
        solid_angles::AbstractVector,
    )
    (; lmax, mmax) = spectrum                   # 1-based
    (; nlat_half) = grid

    # start from the grid's geometric solid angle, order by order
    ΔΩ = NF[solid_angles[j] for _ in 1:mmax, j in 1:nlat_half]
    g₀ = ring_solid_angles(nlons, nlat, nlat_half, solid_angles)

    λ = Array(legendre_polynomials.data)        # (nonzeros, nlat_half), on the CPU for the fit
    inexact = Int[]                             # orders the fit could not make exact

    for m in 1:mmax
        # rings the Legendre shortcut keeps at this order
        rings = [j for j in 1:nlat_half if mmax_truncation[j] + 1 >= m]
        isempty(rings) && continue

        # The reduced (diagonal-only) system is equivalent to the full orthogonality system only
        # where a solution exists, i.e. where at least as many rings contribute as there are
        # degrees. Below that, satisfying the diagonal conditions in a least-squares sense leaves
        # the off-diagonal ones unconstrained and can make the transform worse, so keep the
        # geometric weights instead.
        ndegrees = lmax - m + 1
        if length(rings) < ndegrees
            push!(inexact, m - 1)
            continue
        end

        # Σ_j g_j λ_lm(μ_j)² = 1, solved for the relative change u = (g - g₀)/g₀
        λ² = Float64.(λ[LowerTriangularArrays.get_lm_range(m, lmax - 1), rings]) .^ 2
        u, residual = minimum_relative_change(λ², ones(ndegrees), g₀[rings])

        # An exactly determined system has a unique solution, but on these grids it is an
        # interpolatory quadrature at maximal degree and oscillates the way Newton-Cotes does.
        # Keep the geometric weights rather than trust weights that stray this far.
        if !all(isfinite, u) || !all(r -> QUADRATURE_MIN_RATIO <= r <= QUADRATURE_MAX_RATIO, 1 .+ u)
            push!(inexact, m - 1)
            continue
        end

        # Weights within bounds are still worth using when the solve left a residual — they beat
        # the geometric ones — but the transform is then only approximate at this order, and the
        # bounds alone would not have revealed it. Record it so the warning below is honest.
        residual > QUADRATURE_RTOL && push!(inexact, m - 1)

        for (i, j) in enumerate(rings)
            ΔΩ[m, j] = NF(solid_angles[j] * (1 + u[i]))
        end
    end

    if !isempty(inexact)
        @warn "Quadrature weights could not be made exact for $(length(inexact)) of $mmax orders m " *
            "(m = $(first(inexact)) … $(last(inexact))); the transform stays approximate there. " *
            "The grid is too coarse for this truncation: " *
            "$(nonparametric_type(typeof(grid))) needs nlat_half ≳ $(cld(17 * lmax, 16)) at " *
            "truncation $lmax, but has $nlat_half. Increase the grid resolution, or the dealiasing " *
            "(default $(default_dealiasing(typeof(grid)))), for an exact transform." maxlog = 3
    end

    return ΔΩ
end
