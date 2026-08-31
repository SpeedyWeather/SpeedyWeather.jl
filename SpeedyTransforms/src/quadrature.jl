# Relative singular value below which a direction of the exactness system is treated as noise.
const QUADRATURE_RTOL = 1.0e-8

# Bounds on the optimised weights relative to the equal-area weights. A well-resolved
# configuration stays within [0.9, 1.2]; an under-resolved one overshoots by orders of magnitude
const QUADRATURE_MIN_RATIO = 0.1
const QUADRATURE_MAX_RATIO = 2.0

"""$(TYPEDSIGNATURES)
Whether a grid's geometric solid angles should be replaced by optimised quadrature weights.
Gaussian and Clenshaw-Curtis grids already carry an exact quadrature rule in their latitudes,
the HEALPix grids do not: they are equal-area discretisations whose weights are a Riemann sum,
not a quadrature rule."""
optimize_quadrature(Grid::Type{<:AbstractGrid}) = false
optimize_quadrature(
    ::Type{
        <:Union{
            HEALPixGrid, OctaHEALPixGrid,
            FullHEALPixGrid, FullOctaHEALPixGrid,
        },
    }
) = true

optimize_quadrature(grid::AbstractGrid) = optimize_quadrature(typeof(grid))

"""$(TYPEDSIGNATURES)
Quadrature weights `ΔΩ[m, j]` = sinθ Δθ Δϕ for order `m` (1-based) and northern latitude ring `j`,
as read by the forward (grid to spectral) Legendre transform.

For grids whose latitudes carry an exact quadrature rule every column is the grid's geometric
solid angle, i.e. `ΔΩ[m, j]` does not depend on `m` and this reduces to the previous behaviour.
For the HEALPix grids the weights are optimised per order to make the transform exact.

The transform is exact when, for every order `m`,

```
Σ_j g_j λ_lm(μ_j) λ_l'm(μ_j) = δ_ll'      for all l, l' ≥ m
```

with `λ_lm` the normalised associated Legendre polynomials, `μ_j = sin(lat_j)`, `g_j` the total
solid angle of ring pair `j` (north + south, so `Σ_j g_j = 4π`) and `δ_ll'` the Kronecker delta,
1 for `l == l'` and 0 otherwise. That is, the quadrature has to reproduce the exact orthonormality
of the spherical harmonics: 1 on the diagonal (right normalisation, no gain error) and exactly 0
off it (no leakage between degrees at the same order `m`) — which is what makes analysis ∘
synthesis the identity. That system looks quadratic in the
number of degrees, but it is not: substituting `x = μ²`, the products `λ_lm λ_l'm` with `l + l'`
even span `(1-x)^m` times the polynomials of degree `≤ lmax - m` in `x`, and the diagonal products
`λ_lm²` are already a triangular basis of that space. So the off-diagonal conditions follow from the
diagonal ones and the whole system collapses to one row per degree — *integrate every `|Y_lm|²`
exactly*:

```
Σ_j g_j λ_lm(μ_j)² = 1                    for l = m … lmax
```

which is solved here for the minimum *relative* change from the equal-area weights. Relative,
because the weights span orders of magnitude between polar and equatorial rings, and minimising the
plain distance would let a polar ring swing past zero while barely moving the objective.

A solution exists only if enough rings contribute at order `m`. At `m = 0` that means
`nlat_half ≥ truncation`, and with equality the system is square: its unique solution is an
interpolatory quadrature at maximal degree, which oscillates the way Newton-Cotes does. A margin of
`nlat_half/16` is enough for `HEALPixGrid`, `OctaHEALPixGrid` needs a little more, and the
`nlat_half/9` that `default_dealiasing` provides keeps every weight positive and within 20% of
equal area on both. Orders where no acceptable solution exists keep the geometric weights, leaving
them exactly as accurate as before, and are reported in a warning."""
function quadrature_weights(
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

    # default: the grid's geometric solid angle, independent of order m
    ΔΩ = NF[solid_angles[j] for _ in 1:mmax, j in 1:nlat_half]
    optimize_quadrature(grid) || return ΔΩ

    # total solid angle of ring pair j (north + south), Σ_j g₀ = 4π. The equator, if the grid
    # has one, is a single ring and enters once.
    multiplicity(j) = (nlat - j + 1 == j) ? 1 : 2
    g₀ = [multiplicity(j) * nlons[j] * Float64(solid_angles[j]) for j in 1:nlat_half]

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
        A = λ² * LinearAlgebra.Diagonal(g₀[rings])
        b = 1 .- λ² * g₀[rings]
        F = LinearAlgebra.svd(A)
        keep = 1:count(>(QUADRATURE_RTOL * F.S[1]), F.S)
        u = F.V[:, keep] * ((F.U[:, keep]' * b) ./ F.S[keep])

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
        LinearAlgebra.norm(A * u - b) > QUADRATURE_RTOL * LinearAlgebra.norm(b) && push!(inexact, m - 1)

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
