# Relative singular value below which a direction of the exactness system is treated as noise.
const QUADRATURE_RTOL = 1.0e-8

# Bounds on the optimised weights relative to the equal-area weights. A well-resolved
# configuration stays within [0.9, 1.2]; an under-resolved one overshoots by orders of magnitude
const QUADRATURE_MIN_RATIO = 0.1
const QUADRATURE_MAX_RATIO = 2.0

# How far above 1 `σmax` may sit for `ContractiveQuadrature` to consider an order already
# non-expansive and leave its weights untouched. This is the tolerance the `‖A‖ ≤ 1` guarantee is
# stated to, so it sits at roundoff rather than at `QUADRATURE_RTOL`, which is a fitting tolerance
# and orders of magnitude too loose for it.
const QUADRATURE_NONEXPANSIVE_TOL = 1.0e-12

"""Quadrature scheme used by the forward (grid to spectral) transform to weight each latitude
ring. 

Implementations are `EqualAreaQuadrature`, `RingQuadrature`, `PerOrderQuadrature` and
`ContractiveQuadrature`."""
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
Per-order weights fitted so that the analysis (grid to spectral) operator is non-expansive and as
accurate as any non-expansive weighting can be. See `quadrature_weights` for the derivation.

`PerOrderQuadrature` makes `analysis ∘ synthesis` the identity, which pins down the transform on
the *retained band* but leaves it free above the truncation — where a nonlinear model puts power
every step. This scheme targets the other operator instead: it caps `‖analysis‖ ≤ 1` in the
physical (solid-angle) norm on grid space, so no grid field of any wavenumber, in band or aliased,
can come back with more energy than it had. Because hyperdiffusion is itself a contraction, that
bound survives it, `‖D ∘ analysis‖ ≤ 1`, without the quadrature needing to know anything about the
model's diffusion settings.

The price is that the transform is no longer exact, and that price is a bound rather than a
shortcoming of the fit. Exactness makes analysis a left inverse of synthesis, and every left
inverse obeys `‖A‖ ≥ 1/σmin(S)`, set by the grid's geometric ring areas alone and not by the
weights being chosen. On the HEALPix grids `σmin(S) ≈ 0.95` — the equal-area quadrature
under-integrates some `|Y_lm|²` by ~9.5%, so synthesis loses that energy and an exact analysis has
to amplify the same direction to put it back, which it then also does to any aliased content lying
along it. Exactness therefore forces `‖A‖ ≳ 1.05` at the default dealiasing, and `PerOrderQuadrature`
sits within 0.5% of that floor. Dually, `‖A‖ ≤ 1` forces `‖A S - I‖₂ ≥ 1 - σmin(S)`, and this scheme
fits to within 3% of *that* floor. Grids with an exact quadrature rule have `σmin(S) = 1`, where the
two requirements coincide and this scheme is a no-op.

Within that limit it is a strict improvement on the equal-area weights rather than a rescaling of
them: at the default dealiasing the round-trip error is 1.85–1.92× *below* `EqualAreaQuadrature` on
both HEALPix grids at T32/T64/T128, the round-trip gain is below 1 everywhere (equal area exceeds it
on `OctaHEALPixGrid`), and the residual damping is 4·10⁻⁵ … 1·10⁻⁶ rather than equal area's
one-signed error. Orders too coarsely resolved for the fit fall back to a uniform rescaling that is
still non-expansive, so the bound never depends on the fit succeeding."""
struct ContractiveQuadrature <: AbstractQuadrature end

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
`EqualAreaQuadrature`, `RingQuadrature`, `PerOrderQuadrature` or `ContractiveQuadrature`."""
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

"""$(TYPEDSIGNATURES)
Per-order weights that make analysis non-expansive with the smallest round-trip error any such
weighting can have, see `ContractiveQuadrature`.

Write `Λ[l, j] = λ_lm(μ_j)` for one order `m`, `g⁰_j` for the total solid angle of ring pair `j`
and `g_j` for the weight actually used. Analysis maps a ring's Fourier amplitudes to spectral
coefficients as `x_l = Σ_j g_j Λ[l, j] a_j`, while the physical (real L2) norm of the same ring
data is `Σ_j g⁰_j |a_j|²`. So in the physical metric the analysis and synthesis operators are

```
A = Λ diag(g) diag(g⁰)^(-1/2) ,      S = diag(g⁰)^(1/2) Λᵀ
```

and `A` can only ever remove energy if `σmax(A) ≤ 1`. On grids whose latitudes carry an exact
quadrature rule `g = g⁰` makes `A` the adjoint of `S`, so `σmax(A) = 1` identically and there is
nothing to do. On the HEALPix grids the equal-area weights are not a quadrature rule, the Gram
matrix `G = Λ diag(g⁰) Λᵀ` has eigenvalues spread around 1 (at T128, `nlat_half = 144`: 0.905 to
1.0029), and `σmax(A) = √λmax > 1`.

# What the optimum is

Being non-expansive costs accuracy, and exactly how much is a bound rather than a matter of the
construction. For any `A` at all with `‖A‖ ≤ 1`, dense ones included,

```
‖A S - I‖₂ ≥ 1 - σmin(S) ,           σmin(S)² = λmin(G)
```

because in the SVD `Λ diag(g⁰)^(1/2) = U Σ Vᵀ` the matrix `Y = Uᵀ A V` obeys `‖Y‖ ≤ 1`, so
`‖A S - I‖₂ = ‖Y Σ - I‖₂ ≥ max_i |Y_ii σ_i - 1| ≥ 1 - σmin`. Equality holds at `Y = I`, i.e. at the
orthogonal polar factor `A = U Vᵀ`, whose singular values are all exactly 1 and whose round trip is

```
A S = G^(1/2)
```

That is the target fitted here. It is the mirror image of the bound in `PerOrderQuadrature`'s
direction: exactness forces `‖A‖ ≥ 1/σmin(S)`, non-expansiveness forces `‖A S - I‖₂ ≥ 1 - σmin(S)`,
and the two meet only where `σmin(S) = 1`, i.e. on grids that already carry a quadrature rule.

# How it is fitted

The polar factor itself is a dense matrix per order and parity, too expensive to apply on every
transform. But the map `g ↦ A S` factors through one number per degree — the argument in
`PerOrderQuadrature` — so matching the `n_m` diagonal entries of the target pins down the whole
achievable round trip:

```
Σ_j g_j λ_lm(μ_j)² = (G^(1/2))_ll          for l = m … lmax
```

solved, as there, for the minimum relative change from the equal-area weights. Measured, this
reduced system lands within 3% of the `1 - σmin` floor and agrees with the full `O(n_m²)`-row
system to five digits at a fraction of the cost.

The fit is then rescaled by `1/max(1, σmax(A))`, so `‖A‖ ≤ 1` holds *by construction* rather than by
the fit having succeeded.

# When the fit is declined

Like `PerOrderQuadrature`, the fit needs at least as many contributing rings as degrees, and below
that it returns oscillating or negative weights. Unlike there, the cheap guards are not enough to
detect it: the reduced system only implies the off-diagonal conditions where the target is
*reachable*, so below the sanctioned dealiasing it can match the diagonal to machine precision, keep
every weight well inside `[QUADRATURE_MIN_RATIO, QUADRATURE_MAX_RATIO]`, and still come out several
times worse than equal area. Acceptance is therefore decided on the achieved round trip
`Λ diag(g) Λᵀ`, formed per parity block and compared against the fallback's — which is free, since
rescaling by a scalar takes `G`'s eigenvalues with it.

Declined orders fall back to a plain uniform rescaling `g = g⁰/max(1, σmax)` — the whole of this
scheme before the fit was added. That fallback has no resolution requirement and always exists, so
`‖A‖ ≤ 1` holds at every resolution and dealiasing, and the result is never worse than the uniform
rescaling alone; only the accuracy gain over `EqualAreaQuadrature` is conditional on spectral slack.

The maximum over parities is taken because the transform adds and subtracts the hemispheres, so it
is block diagonal in the parity of `l - m` and the two blocks have different singular values. The
weights array is indexed by order only, so the stricter of the two is applied to both."""
function quadrature_weights(
        ::Type{ContractiveQuadrature},
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

    ΔΩ = NF[solid_angles[j] for _ in 1:mmax, j in 1:nlat_half]
    g₀ = ring_solid_angles(nlons, nlat, nlat_half, solid_angles)
    λ = Array(legendre_polynomials.data)        # (nonzeros, nlat_half), on the CPU for the fit
    unfitted = Int[]                            # orders that fell back to uniform rescaling

    for m in 1:mmax
        # rings the Legendre shortcut keeps at this order
        rings = [j for j in 1:nlat_half if mmax_truncation[j] + 1 >= m]
        isempty(rings) && continue

        Λ = Float64.(λ[LowerTriangularArrays.get_lm_range(m, lmax - 1), rings])
        ndegrees = lmax - m + 1
        parities = [p for p in (0, 1) if length((1 + p):2:ndegrees) > 0]
        B = Λ .* sqrt.(g₀[rings])'              # = Λ diag(g⁰)^(1/2) = Sᵀ, all degrees

        # Diagonal of the polar factor's round trip G^(1/2), and σmax of the unfitted operator, both
        # from one SVD per parity block. From the SVD rather than an eigendecomposition of
        # `block * block'`: `G = U Σ² Uᵀ` so `G^(1/2) = U Σ Uᵀ` either way, but going through the
        # Gram squares the condition number, which is enough to put `σmax` wrong in the 9th digit
        # and silently skip a shrink that was needed.
        target = Vector{Float64}(undef, ndegrees)
        σmax_geometric = 0.0
        geometric_singular_values = Dict{Int, Vector{Float64}}()
        for parity in parities
            block = @view B[(1 + parity):2:end, :]
            factorization = LinearAlgebra.svd(block)
            σ, U = factorization.S, factorization.U
            σmax_geometric = max(σmax_geometric, maximum(σ))
            geometric_singular_values[parity] = σ                # of S, for the fallback's error
            for (i, l) in enumerate((1 + parity):2:ndegrees)
                target[l] = sum(U[i, k]^2 * σ[k] for k in eachindex(σ))
            end
        end

        # An order whose geometric weights already analyse without expanding needs no intervention:
        # this scheme exists to remove the injection, and there is none here. Leaving it untouched
        # keeps every grid that carries a quadrature rule in its latitudes bit-identical to
        # `EqualAreaQuadrature`. The tolerance is the one the `‖A‖ ≤ 1` guarantee is stated to, so
        # it must stay at roundoff — `QUADRATURE_RTOL` is a fitting tolerance, 10⁴ times too loose,
        # and at that width it also swallows HEALPix orders that do need the shrink.
        σmax_geometric <= 1 + QUADRATURE_NONEXPANSIVE_TOL && continue

        λ² = Λ .^ 2

        # The fit needs at least as many contributing rings as degrees, as in PerOrderQuadrature:
        # below that the reduced system stops being equivalent to the full one and the weights
        # oscillate. Fall back rather than trust them.
        accepted = length(rings) >= ndegrees
        local weights
        if accepted
            u, residual = minimum_relative_change(λ², target, g₀[rings])
            accepted = residual <= QUADRATURE_RTOL && all(isfinite, u) &&
                all(r -> QUADRATURE_MIN_RATIO <= r <= QUADRATURE_MAX_RATIO, 1 .+ u)

            if accepted
                # Shrink so that ‖A‖ ≤ 1 holds by construction rather than by the fit succeeding.
                weights = g₀[rings] .* (1 .+ u)
                scaled = Λ .* (weights ./ sqrt.(g₀[rings]))'
                σmax = maximum(LinearAlgebra.opnorm(@view scaled[(1 + p):2:end, :]) for p in parities)
                weights ./= max(1, σmax)

                # The reduced system only implies the off-diagonal conditions where the target is
                # reachable. Below that it matches the diagonal by a route that wrecks the rest, and
                # the weights can stay well inside the bounds above while the transform ends up
                # several times worse. Nothing short of the achieved round trip detects that, so
                # form it and keep the fit only if it beats the fallback it would otherwise get.
                # The fallback's own error is free: rescaling by a scalar takes G's eigenvalues with
                # it, so ‖G/σmax - I‖_F² = Σ_k (λ_k/σmax - 1)².
                shrink_geometric = max(1, σmax_geometric)
                fitted_error = 0.0
                fallback_error = 0.0
                for parity in parities
                    Λp = @view Λ[(1 + parity):2:end, :]
                    fitted_error += sum(abs2, (Λp .* weights') * Λp' - LinearAlgebra.I)
                    fallback_error += sum(abs2, geometric_singular_values[parity] .^ 2 ./ shrink_geometric .- 1)
                end
                accepted = fitted_error < fallback_error
            end
        end

        if !accepted
            # uniform rescaling: always available, no resolution requirement, still contractive
            push!(unfitted, m - 1)
            for j in rings
                ΔΩ[m, j] = NF(solid_angles[j] / σmax_geometric)
            end
            continue
        end

        for (i, j) in enumerate(rings)
            ΔΩ[m, j] = NF(weights[i] / (((nlat - j + 1 == j) ? 1 : 2) * nlons[j]))
        end
    end

    # Declining an order is graceful here — the fallback is still non-expansive, just no more
    # accurate than equal area — so this warns only where resolution is the cause and the user can
    # act on it. A few top orders declining on a well-resolved grid is normal: there the fit simply
    # does not beat the uniform rescaling, and saying "the grid is too coarse" would be wrong.
    required_nlat_half = cld(17 * lmax, 16)
    if !isempty(unfitted) && nlat_half < required_nlat_half
        @warn "Quadrature weights could not be fitted for $(length(unfitted)) of $mmax orders m " *
            "(m = $(first(unfitted)) … $(last(unfitted))); analysis stays non-expansive there but " *
            "no more accurate than EqualAreaQuadrature. The grid is too coarse for this truncation: " *
            "$(nonparametric_type(typeof(grid))) needs nlat_half ≳ $required_nlat_half at " *
            "truncation $lmax, but has $nlat_half. Increase the grid resolution, or the dealiasing " *
            "(default $(default_dealiasing(typeof(grid)))), for the full accuracy gain." maxlog = 3
    end

    return ΔΩ
end
