"""Cubic (16-point) interpolation on ring grids.

`AnvilInterpolator` reads a 4-point stencil — two grid points on each of the two rings bracketing
the target latitude — and averages them bilinearly. That is cheap but strongly damping, which
matters when the interpolation is applied every time step, as in semi-Lagrangian advection.

`CubicInterpolator` reads a 16-point stencil instead: 4 longitude points on each of the 4 rings
`j-1, j, j+1, j+2` around the target. Within a ring the longitude spacing is uniform, so the
longitude weights are the cubic Lagrange polynomials on a uniform stencil. Across rings the
latitudes are not uniformly spaced (Gaussian latitudes, and reduced grids vary), so the latitude
weights are the general Lagrange basis over the four ring latitudes.

The two directions are combined into a single weight per stencil point, `w = w_lat * w_lon`, so
the interpolation kernel is a plain 16-term weighted sum. Weights sum to 1 by construction, and
the scheme reproduces any polynomial of degree ≤ 3 in each direction exactly.

Near the poles fewer than 4 rings are available. The poles are then treated as two extra
"rings" carrying the ring-average value (the same convention `AnvilInterpolator` uses, with flag
indices 0 for the north pole and -1 for the south pole) and the latitude Lagrange basis is built
over however many of the 4 slots are usable, so the weights still sum to 1 and the scheme degrades
gracefully to quadratic and then linear rather than extrapolating."""

export CubicInterpolator

"""Contains the 16-point stencil and its weights for [`CubicInterpolator`](@ref).
The stencil is stored flat with stride `NSTENCIL_CUBIC`: entry `s` of output point `k` lives at
`(k-1)*NSTENCIL_CUBIC + s`. Slot `s = 4(r-1) + p` is longitude point `p ∈ 1:4` on ring slot
`r ∈ 1:4`, where ring slot `r` is grid ring `js + r - 2`.
$(TYPEDFIELDS)"""
@kwdef struct CubicLocator{
        VectorType,
        VectorIntType,
        IntType,
    } <: AbstractLocator

    "number of points to interpolate onto"
    npoints_output::IntType

    "ring indices j such that [j, j+1) contains the point, as for every ring-based locator"
    js::VectorIntType = zeros(Int, npoints_output)

    "distance fractions between rings j and j+1, as for every ring-based locator"
    Δys::VectorType = zero(VectorType(undef, npoints_output))

    "flat indices ij of the 16 stencil points, 0 = north pole, -1 = south pole"
    ijs::VectorIntType = zeros(Int, NSTENCIL_CUBIC * npoints_output)

    "combined latitude*longitude weight of each of the 16 stencil points, summing to 1"
    weights::VectorType = zero(VectorType(undef, NSTENCIL_CUBIC * npoints_output))
end

"""4 longitude points on each of 4 rings."""
const NSTENCIL_CUBIC = 16

Adapt.@adapt_structure CubicLocator

function Architectures.on_architecture(arch::AbstractArchitecture, loc::CubicLocator)
    return CubicLocator(
        npoints_output = loc.npoints_output,
        js = on_architecture(arch, loc.js),
        Δys = on_architecture(arch, loc.Δys),
        ijs = on_architecture(arch, loc.ijs),
        weights = on_architecture(arch, loc.weights),
    )
end

function Base.show(io::IO, L::CubicLocator)
    println(io, styled"{warning:CubicLocator} ($NSTENCIL_CUBIC-point stencil)")
    print(io, styled"└ {info:onto}: $(L.npoints_output) points")
    return nothing
end

"""Interpolator for a cubic 16-point stencil, see [`CubicLocator`](@ref).

NF is the number format used to calculate the interpolation, which can be different from the
input data and/or the interpolated data on the new grid."""
struct CubicInterpolator{NF, Geometry, Locator} <: AbstractInterpolator
    geometry::Geometry
    locator::Locator
end

function Adapt.adapt_structure(to, I::CubicInterpolator{NF}) where {NF}
    geometry = adapt_structure(to, I.geometry)
    locator = adapt_structure(to, I.locator)
    return CubicInterpolator{eltype(NF), typeof(geometry), typeof(locator)}(geometry, locator)
end

Base.eltype(::CubicInterpolator{NF}) where {NF} = NF
grid_type(I::CubicInterpolator) = typeof(I.geometry.grid)
nonparametric_type(::Type{<:CubicInterpolator}) = CubicInterpolator

# the hooks that make the generic machinery pick the right pieces
Locator(::Type{<:CubicInterpolator}) = CubicLocator

function Base.show(io::IO, L::CubicInterpolator{NF}) where {NF}
    NF_str = "{$NF}"
    println(io, styled"{warning:CubicInterpolator}{note:$NF_str} for $(L.geometry.grid)")
    print(io, styled"└ {info:onto}: $(L.locator.npoints_output) points")
    return nothing
end

# generator with NF default based on geometry and locator, mirroring AnvilInterpolator
function CubicInterpolator(
        geometry::AbstractGridGeometry,
        locator::AbstractLocator;
        NF = DEFAULT_NF
    )
    return CubicInterpolator{NF, typeof(geometry), typeof(locator)}(geometry, locator)
end

"""$(TYPEDSIGNATURES)
Cubic Lagrange weights on a *uniform* 4-point stencil at offsets -1, 0, 1, 2, evaluated at
fractional position `Δ ∈ [0, 1)` between stencil points 0 and 1. Used along a ring, where the
longitude spacing `360/nlon` is uniform by construction. Exact for cubics, and the four weights
sum to 1."""
@inline function cubic_lagrange_weights(Δ::NF) where {NF}
    Δm1 = Δ - 1
    Δm2 = Δ - 2
    Δp1 = Δ + 1
    w0 = -Δ * Δm1 * Δm2 / 6         # at offset -1
    w1 = Δp1 * Δm1 * Δm2 / 2        # at offset  0
    w2 = -Δp1 * Δ * Δm2 / 2         # at offset  1
    w3 = Δp1 * Δ * Δm1 / 6          # at offset  2
    return w0, w1, w2, w3
end

"""$(TYPEDSIGNATURES)
All four latitude Lagrange weights at `y` over nodes `y1..y4`, with `v1..v4` flagging which nodes
take part. Invalid nodes are skipped in both numerator and denominator (their factors are replaced
by 1 and their own weight set to 0), so the weights over the remaining nodes still sum to 1 and the
interpolation degrades cubic → quadratic → linear near the poles instead of extrapolating.

Computed in one pass: 4 divisions and a handful of `ifelse`, rather than calling a per-node helper
that re-selects its node with nested `ifelse` chains. `find_grid_indices!` runs over every grid
point every time step, so this is the hot path."""
@inline function lagrange_weights_4(y, y1, y2, y3, y4, v1::Bool, v2::Bool, v3::Bool, v4::Bool)
    one_ = one(y)

    # distances to the target, neutralised to 1 for nodes that do not take part
    d1 = ifelse(v1, y - y1, one_)
    d2 = ifelse(v2, y - y2, one_)
    d3 = ifelse(v3, y - y3, one_)
    d4 = ifelse(v4, y - y4, one_)

    # node separations, likewise neutralised (never zero for valid pairs: ring latitudes differ)
    a12 = ifelse(v2, y1 - y2, one_); a13 = ifelse(v3, y1 - y3, one_); a14 = ifelse(v4, y1 - y4, one_)
    a21 = ifelse(v1, y2 - y1, one_); a23 = ifelse(v3, y2 - y3, one_); a24 = ifelse(v4, y2 - y4, one_)
    a31 = ifelse(v1, y3 - y1, one_); a32 = ifelse(v2, y3 - y2, one_); a34 = ifelse(v4, y3 - y4, one_)
    a41 = ifelse(v1, y4 - y1, one_); a42 = ifelse(v2, y4 - y2, one_); a43 = ifelse(v3, y4 - y3, one_)

    w1 = ifelse(v1, (d2 * d3 * d4) / (a12 * a13 * a14), zero(y))
    w2 = ifelse(v2, (d1 * d3 * d4) / (a21 * a23 * a24), zero(y))
    w3 = ifelse(v3, (d1 * d2 * d4) / (a31 * a32 * a34), zero(y))
    w4 = ifelse(v4, (d1 * d2 * d3) / (a41 * a42 * a43), zero(y))
    return w1, w2, w3, w4
end

"""$(TYPEDSIGNATURES)
Write one ring slot's 4 stencil entries: indices into `ijs` and weights into `weights` at offset
`off`. `jr` is the ring number, `wlat` its latitude weight. Rings outside `[1, nlat]` are the poles
and carry a single value, flagged 0 (north) / -1 (south) with the whole latitude weight on the
first entry.

The longitude wrap is done with `ifelse` rather than `mod`, which would be an integer division per
stencil point (16 per grid point per time step). `i0 ∈ [0, nlon)` so a single conditional add or
subtract is enough, given `nlon >= 4`."""
@inline function write_ring_stencil!(
        ijs, weights, off, jr, wlat, λ,
        lon_offsets, nlons, ring_starts, nlat,
    )
    NF = eltype(weights)

    if (jr < 1) | (jr > nlat)
        pole = ifelse(jr < 1, 0, -1)
        ijs[off + 1] = pole
        ijs[off + 2] = pole
        ijs[off + 3] = pole
        ijs[off + 4] = pole
        weights[off + 1] = wlat
        weights[off + 2] = zero(NF)
        weights[off + 3] = zero(NF)
        weights[off + 4] = zero(NF)
    else
        nlon = nlons[jr]
        start = ring_starts[jr]
        i_a, _, Δ = find_lon_indices(λ, lon_offsets[jr], nlon)
        c1, c2, c3, c4 = cubic_lagrange_weights(NF(Δ))

        i0 = i_a - 1                     # 0-based in-ring index of stencil offset 0, ∈ [0, nlon)
        im = i0 - 1; im = ifelse(im < 0, im + nlon, im)
        ip = i0 + 1; ip = ifelse(ip >= nlon, ip - nlon, ip)
        iq = i0 + 2; iq = ifelse(iq >= nlon, iq - nlon, iq)

        ijs[off + 1] = start + im
        ijs[off + 2] = start + i0
        ijs[off + 3] = start + ip
        ijs[off + 4] = start + iq
        weights[off + 1] = wlat * c1
        weights[off + 2] = wlat * c2
        weights[off + 3] = wlat * c3
        weights[off + 4] = wlat * c4
    end
    return nothing
end

function find_grid_indices!(
        locator::CubicLocator,      # update indices and weights
        geometry::AbstractGridGeometry,
        λs::AbstractArray,           # based on new longitudes λ
        architecture::AbstractArchitecture = architecture(λs)
    )

    (; js, Δys, ijs, weights) = locator
    (; nlons, ring_starts, lon_offsets, nlat, latd) = geometry

    λs_converted = convert.(eltype(lon_offsets), λs)

    launch!(
        architecture,
        LinearWorkOrder,
        size(js),
        find_cubic_indices_kernel!,
        ijs, weights,
        js, Δys,
        λs_converted,
        lon_offsets,
        nlons,
        ring_starts,
        latd,
        nlat,
    )
    return nothing
end

@kernel inbounds = true function find_cubic_indices_kernel!(
        ijs, weights,
        @Const(js), @Const(Δys),
        @Const(λs),
        @Const(lon_offsets),
        @Const(nlons),
        @Const(ring_starts),
        @Const(latd),           # ring latitudes INCLUDING both poles, length nlat+2
        nlat,
    )
    k = @index(Global, Linear)

    j = js[k]           # θ ∈ [latd[j], latd[j+1]) in ring numbering, 0 = north of ring 1
    λ = λs[k]

    # Reconstruct the target latitude from the ring fraction. `latd` includes the poles, so ring r
    # has latitude latd[r+1], and find_rings! defined Δy = (latd[j] - θ)/(latd[j] - latd[j+1]).
    θ = latd[j + 1] + Δys[k] * (latd[j + 2] - latd[j + 1])

    # THE FOUR RING SLOTS j-1, j, j+1, j+2. Ring 0 is the north pole and ring nlat+1 the south
    # pole, both carrying the ring-average value. With j ∈ [0, nlat] the middle two slots are
    # always in [0, nlat+1] and only the outer two can fall off the end.
    v1 = j >= 1
    v4 = j <= nlat - 1

    # latitudes of the four slots, clamped so the reads stay in bounds (invalid slots carry
    # weight 0 anyway, but the index still has to be legal)
    y1 = latd[clamp(j, 1, nlat + 2)]
    y2 = latd[j + 1]
    y3 = latd[j + 2]
    y4 = latd[clamp(j + 3, 1, nlat + 2)]

    w1, w2, w3, w4 = lagrange_weights_4(θ, y1, y2, y3, y4, v1, true, true, v4)

    base = (k - 1) * NSTENCIL_CUBIC
    write_ring_stencil!(ijs, weights, base, j - 1, w1, λ, lon_offsets, nlons, ring_starts, nlat)
    write_ring_stencil!(ijs, weights, base + 4, j, w2, λ, lon_offsets, nlons, ring_starts, nlat)
    write_ring_stencil!(ijs, weights, base + 8, j + 1, w3, λ, lon_offsets, nlons, ring_starts, nlat)
    write_ring_stencil!(ijs, weights, base + 12, j + 2, w4, λ, lon_offsets, nlons, ring_starts, nlat)
end

# the actual interpolation: a plain 16-term weighted sum
function _interpolate!(
        Aout,                               # Out: interpolated values
        A,                                  # gridded values to interpolate from
        locator::CubicLocator,
        geometry::GridGeometry,
        architecture::AbstractArchitecture
    )
    (; npoints_output, ijs, weights) = locator
    (; npoints) = geometry
    (; rings) = geometry.grid # CPU version even on GPU

    @boundscheck size(Aout, 1) == npoints_output || throw(DimensionMismatchArray(Aout, ijs))
    @boundscheck length(A) == npoints ||
        throw(DimensionMismatch("Interpolator ($npoints points) mismatches input grid ($(length(A)) points)."))
    @boundscheck extrema_in(ijs, -1, npoints) || throw(BoundsError)

    A_northpole, A_southpole = average_on_poles(A, rings)

    launch!(
        architecture,
        LinearWorkOrder,
        (npoints_output,),
        _interpolate_cubic_kernel!,
        Aout, A, ijs, weights, A_northpole, A_southpole,
    )

    return Aout
end

@kernel inbounds = true function _interpolate_cubic_kernel!(
        Aout, @Const(A), @Const(ijs), @Const(weights), A_northpole, A_southpole,
    )
    k = @index(Global, Linear)

    base = (k - 1) * NSTENCIL_CUBIC
    sum = zero(eltype(Aout))

    for s in 1:NSTENCIL_CUBIC
        ij = ijs[base + s]
        # 0 and -1 are the pole flags, see CubicLocator; clamp keeps the read in bounds so the
        # branch does not have to guard the indexing itself
        value = ifelse(
            ij == 0, A_northpole,
            ifelse(ij == -1, A_southpole, A[max(ij, 1)])
        )
        sum += weights[base + s] * value
    end

    Aout[k] = sum
end
