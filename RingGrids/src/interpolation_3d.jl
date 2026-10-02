"""Supertype for markers describing how a field is staggered vertically relative to the
σ coordinates it is interpolated onto/from. Subtypes bundle the σ-vector together with
whatever boundary condition applies at the missing end of the vertical bracket, so kernels
never have to hardcode a boundary value."""
abstract type AbstractVerticalStaggering end

"""Fields located on full vertical levels (cell centers): u, v, T, etc. Every bracket
index maps directly onto a stored data layer, so there is no boundary condition to
apply. $(TYPEDFIELDS)"""
struct SigmaCenter{V} <: AbstractVerticalStaggering
    "σ_levels_full: one σ value per stored data layer"
    sigma::V
end

Adapt.@adapt_structure SigmaCenter

"""Fields located on half vertical levels (upper cell face, k-½). Layer `k` in the data
stores the value at the face above full-level `k` (σ at k-½); the bottom face (σ=1,
k=nlayers+½) is not stored and is replaced by `bottom_boundary_condition`, zero by
default as in `SigmaFaceAbove(σ)`. $(TYPEDFIELDS)"""
struct SigmaFaceAbove{V, T} <: AbstractVerticalStaggering
    "σ_levels_half: one σ value per face, length nlayers+1"
    sigma::V
    "Value substituted for the unstored bottom face (σ=1)"
    bottom_boundary_condition::T
end

SigmaFaceAbove(σ::AbstractVector) = SigmaFaceAbove(σ, zero(eltype(σ)))
Adapt.@adapt_structure SigmaFaceAbove

"""Fields located on half vertical levels (lower cell face, k+½), e.g. vertical velocity
w. Layer `k` in the data stores the value at the face below full-level `k` (σ at k+½);
the top face (σ=0, k=½) is not stored and is replaced by `top_boundary_condition`,
zero by default as in `SigmaFaceBelow(σ)` (kinematic boundary condition w(σ=0) = 0).
$(TYPEDFIELDS)"""
struct SigmaFaceBelow{V, T} <: AbstractVerticalStaggering
    "σ_levels_half: one σ value per face, length nlayers+1"
    sigma::V
    "Value substituted for the unstored top face (σ=0)"
    top_boundary_condition::T
end

SigmaFaceBelow(σ::AbstractVector) = SigmaFaceBelow(σ, zero(eltype(σ)))
Adapt.@adapt_structure SigmaFaceBelow

# data column index = bracket index k + shift(staggering); out-of-[1, nlayers] means
# the value isn't stored and boundary(staggering) is substituted instead.
@inline shift(::SigmaCenter) = 0
@inline shift(::SigmaFaceAbove) = 0
@inline shift(::SigmaFaceBelow) = -1

@inline boundary(s::AbstractVerticalStaggering) = zero(eltype(s.sigma))
@inline boundary(s::SigmaFaceAbove) = s.bottom_boundary_condition
@inline boundary(s::SigmaFaceBelow) = s.top_boundary_condition

# number of σ levels a staggering needs for a field with nlayers stored layers
@inline nlevels(::SigmaCenter, nlayers) = nlayers
@inline nlevels(::Union{SigmaFaceAbove, SigmaFaceBelow}, nlayers) = nlayers + 1

# Zero-gradient (von Neumann) BC: σ outside [σ_levels[1], σ_levels[end]]
# returns the boundary value unchanged (gradient assumed zero at boundary).
@inline function find_vertical_bracket(σ, σ_levels)
    n = length(σ_levels)
    n == 1 && return 1, 1, zero(σ)      # single level: field is constant in the vertical
    k_lo = 1
    for k in 1:(n - 1)
        σ_levels[k] <= σ && (k_lo = k)
    end
    k_hi = k_lo + 1
    Δσ = σ_levels[k_hi] - σ_levels[k_lo]
    α = clamp((σ - σ_levels[k_lo]) / Δσ, zero(σ), one(σ))
    return k_lo, k_hi, α
end

# Compute north and south pole ring averages per vertical layer on device.
@kernel inbounds = true function _compute_pole_averages_kernel!(
        north_pole_average, south_pole_average, A_data, geometry
    )
    (; ring_starts, nlons, nlat) = geometry
    k = @index(Global, Linear)
    n_north = nlons[1]
    rs_north = ring_starts[1]
    north_sum = zero(eltype(north_pole_average))
    for i in 0:(n_north - 1)
        north_sum += A_data[rs_north + i, k]
    end
    north_pole_average[k] = north_sum / n_north
    n_south = nlons[nlat]
    rs_south = ring_starts[nlat]
    south_sum = zero(eltype(south_pole_average))
    for i in 0:(n_south - 1)
        south_sum += A_data[rs_south + i, k]
    end
    south_pole_average[k] = south_sum / n_south
end

# anvil interpolation of layer k of A_data onto point i, pole values from the locator
Base.@propagate_inbounds function anvil_average(A_data, k, i, locator)
    (; ij_as, ij_bs, ij_cs, ij_ds, Δabs, Δcds, Δys, north_pole_average, south_pole_average) = locator
    a, b = ij_as[i] == 0 ? (north_pole_average[k], north_pole_average[k]) : (A_data[ij_as[i], k], A_data[ij_bs[i], k])
    c, d = ij_cs[i] == -1 ? (south_pole_average[k], south_pole_average[k]) : (A_data[ij_cs[i], k], A_data[ij_ds[i], k])
    return anvil_average(a, b, c, d, Δabs[i], Δcds[i], Δys[i])
end

# Vertically-blended anvil interpolation kernel, generic over the vertical staggering
# (full level / face above / face below): staggering determines the data-column shift
# and the boundary value substituted when a bracket index falls outside stored data.
@kernel inbounds = true function _interpolate_3D_kernel!(
        Aout, A_data, locator, positions, staggering,
    )
    i = @index(Global, Linear)
    k_lo, k_hi, α = find_vertical_bracket(positions[i].σ, staggering.sigma)

    s = shift(staggering)
    k_data_lo = k_lo + s
    k_data_hi = k_hi + s
    nlayers = size(A_data, 2)

    val_lo = 1 <= k_data_lo <= nlayers ? anvil_average(A_data, k_data_lo, i, locator) : boundary(staggering)
    val_hi = 1 <= k_data_hi <= nlayers ? anvil_average(A_data, k_data_hi, i, locator) : boundary(staggering)

    Aout[i] = val_lo + (val_hi - val_lo) * α
end

"""$(TYPEDSIGNATURES)
Vertically blended anvil interpolation at each particle's σ coordinate for a field with
vertical staggering `staggering` (`SigmaCenter`, `SigmaFaceAbove`, or `SigmaFaceBelow`).
Pole ring averages are stored in `locator.north_pole_average` and
`locator.south_pole_average` (pre-allocated via the 3D `AnvilLocator` constructor) and
reused across calls without CPU allocations."""
function interpolate_3D!(Aout, A, locator, geometry, positions, staggering::AbstractVerticalStaggering)
    (; npoints_output, north_pole_average, south_pole_average) = locator
    A_data = A.data
    nlayers = size(A_data, 2)

    @boundscheck length(Aout) == npoints_output || throw(
        DimensionMismatch(
            "Aout has length $(length(Aout)) but the locator has $npoints_output points."
        )
    )
    @boundscheck length(positions) == npoints_output || throw(
        DimensionMismatch(
            "positions has length $(length(positions)) but the locator has $npoints_output points."
        )
    )
    @boundscheck length(north_pole_average) >= nlayers || throw(
        DimensionMismatch(
            "Locator has pole buffers for $(length(north_pole_average)) layers but the field has $nlayers, " *
                "create it with AnvilLocator(NF, npoints, nlayers)."
        )
    )
    @boundscheck length(staggering.sigma) == nlevels(staggering, nlayers) || throw(
        DimensionMismatch(
            "$(nameof(typeof(staggering))) needs $(nlevels(staggering, nlayers)) σ levels for a field " *
                "with $nlayers layers, got $(length(staggering.sigma))."
        )
    )

    arch = architecture(Aout)
    launch!(
        arch, LinearWorkOrder, (nlayers,), _compute_pole_averages_kernel!,
        north_pole_average, south_pole_average, A_data, geometry
    )
    launch!(
        arch, LinearWorkOrder, (npoints_output,),
        _interpolate_3D_kernel!,
        Aout, A_data, locator, positions, staggering,
    )
    return Aout
end

"""$(TYPEDSIGNATURES)
Convenience method: a bare σ-levels vector is interpreted as full-level (`SigmaCenter`)
staggering, e.g. for u, v, temperature. Pass `σ_levels_full` from `model.geometry`."""
interpolate_3D!(Aout, A, locator, geometry, positions, σ::AbstractVector) =
    interpolate_3D!(Aout, A, locator, geometry, positions, SigmaCenter(σ))

"""$(TYPEDSIGNATURES)
`interpolate!` with vertical positions: a 3D interpolation that also blends vertically at
each point's `σ`, forwarding to [`interpolate_3D!`](@ref). `vertical` is an
`AbstractVerticalStaggering` or a σ-levels vector (interpreted as `SigmaCenter`). Without
`positions` and `vertical`, `interpolate!` is a horizontal interpolation via `interpolate_2D!`."""
interpolate!(
    Aout::AbstractVector,
    A::AbstractField,
    locator::AbstractLocator,
    geometry::AbstractGridGeometry,
    positions::AbstractVector,
    vertical::Union{AbstractVerticalStaggering, AbstractVector},
) = interpolate_3D!(Aout, A, locator, geometry, positions, vertical)
