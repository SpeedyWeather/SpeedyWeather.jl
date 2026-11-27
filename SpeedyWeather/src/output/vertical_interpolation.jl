# VERTICAL INTERPOLATION METHODS

"""Supertype of methods to interpolate between two model layers, used by
[`interpolate_pressure_layers!`](@ref) to interpolate from the model's vertical
layers onto pressure layers."""
abstract type AbstractVerticalInterpolation <: AbstractModelComponent end

"""Interpolate linearly in pressure between two model layers."""
struct LinearInPressure <: AbstractVerticalInterpolation end

"""Interpolate linearly in the logarithm of pressure between two model layers. Default,
as most variables vary more linearly with log(p) than with p."""
struct LinearInLogPressure <: AbstractVerticalInterpolation end

"""$(TYPEDSIGNATURES)
Weight `w` of the lower layer (at the higher pressure `p₂`) when interpolating onto
pressure `p` with `p₁ <= p <= p₂`, so that the interpolated value is `(1 - w)*ξ₁ + w*ξ₂`
with `ξ₁, ξ₂` the values at `p₁, p₂`. Returns 0 at `p = p₁` and 1 at `p = p₂`."""
@inline interpolation_weight(p, p₁, p₂, ::LinearInPressure) = (p - p₁) / (p₂ - p₁)
@inline interpolation_weight(p, p₁, p₂, ::LinearInLogPressure) = log(p / p₁) / log(p₂ / p₁)

# VERTICAL EXTRAPOLATION METHODS

"""Supertype of methods to extrapolate above the top-most or below the bottom-most full
model layer, used by [`interpolate_pressure_layers!`](@ref). Subtypes extend
[`extrapolate_below`](@ref) and optionally [`extrapolate_above`](@ref)."""
abstract type AbstractVerticalExtrapolation <: AbstractModelComponent end

"""Hold the outer-most model layer constant beyond the model layers."""
struct ConstantExtrapolation <: AbstractVerticalExtrapolation end

"""Extrapolate dry-adiabatically below the lowest model layer, `T(p) = T*(p/p_bottom)^κ`,
to obtain e.g. a 1000 hPa temperature below a model layer that sits above it.
Same adiabat as used for the mean sea-level pressure output. Fields are $(TYPEDFIELDS)"""
@kwdef struct DryAdiabaticExtrapolation{NF} <: AbstractVerticalExtrapolation
    "[DERIVED] adiabatic exponent κ = R_dry/cₚ, set from the model's atmosphere at initialize!"
    κ::NF = 2 / 7
end

"""Mask everything below the surface (`p > pₛ`) with `missing_value`, and use
`above_surface` between the lowest model layer and the surface, where a pressure layer
is below the lowest model layer but still above ground. Fields are $(TYPEDFIELDS)"""
@kwdef struct SubsurfaceMask{E, NF} <: AbstractVerticalExtrapolation
    "[OPTION] extrapolation between the lowest model layer and the surface"
    above_surface::E = ConstantExtrapolation()

    "[OPTION] value to mask subsurface points with, ideally `variable.missing_value`"
    missing_value::NF = NaN
end

"""$(TYPEDSIGNATURES)
Value at pressure `p` above the top-most full model layer at `p_top` with value `ξ_top`.

TODO: currently the top-most layer is held constant for every extrapolation method. To
change this, e.g. to mask with a missing value or to continue isothermally, extend this
method for a new `<: AbstractVerticalExtrapolation`."""
@inline extrapolate_above(ξ_top, p, p_top, ::AbstractVerticalExtrapolation) = ξ_top

"""$(TYPEDSIGNATURES)
Value at pressure `p` below the bottom-most full model layer at `p_bottom` with value
`ξ_bottom`, given the surface pressure `pₛ`. Note that `p_bottom < pₛ`, so `p` may well
be above ground; `pₛ` is passed on to distinguish the two."""
@inline extrapolate_below(ξ_bottom, p, p_bottom, pₛ, ::ConstantExtrapolation) = ξ_bottom

# κ is converted to the number format of the field so that the exponentiation stays in it
@inline extrapolate_below(ξ_bottom, p, p_bottom, pₛ, E::DryAdiabaticExtrapolation) =
    ξ_bottom * (p / p_bottom)^convert(typeof(ξ_bottom), E.κ)

@inline extrapolate_below(ξ_bottom, p, p_bottom, pₛ, E::SubsurfaceMask) =
    p > pₛ ? convert(typeof(ξ_bottom), E.missing_value) :
    extrapolate_below(ξ_bottom, p, p_bottom, pₛ, E.above_surface)

# INTERPOLATION FROM MODEL LAYERS ONTO PRESSURE LAYERS

"""$(TYPEDSIGNATURES)
Interpolate `in_field` on the model's vertical layers onto the pressure layers `p` [Pa],
writing into the externally allocated `out_field`. Both are fields on the same horizontal
grid, `in_field` has `nlayers` vertical layers, `out_field` has `length(p)` of them.

The vertical `coordinates` (`SigmaCoordinates` or `SigmaPressureCoordinates`, from
`model.geometry.vertical_coordinates`) together with the `surface_pressure` [Pa] define
the pressure of every model layer in every column, so this works for terrain-following
and hybrid sigma-pressure coordinates alike. Interpolation happens in pressure, the
`interpolation` method chooses linear in `p` or in `log(p)`. Outside the model layers
`extrapolation` decides, see [`extrapolate_above`](@ref), [`extrapolate_below`](@ref).

Columns are independent, so this launches a kernel over the horizontal grid points and
the pressure layers and runs on CPU and GPU. `p` therefore has to be on the same
architecture as the fields. `p` does not have to be sorted (every layer is searched for
independently) but normally would be. If `in_field` has a time step dimension, pass in
the step to interpolate, e.g. via `get_prognostic_step`."""
function interpolate_pressure_layers!(
        out_field::AbstractField,                       # OUTPUT: (lonlat, npressure)
        in_field::AbstractField,                        # INPUT: (lonlat, nlayers)
        surface_pressure::AbstractField,                # INPUT: (lonlat), [Pa]
        p::AbstractVector,                              # pressure layers [Pa]
        coordinates::AbstractVerticalCoordinates,       # of the model
        interpolation::AbstractVerticalInterpolation = LinearInLogPressure(),
        extrapolation::AbstractVerticalExtrapolation = ConstantExtrapolation(),
    )
    grids_match(out_field, in_field, surface_pressure) ||
        throw(DimensionMismatch(out_field, in_field))
    size(out_field, 2) == length(p) || throw(
        DimensionMismatch(
            "$(size(out_field, 2)) pressure layers allocated in the output field, " *
                "but $(length(p)) pressure layers requested."
        )
    )

    nlayers = size(in_field, 2)
    launch!(
        architecture(in_field), RingGridWorkOrder, size(out_field),
        interpolate_pressure_layers_kernel!,
        out_field, in_field, surface_pressure, p, coordinates, nlayers,
        interpolation, extrapolation,
    )
    return out_field
end

@kernel inbounds = true function interpolate_pressure_layers_kernel!(
        out_field, in_field, surface_pressure, p, coordinates, nlayers,
        interpolation, extrapolation,
    )
    ij, k = @index(Global, NTuple)
    out_field[ij, k] = interpolate_pressure_layer(
        ij, p[k], in_field, surface_pressure[ij], coordinates, nlayers,
        interpolation, extrapolation,
    )
end

"""$(TYPEDSIGNATURES)
Value of `in_field` in column `ij` interpolated onto pressure `p` [Pa], given the surface
pressure `pₛ` [Pa] in that column and the vertical `coordinates` of the model.
Column-local, called from [`interpolate_pressure_layers!`](@ref) for every grid point and
pressure layer."""
@propagate_inbounds function interpolate_pressure_layer(
        ij, p, in_field, pₛ, coordinates, nlayers, interpolation, extrapolation,
    )
    # pressure of the top-most and bottom-most full layer in this column. `pressure`
    # dispatches over the vertical coordinates, σ[k]*pₛ for SigmaCoordinates and
    # A[k]*p_ref + B[k]*pₛ for SigmaPressureCoordinates, so both are covered here
    p_top = pressure(1, pₛ, coordinates)
    p_bottom = pressure(nlayers, pₛ, coordinates)

    p <= p_top && return extrapolate_above(in_field[ij, 1], p, p_top, extrapolation)
    p >= p_bottom && return extrapolate_below(in_field[ij, nlayers], p, p_bottom, pₛ, extrapolation)

    # find k with p_full(k-1) < p < p_full(k), starting at k=2 as p > p_top is guaranteed.
    # terminates at k = nlayers at the latest as p < p_bottom is guaranteed too
    k = 2
    p₂ = pressure(k, pₛ, coordinates)
    while p₂ < p
        k += 1
        p₂ = pressure(k, pₛ, coordinates)
    end
    p₁ = pressure(k - 1, pₛ, coordinates)

    w = interpolation_weight(p, p₁, p₂, interpolation)
    return (1 - w) * in_field[ij, k - 1] + w * in_field[ij, k]
end
