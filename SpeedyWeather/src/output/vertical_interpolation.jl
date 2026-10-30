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

# INTERPOLATION FROM MODEL LAYERS ONTO ISENTROPIC LAYERS

"""$(TYPEDSIGNATURES)
Dry potential temperature θ = T (p₀/p)^κ [K] of temperature `T` [K] at pressure `p` [Pa]."""
@inline dry_potential_temperature(T, p, κ, p₀) = T * (p₀ / p)^κ

"""$(TYPEDSIGNATURES)
Interpolate `in_field` on the model's vertical layers onto the isentropic layers `θ` [K]
(dry potential temperature), writing into the externally allocated `out_field`. All are
fields on the same horizontal grid, `in_field` and `temperature` have `nlayers` vertical
layers, `out_field` has `length(θ)` of them.

The dry potential temperature θ(k) = T(k) (p₀/p(k))^κ on every model layer is computed
from the absolute `temperature` [K], the `surface_pressure` [Pa] and the vertical
`coordinates` (`SigmaCoordinates` or `SigmaPressureCoordinates`) with `κ = R_dry/cₚ` and
the reference pressure `p₀` [Pa]. Between two layers θ is linear in ln(p), which gives the
pressure p* of the isentrope in closed form; every variable is then interpolated linearly
in ln(p) at p*, which is the same weight w = (θ - θ₁)/(θ₂ - θ₁). Each column is searched
from the top, so the first crossing from the top is used where θ is not monotonic.

Outside the column's θ range, p* is extrapolated with θ linear in ln(p) from the outermost
two layers and `extrapolation` decides at p*, see [`extrapolate_above`](@ref),
[`extrapolate_below`](@ref), as for pressure layers. If the outermost two layers are not
stably stratified, the Lorenz condition is used: p* = pₛ below, p* = p_top above. An
isentrope below ground (extrapolated p* > pₛ) also follows the Lorenz condition, p* = pₛ,
except that a `SubsurfaceMask` masks it.

Columns are independent, so this launches a kernel over the horizontal grid points and
the isentropic layers and runs on CPU and GPU. `θ` therefore has to be on the same
architecture as the fields. Requires at least 2 model layers."""
function interpolate_isentropic_layers!(
        out_field::AbstractField,                       # OUTPUT: (lonlat, nθ)
        in_field::AbstractField,                        # INPUT: (lonlat, nlayers)
        temperature::AbstractField,                     # INPUT: (lonlat, nlayers), [K]
        surface_pressure::AbstractField,                # INPUT: (lonlat), [Pa]
        θ::AbstractVector,                              # isentropic layers [K]
        coordinates::AbstractVerticalCoordinates,       # of the model
        κ::Real,                                        # R_dry/cₚ
        p₀::Real,                                       # reference pressure [Pa]
        extrapolation::AbstractVerticalExtrapolation = ConstantExtrapolation(),
    )
    grids_match(out_field, in_field, temperature, surface_pressure) ||
        throw(DimensionMismatch(out_field, in_field))
    size(out_field, 2) == length(θ) || throw(
        DimensionMismatch(
            "$(size(out_field, 2)) isentropic layers allocated in the output field, " *
                "but $(length(θ)) isentropic layers requested."
        )
    )
    nlayers = size(in_field, 2)
    size(temperature, 2) == nlayers || throw(
        DimensionMismatch(
            "Temperature has $(size(temperature, 2)) layers but the field to interpolate $nlayers."
        )
    )
    nlayers >= 2 || throw(ArgumentError("Interpolation onto isentropic layers requires at least 2 model layers."))

    # κ, p₀ in the number format of the temperature, so that the kernel stays in it
    T = eltype(temperature)
    launch!(
        architecture(in_field), RingGridWorkOrder, size(out_field),
        interpolate_isentropic_layers_kernel!,
        out_field, in_field, temperature, surface_pressure, θ, coordinates, nlayers,
        convert(T, κ), convert(T, p₀), extrapolation,
    )
    return out_field
end

@kernel inbounds = true function interpolate_isentropic_layers_kernel!(
        out_field, in_field, temperature, surface_pressure, θ, coordinates, nlayers,
        κ, p₀, extrapolation,
    )
    ij, k = @index(Global, NTuple)
    out_field[ij, k] = interpolate_isentropic_layer(
        ij, θ[k], in_field, temperature, surface_pressure[ij], coordinates, nlayers,
        κ, p₀, extrapolation,
    )
end

"""$(TYPEDSIGNATURES)
Value of `in_field` in column `ij` interpolated onto the isentrope `θ` [K], given the
`temperature` [K], the surface pressure `pₛ` [Pa] in that column and the vertical
`coordinates` of the model. Column-local, called from
[`interpolate_isentropic_layers!`](@ref) for every grid point and isentropic layer."""
@propagate_inbounds function interpolate_isentropic_layer(
        ij, θ, in_field, temperature, pₛ, coordinates, nlayers, κ, p₀, extrapolation,
    )
    # search from the top for the first pair of layers k-1, k with θ between θ(k-1) and θ(k)
    p₁ = pressure(1, pₛ, coordinates)
    θ₁ = dry_potential_temperature(temperature[ij, 1], p₁, κ, p₀)
    for k in 2:nlayers
        p₂ = pressure(k, pₛ, coordinates)
        θ₂ = dry_potential_temperature(temperature[ij, k], p₂, κ, p₀)
        if min(θ₁, θ₂) <= θ <= max(θ₁, θ₂)
            # θ linear in ln(p): ln(p*) = ln(p₁) + w*(ln(p₂) - ln(p₁)), variables linear in
            # ln(p) at p* have the same weight w. θ₁ == θ₂ (== θ) takes layer k-1, avoids 0/0
            w = ifelse(θ₁ == θ₂, zero(θ₁), (θ - θ₁) / (θ₂ - θ₁))
            return (1 - w) * in_field[ij, k - 1] + w * in_field[ij, k]
        end
        p₁, θ₁ = p₂, θ₂
    end

    # no crossing: θ is outside the column's θ range, above the top or below the lowest layer
    p_top = pressure(1, pₛ, coordinates)
    θ_top = dry_potential_temperature(temperature[ij, 1], p_top, κ, p₀)

    if θ > θ_top    # above the top-most layer
        p₂ = pressure(2, pₛ, coordinates)
        θ₂ = dry_potential_temperature(temperature[ij, 2], p₂, κ, p₀)
        # extrapolate p* with θ linear in ln(p) from the top two layers if stably stratified
        # (θ decreasing with p), otherwise the Lorenz condition p* = p_top
        p = θ_top > θ₂ ? exp(log(p_top) + (θ - θ_top) * log(p_top / p₂) / (θ_top - θ₂)) : p_top
        return extrapolate_above(in_field[ij, 1], p, p_top, extrapolation)
    else            # below the lowest layer
        p_bottom = pressure(nlayers, pₛ, coordinates)
        θ_bottom = dry_potential_temperature(temperature[ij, nlayers], p_bottom, κ, p₀)
        p₁ = pressure(nlayers - 1, pₛ, coordinates)
        θ₁ = dry_potential_temperature(temperature[ij, nlayers - 1], p₁, κ, p₀)
        # extrapolate p* with θ linear in ln(p) from the lowest two layers if stably
        # stratified (θ decreasing with p), otherwise the Lorenz condition p* = pₛ
        p = θ₁ > θ_bottom ? exp(log(p_bottom) + (θ - θ_bottom) * log(p_bottom / p₁) / (θ_bottom - θ₁)) : pₛ
        return extrapolate_isentrope_below(in_field[ij, nlayers], p, p_bottom, pₛ, extrapolation)
    end
end

"""$(TYPEDSIGNATURES)
Value at the pressure `p` [Pa] of an isentrope below the bottom-most full model layer at
`p_bottom` with value `ξ_bottom`, given the surface pressure `pₛ`. An isentrope below
ground (`p > pₛ`) follows the Lorenz condition and is evaluated at the surface, `p = pₛ`,
as an extrapolated `p` can be arbitrarily large for a weakly stratified lowest layer pair.
Between the lowest layer and the surface `extrapolate_below` applies at `p`."""
@inline extrapolate_isentrope_below(ξ_bottom, p, p_bottom, pₛ, E::AbstractVerticalExtrapolation) =
    extrapolate_below(ξ_bottom, min(p, pₛ), p_bottom, pₛ, E)

# a SubsurfaceMask still masks where the extrapolated p is below ground
@inline extrapolate_isentrope_below(ξ_bottom, p, p_bottom, pₛ, E::SubsurfaceMask) =
    p > pₛ ? convert(typeof(ξ_bottom), E.missing_value) :
    extrapolate_below(ξ_bottom, p, p_bottom, pₛ, E.above_surface)
