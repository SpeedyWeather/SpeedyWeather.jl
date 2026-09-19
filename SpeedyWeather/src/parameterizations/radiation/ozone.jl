abstract type AbstractOzone <: AbstractModelComponent end

export NoOzone

"""No ozone, the shortwave radiation is not absorbed in the stratosphere by ozone."""
struct NoOzone <: AbstractOzone end
Adapt.@adapt_structure NoOzone
NoOzone(::SpectralGrid) = NoOzone()
initialize!(::NoOzone, ::AbstractModel) = nothing
ozone_absorption!(vars, ::NoOzone, model) = nothing

"""$(TYPEDSIGNATURES)
Fraction of the top-of-atmosphere shortwave flux absorbed by ozone in layer `k`, zero for `NoOzone`."""
@propagate_inbounds ozone_absorption(ij, k, vars, ::NoOzone, model) =
    zero(eltype(vars.parameterizations.cos_zenith))

export SeasonalOzone

"""Ozone absorption of shortwave radiation following Fortran SPEEDY (Molteni, 2003).
A fraction of the top-of-atmosphere (TOA) shortwave flux is absorbed in two stratospheric
layers: The upper stratosphere absorbs `upper_fraction * absorption` everywhere, the lower
stratosphere absorbs a fraction that depends on latitude ϕ and season

    lower_fraction * absorption * (1 + a_s * max(0, -δ/|δₘₐₓ|) * sin(ϕ) + a_ϕ * (3/2 sin²(ϕ) - 1/2))

with `a_s = seasonal_amplitude` and `a_ϕ = latitudinal_amplitude`. The seasonal cycle follows
the solar declination δ of the solar zenith angle (`model.solar_zenith`) normalized by the
planet's axial tilt δₘₐₓ, so it is synchronized with the orbit (`length_of_year`, `equinox`,
`axial_tilt`, `seasonal_cycle`) of the planet. For Earth's sinusoidal declination
`-δ/δₘₐₓ = cos(α)` with α the angle of the year from the northern winter solstice as in SPEEDY.
More ozone absorption therefore in high latitudes and in the northern hemisphere during
its winter half-year. The absorbed flux is further multiplied
by the zenith angle correction factor of the shortwave transmissivity. The absorption is
distributed on the model layers following their overlap with the σ-intervals
[0, `σ_upper`] and [`σ_upper`, `σ_lower`]. Fields are
$(TYPEDFIELDS)"""
@parameterized @kwdef struct SeasonalOzone{NF} <: AbstractOzone
    "[OPTION] Reference fraction of TOA shortwave flux absorbed by ozone (SPEEDY epssw) [1]"
    @param absorption::NF = 0.02 (bounds = 0 .. 1,)

    "[OPTION] Fraction of `absorption` in the upper stratosphere [1]"
    @param upper_fraction::NF = 0.5 (bounds = Nonnegative,)

    "[OPTION] Fraction of `absorption` in the lower stratosphere (at the equator, excluding latitudinal and seasonal terms) [1]"
    @param lower_fraction::NF = 0.4 (bounds = Nonnegative,)

    "[OPTION] Amplitude of hemispheric seasonal asymmetry in lower stratospheric ozone (SPEEDY coz1) [1]"
    @param seasonal_amplitude::NF = 1 (bounds = Nonnegative,)

    "[OPTION] Amplitude of the equator-to-pole contrast of lower stratospheric ozone (SPEEDY coz2) [1]"
    @param latitudinal_amplitude::NF = 1.8 (bounds = Nonnegative,)

    "[OPTION] Lower σ boundary of the upper stratosphere [1]"
    σ_upper::NF = 0.05

    "[OPTION] Lower σ boundary of the lower stratosphere [1]"
    σ_lower::NF = 0.14
end

Adapt.@adapt_structure SeasonalOzone
SeasonalOzone(SG::SpectralGrid; kwargs...) = SeasonalOzone{SG.NF}(; kwargs...)
initialize!(::SeasonalOzone, ::AbstractModel) = nothing

variables(::SeasonalOzone) = (
    ParameterizationVariable(:ozone_absorption_lower, Grid2D(), desc = "Fraction of TOA shortwave absorbed by ozone in lower stratosphere", units = "1"),
)

"""$(TYPEDSIGNATURES)
Update the latitude- and season-dependent ozone absorption in the lower stratosphere,
with the seasonal cycle synchronized to the solar declination of `model.solar_zenith`."""
function ozone_absorption!(vars, ozone::SeasonalOzone, model)
    field = vars.parameterizations.ozone_absorption_lower
    NF = eltype(field)

    # seasonal cycle synchronized with the solar declination δ, normalized by the axial tilt,
    # 0 for the northern summer half-year, 1 at the northern winter solstice
    δ = solar_declination(NF, model.solar_zenith, vars.prognostic.clock.orbit_time)
    δₘₐₓ = abs(deg2rad(convert(NF, model.planet.axial_tilt)))
    northern_winter = δₘₐₓ > 0 ? clamp(-δ / δₘₐₓ, 0, 1) : zero(NF)
    seasonal_term = ozone.seasonal_amplitude * northern_winter
    lower = ozone.absorption * ozone.lower_fraction

    launch!(
        architecture(field), LinearWorkOrder, size(field), ozone_absorption_kernel!,
        field, lower, seasonal_term, ozone.latitudinal_amplitude, model.geometry.sinlat, field.grid.whichring
    )
    return nothing
end

@kernel inbounds = true function ozone_absorption_kernel!(field, lower, seasonal_term, latitudinal_amplitude, sinlat, whichring)
    ij = @index(Global, Linear)
    sinϕ = sinlat[whichring[ij]]
    P₂ = (3 * sinϕ^2 - 1) / 2       # second Legendre polynomial, equator-to-pole contrast
    field[ij] = lower * max(0, 1 + seasonal_term * sinϕ + latitudinal_amplitude * P₂)
end

"""$(TYPEDSIGNATURES)
Fraction of the top-of-atmosphere shortwave flux absorbed by ozone in layer `k`
(before zenith angle correction), distributing the upper and lower stratospheric
absorption by the overlap of layer `k` with the respective σ-intervals."""
@propagate_inbounds function ozone_absorption(ij, k, vars, ozone::SeasonalOzone, model)
    σ_half = model.geometry.σ_levels_half
    σ_top, σ_bottom = σ_half[k], σ_half[k + 1]
    (; σ_upper, σ_lower) = ozone

    # fractions of upper [0, σ_upper] and lower [σ_upper, σ_lower] stratosphere in layer k
    upper_overlap = max(0, min(σ_bottom, σ_upper) - σ_top) / σ_upper
    lower_overlap = max(0, min(σ_bottom, σ_lower) - max(σ_top, σ_upper)) / (σ_lower - σ_upper)

    upper = ozone.absorption * ozone.upper_fraction
    lower = vars.parameterizations.ozone_absorption_lower[ij]
    return upper * upper_overlap + lower * lower_overlap
end
