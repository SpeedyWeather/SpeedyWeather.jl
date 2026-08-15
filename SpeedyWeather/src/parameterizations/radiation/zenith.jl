abstract type AbstractSolarDeclination <: AbstractModelComponent end
abstract type AbstractSolarTimeCorrection <: AbstractModelComponent end
abstract type AbstractZenith <: AbstractModelComponent end

# set these as constants as hour/year angle need them but
# they should not be seen as variable as we 
const LENGTH_OF_DAY = Second(EARTH_DAY).value
const LENGTH_OF_YEAR = Second(EARTH_YEAR).value

"""Coefficients to calculate the solar declination angle δ [radians] based on a simple
sine function, with Earth's axial tilt as amplitude, equinox as phase shift.
$(TYPEDFIELDS)"""
@kwdef struct SinSolarDeclination{NF, DT} <: AbstractSolarDeclination
    axial_tilt::NF = AXIAL_TILT
    equinox::DT = EARTH_EQUINOX
end

SinSolarDeclination(SG::SpectralGrid; kwargs...) = SinSolarDeclination{SG.NF}(; kwargs...)
SinSolarDeclination(P::AbstractPlanet) = SinSolarDeclination(P.axial_tilt, P.equinox)

"""
$(TYPEDSIGNATURES)
SinSolarDeclination functor, computing the solar declination angle of
angular fraction of year g [radians] using the coefficients of the
SinSolarDeclination struct. Uses LENGTH_OF_DAY and LENGTH_OF_YEAR as
Earth's constants and not from planet as Dates functionality here
is relative to the Earth calendar. The year angle `g` should reflect
slower/faster seasons."""
function (S::SinSolarDeclination)(g::NF) where {NF}
    axial_tilt = deg2rad(S.axial_tilt)
    equinox = LENGTH_OF_DAY * Dates.dayofyear(S.equinox) / LENGTH_OF_YEAR
    return axial_tilt * sin(g - 2 * (π * convert(NF, equinox)))
end

"""Coefficients to calculate the solar declination angle δ from

    δ = 0.006918 - 0.399912*cos(g)  + 0.070257*sin(g)
                 - 0.006758*cos(2g) + 0.000907*sin(2g)
                 - 0.002697*cos(3g) + 0.001480*sin(3g)

with g the angular fraction of the year in radians. Following Spencer 1971,
Fourier series representation of the position of the sun. Search 2(5):172.

Note that this Declination type does not support an axial tilt different
from Earth's 23.4 degree. Use SinSolarDeclination instead.
$(TYPEDFIELDS)"""
@kwdef struct SolarDeclination{NF} <: AbstractSolarDeclination
    a::NF = 0.006918      # the offset +
    s1::NF = 0.070257     # s1*sin(g) +
    c1::NF = -0.399912    # c1*cos(g) +
    s2::NF = 0.000907     # s2*sin(2g) +
    c2::NF = -0.006758    # c2*cos(2g) +
    s3::NF = 0.00148      # s3*sin(3g) +
    c3::NF = -0.002697    # c3*cos(3g)
end

"""Generator function pulling the number format NF from a SpectralGrid."""
SolarDeclination(SG::SpectralGrid; kwargs...) = SolarDeclination{SG.NF}(; kwargs...)
SolarDeclination(P::AbstractPlanet; kwargs...) = SolarDeclination{typeof(P.axial_tilt)}(; kwargs...)

"""
$(TYPEDSIGNATURES)
SolarDeclination functor, computing the solar declination angle of
angular fraction of year g [radians] using the coefficients of the
SolarDeclination struct."""
function (SD::SolarDeclination)(g)
    (; a, s1, s2, s3, c1, c2, c3) = SD
    sin1g, cos1g = sincos(g)
    sin2g, cos2g = sincos(2g)
    sin3g, cos3g = sincos(3g)
    return a + s1 * sin1g + c1 * cos1g + s2 * sin2g + c2 * cos2g + s3 * sin3g + c3 * cos3g
end

"""Coefficients for the solar time correction (also called
Equation of time) which adjusts the solar hour to an oscillation
of sunrise/set by about +-16min throughout the year."""
@kwdef struct SolarTimeCorrection{NF <: AbstractFloat} <: AbstractSolarTimeCorrection
    a::NF = 0.004297      # the offset +
    s1::NF = -1.837877    # s1*sin(g) +
    c1::NF = 0.107029     # c1*cos(g) +
    s2::NF = -2.340475    # s2*sin(2g) +
    c2::NF = -0.837378    # c2*cos(2g)
end

Adapt.@adapt_structure SolarTimeCorrection
SolarTimeCorrection(SG::SpectralGrid; kwargs...) = SolarTimeCorrection{SG.NF}(; kwargs...)

"""
$(TYPEDSIGNATURES)
Functor that returns the time correction for a angular
fraction of the year g [radians], so that g=0 for Jan-01 and g=2π for Dec-31."""
function (STC::SolarTimeCorrection)(g)
    (; a, s1, s2, c1, c2) = STC
    sin1g, cos1g = sincos(g)
    sin2g, cos2g = sincos(2g)
    return deg2rad(a + s1 * sin1g + c1 * cos1g + s2 * sin2g + c2 * cos2g)
end

struct NoTimeCorrection <: AbstractSolarTimeCorrection end
Adapt.@adapt_structure NoTimeCorrection
NoTimeCorrection(::SpectralGrid; kwargs...) = NoTimeCorrection()
(STC::NoTimeCorrection)(g) = zero(g)

"""$(TYPEDSIGNATURES)
Chooses from SolarZenith (daily and seasonal cycle) or SolarZenithSeason
given the parameters in model.planet. In both cases the seasonal cycle can
be disabled, calculating the solar declination from the initial time
instead of current (orbit) time."""
function WhichZenith(
    SG::SpectralGrid,
    planet::AbstractPlanet;
    solar_declination = SinSolarDeclination(planet), 
    time_correction = SolarTimeCorrection(SG),
)
    Z = planet.daily_cycle ? SolarZenith : SolarZenithSeason
    return Z(SG; planet.seasonal_cycle, solar_declination, time_correction)
end

# function barrier
function parameterization!(vars::Variables, zenith::AbstractZenith, model::PrimitiveEquation)
    (; rotation_time, orbit_time) = vars.prognostic.clock
    (; geometry) = model
    (; cos_zenith) = vars.parameterizations
    cos_zenith!(cos_zenith, zenith, rotation_time, orbit_time, geometry)
    return nothing
end

# as a global parameterization define the column parameterization as doing nothing
parameterization!(ij, vars::Variables, zenith::AbstractZenith, model) = nothing

export SolarZenith

"""Solar zenith angle varying with daily and seasonal cycle.
$(TYPEDFIELDS)"""
@kwdef struct SolarZenith{SD, STC, RV} <: AbstractZenith
    "[OPTION] Seasonal cycle? Otherwise daily cycle only."
    seasonal_cycle::Bool = true

    "[OPTION] Calculate seasonal cycle from this solar declination function."
    solar_declination::SD
    
    "[DERIVED] Use this time, set via `initialize!(model, time=...), to hold seasonal cycle fixed."
    initial_time::RV = Ref(DEFAULT_DATE)
    
    "[OPTION] Time correction for seasonal wobble of sunset/sunrise times"
    time_correction::STC
end

Adapt.@adapt_structure SolarZenith
SolarZenith(SG::SpectralGrid; kwargs...) =
    SolarZenith(; solar_declination = SolarDeclination(SG), time_correction = SolarTimeCorrection(SG), kwargs...)
Base.show(io::IO, Z::AbstractZenith) = show(io, Z, values=false)

function variables(::AbstractZenith)
    return (
        ParameterizationVariable(:cos_zenith, Grid2D(), desc = "Cosine of solar zenith angle", units = "1"),
    )
end

function initialize!(
        S::AbstractZenith,
        initial_time::DateTime,
        model::AbstractModel
    )
    return S.initial_time[] = initial_time     # to fix the season if no seasonal cycle
end

"""$(TYPEDSIGNATURES)
Fraction of year as angle in radians [0...2π]. Always calculated relative to the
Earth calendar, so LENGTH_OF_DAY, LENGTH_OF_YEAR are used as constants.
This is because `Dates` functions assume an Earth calendar.
Length of the seasonal cycle is instead controlled via a faster orbit time,
see `Earth`."""
function year_angle(::Type{T}, time::DateTime) where {T}
    year2rad = T(2π) / LENGTH_OF_YEAR
    sec_of_day = Dates.second(Dates.Time(time).instant)     # use secondofday?
    return year2rad * (Dates.dayofyear(time) * LENGTH_OF_DAY + sec_of_day)
end

"""$(TYPEDSIGNATURES)
Fraction of day in `time` as angle in radians [0...2π], noon to noon, at longitude `λ`.
Always calculated relative to the Earth calendar, so LENGTH_OF_DAY is used as constant.
This is because `Dates` functions assume an Earth date. Length of the daily cycle is instead
controlled via a faster rotation time, see `Earth`."""
function solar_hour_angle(::Type{T}, time::DateTime, λ=0) where {T}
    day2rad = T(2π) / LENGTH_OF_DAY
    noon_in_sec = LENGTH_OF_DAY ÷ 2
    sec_of_day = Dates.second(Dates.Time(time).instant)     # use secondofday?
    return (sec_of_day - noon_in_sec) * day2rad + convert(T, λ)
end

"""
$(TYPEDSIGNATURES)
Calculate cos of solar zenith angle with a daily cycle
at rotation_time `rotation_time` and orbit_time `orbit_time`.
Seasonal cycle or time correction may be disabled,
depending on parameters in SolarZenith."""
function cos_zenith!(
        cos_zenith::AbstractField{NF},
        S::SolarZenith,
        rotation_time::DateTime,
        orbit_time::DateTime,
        geometry::AbstractGeometry,
    ) where {NF}

    (; sinlat, coslat, lons) = geometry
    @boundscheck geometry.spectral_grid.grid == cos_zenith.grid ||
        throw(DimensionMismatch(geometry.spectral_grid.grid, cos_zenith.grid))

    # g: angular fraction of year [0...2π] for Jan-01 to Dec-31
    time_of_year = S.seasonal_cycle ? orbit_time : S.initial_time[]
    g = year_angle(NF, time_of_year)

    # time correction [radians] due to the equation of time (sunrise/set oscillation)
    tc = S.time_correction(g)

    # solar hour angle at 0˚E (longitude offset added later)
    λ = 0
    solar_hour_angle_0E = solar_hour_angle(NF, rotation_time, λ) + tc

    # solar declination angle [radians] changing from tropic of cancer to capricorn
    # throughout the year measured by g [radians]
    δ = S.solar_declination(g)
    sinδ, cosδ = sincos(δ)

    # Launch kernel for solar zenith calculation
    launch!(
        architecture(cos_zenith), LinearWorkOrder, size(cos_zenith), solar_zenith_kernel!,
        cos_zenith, solar_hour_angle_0E, sinδ, cosδ, sinlat, coslat, lons, cos_zenith.grid.whichring
    )
    return nothing
end

# Kernel for solar zenith calculation with daily cycle
@kernel inbounds = true function solar_zenith_kernel!(
        cos_zenith,
        solar_hour_angle_0E, sinδ, cosδ, sinlat, coslat, lons, whichring
    )

    ij = @index(Global, Linear)
    j = whichring[ij]
    
    sinδsinϕ = sinδ * sinlat[j]
    cosδcosϕ = cosδ * coslat[j]
    h = solar_hour_angle_0E + lons[ij]      # solar hour angle at longitude λ in radians
    cos_zenith[ij] = max(0, sinδsinϕ + cosδcosϕ * cos(h))
end

export SolarZenithSeason

"""Solar zenith angle varying with seasonal cycle only.
$(TYPEDFIELDS)"""
@kwdef struct SolarZenithSeason{SD, RV} <: AbstractZenith
    "[OPTION] Seasonal cycle? Otherwise daily cycle only."
    seasonal_cycle::Bool = true

    "[OPTION] Calculate seasonal cycle from this solar declination function."
    solar_declination::SD
    
    "[DERIVED] Use this time, set via `initialize!(model, time=...), to hold seasonal cycle fixed."
    initial_time::RV = Ref(DEFAULT_DATE)
end

# constructor, add unused `time_correction` kwarg to mirror constructor for SolarZenith
SolarZenithSeason(SG::SpectralGrid; solar_declination = SolarDeclination(SG), time_correction = nothing, kwargs...) =
    SolarZenithSeason(; solar_declination, kwargs...)

Adapt.@adapt_structure SolarZenithSeason

"""
$(TYPEDSIGNATURES)
Calculate cos of solar zenith angle as daily average
at rotation_time `rotation_time` and orbit_time `orbit_time`.
Seasonal cycle or time correction may be disabled,
depending on parameters in SolarZenithSeason."""
function cos_zenith!(
        cos_zenith::AbstractField{NF},
        S::SolarZenithSeason,
        rotation_time::DateTime,        # not used, but to keep the same function barrier as SolarZenith
        orbit_time::DateTime,
        geometry::AbstractGeometry,
    ) where {NF}

    (; sinlat, coslat, lat) = geometry
    @boundscheck geometry.spectral_grid.grid == cos_zenith.grid ||
        throw(DimensionMismatch(geometry.spectral_grid.grid, cos_zenith.grid))

    # g: angular fraction of year [0...2π] for Jan-01 to Dec-31
    time_of_year = S.seasonal_cycle ? orbit_time : S.initial_time[]
    g = year_angle(NF, time_of_year)

    # solar declination angle [radians] changing from tropic of cancer to capricorn
    # throughout the year measured by g [radians]
    δ = S.solar_declination(g)
    sinδ, cosδ = sincos(δ)

    # Launch kernel for seasonal solar zenith calculation
    launch!(
        architecture(cos_zenith), LinearWorkOrder, size(cos_zenith), solar_zenith_season_kernel!,
        cos_zenith, δ, sinδ, cosδ, sinlat, coslat, lat, cos_zenith.grid.whichring
    )
    return nothing
end

# Kernel for seasonal solar zenith calculation (daily average)
@kernel inbounds = true function solar_zenith_season_kernel!(
        cos_zenith,
        δ, sinδ, cosδ, sinlat, coslat, lat, whichring
    )

    ij = @index(Global, Linear)
    j = whichring[ij]
    NF = eltype(cos_zenith)

    ϕ = lat[j]
    h₀ = NF(ifelse(2*(abs(δ) + abs(ϕ)) < π,     # polar day/night?
        acos(clamp(-tan(ϕ) * tan(δ), -one(NF), one(NF))),  # length of day
        ifelse(ϕ * δ > 0, π, 0)))               # polar day if signs are equal, otherwise polar night

    sinϕ, cosϕ = sinlat[j], coslat[j]
    cos_zenith_j = h₀ * sinδ * sinϕ + cosδ * cosϕ * sin(h₀)
    cos_zenith[ij] = cos_zenith_j / π
end