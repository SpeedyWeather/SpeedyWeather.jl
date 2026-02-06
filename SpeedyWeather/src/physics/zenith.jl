abstract type AbstractSolarDeclination end
abstract type AbstractSolarTimeCorrection end
abstract type AbstractZenith end

"""Coefficients to calculate the solar declination angle δ [radians] based on a simple
sine function, with Earth's axial tilt as amplitude, equinox as phase shift.
$(TYPEDFIELDS)"""
struct SinSolarDeclination{P} <: AbstractSolarDeclination
    planet::P
end

"""
$(TYPEDSIGNATURES)
SinSolarDeclination functor, computing the solar declination angle of
angular fraction of year g [radians] using the coefficients of the
SinSolarDeclination struct."""
function (S::SinSolarDeclination)(g::NF) where {NF}
    planet = S.planet
    axial_tilt = deg2rad(planet.axial_tilt)
    equinox = planet.length_of_day.value * Dates.dayofyear(planet.equinox) / planet.length_of_year.value
    return axial_tilt * sin(g - 2 * (π * convert(NF, equinox)))
end

"""Coefficients to calculate the solar declination angle δ from

    δ = 0.006918 - 0.399912*cos(g)  + 0.070257*sin(g)
                 - 0.006758*cos(2g) + 0.000907*sin(2g)
                 - 0.002697*cos(3g) + 0.001480*sin(3g)

with g the angular fraction of the year in radians. Following Spencer 1971,
Fourier series representation of the position of the sun. Search 2(5):172.
$(TYPEDFIELDS)"""
@parameterized Base.@kwdef struct SolarDeclination{NF <: AbstractFloat} <: AbstractSolarDeclination
    @param a::NF = 0.006918      # the offset +
    @param s1::NF = 0.070257     # s1*sin(g) +
    @param c1::NF = -0.399912    # c1*cos(g) +
    @param s2::NF = 0.000907     # s2*sin(2g) +
    @param c2::NF = -0.006758    # c2*cos(2g) +
    @param s3::NF = 0.00148      # s3*sin(3g) +
    @param c3::NF = -0.002697    # c3*cos(3g)
end

"""Generator function pulling the number format NF from a SpectralGrid."""
SolarDeclination(SG::SpectralGrid; kwargs...) = SolarDeclination{SG.NF}(; kwargs...)

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

function Base.show(io::IO, L::AbstractSolarDeclination)
    println(io, "$(typeof(L)) <: AbstractSolarDeclination")
    keys = propertynames(L)
    return print_fields(io, L, keys)
end

"""Coefficients for the solar time correction (also called
Equation of time) which adjusts the solar hour to an oscillation
of sunrise/set by about +-16min throughout the year."""
@parameterized Base.@kwdef struct SolarTimeCorrection{NF <: AbstractFloat} <: AbstractSolarTimeCorrection
    @param a::NF = 0.004297      # the offset +
    @param s1::NF = -1.837877    # s1*sin(g) +
    @param c1::NF = 0.107029     # c1*cos(g) +
    @param s2::NF = -2.340475    # s2*sin(2g) +
    @param c2::NF = -0.837378    # c2*cos(2g)
end

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

function Base.show(io::IO, L::AbstractSolarTimeCorrection)
    println(io, "$(typeof(L)) <: AbstractSolarTimeCorrection")
    keys = propertynames(L)
    return print_fields(io, L, keys)
end

"""
$(TYPEDSIGNATURES)
Chooses from SolarZenith (daily and seasonal cycle) or SolarZenithSeason
given the parameters in model.planet. In both cases the seasonal cycle can
be disabled, calculating the solar declination from the initial time
instead of current time."""
function WhichZenith(SG::SpectralGrid, P::AbstractPlanet; kwargs...)
    (; NF) = SG
    (; daily_cycle, seasonal_cycle, length_of_day, length_of_year) = P
    solar_declination = SinSolarDeclination(P)

    if daily_cycle
        return SolarZenith{NF, typeof(solar_declination)}(;
            length_of_day, length_of_year, solar_declination, seasonal_cycle, kwargs...
        )

    else
        return SolarZenithSeason{NF, typeof(solar_declination)}(;
            length_of_day, length_of_year, solar_declination, seasonal_cycle, kwargs...
        )
    end
end

# function barrier
function cos_zenith!(diagn::DiagnosticVariables, time::DateTime, model::PrimitiveEquation)
    (; solar_zenith, geometry) = model
    (; cos_zenith) = diagn.physics
    rotation_time = time
    orbit_time = time
    return cos_zenith!(cos_zenith, solar_zenith, rotation_time, orbit_time, geometry)
end

function Base.show(io::IO, L::AbstractZenith)
    println(io, "$(typeof(L)) <: AbstractZenith")
    keys = propertynames(L)
    return print_fields(io, L, keys)
end

export SolarZenith

"""Solar zenith angle varying with daily and seasonal cycle.
$(TYPEDFIELDS)"""
@parameterized @kwdef struct SolarZenith{NF <: AbstractFloat, SD <: AbstractSolarDeclination} <: AbstractZenith
    # OPTIONS
    length_of_day::Second = Hour(24)
    length_of_year::Second = Day(365.25)
    equation_of_time::Bool = true
    seasonal_cycle::Bool = true

    # COEFFICIENTS
    @param solar_declination::SD = SinSolarDeclination(Earth{NF}()) (group = :solar_declination,)
    @param time_correction::SolarTimeCorrection{NF} = SolarTimeCorrection{NF}() (group = :time_correction,)

    initial_time::Base.RefValue{DateTime} = Ref(DEFAULT_DATE)
end

SolarZenith(SG::SpectralGrid; kwargs...) = SolarZenith{SG.NF, SinSolarDeclination{Earth{SG.NF}}}(; kwargs...)

function initialize!(
        S::AbstractZenith,
        initial_time::DateTime,
        model::AbstractModel
    )
    return S.initial_time[] = initial_time     # to fix the season if no seasonal cycle
end

"""
$(TYPEDSIGNATURES)
Fraction of year as angle in radians [0...2π].
TODO: Takes length of day/year as argument, but calls to Dates.Time(), Dates.dayofyear()
currently have these hardcoded."""
function year_angle(::Type{T}, time::DateTime, length_of_day::Second, length_of_year::Second) where {T}
    year2rad = convert(T, 2π / length_of_year.value)
    sec_of_day = Dates.second(Dates.Time(time).instant)
    return year2rad * (Dates.dayofyear(time) * length_of_day.value + sec_of_day)
end

"""
$(TYPEDSIGNATURES)
Fraction of day as angle in radians [0...2π].
TODO: Takes length of day as argument, but a call to Dates.Time()
currently have this hardcoded anyway."""
function solar_hour_angle(
        ::Type{T},
        time::DateTime,
        λ,                      # longitude in radians
        length_of_day::Second
    ) where {T}
    day2rad = convert(T, 2π / length_of_day.value)
    noon_in_sec = length_of_day.value ÷ 2
    sec_of_day = Dates.second(Dates.Time(time).instant)
    return (sec_of_day - noon_in_sec) * day2rad + convert(T, λ)
end

"""
$(TYPEDSIGNATURES)
Calculate cos of solar zenith angle with a daily cycle
at rotation_time `rotation_time` and orbit_time `orbit_time`.
Seasonal cycle or time correction may be disabled,
depending on parameters in SolarZenith."""
function cos_zenith!(
        cos_zenith::AbstractField,
        S::SolarZenith,
        rotation_time::DateTime,
        orbit_time::DateTime,
        geometry::AbstractGeometry,
    )
    NF = eltype(cos_zenith)
    (; sinlat, coslat, lons) = geometry
    (; length_of_day, length_of_year) = S
    @boundscheck geometry.spectral_grid.grid == cos_zenith.grid ||
        throw(DimensionMismatch(geometry.spectral_grid.grid, cos_zenith.grid))

    # g: angular fraction of year [0...2π] for Jan-01 to Dec-31
    time_of_year = S.seasonal_cycle ? orbit_time : S.initial_time[]
    g = year_angle(NF, time_of_year, length_of_day, length_of_year)

    # time correction [radians] due to the equation of time (sunrise/set oscillation)
    tc = S.equation_of_time ? S.time_correction(g) : zero(NF)

    # solar hour angle at 0˚E (longtiude offset added later)
    λ = 0
    solar_hour_angle_0E = solar_hour_angle(NF, rotation_time, λ, length_of_day) + tc

    # solar declination angle [radians] changing from tropic of cancer to capricorn
    # throughout the year measured by g [radians]
    δ = S.solar_declination(g)
    sinδ, cosδ = sincos(δ)

    # Launch kernel for solar zenith calculation
    return launch!(
        architecture(cos_zenith), LinearWorkOrder, size(cos_zenith), solar_zenith_kernel!,
        cos_zenith, solar_hour_angle_0E, sinδ, cosδ, sinlat, coslat, lons, cos_zenith.grid.whichring
    )
end

# Kernel for solar zenith calculation with daily cycle
@kernel inbounds = true function solar_zenith_kernel!(
        cos_zenith,
        @Const(solar_hour_angle_0E), @Const(sinδ), @Const(cosδ), @Const(sinlat), @Const(coslat), @Const(lons), @Const(whichring)
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
@parameterized @kwdef struct SolarZenithSeason{NF <: AbstractFloat, SD <: AbstractSolarDeclination} <: AbstractZenith
    # OPTIONS
    length_of_day::Second = Hour(24)
    length_of_year::Second = Day(365.25)
    seasonal_cycle::Bool = true

    # COEFFICIENTS
    @param solar_declination::SD = SinSolarDeclination(Earth{NF}()) (group = :solar_declination,)

    initial_time::Base.RefValue{DateTime} = Ref(DEFAULT_DATE)
end

SolarZenithSeason(SG::SpectralGrid; kwargs...) = SolarZenithSeason{SG.NF, SinSolarDeclination{Earth{SG.NF}}}(; kwargs...)

"""
$(TYPEDSIGNATURES)
Calculate cos of solar zenith angle as daily average
at rotation_time `rotation_time` and orbit_time `orbit_time`.
Seasonal cycle or time correction may be disabled,
depending on parameters in SolarZenithSeason."""
function cos_zenith!(
        cos_zenith::AbstractField,
        S::SolarZenithSeason,
        rotation_time::DateTime,
        orbit_time::DateTime,
        geometry::AbstractGeometry,
    )
    NF = eltype(cos_zenith)
    (; sinlat, coslat, lat) = geometry
    (; length_of_day, length_of_year) = S
    @boundscheck geometry.spectral_grid.grid == cos_zenith.grid ||
        throw(DimensionMismatch(geometry.spectral_grid.grid, cos_zenith.grid))

    # g: angular fraction of year [0...2π] for Jan-01 to Dec-31
    time_of_year = S.seasonal_cycle ? orbit_time : S.initial_time[]
    g = year_angle(NF, time_of_year, length_of_day, length_of_year)

    # solar declination angle [radians] changing from tropic of cancer to capricorn
    # throughout the year measured by g [radians]
    δ = S.solar_declination(g)
    sinδ, cosδ = sincos(δ)

    # Launch kernel for seasonal solar zenith calculation
    return launch!(
        architecture(cos_zenith), LinearWorkOrder, size(cos_zenith), solar_zenith_season_kernel!,
        cos_zenith, δ, sinδ, cosδ, sinlat, coslat, lat, cos_zenith.grid.whichring
    )
end

# Kernel for seasonal solar zenith calculation (daily average)
@kernel inbounds = true function solar_zenith_season_kernel!(
        cos_zenith,
        @Const(δ), @Const(sinδ), @Const(cosδ), @Const(sinlat), @Const(coslat), @Const(lat), @Const(whichring)
    )

    ij = @index(Global, Linear)
    j = whichring[ij]

    NF = eltype(cos_zenith)         # force type stability
    local h₀::NF                    # hour angle sunrise to sunset
    local cos_zenith_j::NF          # at latitude j

    ϕ = lat[j]
    h₀ = abs(δ) + abs(ϕ) < π / 2 ?    # polar day/night?
        acos(-tan(ϕ) * tan(δ)) :   # if not: calculate length of day
        ϕ * δ > 0 ? π : 0            # polar day if signs are equal, otherwise polar night

    sinϕ, cosϕ = sinlat[j], coslat[j]
    cos_zenith_j = h₀ * sinδ * sinϕ + cosδ * cosϕ * sin(h₀)
    cos_zenith_j /= π

    cos_zenith[ij] = cos_zenith_j
end

function variables(::AbstractZenith)
    return (
        DiagnosticVariable(name = :cos_zenith, dims = Grid2D(), desc = "Cosine of solar zenith angle", units = "1"),
    )
end
