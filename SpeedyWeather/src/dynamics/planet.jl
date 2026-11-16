abstract type AbstractPlanet <: AbstractModelComponent end

# EARTH DEFAULTS
const DEFAULT_RADIUS = 6.371e6          # radius of planet [m]
const DEFAULT_ROTATION = 7.29e-5        # angular frequency of Earth's rotation [1/s]
const DEFAULT_GRAVITY = 9.81            # gravitational acceleration [m/s²]
const EARTH_DAY = Hour(24)              # length of day
const EARTH_YEAR = Day(365) + Hour(6)   # length of year
const EARTH_EQUINOX = DateTime(2000, 3, 20)
const AXIAL_TILT = 23.4

export AbstractOrbit, CircularOrbit, KeplerOrbit

"""Abstract supertype for the orbit of a planet around its star, used by the solar zenith
calculation to obtain the solar declination and the Sun–planet distance factor."""
abstract type AbstractOrbit end

"""Circular orbit: the solar declination follows the solar zenith's `solar_declination`
component and the Sun–planet distance is constant, so the distance factor (a/r)² is 1."""
struct CircularOrbit <: AbstractOrbit end

"""$(TYPEDSIGNATURES)
Elliptical (Kepler) orbit. The calendar is fixed by the planet's vernal `equinox`, as in the
PMIP protocol, so that precession moves perihelion through the year. Solar declination and the
distance factor (a/r)² are computed from the orbit, the solar zenith's `solar_declination`
component is then not used. Fields are
$(TYPEDFIELDS)"""
@kwdef struct KeplerOrbit{NF} <: AbstractOrbit
    "eccentricity of the orbit [1], 0 ≤ e < 1"
    eccentricity::NF = 0.0167

    "longitude of perihelion [˚]: true solar longitude at perihelion measured from the vernal equinox, ≈ 283˚ for present-day Earth (perihelion early January)"
    longitude_of_perihelion::NF = 282.9
end

export Earth

"""$(TYPEDSIGNATURES)
Create a struct `Earth<:AbstractPlanet`, with the following physical/orbital
characteristics. Note that `radius` is not part of it as this should be chosen
in `SpectralGrid`. Keyword arguments are
$(TYPEDFIELDS)"""
@parameterized @kwdef struct Earth{NF, S, DT, B, O} <: AbstractPlanet

    "Earth's radius [m]"
    @param radius::NF = DEFAULT_RADIUS

    "angular frequency of Earth's rotation [rad/s]"
    @param rotation::NF = DEFAULT_ROTATION

    "gravitational acceleration [m/s^2]"
    @param gravity::NF = DEFAULT_GRAVITY (bounds = Nonnegative,)

    "switch on/off daily cycle"
    daily_cycle::B = true

    "Seconds in a daily rotation"
    length_of_day::S = Second(EARTH_DAY)

    "switch on/off seasonal cycle"
    seasonal_cycle::B = true

    "Seconds in an orbit around the sun"
    length_of_year::S = Second(EARTH_YEAR)

    "time of spring equinox (year irrelevant)"
    equinox::DT = EARTH_EQUINOX

    "angle [˚] rotation axis tilt wrt to orbit"
    @param axial_tilt::NF = AXIAL_TILT (bounds = -90 .. 90,)

    "Total solar irradiance at the distance of 1 AU [W/m²]"
    @param solar_constant::NF = 1365.0 (bounds = Nonnegative,)

    "orbit around the star, `CircularOrbit()` or `KeplerOrbit(; eccentricity, longitude_of_perihelion)`"
    orbit::O = CircularOrbit()
end

Adapt.@adapt_structure Earth

Earth(SG::SpectralGrid; kwargs...) = Earth(SG.NF; kwargs...)
Earth(::Type{NF}; orbit = CircularOrbit(), kwargs...) where {NF} =
    Earth{NF, Dates.Second, Dates.DateTime, Bool, typeof(orbit)}(; orbit, kwargs...)
