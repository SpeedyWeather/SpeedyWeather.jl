abstract type AbstractPlanet <: AbstractModelComponent end

# EARTH DEFAULTS
const DEFAULT_RADIUS = 6.371e6          # radius of planet [m]
const DEFAULT_ROTATION = 7.29e-5        # angular frequency of Earth's rotation [1/s]
const DEFAULT_GRAVITY = 9.81            # gravitational acceleration [m/s²]
const EARTH_DAY = Hour(24)              # length of day
const EARTH_YEAR = Day(365) + Hour(6)   # length of year
const EARTH_EQUINOX = DateTime(2000, 3, 20)
const AXIAL_TILT = 23.4
const EARTH_ECCENTRICITY = 0.0167
const EARTH_LONGITUDE_OF_PERIHELION = 282.9     # [˚] from the vernal equinox, perihelion early January

export AbstractOrbit

"""Abstract supertype for the orbit of a planet around its star. The orbit defines the
variables needed for its computation (`variables`), the solar zenith executes that
computation once per time step before the zenith angle (`orbit_position!`), and the
shortwave radiation reads the Sun-planet distance factor (`solar_distance_factor`)."""
abstract type AbstractOrbit <: AbstractModelComponent end

export CircularOrbit

"""Circular orbit: the solar declination follows the solar zenith's `solar_declination`
component and the Sun-planet distance is constant, so the distance factor (a/r)² is 1
and no variable is needed for it."""
struct CircularOrbit <: AbstractOrbit end

"""$(TYPEDSIGNATURES)
Mean anomaly [radians] at the vernal equinox of an orbit with eccentricity `e` and longitude
of perihelion `longitude_of_perihelion` [˚], where the true anomaly is minus the longitude of
perihelion. Throws a `DomainError` unless 0 ≤ e < 1 (e = 1 is a parabola, not an orbit)."""
function mean_anomaly_at_equinox(e, longitude_of_perihelion)
    0 <= e < 1 || throw(DomainError(e, "eccentricity has to be 0 ≤ e < 1"))
    ϖ = deg2rad(longitude_of_perihelion)
    E = 2 * atan(sqrt(1 - e) * sin(-ϖ / 2), sqrt(1 + e) * cos(-ϖ / 2))    # eccentric anomaly
    return E - e * sin(E)
end

export KeplerOrbit

"""$(TYPEDSIGNATURES)
Elliptical (Kepler) orbit. `KeplerOrbit` can represent any elliptical orbit but its defaults
are chosen for present-day Earth: note that `KeplerOrbit()` is almost identical to
`CircularOrbit()` in the solar declination but not in insolation, as the distance factor
(a/r)² varies by about ±3.4% over the year, roughly ±46 W/m² at the top of the atmosphere.
The calendar is fixed by the planet's vernal `equinox`, as in the PMIP protocol, so that
precession moves perihelion through the year. Solar declination and the distance factor are
computed from the orbit, the solar zenith's `solar_declination` component is then not used.
The equation of time (`SolarTimeCorrection`) remains that of present-day Earth and is not
derived from the orbit, so for strongly eccentric orbits the daily cycle is slightly off.
Fields are
$(TYPEDFIELDS)"""
@parameterized @kwdef struct KeplerOrbit{NF} <: AbstractOrbit
    "[OPTION] Eccentricity of the orbit [1], 0 ≤ e < 1"
    @param eccentricity::NF = EARTH_ECCENTRICITY (bounds = 0 .. 1,)

    "[OPTION] Longitude of perihelion [˚]: true solar longitude at perihelion measured from the vernal equinox"
    @param longitude_of_perihelion::NF = EARTH_LONGITUDE_OF_PERIHELION

    "[DERIVED] Mean anomaly at the vernal equinox [radians]"
    mean_anomaly_equinox::NF = mean_anomaly_at_equinox(eccentricity, longitude_of_perihelion)
end

KeplerOrbit(SG::SpectralGrid; kwargs...) = KeplerOrbit{SG.NF}(; kwargs...)

# the orbit defines the variables needed for its computation
variables(::AbstractOrbit) = ()
variables(::KeplerOrbit) = (
    ParameterizationVariable(:solar_distance_factor, ScalarDim(1), desc = "Sun-planet distance factor (a/r)²", units = "1"),
)

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

# the planet's variables are those its orbit requires
variables(planet::Earth) = variables(planet.orbit)
