const DEFAULT_RADIUS = 6371000  # Earth mean radius in meters, defined as integer for more type stability 

"""
    abstract type AbstractSphericalDistance end

Super type of formulas to calculate the spherical distance or great-circle distance.
To define a `NewFormula`, define `struct NewFormula <: AbstractSphericalDistance end`
and the actual calculation as a functor

    function NewFormula(lonlat1::Tuple, lonlat2::Tuple; radius=DEFAULT_RADIUS, kwargs...)

assuming inputs in degrees and returning the distance in meters (or radians for radius=1).
Then use the general interface `spherical_distance(NewFormula, args...; kwargs...)`"""
abstract type AbstractSphericalDistance end
struct Haversine <: AbstractSphericalDistance end

"""$(TYPEDSIGNATURES) 
Haversine formula calculating the great-circle or spherical distance (in meters) on the sphere
between two tuples of  longitude-latitude points in degrees ˚E, ˚N. Use keyword argument `radius`
to change the radius of the sphere (default 6371e3 meters, Earth's radius), use `radius=1`
to return the central angle in radians."""
function Haversine(             # functor for Haversine struct
    lonlat1::Tuple,             # point 1 in spherical coordinates
    lonlat2::Tuple;             # point 2
    radius = DEFAULT_RADIUS,    # radius of the sphere [m]
)
    lon1, lat1 = lonlat1
    lon2, lat2 = lonlat2    

    φ1 = deg2rad(lat1)
    φ2 = deg2rad(lat2)

    Δφ = deg2rad(lat2 - lat1)
    Δλ = deg2rad(lon2 - lon1)

    # Haversine formula, see https://en.wikipedia.org/wiki/Haversine_formula
    a = sin(Δφ / 2)^2 + cos(φ1) * cos(φ2) * sin(Δλ / 2)^2
    c = 2 * atan(sqrt(a), sqrt(1 - a))
    return c*radius             # in meters
end

# allow for any <:AbstractSphericalDistance also non-tupled arguments lon1, lat1, lon2, lat2
(F::Type{<:AbstractSphericalDistance})(lon1, lat1, lon2, lat2; kwargs...) = F((lon1, lat1), (lon2, lat2); kwargs...) 

"""$(TYPEDSIGNATURES)
Spherical distance, or great-circle distance, between two points `lonlat1` and `lonlat2`
using the `Formula` (default `Haversine`)."""
spherical_distance(Formula::Type{<:AbstractSphericalDistance}, args...; kwargs...) =
    Formula(args...; kwargs...)

# define Haversine as default
spherical_distance(args...; kwargs...) = spherical_distance(Haversine, args...; kwargs...)


