const EARTH_RADIUS = 6371e3  # Earth mean radius in meters.

"""$(TYPEDSIGNATURES) 
Haversine formula calculating the great-circle distance between a pair of 
longitude-latitude points. Assumes the Earth is a perfect sphere with radius 
6371 km."""
function haversine((lon1, lat1), (lon2, lat2))
    phi_1 = deg2rad(lat1)
    phi_2 = deg2rad(lat2)

    delta_phi = deg2rad(lat2 - lat1)
    delta_lambda = deg2rad(lon2 - lon1)

    a = sin(delta_phi / 2)^2 + cos(phi_1) * cos(phi_2) * sin(delta_lambda / 2)^2
    c = 2 * atan(sqrt(a), sqrt(1 - a))
    return EARTH_RADIUS * c # in meters
end