"""$(TYPEDSIGNATURES)
The Gaussian weights for a Gaussian grid (full or octahedral) of size nlat_half.
Gaussian weights are of length `nlat`, i.e. a vector for every latitude ring, pole to pole.
`sum(gaussian_weights(nlat_half))` is always `2` as int_0^π sin(x) dx = 2 (colatitudes),
or equivalently int_-pi/2^pi/2 cos(x) dx (latitudes).

Integration (and therefore the spectral transform) is _exact_ (only rounding errors)
when using Gaussian grids provided that nlat >= 3(T + 1)/2, meaning that a grid resolution
of at least 96x48 (nlon x nlat) is sufficient for an exact transform with a T=31 spectral
truncation."""
gaussian_weights(nlat_half::Integer) = FastGaussQuadrature.gausslegendre(2nlat_half)[2]


"""$(TYPEDSIGNATURES)
The Clenshaw-Curtis weights for a Clenshaw grid (full or octahedral) of size nlat_half.
Clenshaw-Curtis weights are of length `nlat`, i.e. a vector for every latitude ring, pole to pole.
`sum(clenshaw_curtis_weights(nlat_half))` is always `2` as int_0^π sin(x) dx = 2 (colatitudes),
or equivalently int_-pi/2^pi/2 cos(x) dx (latitudes).

Integration (and therefore the spectral transform) is _exact_ (only rounding errors)
when using Clenshaw grids provided that nlat >= 2(T + 1), meaning that a grid resolution
of at least 128x64 (nlon x nlat) is sufficient for an exact transform with a T=31 spectral
truncation."""
function clenshaw_curtis_weights(nlat_half::Integer)
    nlat = get_nlat(FullClenshawArray, nlat_half)
    θs = get_colat(FullClenshawGrid, nlat_half)
    return [4sin(θj)/(nlat+1)*sum([sin(p*θj)/p for p in 1:2:nlat]) for θj in θs]
end

"""$(TYPEDSIGNATURES)
The equal-area weights used for the HEALPix grids (original or OctaHEALPix) of size nlat_half.
The weights are of length `nlat`, i.e. a vector for every latitude ring, pole to pole.
`sum(equal_area_weights(nlat_half))` is always `2` as int_0^π sin(x) dx = 2 (colatitudes),
or equivalently int_-pi/2^pi/2 cos(x) dx (latitudes). Integration (and therefore the
spectral transform) is not exact with these grids but errors reduce for higher resolution."""
function equal_area_weights(Grid::Type{<:AbstractGrid}, nlat_half::Integer)
    nlat = get_nlat(Grid, nlat_half)
    npoints = get_npoints(Grid, nlat_half)
    weights = zeros(nlat)
    for j in 1:nlat
        nlon = get_nlon_per_ring(Grid, nlat_half, j)
        weights[j] = 2nlon/npoints
    end
    return weights
end

# SOLID ANGLES ΔΩ = sinθ Δθ Δϕ
get_solid_angles(Grid::Type{<:AbstractGrid}, nlat_half::Integer) =
    get_quadrature_weights(Grid, nlat_half) .* (2π./get_nlons(Grid, nlat_half))
get_solid_angles(Grid::Type{<:Union{HEALPixGrid, OctaHEALPixGrid}}, nlat_half::Integer) =
    4π/get_npoints(Grid, nlat_half)*ones(get_nlat(Grid, nlat_half))