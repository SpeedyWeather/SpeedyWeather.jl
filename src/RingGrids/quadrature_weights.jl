# QUADRATURE WEIGHTS (EXACT)
# gaussian_weights are exact for Gaussian latitudes when nlat > (2T+1)/2
# clenshaw_curtis_weights are exact for equi-angle latitudes when nlat > 2T+1
gaussian_weights(nlat_half::Integer) = FastGaussQuadrature.gausslegendre(2nlat_half)[2][1:nlat_half]
function clenshaw_curtis_weights(nlat_half::Integer)
    nlat = get_nlat(FullClenshawArray, nlat_half)
    θs = get_colat(FullClenshawGrid, nlat_half)
    return [4sin(θj)/(nlat+1)*sum([sin(p*θj)/p for p in 1:2:nlat]) for θj in θs[1:nlat_half]]
end

# QUADRATURE WEIGHTS (INEXACT), for HEALPix full grids
function equal_area_weights(Grid::Type{<:AbstractGridArray}, nlat_half::Integer)
    weights = zeros(nlat_half)
    for j in 1:nlat_half
        nlon = get_nlon_per_ring(Grid, nlat_half, j)
        weights[j] = 2nlon/get_npoints2D(Grid, nlat_half)
    end
    return weights
end

healpix_weights(nlat_half::Integer) = equal_area_weights(HEALPixArray, nlat_half)
octahealpix_weights(nlat_half::Integer) = equal_area_weights(OctaHEALPixArray, nlat_half)

# SOLID ANGLES ΔΩ = sinθ Δθ Δϕ
get_solid_angles(Grid::Type{<:AbstractGridArray}, nlat_half::Integer) = 
    get_quadrature_weights(Grid, nlat_half) .* (2π./get_nlons(Grid, nlat_half))
get_solid_angles(Grid::Type{<:Union{HEALPixArray, OctaHEALPixArray}}, nlat_half::Integer) =
    4π/get_npoints(Grid, nlat_half)*ones(nlat_half)
