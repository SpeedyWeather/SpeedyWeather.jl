# QUADRATURE WEIGHTS (EXACT)
# gaussian_weights are exact for Gaussian latitudes when nlat > (2T+1)/2
# clenshaw_curtis_weights are exact for equi-angle latitudes when nlat > 2T+1
gaussian_weights(nlat_half::Integer) = FastGaussQuadrature.gausslegendre(2nlat_half)[2][1:nlat_half]
function clenshaw_curtis_weights(nlat_half::Integer)
    nlat = 2nlat_half - 1
    θs = get_colat(FullClenshawGrid,nlat_half)
    return [4sin(θj)/(nlat+1)*sum([sin(p*θj)/p for p in 1:2:nlat]) for θj in θs[1:nlat_half]]
end

# QUADRATURE WEIGHTS (INEXACT), for HEALPix full grids
healpix_weights(nlat_half::Integer) =
    [nlon_healpix(nlat_half,j) for j in 1:nlat_half].*(2/npoints_healpix(nlat_half))
healpix4_weights(nlat_half::Integer) =
    [nlon_healpix4(nlat_half,j) for j in 1:nlat_half].*(2/npoints_healpix4(nlat_half))

# SOLID ANGLES ΔΩ = sinθ Δθ Δϕ
get_solid_angles(Grid::Type{<:AbstractGrid},nlat_half::Integer) = 
    get_quadrature_weights(Grid,nlat_half) .* (2π./get_nlons(Grid,nlat_half))
get_solid_angles(Grid::Type{<:Union{HEALPixGrid,HEALPix4Grid}},nlat_half::Integer) =
    4π/get_npoints(Grid,nlat_half)*ones(nlat_half)
