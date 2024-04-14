struct FullHEALPixArray{T, N, ArrayType <: AbstractArray{T, N}} <: AbstractFullGridArray{T, N, ArrayType}
    data::ArrayType                 # data array, ring by ring, north to south
    nlat_half::Int                  # number of latitudes on one hemisphere
    rings::Vector{UnitRange{Int}}   # TODO make same array type as data?

    FullHEALPixArray(data::A, nlat_half, rings) where {A <: AbstractArray{T, N}} where {T, N} =
        check_inputs(data, nlat_half, rings, FullHEALPixArray) ?
        new{T, N, A}(data, nlat_half, rings) :
        error_message(data, nlat_half, rings, FullHEALPixArray, T, N, A)
end

# TYPES
const FullHEALPixGrid{T} = FullHEALPixArray{T, 1, Vector{T}}
nonparametric_type(::Type{<:FullHEALPixArray}) = FullHEALPixArray
horizontal_grid_type(::Type{<:FullHEALPixArray}) = FullHEALPixGrid

# SIZE
nlat_odd(::Type{<:FullHEALPixArray}) = true
get_npoints2D(::Type{<:FullHEALPixArray}, nlat_half::Integer) = 4nlat_half * (2nlat_half-1)
get_nlat_half(::Type{<:FullHEALPixArray}, npoints2D::Integer) = round(Int, 1/4 + sqrt(1/16 + npoints2D/8))
get_nlon(::Type{<:FullHEALPixArray}, nlat_half::Integer) = 4nlat_half

## COORDINATES
get_colat(::Type{<:FullHEALPixArray}, nlat_half::Integer) = get_colat(HEALPixGrid, nlat_half)
get_lon(::Type{<:FullHEALPixArray}, nlat_half::Integer) = get_lon(FullGaussianArray, nlat_half)

# QUADRATURE
get_quadrature_weights(::Type{<:FullHEALPixArray}, nlat_half::Integer) = healpix_weights(nlat_half)


"""
    G = FullHEALPixGrid{T}

A full HEALPix grid is a regular latitude-longitude grid that uses `nlat` latitudes from the HEALPix grid,
and the same `nlon` longitudes for every latitude ring. The grid points are closer in zonal direction
around the poles. The values of all grid points are stored in a vector field `v` that unravels
the data 0 to 360˚, then ring by ring, which are sorted north to south."""
FullHEALPixGrid