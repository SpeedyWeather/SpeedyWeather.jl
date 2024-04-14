struct FullClenshawArray{T, N, ArrayType <: AbstractArray{T, N}} <: AbstractFullGridArray{T, N, ArrayType}
    data::ArrayType                 # data array, ring by ring, north to south
    nlat_half::Int                  # number of latitudes on one hemisphere
    rings::Vector{UnitRange{Int}}   # TODO make same array type as data?

    FullClenshawArray(data::A, nlat_half, rings) where {A <: AbstractArray{T, N}} where {T, N} =
        check_inputs(data, nlat_half, rings, FullClenshawArray) ?
        new{T, N, A}(data, nlat_half, rings) :
        error_message(data, nlat_half, rings, FullClenshawArray, T, N, A)
end

# TYPES
const FullClenshawGrid{T} = FullClenshawArray{T, 1, Vector{T}}
nonparametric_type(::Type{<:FullClenshawArray}) = FullClenshawArray
horizontal_grid_type(::Type{<:FullClenshawArray}) = FullClenshawGrid

# SIZE
nlat_odd(::Type{<:FullClenshawArray}) = true
get_npoints2D(::Type{<:FullClenshawArray}, nlat_half::Integer) = 8 * nlat_half^2 - 4nlat_half
get_nlat_half(::Type{<:FullClenshawArray}, npoints2D::Integer) = round(Int, 1/4 + sqrt(1/16 + npoints2D/8))
get_nlon(::Type{<:FullClenshawArray}, nlat_half::Integer) = 4nlat_half

## COORDINATES
get_colat(::Type{<:FullClenshawArray}, nlat_half::Integer) = [j/(2nlat_half)*π for j in 1:2nlat_half-1]
get_lon(::Type{<:FullClenshawArray}, nlat_half::Integer) = get_lon(FullGaussianArray, nlat_half)

# QUADRATURE
get_quadrature_weights(::Type{<:FullClenshawArray}, nlat_half::Integer) = clenshaw_curtis_weights(nlat_half)


"""
    G = FullClenshawGrid{T}

A FullClenshawGrid is a regular latitude-longitude grid with an odd number of `nlat` equi-spaced
latitudes, the central latitude ring is on the Equator. The same `nlon` longitudes for every latitude ring.
The grid points are closer in zonal direction around the poles. The values of all grid points are stored
in a vector field `data` that unravels the data 0 to 360˚, then ring by ring, which are sorted north to south."""
FullClenshawGrid