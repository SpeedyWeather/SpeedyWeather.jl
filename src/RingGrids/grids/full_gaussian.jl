"""A `FullGaussianArray` is an array of full grids, subtyping `AbstractFullGridArray`, that use
Gaussian latitudes for each ring. First dimension of the underlying `N`-dimensional `data`
represents the horizontal dimension, in ring order (0 to 360˚E, then north to south),
other dimensions are used for the vertical and/or time or other dimensions.
The resolution parameter of the horizontal grid is `nlat_half` (number of latitude rings
on one hemisphere, Equator included) and the ring indices are precomputed in `rings`.
Fields are
$(TYPEDFIELDS)"""
struct FullGaussianArray{T, N, ArrayType <: AbstractArray{T, N}} <: AbstractFullGridArray{T, N, ArrayType}
    "Data array, west to east, ring by ring, north to south."
    data::ArrayType
    
    "Number of latitudes on one hemisphere"
    nlat_half::Int
    
    "Precomputed ring indices, ranging from first to last grid point on every ring."
    rings::Vector{UnitRange{Int}}   # TODO make same array type as data?

    FullGaussianArray(data::A, nlat_half, rings) where {A <: AbstractArray{T, N}} where {T, N} =
        check_inputs(data, nlat_half, rings, FullGaussianArray) ?
        new{T, N, A}(data, nlat_half, rings) :
        error_message(data, nlat_half, rings, FullGaussianArray, T, N, A)
end

# TYPES
const FullGaussianGrid{T} = FullGaussianArray{T, 1, Vector{T}}
nonparametric_type(::Type{<:FullGaussianArray}) = FullGaussianArray     # identity for full grids
horizontal_grid_type(::Type{<:FullGaussianArray}) = FullGaussianGrid

"""A `FullGaussianArray` but constrained to `N=1` dimensions (horizontal only) and data is a `Vector{T}`."""
FullGaussianGrid

# SIZE
nlat_odd(::Type{<:FullGaussianArray}) = false   # Gaussian latitudes always even
get_npoints2D(::Type{<:FullGaussianArray}, nlat_half::Integer) = 8 * nlat_half^2
get_nlat_half(::Type{<:FullGaussianArray}, npoints2D::Integer) = round(Int, sqrt(npoints2D/8))
get_nlon(::Type{<:FullGaussianArray}, nlat_half::Integer) = 4nlat_half

## COORDINATES
function get_latd(::Type{<:FullGaussianArray}, nlat_half::Integer)
    return acosd.(FastGaussQuadrature.gausslegendre(2nlat_half)[1]) .- 90
end

function get_lond(::Type{<:FullGaussianArray}, nlat_half::Integer)
    nlat_half == 0 && return Float64[]      # necessary to avoid error from /0 below
    nlon = get_nlon(FullGaussianArray, nlat_half)
    return collect(range(0, 360 - 180/nlon, step=360/nlon))
end

# QUADRATURE
get_quadrature_weights(::Type{<:FullGaussianArray}, nlat_half::Integer) = gaussian_weights(nlat_half)