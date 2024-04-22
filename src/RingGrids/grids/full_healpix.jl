"""A `FullHEALPixArray` is an array of full grids, subtyping `AbstractFullGridArray`, that use
HEALPix latitudes for each ring. This type primarily equists to interpolate data from
the (reduced) HEALPixGrid onto a full grid for output.

First dimension of the underlying `N`-dimensional `data` represents the horizontal dimension,
in ring order (0 to 360˚E, then north to south), other dimensions are used for the vertical
and/or time or other dimensions. The resolution parameter of the horizontal grid is
`nlat_half` (number of latitude rings on one hemisphere, Equator included) and the ring indices
are precomputed in `rings`. Fields are
$(TYPEDFIELDS)"""
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

"""A `FullHEALPixArray` but constrained to `N=1` dimensions (horizontal only) and data is a `Vector{T}`."""
FullHEALPixGrid

# SIZE
nlat_odd(::Type{<:FullHEALPixArray}) = true
get_npoints2D(::Type{<:FullHEALPixArray}, nlat_half::Integer) = 4nlat_half * (2nlat_half-1)
get_nlat_half(::Type{<:FullHEALPixArray}, npoints2D::Integer) = round(Int, 1/4 + sqrt(1/16 + npoints2D/8))
get_nlon(::Type{<:FullHEALPixArray}, nlat_half::Integer) = 4nlat_half

## COORDINATES
get_colat(::Type{<:FullHEALPixArray}, nlat_half::Integer) = get_colat(HEALPixGrid, nlat_half)
get_lon(::Type{<:FullHEALPixArray}, nlat_half::Integer) = get_lon(FullGaussianArray, nlat_half)

# QUADRATURE
get_quadrature_weights(::Type{<:FullHEALPixArray}, nlat_half::Integer) = equal_area_weights(FullHEALPixArray, nlat_half)