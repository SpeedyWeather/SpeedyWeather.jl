"""A `FullClenshawArray` is an array of full grid, subtyping `AbstractFullGridArray`, that use
equidistant latitudes for each ring (a regular lon-lat grid). These require the
Clenshaw-Curtis quadrature in the spectral transform, hence the name. One ring is on the equator,
total number of rings is odd, no rings on the north or south pole. First dimension of the
underlying `N`-dimensional `data` represents the horizontal dimension, in ring order
(0 to 360˚E, then north to south), other dimensions are used for the vertical and/or
time or other dimensions. The resolution parameter of the horizontal grid is `nlat_half`
(number of latitude rings on one hemisphere, Equator included) and the ring indices are
precomputed in `rings`. Fields are
$(TYPEDFIELDS)"""
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

"""A `FullClenshawArray` but constrained to `N=1` dimensions (horizontal only) and data is a `Vector{T}`."""
FullClenshawGrid

# SIZE
nlat_odd(::Type{<:FullClenshawArray}) = true
get_npoints2D(::Type{<:FullClenshawArray}, nlat_half::Integer) = 8 * nlat_half^2 - 4nlat_half
get_nlat_half(::Type{<:FullClenshawArray}, npoints2D::Integer) = round(Int, 1/4 + sqrt(1/16 + npoints2D/8))
get_nlon(::Type{<:FullClenshawArray}, nlat_half::Integer) = 4nlat_half

## COORDINATES
get_latd(::Type{<:FullClenshawArray}, nlat_half::Integer) = [90 - 90j/nlat_half for j in 1:2nlat_half-1]
get_lond(::Type{<:FullClenshawArray}, nlat_half::Integer) = get_lond(FullGaussianArray, nlat_half)

# QUADRATURE
get_quadrature_weights(::Type{<:FullClenshawArray}, nlat_half::Integer) = clenshaw_curtis_weights(nlat_half)