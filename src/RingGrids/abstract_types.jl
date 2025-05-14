# """Abstract supertype for all arrays of ring grids, representing `N`-dimensional
# data on the sphere in two dimensions (but unravelled into a vector in the first dimension,
# the actual "ring grid") plus additional `N-1` dimensions for the vertical and/or time etc.
# Parameter `T` is the `eltype` of the underlying data, held as in the array type `ArrayType`
# (Julia's `Array` for CPU or others for GPU).

# Ring grids have several consecuitive grid points share the same latitude (= a ring),
# grid points on a given ring are equidistant. Grid points are ordered 0 to 360˚E,
# starting around the north pole, ring by ring to the south pole. """
# # abstract type AbstractGridArray{T, N, ArrayType <: AbstractArray{T, N}} <: AbstractArray{T, N} end

# """Abstract supertype for all ring grids, representing 2-dimensional data on the
# sphere unravelled into a Julia `Vector`. Subtype of `AbstractGridArray` with
# `N=1` and `ArrayType=Vector{T}` of `eltype T`."""
# # const AbstractGrid{T} = AbstractGridArray{T, 1, Vector{T}}

# """An `AbstractFullGrid` is a horizontal grid with a constant number of longitude
# points across latitude rings. Different latitudes can be used, Gaussian latitudes,
# equi-angle latitudes (also called Clenshaw from Clenshaw-Curtis quadrature), or others."""
# """Subtype of `AbstractGridArray` for all N-dimensional arrays of ring grids that have the
# same number of longitude points on every ring. As such these (horizontal) grids are representable
# as a matrix, with denser grid points towards the poles."""

# """Subtype of `AbstractGridArray` for arrays of rings grids that have a reduced number
# of longitude points towards the poles, i.e. they are not "full", see `AbstractFullGridArray`.
# Data on these grids cannot be represented as matrix and has to be unravelled into a vector,
# ordered 0 to 360˚E then north to south, ring by ring. Examples for reduced grids are 
# the octahedral Gaussian or Clenshaw grids, or the HEALPix grid."""

abstract type AbstractGrid end
abstract type AbstractFullGrid <: AbstractGrid end
abstract type AbstractReducedGrid <: AbstractGrid end

abstract type AbstractField{T, N, ArrayType, Grid} <: AbstractArray{T, N} end