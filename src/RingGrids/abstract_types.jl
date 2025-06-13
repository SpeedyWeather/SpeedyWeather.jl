"""Abstract supertype for all grids in RingGrids.jl, representing a discretization (particularly a tessellation or tiling)
of the sphere. A "grid" does not contain any data only its resolution is defined by `nlat_half`
(number of latitude rings on one hemisphere including the equator). A grid does not contain the coordinates of the grid points,
vertices, latitudes, longtiudes etc but they can be recomputed from the grid object (or its type and resolution) at any time.

A grid is regarded as 2D but a field (data on a grid) can have N additional dimensions, e.g. for vertical levels or time.
In that sense, a grid is does not have a number format / eltype. This is a property of a `Field` which can be different for
every field even when using the same grid. Points on the grid (cell centres) are unravalled into a vector ordered 0 to 360˚E,
starting at the north pole, then going ring by ring to the south pole. This way both full (those representable as a matrix)
and reduced grids (not representable as a matrix, with fewer longitude points towards the poles) can be represented.

A grid has a parameter `Architecture` (the type of the field `architecture`) that can be used to store information
about the architecture the grid is on, e.g. CPU or GPU.

Furthermore all grids have a `rings` to allow ring-by-ring indexing. Other precomputed indices can be added for specific grids."""
abstract type AbstractGrid{Architecture} end

"""Abstract supertype for all full grids, representing a horizontal grid with a constant number of longitude
points across latitude rings. Different latitudes can be used, Gaussian latitudes,
equi-angle latitudes (also called Clenshaw from Clenshaw-Curtis quadrature), or others."""
abstract type AbstractFullGrid{Architecture} <: AbstractGrid{Architecture} end

"""Abstract supertype for all reduced grids, representing arrays of ring grids that have a reduced number
of longitude points towards the poles, i.e. they are not "full", see `AbstractFullGrid`.
Data on these grids (a `Field`) cannot be represented as matrix and has to be unravelled into a vector,
ordered 0 to 360˚E then north to south, ring by ring. Examples for reduced grids are
the octahedral Gaussian or Clenshaw grids, or the HEALPix grid."""
abstract type AbstractReducedGrid{Architecture} <: AbstractGrid{Architecture} end

"""Abstract supertype for all fields, i.e. data on a grid. A field is an `AbstractArray` with a number format
`T`, a number of dimensions `N`, an `ArrayType` (e.g. `Vector`, `Matrix`, `Tensor`) and a grid type `Grid`.
Fields can be full or reduced, 2D, 3D, or 4D, depending on the grid and the number of dimensions.
Fields are used to store data on the grid, e.g. temperature, pressure, or any other variable.
Every `Field` has `data` the grid-free array of the data and `grid` which contains information about the grid
the field is defined on."""
abstract type AbstractField{T, N, ArrayType, Grid} <: AbstractArray{T, N} end

"""Abstract supertype for all fields on full grids, i.e. grids with a constant number of longitude points across latitude rings."""
const AbstractFullField = AbstractField{T, N, ArrayType, Grid} where {T, N, ArrayType, Grid<:AbstractFullGrid}

"""Abstract supertype for all fields on reduced grids, i.e. grids with a reduced number of longitude points towards the poles."""
const AbstractReducedField = AbstractField{T, N, ArrayType, Grid} where {T, N, ArrayType, Grid<:AbstractReducedGrid}

"""Abstract supertype for all 2D fields, i.e. fields with the horizontal dimensions only. Note that this is a `<:AbstractVector`
as the horizontal dimensions are unravelled into a vector for all grids to be conistent with the reduced grids that cannot be
represented as a matrix."""
const AbstractField2D = AbstractField{T, 1} where T

"""Abstract supertype for all 3D fields, i.e. fields with horizontal and one vertical (or time etc) dimension."""
const AbstractField3D = AbstractField{T, 2} where T

"""Abstract supertype for all 4D fields, i.e. fields with horizontal and (in most cases) a vertical and a time dimensions,
though these additional dimensions are arbitrary."""
const AbstractField4D = AbstractField{T, 3} where T

const AbstractFullField2D = AbstractFullField{T, 1} where T
const AbstractFullField3D = AbstractFullField{T, 2} where T
const AbstractFullField4D = AbstractFullField{T, 3} where T
const AbstractReducedField2D = AbstractReducedField{T, 1} where T
const AbstractReducedField3D = AbstractReducedField{T, 2} where T
const AbstractReducedField4D = AbstractReducedField{T, 3} where T