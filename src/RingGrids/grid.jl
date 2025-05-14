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

isreduced(::Type{<:AbstractGrid}) = true
isreduced(::Type{<:AbstractFullGrid}) = false
isreduced(grid) = isreduced(typeof(grid))
isfull(x) = ~isreduced(x)

function Base.show(io::IO, grid::AbstractGrid)
    Grid_ = nonparametric_type(grid)
    nlat_half = get_nlat_half(grid)
    nlat = get_nlat(grid)
    npoints = get_npoints(grid)
    full_or_reduced = isfull(grid) ? "full" : "reduced"
    average_resolution = Printf.@sprintf("%.3f˚", sqrt(4π/npoints)*360/2π)

    println(io, "$nlat-ring $Grid_")
    println(io, "├ nlat_half=$nlat_half ($npoints points, $full_or_reduced)")
    print(io,   "└ Resolution: $average_resolution (average)")
end

## TYPES
"""$(TYPEDSIGNATURES) For any instance of `AbstractGrid` type its n-dimensional type
(*Grid{T, N, ...} returns *Array) but without any parameters `{T, N, ArrayType}`"""
nonparametric_type(grid::AbstractGrid) = nonparametric_type(typeof(grid))

# also needed for other array types, defined in extensions
nonparametric_type(::Type{<:Array}) = Array

# """$(TYPEDSIGNATURES) Full grid array type for `grid`. Always returns the N-dimensional `*Array`
# not the two-dimensional (`N=1`) `*Grid`. For reduced grids the corresponding full grid that
# share the same latitudes."""
# full_array_type(grid::AbstractGridArray) = full_array_type(typeof(grid))

# """$(TYPEDSIGNATURES) Full (horizontal) grid type for `grid`. Always returns the two-dimensional
# (`N=1`) `*Grid` type. For reduced grids the corresponding full grid that share the same latitudes."""
# full_grid_type(grid::AbstractGridArray) = horizontal_grid_type(full_array_type(grid))
# full_grid_type(Grid::Type{<:AbstractGridArray}) = horizontal_grid_type(full_array_type(Grid))

# """$(TYPEDSIGNATURES) The two-dimensional (`N=1`) `*Grid` for `grid`, which can be an N-dimensional
# `*GridArray`."""
# horizontal_grid_type(grid::AbstractGridArray) = horizontal_grid_type(typeof(grid))

# ## SIZE
# Base.length(G::AbstractGridArray) = length(G.data)  # total number of grid points
# Base.size(G::AbstractGridArray) = size(G.data)      # size of underlying data (horizontal is unravelled)

# """$(TYPEDSIGNATURES) Size of underlying data array plus precomputed ring indices."""
# Base.sizeof(G::AbstractGridArray) = sizeof(G.data) + sizeof(G.rings)

"""$(TYPEDSIGNATURES) Resolution paraemeters `nlat_half` of a `grid`.
Number of latitude rings on one hemisphere, Equator included."""
get_nlat_half(grid::AbstractGrid) = grid.nlat_half

"""$(TYPEDSIGNATURES) True for a `grid` that has an odd number of
latitude rings `nlat` (both hemispheres)."""
nlat_odd(grid::AbstractGrid) = nlat_odd(typeof(grid))

# get total number of latitude rings, *(nlat_half > 0) to return 0 for nlat_half = 0
get_nlat(Grid::Type{<:AbstractGrid}, nlat_half::Integer) = 2nlat_half - nlat_odd(Grid)*(nlat_half > 0)

"""$(TYPEDSIGNATURES) Get number of latitude rings, pole to pole."""
get_nlat(grid::AbstractGrid) = get_nlat(typeof(grid), grid.nlat_half)

"""$(TYPEDSIGNATURES) Total number of grid points in all dimensions of `grid`.
Equivalent to length of the underlying data array."""
get_npoints(grid::Grid) where {Grid<:AbstractGrid} = get_npoints(Grid, grid.nlat_half)
get_npoints(G::Type{<:AbstractGrid}, nlat_half::Integer, k::Integer...) = prod(k) * get_npoints(G, nlat_half)

"""$(TYPEDSIGNATURES) Size of the matrix of the horizontal grid if representable as such (not all grids)."""
matrix_size(grid::Grid) where {Grid<:AbstractGrid} = matrix_size(Grid, get_nlat_half(grid))

# CONSTRUCTORS
"""$(TYPEDSIGNATURES) Create a new `grid` of type `Grid` with resolution parameter `nlat_half`.
`architecture` is the device type (CPU/GPU). Precomputes the ring indices `rings`."""
function (::Type{Grid})(nlat_half::Integer, architecture=DEFAULT_ARCHITECTURE) where {Grid<:AbstractGrid}
    Grid_ = nonparametric_type(Grid)    # strip away parameters of type, obtain from arguments
    rings = eachring(Grid, nlat_half)   # precompute indices to access the variable-length rings
    return Grid_(nlat_half, architecture, rings)
end

## COORDINATES

"""$(TYPEDSIGNATURES) Longitudes (degrees, 0-360˚E), latitudes (degrees, 90˚N to -90˚N) for
every (horizontal) grid point in `grid` in ring order (0-360˚E then north to south)."""
get_londlatds(grid::Grid) where {Grid<:AbstractGrid} = get_londlatds(Grid, grid.nlat_half)

"""$(TYPEDSIGNATURES) Longitudes (radians, 0-2π), latitudes (degrees, π/2 to -π/2) for
every (horizontal) grid point in `grid` in ring order (0-360˚E then north to south)."""
get_lonlats(grid::Grid) where {Grid<:AbstractGrid} = get_lonlats(Grid, grid.nlat_half)

"""$(TYPEDSIGNATURES) Longitudes (radians, 0-2π), colatitudes (degrees, 0 to π) for
every (horizontal) grid point in `grid` in ring order (0-360˚E then north to south)."""
get_loncolats(grid::Grid) where {Grid<:AbstractGrid} = get_loncolats(Grid, grid.nlat_half)

# radians
function get_lonlats(Grid::Type{<:AbstractGrid}, nlat_half::Integer)
    londs, latds = get_londlatds(Grid, nlat_half)  # longitudes, latitudes in degrees
    return londs * π/180, latds * π/180             # to radians
end

# radians and colatitudes
function get_loncolats(Grid::Type{<:AbstractGrid}, nlat_half::Integer)
    lons, lats = get_lonlats(Grid, nlat_half)
    colats = π/2 .- lats
    return lons, colats
end

"""$(TYPEDSIGNATURES) Latitude (radians) for each ring in `grid`, north to south."""
get_lat(grid::Grid) where {Grid<:AbstractGrid} = get_lat(Grid, grid.nlat_half)

"""$(TYPEDSIGNATURES) Latitude (degrees) for each ring in `grid`, north to south."""
get_latd(grid::Grid) where {Grid<:AbstractGrid} = get_latd(Grid, grid.nlat_half)

"""$(TYPEDSIGNATURES) Colatitudes (radians) for each ring in `grid`, north to south."""
get_colat(grid::Grid) where {Grid<:AbstractGrid} = get_colat(Grid, grid.nlat_half)

"""$(TYPEDSIGNATURES) Longitude (degrees). Full grids only."""
get_lond(grid::Grid) where {Grid<:AbstractGrid} = get_lond(Grid, grid.nlat_half)

"""$(TYPEDSIGNATURES) Longitude (radians). Full grids only."""
get_lon(grid::Grid) where {Grid<:AbstractGrid} = get_lon(Grid, grid.nlat_half)

get_nlon_max(grid::Grid) where {Grid<:AbstractGrid} = get_nlon_max(Grid, grid.nlat_half)

get_lat(Grid::Type{<:AbstractGrid}, nlat_half::Integer) = get_latd(Grid, nlat_half) * (π/180)
get_colat(Grid::Type{<:AbstractGrid}, nlat_half::Integer) = π/2 .- get_lat(Grid, nlat_half)

"""$(TYPEDSIGNATURES)
Returns a vector `nlons` for the number of longitude points per latitude ring, north to south.
Provide grid `Grid` and its resolution parameter `nlat_half`. For keyword argument
`both_hemispheres=false` only the northern hemisphere (incl Equator) is returned."""
function get_nlons(Grid::Type{<:AbstractGrid}, nlat_half::Integer; both_hemispheres::Bool=true)
    n = both_hemispheres ? get_nlat(Grid, nlat_half) : nlat_half
    return [get_nlon_per_ring(Grid, nlat_half, j) for j in 1:n]
end

"""$(TYPEDSIGNATURES)
Number of longitude points per latitude ring `j`."""
function get_nlon_per_ring(grid::AbstractGrid, j::Integer)
    return get_nlon_per_ring(typeof(grid), grid.nlat_half, j)
end

## ITERATOR

"""$(TYPEDSIGNATURES)
Vector{UnitRange} `rings` to loop over every ring of `grid`
and then each grid point per ring. To be used like

    rings = eachring(grid)
    for ring in rings
        for ij in ring
            field[ij]

Accesses precomputed `grid.rings`."""
@inline eachring(grid::AbstractGrid) = grid.rings

"""$(TYPEDSIGNATURES)
Computes the ring indices `i0:i1` for start and end of every longitudinal point
on a given ring `j` of `Grid` at resolution `nlat_half`. Used to loop
over rings of a grid. These indices are also precomputed in every `grid.rings`."""
function eachring(Grid::Type{<:AbstractGrid}, nlat_half::Integer)
    rings = Vector{UnitRange{Int}}(undef, get_nlat(Grid, nlat_half))    # allocate
    each_index_in_ring!(rings, Grid, nlat_half)                         # calculate iteratively
    return rings
end

"""$(TYPEDSIGNATURES) Same as `eachring(grid)` but performs a bounds check to assess
that all `grids` according to `grids_match` (non-parametric grid type, nlat_half and length)."""
function eachring(grid1::AbstractGrid, grids::AbstractGrid...)
    grids_match(grid1, grids...) || throw(DimensionMismatch(grid1, grids...))
    return eachring(grid1)
end

function Base.DimensionMismatch(grid1::AbstractGrid, grids::AbstractGrid...)
    s = "Grids do not match; $(size(grid1)) $(nonparametric_type(grid1))"
    for grid in grids
        s *= ", $(size(grid))-$(nonparametric_type(grid))"
    end
    return DimensionMismatch(s)
end

# equality and comparison, somehow needed as not covered by broadcasting
Base.:(==)(G1::AbstractGrid, G2::AbstractGrid) = grids_match(G1, G2)

"""$(TYPEDSIGNATURES) True if both `A` and `B` are of the same nonparametric grid type
(e.g. OctahedralGaussianArray, regardless type parameter `T` or underyling array type `ArrayType`)
and of same resolution (`nlat_half`) and total grid points (`length`). Sizes of `(4,)` and `(4,1)`
would match for example, but `(8,1)` and `(4,2)` would not (`nlat_half` not identical)."""
function grids_match(A::AbstractGrid, B::AbstractGrid)
    horizontal_match = get_nlat_half(A) == get_nlat_half(B)
    type_match = grids_match(typeof(A), typeof(B))
    return horizontal_match && type_match
end

"""$(TYPEDSIGNATURES) True if both `A` and `B` are of the same nonparametric grid type."""
function grids_match(A::Type{<:AbstractGrid}, B::Type{<:AbstractGrid})
    return nonparametric_type(A) == nonparametric_type(B)
end

"""$(TYPEDSIGNATURES) True if all grids `A, B, C, ...` provided as arguments
match according to `grids_match` wrt to `A` (and therefore all)."""
function grids_match(A::AbstractGrid, Bs::AbstractGrid...)
    match = true    # single grid A always matches itself
    for B in Bs     # check for all matching respectively with A
        match &= grids_match(A, B)
    end
    return match
end

"""$(TYPEDSIGNATURES) UnitRange to access data on grid `grid` on ring `j`."""
each_index_in_ring(grid::AbstractGrid, j::Integer) = grid.rings[j]

""" $(TYPEDSIGNATURES) UnitRange to access each horizontal grid point on grid `grid`.
For a `NxM` (`N` horizontal grid points, `M` vertical layers) `OneTo(N)` is returned."""
eachgridpoint(grid::AbstractGrid) = Base.OneTo(get_npoints(grid))

""" $(TYPEDSIGNATURES) Like `eachgridpoint(::AbstractGridArray)` but checks for
equal size between input arguments first."""
function eachgridpoint(grid1::AbstractGrid, grids::AbstractGrid...)
    grids_match(grid1, grids...) || throw(DimensionMismatch(grid1, grids...))
    return eachgridpoint(grid1)
end

"""$(TYPEDSIGNATURES) Obtain ring index `j` from gridpoint `ij` and `rings`
describing rind indices as obtained from `eachring(::Grid)`"""
function whichring(ij::Integer, rings::AbstractVector)
    @boundscheck 0 < ij <= rings[end][end] || throw(BoundsError)
    j = 1
    @inbounds while ij > rings[j][end]
        j += 1
    end
    return j
end