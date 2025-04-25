"""Abstract supertype for all arrays of ring grids, representing `N`-dimensional
data on the sphere in two dimensions (but unravelled into a vector in the first dimension,
the actual "ring grid") plus additional `N-1` dimensions for the vertical and/or time etc.
Parameter `T` is the `eltype` of the underlying data, held as in the array type `ArrayType`
(Julia's `Array` for CPU or others for GPU).

Ring grids have several consecuitive grid points share the same latitude (= a ring),
grid points on a given ring are equidistant. Grid points are ordered 0 to 360˚E,
starting around the north pole, ring by ring to the south pole. """
abstract type AbstractGridArray{T, N, ArrayType <: AbstractArray{T, N}} <: AbstractArray{T, N} end

"""Abstract supertype for all ring grids, representing 2-dimensional data on the
sphere unravelled into a Julia `Vector`. Subtype of `AbstractGridArray` with
`N=1` and `ArrayType=Vector{T}` of `eltype T`."""
const AbstractGrid{T} = AbstractGridArray{T, 1, Vector{T}}

## TYPES
"""$(TYPEDSIGNATURES) For any instance of `AbstractGridArray` type its n-dimensional type
(*Grid{T, N, ...} returns *Array) but without any parameters `{T, N, ArrayType}`"""
nonparametric_type(grid::AbstractGridArray) = nonparametric_type(typeof(grid))

# also needed for other array types, defined in extensions
nonparametric_type(::Type{<:Array}) = Array

# needed for unalias
@inline Base.dataids(grid::AbstractGridArray) = Base.dataids(grid.data)

"""$(TYPEDSIGNATURES) Full grid array type for `grid`. Always returns the N-dimensional `*Array`
not the two-dimensional (`N=1`) `*Grid`. For reduced grids the corresponding full grid that
share the same latitudes."""
full_array_type(grid::AbstractGridArray) = full_array_type(typeof(grid))

"""$(TYPEDSIGNATURES) Full (horizontal) grid type for `grid`. Always returns the two-dimensional
(`N=1`) `*Grid` type. For reduced grids the corresponding full grid that share the same latitudes."""
full_grid_type(grid::AbstractGridArray) = horizontal_grid_type(full_array_type(grid))
full_grid_type(Grid::Type{<:AbstractGridArray}) = horizontal_grid_type(full_array_type(Grid))

"""$(TYPEDSIGNATURES) The two-dimensional (`N=1`) `*Grid` for `grid`, which can be an N-dimensional
`*GridArray`."""
horizontal_grid_type(grid::AbstractGridArray) = horizontal_grid_type(typeof(grid))

## SIZE
Base.length(G::AbstractGridArray) = length(G.data)  # total number of grid points
Base.size(G::AbstractGridArray) = size(G.data)      # size of underlying data (horizontal is unravelled)

"""$(TYPEDSIGNATURES) Size of underlying data array plus precomputed ring indices."""
Base.sizeof(G::AbstractGridArray) = sizeof(G.data) + sizeof(G.rings)

"""$(TYPEDSIGNATURES) Resolution paraemeters `nlat_half` of a `grid`.
Number of latitude rings on one hemisphere, Equator included."""
get_nlat_half(grid::AbstractGridArray) = grid.nlat_half

"""$(TYPEDSIGNATURES) True for a `grid` that has an odd number of
latitude rings `nlat` (both hemispheres)."""
nlat_odd(grid::AbstractGridArray) = nlat_odd(typeof(grid))

# get total number of latitude rings, *(nlat_half > 0) to return 0 for nlat_half = 0
get_nlat(Grid::Type{<:AbstractGridArray}, nlat_half::Integer) = 2nlat_half - nlat_odd(Grid)*(nlat_half > 0)

"""$(TYPEDSIGNATURES) Get number of latitude rings, pole to pole."""
get_nlat(grid::Grid) where {Grid<:AbstractGridArray} = get_nlat(Grid, grid.nlat_half)

"""$(TYPEDSIGNATURES) Total number of grid points in all dimensions of `grid`.
Equivalent to length of the underlying data array."""
get_npoints(grid::Grid) where {Grid<:AbstractGridArray} = get_npoints(Grid, grid.nlat_half)
get_npoints(G::Type{<:AbstractGridArray}, nlat_half::Integer, k::Integer...) = prod(k) * get_npoints2D(G, nlat_half)

"""$(TYPEDSIGNATURES) Number of grid points in the horizontal dimension only,
even if `grid` is N-dimensional."""
get_npoints2D(grid::Grid) where {Grid<:AbstractGridArray} = get_npoints2D(Grid, grid.nlat_half)

"""$(TYPEDSIGNATURES) Size of the matrix of the horizontal grid if representable as such (not all grids)."""
matrix_size(grid::Grid) where {Grid<:AbstractGridArray} = matrix_size(Grid, grid.nlat_half)

## INDEXING
# simply propagate all indices forward
Base.@propagate_inbounds Base.getindex(G::AbstractGridArray, ijk...) = getindex(G.data, ijk...)

@inline function Base.getindex(
    G::GridArray,
    col::Colon,
    k...,
) where {GridArray<:AbstractGridArray}
    GridArray_ = nonparametric_type(GridArray)  # obtain parameters from G.data
    return GridArray_(getindex(G.data, col, k...), G.nlat_half, G.rings)
end

# simply propagate all indices forward
Base.@propagate_inbounds Base.setindex!(G::AbstractGridArray, x, ijk...) = setindex!(G.data, x, ijk...)
Base.fill!(G::AbstractGridArray, x) = fill!(G.data, x)

## CONSTRUCTORS
"""$(TYPEDSIGNATURES) True for `data`, `nlat_half` and `rings` that all match in size
to construct a grid of type `Grid`."""
function check_inputs(data, nlat_half, rings, Grid)
    check = true
    check &= size(data, 1) == get_npoints2D(Grid, nlat_half)    # test number of 2D grid points
    check &= length(rings) == get_nlat(Grid, nlat_half)         # test number of rings == nlat
    # TODO also check that rings map to all and only valid grid points?
    return check
end

function error_message(data, nlat_half, rings, G, T, N, A)
    nlat = get_nlat(G, nlat_half)
    nrings = length(rings)
    if nlat != nrings
        return error("$nrings-element ring indices "*
            "cannot be used to create a $nlat-ring $G{$T, $N, $A}.")
    else
        return error("$(summary(data)) cannot be used to create a $nlat-ring $G{$T, $N, $A}")
    end
end

# if no rings are provided, calculate them
function (::Type{Grid})(data::AbstractArray, nlat_half::Integer) where {Grid<:AbstractGridArray}
    GridArray_ = nonparametric_type(Grid)   # strip away parameters of type, obtain parameters from data
    rings = eachring(Grid, nlat_half)       # precompute indices to access the variable-length rings
    return GridArray_(data, nlat_half, rings)
end

# if no nlat_half provided calculate it
(::Type{Grid})(M::AbstractArray; input_as=Vector) where Grid<:AbstractGridArray = Grid(M, input_as)

function (::Type{Grid})(data::AbstractArray, input_as::Type{Vector}) where Grid<:AbstractGridArray
    npoints2D = size(data, 1)                   # from 1st dim of data
    nlat_half = get_nlat_half(Grid, npoints2D)  # get nlat_half of Grid
    return Grid(data, nlat_half)
end

function (::Type{Grid})(data::AbstractArray, input_as::Type{Matrix})  where Grid<:AbstractGridArray
    error("Only full grids can be created from matrix input")
end

for f in (:zeros, :ones, :rand, :randn)
    @eval begin
        # general version with ArrayType(zeros(...)) conversion
        function Base.$f(
            ::Type{Grid},
            nlat_half::Integer,
            k::Integer...,
        ) where {Grid<:AbstractGridArray{T, N, ArrayType}} where {T, N, ArrayType}
            return Grid(ArrayType($f(T, get_npoints2D(Grid, nlat_half), k...)), nlat_half)
        end

        # CPU version with zeros(T, ...) producing Array
        function Base.$f(
            ::Type{Grid},
            nlat_half::Integer,
            k::Integer...,
        ) where {Grid<:AbstractGridArray{T}} where T
            return Grid($f(T, get_npoints2D(Grid, nlat_half), k...), nlat_half)
        end

        # use Float64 if no type provided
        function Base.$f(
            ::Type{Grid},
            nlat_half::Integer,
            k::Integer...,
        ) where {Grid<:AbstractGridArray}
            return $f(Grid{Float64}, nlat_half, k...)
        end
    end
end

# zero element of an AbstractGridArray instance grid by creating new zero(grid.data)
Base.zero(grid::Grid) where {Grid<:AbstractGridArray} =
    nonparametric_type(Grid)(zero(grid.data), grid.nlat_half, grid.rings)

# similar data but everything else identical
function Base.similar(grid::Grid) where {Grid<:AbstractGridArray}
    return nonparametric_type(Grid)(similar(grid.data), grid.nlat_half, grid.rings)
end

# data with new type T but everything else identical
function Base.similar(grid::Grid, ::Type{T}) where {Grid<:AbstractGridArray, T}
    return nonparametric_type(Grid)(similar(grid.data, T), grid.nlat_half, grid.rings)
end

# data with same type T but new size
function Base.similar(
    grid::Grid,
    nlat_half::Integer,
    k::Integer...
) where {Grid<:AbstractGridArray{T, N, ArrayType}} where {T, N, ArrayType}
    similar_data = similar(grid.data, get_npoints2D(Grid, nlat_half), k...)
    return nonparametric_type(Grid)(similar_data, nlat_half)
end

# data with new type T and new size
function Base.similar(
    grid::Grid,
    ::Type{Tnew},
    nlat_half::Integer,
    k::Integer...
) where {Grid<:AbstractGridArray, Tnew}
    similar_data = similar(grid.data, Tnew, get_npoints2D(Grid, nlat_half), k...)
    return nonparametric_type(Grid)(similar_data, nlat_half)
end

# general version with ArrayType{T, N}(undef, ...) generator
function (::Type{Grid})(
    ::UndefInitializer,
    nlat_half::Integer,
    k::Integer...,
) where {Grid<:AbstractGridArray{T, N, ArrayType}} where {T, N, ArrayType}
    ArrayType_ = nonparametric_type(ArrayType)
    return Grid(ArrayType_{T, N}(undef, get_npoints2D(Grid, nlat_half), k...), nlat_half)
end

# CPU version with Array{T, N}(undef, ...) generator
function (::Type{Grid})(
    ::UndefInitializer,
    nlat_half::Integer,
    k::Integer...,
) where {Grid<:AbstractGridArray{T}} where T
    return Grid(Array{T}(undef, get_npoints2D(Grid, nlat_half), k...), nlat_half)
end

# use Float64 if no type provided
function (::Type{Grid})(
    ::UndefInitializer,
    nlat_half::Integer,
    k::Integer...,
) where {Grid<:AbstractGridArray}
    return Grid(Array{Float64}(undef, get_npoints2D(Grid, nlat_half), k...), nlat_half)
end

function Base.convert(
    ::Type{Grid},
    grid::AbstractGridArray,
) where {Grid<:AbstractGridArray{T, N, ArrayType}} where {T, N, ArrayType}
    return Grid(ArrayType(grid.data))
end

## COORDINATES

"""$(TYPEDSIGNATURES) Longitudes (degrees, 0-360˚E), latitudes (degrees, 90˚N to -90˚N) for
every (horizontal) grid point in `grid` in ring order (0-360˚E then north to south)."""
get_londlatds(grid::Grid) where {Grid<:AbstractGridArray} = get_londlatds(Grid, grid.nlat_half)

"""$(TYPEDSIGNATURES) Longitudes (radians, 0-2π), latitudes (degrees, π/2 to -π/2) for
every (horizontal) grid point in `grid` in ring order (0-360˚E then north to south)."""
get_lonlats(grid::Grid) where {Grid<:AbstractGridArray} = get_lonlats(Grid, grid.nlat_half)

"""$(TYPEDSIGNATURES) Longitudes (radians, 0-2π), colatitudes (degrees, 0 to π) for
every (horizontal) grid point in `grid` in ring order (0-360˚E then north to south)."""
get_loncolats(grid::Grid) where {Grid<:AbstractGridArray} = get_loncolats(Grid, grid.nlat_half)

# radians
function get_lonlats(Grid::Type{<:AbstractGridArray}, nlat_half::Integer)
    londs, latds = get_londlatds(Grid, nlat_half)  # longitudes, latitudes in degrees
    return londs * π/180, latds * π/180             # to radians
end

# radians and colatitudes
function get_loncolats(Grid::Type{<:AbstractGridArray}, nlat_half::Integer)
    lons, lats = get_lonlats(Grid, nlat_half)
    colats = π/2 .- lats
    return lons, colats
end

"""$(TYPEDSIGNATURES) Latitude (radians) for each ring in `grid`, north to south."""
get_lat(grid::Grid) where {Grid<:AbstractGridArray} = get_lat(Grid, grid.nlat_half)

"""$(TYPEDSIGNATURES) Latitude (degrees) for each ring in `grid`, north to south."""
get_latd(grid::Grid) where {Grid<:AbstractGridArray} = get_latd(Grid, grid.nlat_half)

"""$(TYPEDSIGNATURES) Colatitudes (radians) for each ring in `grid`, north to south."""
get_colat(grid::Grid) where {Grid<:AbstractGridArray} = get_colat(Grid, grid.nlat_half)

"""$(TYPEDSIGNATURES) Longitude (degrees). Full grids only."""
get_lond(grid::Grid) where {Grid<:AbstractGridArray} = get_lond(Grid, grid.nlat_half)

"""$(TYPEDSIGNATURES) Longitude (radians). Full grids only."""
get_lon(grid::Grid) where {Grid<:AbstractGridArray} = get_lon(Grid, grid.nlat_half)

get_nlon_max(grid::Grid) where {Grid<:AbstractGridArray} = get_nlon_max(Grid, grid.nlat_half)

get_lat(Grid::Type{<:AbstractGridArray}, nlat_half::Integer) = get_latd(Grid, nlat_half) * (π/180)
get_colat(Grid::Type{<:AbstractGridArray}, nlat_half::Integer) = π/2 .- get_lat(Grid, nlat_half)

"""
$(TYPEDSIGNATURES)
Returns a vector `nlons` for the number of longitude points per latitude ring, north to south.
Provide grid `Grid` and its resolution parameter `nlat_half`. For keyword argument
`both_hemispheres=false` only the northern hemisphere (incl Equator) is returned."""
function get_nlons(Grid::Type{<:AbstractGridArray}, nlat_half::Integer; both_hemispheres::Bool=true)
    n = both_hemispheres ? get_nlat(Grid, nlat_half) : nlat_half
    return [get_nlon_per_ring(Grid, nlat_half, j) for j in 1:n]
end

"""$(TYPEDSIGNATURES)
Number of longitude points per latitude ring `j`."""
function get_nlon_per_ring(grid::AbstractGridArray, j::Integer)
    return get_nlon_per_ring(typeof(grid), grid.nlat_half, j)
end

## ITERATORS
"""
$(TYPEDSIGNATURES)
CartesianIndices for the 2nd to last dimension of an AbstractGridArray,
to be used like

for k in eachgrid(grid)
    for ring in eachring(grid)
        for ij in ring
            grid[ij, k]"""
@inline eachgrid(grid::AbstractGridArray) = CartesianIndices(size(grid)[2:end])

# several arguments to check for matching grids
function eachgrid(grid1::AbstractGridArray, grids::AbstractGridArray...; kwargs...)
    grids_match(grid1, grids...; kwargs...) || throw(DimensionMismatch(grid1, grids...))
    return eachgrid(grid1)
end

"""
$(TYPEDSIGNATURES)
Vector{UnitRange} `rings` to loop over every ring of grid `grid`
and then each grid point per ring. To be used like

    rings = eachring(grid)
    for ring in rings
        for ij in ring
            grid[ij]

Accesses precomputed `grid.rings`."""
@inline eachring(grid::AbstractGridArray) = grid.rings

"""$(TYPEDSIGNATURES)
Computes the ring indices `i0:i1` for start and end of every longitudinal point
on a given ring `j` of `Grid` at resolution `nlat_half`. Used to loop
over rings of a grid. These indices are also precomputed in every `grid.rings`."""
function eachring(Grid::Type{<:AbstractGridArray}, nlat_half::Integer)
    rings = Vector{UnitRange{Int}}(undef, get_nlat(Grid, nlat_half))    # allocate
    each_index_in_ring!(rings, Grid, nlat_half)                         # calculate iteratively
    return rings
end

"""$(TYPEDSIGNATURES) Same as `eachring(grid)` but performs a bounds check to assess
that all `grids` according to `grids_match` (non-parametric grid type, nlat_half and length)."""
function eachring(grid1::AbstractGridArray, grids::AbstractGridArray...)
    # for eachring grids only need to match in the horizontal, but can different vertical (or other) dimensions
    grids_match(grid1, grids...; horizontal_only = true) || throw(DimensionMismatch(grid1, grids...))
    return eachring(grid1)
end

function Base.DimensionMismatch(grid1::AbstractGridArray, grids::AbstractGridArray...)
    s = "AbstractGridArrays do not match; $(size(grid1)) $(nonparametric_type(grid1))"
    for grid in grids
        s *= ", $(size(grid))-$(nonparametric_type(grid))"
    end
    return DimensionMismatch(s)
end

# equality and comparison, somehow needed as not covered by broadcasting
Base.:(==)(G1::AbstractGridArray, G2::AbstractGridArray) = grids_match(G1, G2) && G1.data == G2.data
Base.all(G::AbstractGridArray) = all(G.data)
Base.any(G::AbstractGridArray) = any(G.data)

"""$(TYPEDSIGNATURES) True if both `A` and `B` are of the same nonparametric grid type
(e.g. OctahedralGaussianArray, regardless type parameter `T` or underyling array type `ArrayType`)
and of same resolution (`nlat_half`) and total grid points (`length`). Sizes of `(4,)` and `(4,1)`
would match for example, but `(8,1)` and `(4,2)` would not (`nlat_half` not identical)."""
function grids_match(
    A::AbstractGridArray,
    B::AbstractGridArray;
    horizontal_only::Bool = false,
    vertical_only::Bool = false,
)
    @assert ~(horizontal_only && vertical_only) "Conflicting options: horizontal_only = $horizontal_ony and vertical_only = $vertical_only"

    horizontal_match = get_nlat_half(A) == get_nlat_half(B)
    vertical_match = size(A)[2:end] == size(B)[2:end]
    type_match = grids_match(typeof(A), typeof(B))

    if horizontal_only
        # type also has to match as two different grid types can have the same nlat_half
        return horizontal_match && type_match
    elseif vertical_only
        return vertical_match
    else
        return horizontal_match && vertical_match && type_match
    end
end

# eltypes can be different and also array types of underlying data
function grids_match(A::Type{<:AbstractGridArray}, B::Type{<:AbstractGridArray})
    return nonparametric_type(A) == nonparametric_type(B)
end

"""$(TYPEDSIGNATURES) True if all grids `A, B, C, ...` provided as arguments
match according to `grids_match` wrt to `A` (and therefore all)."""
function grids_match(A::AbstractGridArray, B::AbstractGridArray...; kwargs...)
    match = true    # single grid A always matches itself
    for Bi in B     # check for all matching respectively with A
        match &= grids_match(A, Bi; kwargs...)
    end
    return match
end

"""$(TYPEDSIGNATURES) UnitRange to access data on grid `grid` on ring `j`."""
function each_index_in_ring(grid::Grid, j::Integer) where {Grid<:AbstractGridArray}
    return grid.rings[j]
end

""" $(TYPEDSIGNATURES) UnitRange to access each horizontal grid point on grid `grid`.
For a `NxM` (`N` horizontal grid points, `M` vertical layers) `OneTo(N)` is returned."""
eachgridpoint(grid::AbstractGridArray) = Base.OneTo(get_npoints(grid))

""" $(TYPEDSIGNATURES) Like `eachgridpoint(::AbstractGridArray)` but checks for
equal size between input arguments first."""
function eachgridpoint(grid1::Grid, grids::Grid...) where {Grid<:AbstractGridArray}
    n = length(grid1)
    # TODO check only nonparametric_type and nlat_half identical!
    Base._all_match_first(X->length(X), n, grid1, grids...) || throw(BoundsError)
    return eachgridpoint(grid1)
end

"""$(TYPEDSIGNATURES) Obtain ring index `j` from gridpoint `ij` and `rings`
describing rind indices as obtained from `eachring(::Grid)`"""
function whichring(ij::Integer, rings::Vector{UnitRange{Int}})
    @boundscheck 0 < ij <= rings[end][end] || throw(BoundsError)
    j = 1
    @inbounds while ij > rings[j][end]
        j += 1
    end
    return j
end

## BROADCASTING
# following https://docs.julialang.org/en/v1/manual/interfaces/#man-interfaces-broadcasting
import Base.Broadcast: BroadcastStyle, Broadcasted, DefaultArrayStyle

# {1} as grids are <:AbstractVector, Grid here is the non-parameteric Grid type!
struct AbstractGridArrayStyle{N, Grid} <: Broadcast.AbstractArrayStyle{N} end

# important to remove Grid{T} parameter T (eltype/number format) here to broadcast
# automatically across the same grid type but with different T
# e.g. FullGaussianGrid{Float32} and FullGaussianGrid{Float64}
Base.BroadcastStyle(::Type{Grid}) where {Grid<:AbstractGridArray{T, N, ArrayType}} where {T, N, ArrayType} =
    AbstractGridArrayStyle{N, nonparametric_type(Grid)}()

# allocation for broadcasting, create a new Grid with undef of type/number format T
function Base.similar(bc::Broadcasted{AbstractGridArrayStyle{N, Grid}}, ::Type{T}) where {N, Grid, T}
    return Grid(Array{T}(undef, size(bc)))
end

# ::Val{0} for broadcasting with 0-dimensional, ::Val{1} for broadcasting with vectors, etc
# when there's a dimension mismatch always choose the larger dimension
AbstractGridArrayStyle{N, Grid}(::Val{N}) where {N, Grid} = AbstractGridArrayStyle{N, Grid}()
AbstractGridArrayStyle{1, Grid}(::Val{2}) where {Grid} = AbstractGridArrayStyle{2, Grid}()
AbstractGridArrayStyle{1, Grid}(::Val{0}) where {Grid} = AbstractGridArrayStyle{1, Grid}()
AbstractGridArrayStyle{2, Grid}(::Val{3}) where {Grid} = AbstractGridArrayStyle{3, Grid}()
AbstractGridArrayStyle{2, Grid}(::Val{1}) where {Grid} = AbstractGridArrayStyle{2, Grid}()
AbstractGridArrayStyle{3, Grid}(::Val{4}) where {Grid} = AbstractGridArrayStyle{4, Grid}()
AbstractGridArrayStyle{3, Grid}(::Val{2}) where {Grid} = AbstractGridArrayStyle{3, Grid}()

## GPU
struct AbstractGPUGridArrayStyle{N, ArrayType, Grid} <: GPUArrays.AbstractGPUArrayStyle{N} end

function Base.BroadcastStyle(
    ::Type{Grid}
) where {Grid<:AbstractGridArray{T, N, ArrayType}} where {T, N, ArrayType <: GPUArrays.AbstractGPUArray}
    return AbstractGPUGridArrayStyle{N, ArrayType, nonparametric_type(Grid)}()
end

# ::Val{0} for broadcasting with 0-dimensional, ::Val{1} for broadcasting with vectors, etc
# when there's a dimension mismatch always choose the larger dimension
AbstractGPUGridArrayStyle{N, ArrayType, Grid}(::Val{N}) where {N, ArrayType, Grid} =
    AbstractGPUGridArrayStyle{N, ArrayType, Grid}()

AbstractGPUGridArrayStyle{1, ArrayType, Grid}(::Val{2}) where {ArrayType, Grid} = AbstractGPUGridArrayStyle{2, ArrayType, Grid}()
AbstractGPUGridArrayStyle{1, ArrayType, Grid}(::Val{0}) where {ArrayType, Grid} = AbstractGPUGridArrayStyle{1, ArrayType, Grid}()
AbstractGPUGridArrayStyle{2, ArrayType, Grid}(::Val{3}) where {ArrayType, Grid} = AbstractGPUGridArrayStyle{3, ArrayType, Grid}()
AbstractGPUGridArrayStyle{2, ArrayType, Grid}(::Val{1}) where {ArrayType, Grid} = AbstractGPUGridArrayStyle{2, ArrayType, Grid}()
AbstractGPUGridArrayStyle{3, ArrayType, Grid}(::Val{4}) where {ArrayType, Grid} = AbstractGPUGridArrayStyle{4, ArrayType, Grid}()
AbstractGPUGridArrayStyle{3, ArrayType, Grid}(::Val{2}) where {ArrayType, Grid} = AbstractGPUGridArrayStyle{3, ArrayType, Grid}()

function KernelAbstractions.get_backend(
    g::Grid
) where {Grid <: AbstractGridArray{T, N, ArrayType}} where {T, N, ArrayType <: GPUArrays.AbstractGPUArray}
    return KernelAbstractions.get_backend(g.data)
end

function Base.similar(
    bc::Broadcasted{AbstractGPUGridArrayStyle{N, ArrayType, Grid}},
    ::Type{T},
) where {N, ArrayType, Grid, T}
    ArrayType_ = nonparametric_type(ArrayType)
    return Grid(ArrayType_{T}(undef, size(bc)))
end

function Adapt.adapt_structure(to, grid::Grid) where {Grid <: AbstractGridArray}
    Grid_ = nonparametric_type(Grid)
    return Grid_(Adapt.adapt(to, grid.data), grid.nlat_half, grid.rings)
end