abstract type AbstractGridArray{T, N, ArrayType <: AbstractArray{T, N}} <: AbstractArray{T, N} end

# Horizontal 2D grids with N=1
const AbstractGrid{T} = AbstractGridArray{T, 1, Vector{T}}

## TYPES
nonparametric_type(grid::AbstractGridArray) = nonparametric_type(typeof(grid))
full_array_type(grid::AbstractGridArray) = full_array_type(typeof(grid))
full_grid_type(grid::AbstractGridArray) = horizontal_grid_type(full_array_type(grid))
horizontal_grid_type(grid::AbstractGridArray) = horizontal_grid_type(typeof(grid))

## SIZE
Base.length(G::AbstractGridArray) = length(G.data)
Base.size(G::AbstractGridArray) = size(G.data)
Base.sizeof(G::AbstractGridArray) = sizeof(G.data)

get_nlat_half(grid::AbstractGridArray) = grid.nlat_half
nlat_odd(grid::AbstractGridArray) = nlat_odd(typeof(grid))

# get total number of latitude rings, *(nlat_half > 0) to return 0 for nlat_half = 0
get_nlat(Grid::Type{<:AbstractGridArray}, nlat_half::Integer) = 2nlat_half - nlat_odd(Grid)*(nlat_half > 0)
get_nlat(grid::Grid) where {Grid<:AbstractGridArray} = get_nlat(Grid, grid.nlat_half)

# get total number of grid points vs horizontal (2D) only
get_npoints(grid::Grid) where {Grid<:AbstractGridArray} = get_npoints(Grid, grid.nlat_half)
get_npoints(G::Type{<:AbstractGridArray}, nlat_half::Integer, k::Integer...) = prod(k) * get_npoints2D(G, nlat_half)
get_npoints2D(grid::Grid) where {Grid<:AbstractGridArray} = get_npoints2D(Grid, grid.nlat_half)

# size of the matrix of the horizontal grid if representable as such (not all grids)
matrix_size(grid::Grid) where {Grid<:AbstractGridArray} = matrix_size(Grid, grid.nlat_half)

## INDEXING
@inline Base.getindex(G::AbstractGridArray, ijk::Integer...) = getindex(G.data,ijk...)
@inline Base.getindex(G::AbstractGridArray, r::AbstractRange, k::Integer...) = getindex(G.data, r, k...)
@inline Base.getindex(G::AbstractGridArray, ij::Integer, k::AbstractRange) = getindex(G.data, ij, k)
@inline Base.getindex(G::AbstractGridArray, I::CartesianIndex) = getindex(G, Tuple(I)...)

@inline function Base.getindex(
    G::GridArray,
    col::Colon,
    k...,
) where {GridArray<:AbstractGridArray{T, N, ArrayType}} where {T, N, ArrayType}
    GridArray_ = nonparametric_type(GridArray)
    return GridArray_{T, 1}(getindex(G.data, col, k...), G.nlat_half)
end

@inline Base.setindex!(G::AbstractGridArray, x, ijk::Integer...) =
    setindex!(G.data, x, ijk...)
@inline Base.setindex!(G::AbstractGridArray, x::AbstractVector, ij::AbstractRange, k::Integer...) =
    setindex!(G.data, x, ij, k...)
@inline Base.setindex!(G::AbstractGridArray, x::AbstractVector, ij::Integer, k::AbstractRange) =
    setindex!(G.data, x, ij, k)

## CONSTRUCTORS
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
    GridArray_ = nonparametric_type(Grid)
    rings = eachring(Grid, nlat_half)
    return GridArray_(data, nlat_half, rings)
end

# if no nlat_half provided calculate it
function (::Type{Grid})(data::AbstractArray) where {Grid<:AbstractGridArray}
    npoints2D = size(data, 1)
    nlat_half = get_nlat_half(Grid, npoints2D)
    return Grid(data, nlat_half)
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
    return Grid(ArrayType{T, N}(undef, get_npoints2D(Grid, nlat_half), k...), nlat_half)
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

## COORDINATES
get_latdlonds(grid::Grid) where {Grid<:AbstractGridArray} = get_latdlonds(Grid, grid.nlat_half)

function get_latdlonds(Grid::Type{<:AbstractGridArray}, nlat_half::Integer)
    colats, lons = get_colatlons(Grid, nlat_half)   # colatitudes, longitudes in radians
    latds = colats                                  # flat copy rename before conversion
    londs = lons
    latds .= π/2 .- colats                          # colatiudes to latitudes in radians
    latds .*= (360/2π)                              # now in degrees 90˚...-90˚
    londs .*= (360/2π)
    return latds, londs
end

function get_latlons(Grid::Type{<:AbstractGridArray}, nlat_half::Integer)
    colats, lons = get_colatlons(Grid, nlat_half)   # colatitudes, longitudes in radians
    lats = colats                                   # flat copy rename before conversion
    lats .= π/2 .- colats                           # colatiudes to latitudes in radians
    return lats, lons
end

get_lat(grid::Grid) where {Grid<:AbstractGridArray} = get_lat(Grid, grid.nlat_half)
get_latd(grid::Grid) where {Grid<:AbstractGridArray} = get_latd(Grid, grid.nlat_half)
get_lond(grid::Grid) where {Grid<:AbstractGridArray} = get_lond(Grid, grid.nlat_half)
get_lon(grid::Grid) where {Grid<:AbstractGridArray} = get_lon(Grid, grid.nlat_half)
get_colat(grid::Grid) where {Grid<:AbstractGridArray} = get_colat(Grid, grid.nlat_half)
get_colatlons(grid::Grid) where {Grid<:AbstractGridArray} = get_colatlons(Grid, grid.nlat_half)
get_nlon_max(grid::Grid) where {Grid<:AbstractGridArray} = get_nlon_max(Grid, grid.nlat_half)

function get_lat(Grid::Type{<:AbstractGridArray}, nlat_half::Integer)
    return π/2 .- get_colat(Grid, nlat_half)
end

function get_latd(Grid::Type{<:AbstractGridArray}, nlat_half::Integer)
    return get_lat(Grid, nlat_half) * (360/2π)
end

# only defined for full grids, empty vector as fallback
get_lon(::Type{<:AbstractGridArray}, nlat_half::Integer) = Float64[]
get_lond(::Type{<:AbstractGridArray}, nlat_half::Integer) = Float64[]

"""
$(TYPEDSIGNATURES)
Returns a vector `nlons` for the number of longitude points per latitude ring, north to south.
Provide grid `Grid` and its resolution parameter `nlat_half`. For both_hemisphere==false only
the northern hemisphere (incl Equator) is returned."""
function get_nlons(Grid::Type{<:AbstractGridArray}, nlat_half::Integer; both_hemispheres::Bool=false)
    n = both_hemispheres ? get_nlat(Grid, nlat_half) : nlat_half
    return [get_nlon_per_ring(Grid, nlat_half, j) for j in 1:n]
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

"""
$(TYPEDSIGNATURES)
Vector{UnitRange} `rings` to loop over every ring of grid `grid`
and then each grid point per ring. To be used like

    rings = eachring(grid)
    for ring in rings
        for ij in ring
            grid[ij]"""
@inline eachring(grid::AbstractGridArray) = grid.rings

function eachring(Grid::Type{<:AbstractGridArray}, nlat_half::Integer)
    rings = Vector{UnitRange{Int}}(undef, get_nlat(Grid, nlat_half))
    each_index_in_ring!(rings, Grid, nlat_half)
    return rings
end

"""
$(TYPEDSIGNATURES)
Same as `eachring(grid)` but performs a bounds check to assess that all grids
in `grids` are of same size."""
function eachring(grid1::Grid, grids::Grid...) where {Grid<:AbstractGridArray}
    @inline 
    n = length(grid1)
    Base._all_match_first(X->length(X), n, grid1, grids...) || throw(BoundsError)
    return eachring(grid1)
end

function grids_match(A::AbstractGridArray, B::AbstractGridArray)
    length(A) == length(B) && return grids_match(typeof(A), typeof(B))
    return false
end

function grids_match(A::Type{<:AbstractGridArray}, B::Type{<:AbstractGridArray})
    return nonparametric_type(A) == nonparametric_type(B)
end

"""
$(TYPEDSIGNATURES)
UnitRange to access data on grid `grid` on ring `j`."""
function each_index_in_ring(grid::Grid, j::Integer) where {Grid<:AbstractGridArray}
    return each_index_in_ring(Grid, j, grid.nlat_half)
end

"""
$(TYPEDSIGNATURES)
UnitRange to access each grid point on grid `grid`."""
eachgridpoint(grid::AbstractGridArray) = Base.OneTo(get_npoints(grid))
function eachgridpoint(grid1::Grid, grids::Grid...) where {Grid<:AbstractGridArray}
    n = length(grid1)
    Base._all_match_first(X->length(X), n, grid1, grids...) || throw(BoundsError)
    return eachgridpoint(grid1)
end

"""
$(TYPEDSIGNATURES)
Obtain ring index j from gridpoint ij and Vector{UnitRange} describing rind indices
as obtained from eachring(::Grid)"""
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
import Base.Broadcast: BroadcastStyle, Broadcasted

# {1} as grids are <:AbstractVector, Grid here is the non-parameteric Grid type!
struct AbstractGridArrayStyle{N, Grid} <: Broadcast.AbstractArrayStyle{N} end

# important to remove Grid{T} parameter T (eltype/number format) here to broadcast
# automatically across the same grid type but with different T
# e.g. FullGaussianGrid{Float32} and FullGaussianGrid{Float64}
Base.BroadcastStyle(::Type{Grid}) where {Grid<:AbstractGridArray{T, N, ArrayType}} where {T, N, ArrayType} =
    AbstractGridArrayStyle{N, nonparametric_type(Grid)}()

# allocation for broadcasting, create a new Grid with undef of type/number format T
function Base.similar(bc::Broadcasted{AbstractGridArrayStyle{N, Grid}}, ::Type{T}) where {N, Grid, T}
    return Grid(Array{T}(undef,size(bc)...))
end

# ::Val{0} for broadcasting with 0-dimensional, ::Val{1} for broadcasting with vectors, etc
AbstractGridArrayStyle{N, Grid}(::Val{M}) where {N, Grid, M} = AbstractGridArrayStyle{N, Grid}()