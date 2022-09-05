"""
    abstract type AbstractGrid{T} <: AbstractVector{T} end

The abstract supertype for all spatial grids on the sphere supported by SpeedyWeather.jl."""
abstract type AbstractGrid{T} <: AbstractVector{T} end

# all AbstractGrids have their grid points stored in a vector field `data`
# propagate length, size, getindex, setindex! for that
Base.length(G::AbstractGrid) = length(G.data)
Base.size(G::AbstractGrid) = size(G.data)
@inline function Base.getindex(G::AbstractGrid,k::Integer)
    @boundscheck 0 < k <= length(G.data) || throw(BoundsError(G,k))
    @inbounds r = G.data[k]
    return r
end

@inline Base.setindex!(G::AbstractGrid,x,k::Integer) = setindex!(G.data,x,k)

# with ranges
@inline Base.getindex(G::AbstractGrid,r::AbstractRange) = G.data[r]
@inline Base.setindex!(G::AbstractGrid,x::AbstractVector,r::AbstractRange) = setindex!(G.data,x,r)


# List all the methods required to be implemented for any Grid types here
# and define the  default behaviour as raising errors
# This would clarify what interface we need to implement when introducing a new Grid

truncation_order(::Type{G}) where {G<:AbstractGrid} =
    error("truncation_order(::Type{G}) not implemented for G==$G")
get_truncation(::Type{G},n::Integer) where {G<:AbstractGrid} =
    error("get_truncation(::Type{G},n::Integer) not implemented for G==$G")
get_resolution(::Type{G},trunc::Integer) where {G<:AbstractGrid} =
    error("get_resolution(::Type{G},trunc::Integer) not implemented for G==$G")
# define nlat_half for all grids (HEALPixGrid is different as it doesn't use nlat_half as resolution parameter)
get_nresolution(G::AbstractGrid) =
    error("get_nresolution(G) not implemented for G==$G")
# define whether there's an odd number of latitude rings for a grid
nlat_odd(::Type{G}) where {G<:AbstractGrid} =
    error("nlat_odd(::Type{G}) not implemented for G==$G")
# return the maxmimum number of longitude points for a grid and its resolution parameter nlat_half/nside
get_nlon(::Type{G},nresolution::Integer) where {G<:AbstractGrid} =
    error("get_nlon(::Type{G},nresolution::Integer) not implemented for G==$G")
# get the number of longitude points at given latitude ring j number 1...J from north to south
get_nlon_per_ring(::Type{G},nresolution::Integer,j::Integer) where {G<:AbstractGrid} =
    error("get_nlon_per_ring(::Type{G},nresolution::Integer,j::Integer) not implemented for G==$G")
# total number of grid points per grid type
get_npoints(::Type{G},nresolution::Integer) where {G<:AbstractGrid} =
    error("get_npoints(::Type{G},nresolution::Integer) not implemented for G==$G")
# colatitude [radians] vectors
get_colat(::Type{G},nresolution::Integer) where {G<:AbstractGrid} =
    error("get_colat(::Type{G},nresolution::Integer) not implemented for G==$G")
# lon [radians] vectors for full grids (empty vectors otherwise)
get_lon(::Type{G},nresolution::Integer) where {G<:AbstractGrid} =
    error("get_lon(::Type{G},nresolution::Integer) not implemented for G==$G")
# get coordinates
get_colatlons(::Type{G},nresolution::Integer) where {G<:AbstractGrid} =
    error("get_colatlons(::Type{G},nresolution::Integer) not implemented for G==$G")
each_index_in_ring(::Type{G},j::Integer,nresolution::Integer) where {G<:AbstractGrid} =
    error("each+index_in_ring(::Type{G},j::Integer,nresolution::Integer) not implemented for G==$G")
get_quadrature_weights(::Type{G},nresolution::Integer) where {G<:AbstractGrid} =
    error("get_quadrature_weights(::Type{G},nresolution::Integer) not implemented for G==$G")

# Define methods that are universally applicable to any G<:AbstractGrid here
# generator functions for grid
Base.zeros(::Type{G},n::Integer) where {G<:AbstractGrid{T}} where T =
    G(zeros(T,get_npoints(G,n)),n)
# use Float64 if not provided
Base.zeros(::Type{G},n::Integer) where {G<:AbstractGrid} = zeros(G{Float64},n)
# zero element of an AbstractGrid instance G by packing a zero(::Vector) into G
Base.zero(g::G) where {G<:AbstractGrid} = G(zero(g.data))

get_truncation(grid::G) where {G<:AbstractGrid} = get_truncation(G,get_nresolution(grid))
get_nlat_half(::G,n::Integer) where {G<:AbstractGrid} = get_nlat_half(G,n) 
nlat_odd(grid::AbstractGrid) = nlat_odd(typeof(grid))
get_nlat_half(::Type{<:AbstractGrid},nlat_half::Integer) = nlat_half
get_nlat(Grid::Type{<:AbstractGrid},n::Integer) = 2get_nlat_half(Grid,n) - nlat_odd(Grid)
get_nlat(grid::Grid) where {Grid<:AbstractGrid} = get_nlat(Grid,get_nresolution(grid))

function each_index_in_ring(grid::G,j::Integer) where {G<:AbstractGrid}
    return each_index_in_ring(G,j,get_nresolution(grid))
end

eachgridpoint(grid::G) where {G<:AbstractGrid} = Base.OneTo(get_npoints(G,get_nresolution(grid)))

function eachring(grid::G) where {G<:AbstractGrid}
    rings = [each_index_in_ring(grid,j) for j in 1:get_nlat(grid)]
    return rings                                            # return Vector{UnitRange}
end

function eachring(grids::Grid...) where {Grid<:AbstractGrid}
    for grid in grids
        @boundscheck length(grid) == length(grids[1]) || throw(BoundsError)
    end
    return eachring(grids[1])
end
